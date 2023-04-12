from ase import Atoms
#from numba import jit
import shutil
import scipy
from ase.neighborlist import NeighborList, natural_cutoffs
from numpy.linalg import norm
from itertools import combinations
import itertools
import numpy as np
from default_func import get_atomic_number
from default_func import get_radicial,get_neighbor_list,get_thermal_correction_from_vaspkit_output,get_gas_thermal_correction_from_vaspkit_output
from ase.io import read,write
from struc_to_graph import atom_to_graph
import networkx as nx
import os
import networkx.algorithms.isomorphism as iso


bond_match = iso.categorical_edge_match('bond','')  #graph要比较的属性，属性的默认值
ads_match = iso.categorical_node_match(['symbol'], [-1, False])

c_form=-8.6789 
o_form=-7.3178 
h_form=-3.4587 
n_form=-4.2458 


def add_h(ad_graph_1,ad_graph_2): #用图论判断加氢
        flag=0
        node_2=[i for i in ad_graph_2.nodes]
        for j in range(1,len(node_2)):
            for k in itertools.combinations(node_2,j):
                sub_g_1=ad_graph_2.subgraph(list(k))
                if nx.is_isomorphic(ad_graph_1,sub_g_1,edge_match=bond_match,node_match=ads_match):
                    tmpnode=[node for node in sub_g_1.nodes]  #加氢前物种的原子序号
                    tmp_gra_1=ad_graph_2.copy()
                    for k1 in tmpnode:
                        tmp_gra_1.remove_node(k1)   #依次给加氢后的物种脱掉 加氢前物种的原子序号

                    ele = [i for i in tmp_gra_1._node.values()]
                    if len(ele)==1:
                        if ele[0]['symbol'] == 'H':
                            flag = 1
                            break
        return flag


def drop_oh(ad_graph_1,ad_graph_2):  #图论判断脱水
    flag=0
    node_1=[i for i in ad_graph_1.nodes]
    for j in range(1,len(node_1)):
        for k in itertools.combinations(node_1,j):
            sub_g_1=ad_graph_1.subgraph(list(k))
            if nx.is_isomorphic(ad_graph_2,sub_g_1,edge_match=bond_match,node_match=ads_match):
                tmpnode=[node for node in sub_g_1.nodes]
                tmp_gra_1=ad_graph_1.copy()
                for k1 in tmpnode:
                    tmp_gra_1.remove_node(k1)
                #if nx.is_connected(tmp_gra_1) and sorted([node[0] for node in tmp_gra_1.nodes])==["H","O"]:
                #        flag=1
                ele = [i for i in tmp_gra_1._node.values()]
                ele_set = set()
                if len(ele)==2:
                    ele_set.add(ele[0]['symbol'])
                    ele_set.add(ele[1]['symbol'])
                else:
                    break         
                if nx.is_connected(tmp_gra_1) and ('H' in ele_set) and ('O' in ele_set):
                    flag = 1        
    return flag


N_species_path = '/share/home/panjj/graph_catalysis/code0907/configuration/N_species'
N_species = []
for file in os.listdir(N_species_path):
    N_atoms = read(N_species_path+'/'+file)
    N_neighbor_list = get_neighbor_list(N_atoms,0.001,1.25)
    full,ad_graph=atom_to_graph(N_atoms,N_neighbor_list,ad_atoms=[])
    tmp=[file[:-5],full]
    N_species.append(tmp)

def judge_coupling_3(ad_graph_1,ad_graph_2):  #判断耦合反应，返回1表示能够耦合，N_part表示耦合反应的氮源
    flag=0
    N_part=0
    N_name=0
    node_2=[i for i in ad_graph_2.nodes]
    for j in range(1,len(node_2)):
        for k in itertools.combinations(node_2,j):
            sub_g_1=ad_graph_2.subgraph(list(k))
            if nx.is_isomorphic(ad_graph_1,sub_g_1,edge_match=bond_match,node_match=ads_match):
                tmpnode=[node for node in sub_g_1.nodes]  #加氢前物种的原子序号
                tmp_gra_1=ad_graph_2.copy()
                for k1 in tmpnode:
                    tmp_gra_1.remove_node(k1)   #依次给加氢后的物种脱掉 加氢前物种的原子序号

                for N_graph in N_species:
                    if iso.is_isomorphic(tmp_gra_1,N_graph[1],edge_match=bond_match,node_match=ads_match) and nx.is_connected(tmp_gra_1):
                        flag = 1
                        N_part = N_graph[1]
                        N_name = N_graph[0]
                        break
    return flag,N_part,N_name


def get_reaction_network2(total_file):
    co_couple_path='/share/home/panjj/graph_catalysis/CN2_product/NEB/couple_energy'
    co_couple_free_energy_dict={}
    for file in os.listdir(co_couple_path):
        co_couple_energy=0
        if os.path.exists(co_couple_path+"/"+file+"/scf/output.relax"):
            for line in open(co_couple_path+"/"+file+"/scf/output.relax"):
                    if "1 F" in line:
                        co_couple_energy=float(line.split()[4])
        else:
            for line in open(co_couple_path+"/"+file+"/scf/output.scf"):
                    if "1 F" in line:
                        co_couple_energy=float(line.split()[4])
        thermal_correction=get_thermal_correction_from_vaspkit_output(co_couple_path+"/"+file)
        co_couple_free_energy=co_couple_energy+thermal_correction
        couple_species=file.split('_')
        couple_reactant_plus_product=couple_species[0]+'_'+couple_species[2]
        co_couple_free_energy_dict[couple_reactant_plus_product]=  co_couple_free_energy

    file_list = []
    reaction_network=nx.DiGraph()
    reaction_count=0
    for i in total_file:
        #print("i",i)
        ad_graph_1=i[1]  #i[0]0是文件名, 3是能量
        ad_graph_edge_1=sorted([edge for edge in ad_graph_1.edges])
        #basis_graph_1=i[2]
        #basis_graph_edge_1=deal_basis_h(sorted([edge for edge in basis_graph_1.edges]))
        for k in total_file:
            if i!=k:
                ad_graph_2=k[1]

                flag=add_h(ad_graph_1,ad_graph_2)
                #print(flag)
                if flag==1:
                    det_e=k[2]-i[2]+3.4587 
                    #print("{}+H>{} det_e={}".format(i[0], k[0],det_e))  
                    tmpedge=[(i[0],k[0],{"deltE":det_e,"type":"add_H"})]
                    file_list.append('{},{},{},add_H'.format(i[0],k[0],det_e))
                    reaction_network.add_edges_from(tmpedge)
                    reaction_count+=1


                if flag==0:
                    flag=drop_oh(ad_graph_1,ad_graph_2)
                    if flag==1:
                        det_e=k[2]-i[2]-10.77647575 
                        #print("{}+H>{}+H2O   det_e={}".format(i[0], k[0],det_e))
                        tmpedge=[(i[0],k[0],{"deltE":det_e,"type":"drop_H2O"})]
                        file_list.append('{},{},{},drop_H2O'.format(i[0],k[0],det_e))
                        reaction_network.add_edges_from(tmpedge)
                        reaction_count+=1

                      
                if flag==0:
                    flag,couple_part,N_name=judge_coupling_3(ad_graph_1,ad_graph_2)
                    if flag==1:
                        couple_rea_and_pro=i[0]+'_'+k[0]
                        if couple_rea_and_pro not in co_couple_free_energy_dict:
                            continue 
                        co_couple_name=i[0]+'_'+N_name

                        gas_fomation_e=0
                        gas_node=[i[2] for i in couple_part.nodes]
                        for node in gas_node:
                            if node=="N":
                                gas_fomation_e+=n_form
                            if node=="O":
                                gas_fomation_e+=o_form
                            if node=="H":
                                gas_fomation_e+=h_form   
                        #det_e=k[2]-i[2]-gas_fomation_e
                        #print("{}+{}>{} couple_ad det_e={}".format(i[0],N_name,k[0],det_e))
                        #print(l[0])
                        
                        det_e1=co_couple_free_energy_dict[couple_rea_and_pro]-i[2]-gas_fomation_e
                        det_e2=k[2]-co_couple_free_energy_dict[couple_rea_and_pro]
                        tmpedge1=[(i[0],co_couple_name,{"deltE":det_e1,"type":"couple with {}".format(N_name)})]
                        file_list.append('{},{},{},couple_with {}'.format(i[0],co_couple_name,det_e1,N_name))
                        tmpedge2=[(co_couple_name,k[0],{"deltE":det_e2,"type":"couple"})]
                        file_list.append('{},{},{},couple'.format(co_couple_name,k[0],det_e2))                        
                        reaction_network.add_edges_from(tmpedge1)
                        reaction_network.add_edges_from(tmpedge2)
                        reaction_count+=1

    print(reaction_count)
    return reaction_network,file_list
