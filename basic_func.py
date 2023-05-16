import  re
import numpy as np
from ase.data import chemical_symbols as sym
from ase.io import read
from ase import Atoms
import os
from ase.neighborlist import NeighborList, natural_cutoffs
from ase import neighborlist
import networkx.algorithms.isomorphism as iso
import networkx as nx
from networkx import Graph, MultiGraph
import copy
from rdkit import Chem


radicial={1:1,6:4,7:4,8:2}  # the number of bonds that can be formed, this can be change in specific reaction
#print(radicial)

def compare_graph(compare_structure,tmp):   #图同构判断
    flag = 1
    neighbor_list = get_neighbor_list(tmp,0.001,1.25)
    mole_graph, chem_ad = atom_to_graph(tmp, neighbor_list, ad_atoms=[])
    bond_match = iso.categorical_edge_match('bond','')
    ads_match = iso.categorical_node_match('symbol', "")
    if compare_structure:
        for each in compare_structure:
            if iso.is_isomorphic(mole_graph,each,edge_match=bond_match,node_match=ads_match):  #同构
                flag = 0
                break
    if flag:
        compare_structure.append(mole_graph)
        return 1
    else:
        return 0 

def atom_to_graph(structure,neigh_list,ad_atoms=[]):
    total_graph=nx.Graph()
    for index,atom in enumerate(structure):
        #print(index)
        node_name=str(index)+":"+atom.symbol
        #print(node_name)
        total_graph.add_node(node_name,symbol=atom.symbol)
        neighbors, offsets = neigh_list.get_neighbors(index)
        for i in neighbors: 
                node_name_1=str(i)+":"+structure[i].symbol
                if node_name!=node_name_1:
                    bond_type="".join(sorted(atom.symbol+structure[i].symbol))
                    total_graph.add_edge(node_name,node_name_1,bond=bond_type)

    ad_node=[]
    if len(ad_atoms)!=0:
        for i in ad_atoms:
            sym=structure[i].symbol
            ad_node.append(str(i)+":"+sym)
    
    ad_graph=nx.subgraph(total_graph,ad_node)
    return total_graph,ad_graph

def get_neighbor_list(atoms,skin,cutoff=1): #返回分子图的邻接矩阵信息
    radii = neighborlist.natural_cutoffs(atoms,mult=cutoff)  #这个radii可以认为设置把，根据需要，未必采用ase的方法  获得近邻原子
    #print(radii)
    neighbor_list = neighborlist.NeighborList(radii,skin,bothways=True)
    neighbor_list.update(atoms)
    return neighbor_list


def get_radii(element): #原子的共价半径
    r_dict = {'C':0.76,'N':0.71,'O':0.66,'H':0.31,"Cu":1.32,"Ni":1.24} 
    return r_dict[element]


def get_thermal_correction_from_vaspkit_output(path): #自由能校正
    correction=0
    if os.path.exists(path+"/gibs"):
            os.chdir(path+"/gibs")
            os.system("echo -e \"501\n298.15\n\" | vaspkit >result.txt")
            if os.path.exists(path+"/gibs/result.txt"):
                zpe,e=0,0
                for line in open(path+"/gibs/result.txt"):
                    if "energy E_ZPE" in line:
                        zpe=float(line.split()[-2])
                    if "Entropy S" in line:
                        s=float(line.split()[-2])*298.15
                correction=zpe-s
    if correction==0:
        print(path)
        print("correction wrong")
        return correction
    else:
        return correction
def get_gas_thermal_correction_from_vaspkit_output(path):
    correction=0
    if os.path.exists(path+"/gibs"):
            os.chdir(path+"/gibs")
            os.system("echo -e \"502\n298.15\n1\n1\n\" | vaspkit >result.txt")
            if os.path.exists(path+"/gibs/result.txt"):
                zpe,e=0,0
                for line in open(path+"/gibs/result.txt"):
                    if "energy E_ZPE" in line:
                        zpe=float(line.split()[-2])
                    if "Entropy S" in line:
                        s=float(line.split()[-2])*298.15
                correction=zpe-s
    if correction==0:
        print("correction wrong")
        return correction
    else:
        return correction


def smile_to_graph(smiles):  #将smiles变成分子图
    structure = Chem.MolFromSmiles(smiles, sanitize=False)
    node_num = 0
    #print(type(structure))
    mol = nx.Graph()
    node_num=0
    for atom in structure.GetAtoms():
        node_num_1=atom.GetIdx()
        node_symbol_1=atom.GetSymbol()
        mol.add_node(node_num_1,symbol=node_symbol_1)
        for x in atom.GetNeighbors():
            node_num_2=x.GetIdx()
            node_symbol_2=x.GetSymbol()
            mol.add_edge(node_num_1,node_num_2)
    #print(nx.get_node_attributes(mol,"symbol"))
    return mol

        
def smile_to_formula(Intermediate):   #把smile变成“分子式”
    N_num = 0 
    N_list = []
    str_C = ""
    structure = Chem.MolFromSmiles(Intermediate, sanitize=False)
    for atom in structure.GetAtoms():           
        if atom.GetSymbol() == "C":
            C_H_num = 0
            C_O_num = 0         
            for x in atom.GetNeighbors():
                if x.GetSymbol() == "H":      #如果C旁边有H
                    C_H_num += 1
                if x.GetSymbol() == "O":      #如果C旁边有O
                    C_O_num += 1
            if C_H_num == 0:
                str_C_H = "C"
            elif C_H_num == 1:
                str_C_H = "CH"
            elif C_H_num == 2:
                str_C_H = "CH2"
            elif C_H_num == 3:
                str_C_H = "CH3"  
            elif C_H_num == 4:
                str_C_H = "CH4"                                    
            if C_O_num == 0:
                    str_C_O = ""
            elif C_O_num == 1:
                for x in atom.GetNeighbors():
                    if x.GetSymbol() == "O": 
                        for y in x.GetNeighbors():
                            if y.GetSymbol() == "H": 
                                str_C_O = "OH" 
                                break
                            else:
                                str_C_O = "O"
            elif C_O_num == 2:
                C_O_H_num = 0
                for x in atom.GetNeighbors():
                    if x.GetSymbol() == "O": 
                        for y in x.GetNeighbors():
                            if y.GetSymbol() == "H": 
                                C_O_H_num += 1 
                if C_O_H_num == 0:
                    str_C_O = "OO"
                elif C_O_H_num == 1:
                    str_C_O = "OOH"
                elif C_O_H_num == 2:
                    str_C_O = "OHOH" 
            str_C = str_C_H + str_C_O 
                    
        if atom.GetSymbol() == "N":
            N_num += 1
            N_H_num = 0
            N_O_num = 0
            for x in atom.GetNeighbors():     
                if x.GetSymbol() == "H":      #如果N旁边有H
                    N_H_num += 1
                if x.GetSymbol() == "O":      #如果N旁边有O
                    N_O_num += 1    
            if N_H_num == 0:
                str_N_H = "N"
            elif N_H_num == 1:
                str_N_H = "NH"
            elif N_H_num == 2:
                str_N_H = "NH2"
            elif N_H_num == 3:
                str_N_H = "NH3"                
            if N_O_num == 0:
                    str_N_O = ""
            elif N_O_num == 1:
                for x in atom.GetNeighbors():
                    if x.GetSymbol() == "O": 
                        for y in x.GetNeighbors():
                            if y.GetSymbol() == "H": 
                                str_N_O = "OH" 
                                break
                            else:
                                str_N_O = "O"
            elif N_O_num == 2:
                N_O_H_num = 0
                for x in atom.GetNeighbors():
                    if x.GetSymbol() == "O": 
                        for y in x.GetNeighbors():
                            if y.GetSymbol() == "H": 
                                N_O_H_num += 1 
                if N_O_H_num == 0:
                    str_N_O = "OO"
                elif N_O_H_num == 1:
                    str_N_O = "OOH"
                elif N_O_H_num == 2:
                    str_N_O = "OHOH" 
            str_N = str_N_H + str_N_O
            N_list.append(str_N)                                                                                                                                                                                                                   
    if len(N_list)==0:
        if str_C == "COO":
            return "CO2"
        else:
            return str_C
    elif len(N_list)==1:
        if str_C == "" and N_list[0] == "NOO":
            return "NO2"
        else:
            return str_C + N_list[0]        
    elif len(N_list)==2:
        N_list.sort(key=len)    
        return str_C + N_list[0] + N_list[1]


#path=r'C:\Users\dell\Desktop\T\test\CONTCAR.vasp'

#smile=["OCO"]
#smile_to_graph_new(smile)
