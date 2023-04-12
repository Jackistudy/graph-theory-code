import os
from re import S
from creat_ad_structure import get_double_bond_atom,get_cross_bond_atom,get_single_bond_atom
from creat_ad_structure import single_ad_structure,double_ad_structure,double_cross_ad_structure
from get_bond_order import get_smx_bond_order
import numpy as np
from ase.io import read,write
import numpy as np
import networkx as nx
from rdkit import Chem
#from package.structure_to_graph import atom_to_graph
from struc_to_graph import atom_to_graph
import networkx.algorithms.isomorphism as iso
from default_func import get_neighbor_list,smile_to_graph,smile_to_formula,get_tri_inter_central
#from ase.data import covalent_radii

obpath='/home/obpath'

slab=read('/home/CONTCAR')
smile_file = "/home/smile_file.txt"  

site=[47,48]

def compare_graph(compare_structure,tmp):   #图同构判断
    flag = 1
    neighbor_list = get_neighbor_list(tmp,0.001,1.25)
    mole_graph, chem_ad = atom_to_graph(tmp, neighbor_list, ad_atoms=[])
    bond_match = iso.categorical_edge_match('bond','')
    ads_match = iso.categorical_node_match('symbol', "")
    if compare_structure:
        for each in compare_structure:
            if iso.is_isomorphic(mole_graph,each,edge_match=bond_match,node_match=ads_match):
                flag = 0
                break
    if flag:
        compare_structure.append(mole_graph)
        return 1
    else:
        return 0            


total_count=0
#smiles = ["NC([H])(O)N([H])O"]

smiles=[]
f = open(smile_file,'r')      #读取txt里面的所有中间产物
for line in f.readlines():
    line=line.strip('\n')
    smiles.append(line)
f.close()


N_num=0
Cu_num=0
bridge_num=0
double_num=0
double_cross_num=0

for file in smiles:
       #保存“分子图”节点序号对应的元素符号，在get_graph生成元素信息
    structure=[]
    qiao_structure=[]
    compare_structure = []
    count=0
    mol = smile_to_graph(file)  #生成“分子图”     
    formula=smile_to_formula(file)  #smile转化为“分子式”，用于文件命名
    #if not os.path.exists(obpath + "/" +formula ):
    #    os.mkdir(obpath + "/" + formula)
    single_site=get_single_bond_atom(mol)
    double_site=get_double_bond_atom(mol)
    double_cross_site=get_cross_bond_atom(mol)
    
    
    for j in single_site:  #位点1
            flag=get_smx_bond_order(mol,mole_bond_atom=[j])
            if flag==1:
                tmp1=single_ad_structure(slab,mol,ad_site=[site[0]],mole_bond_atom=[j])
                #if compare_graph(compare_structure,tmp1):
                structure.append(tmp1) 
                N_num+=1
    for j in single_site:  #位点1
            flag=get_smx_bond_order(mol,mole_bond_atom=[j])
            if flag==1:
                tmp2=single_ad_structure(slab,mol,ad_site=[site[1]],mole_bond_atom=[j])   
                #if compare_graph(compare_structure,tmp2):
                structure.append(tmp2)
                Cu_num+=1  
                          
    for j in single_site: #桥位
            flag=get_smx_bond_order(mol,mole_bond_atom=[j],bond_to_basis=2)
            if flag==1:
                tmp=single_ad_structure(slab,mol,ad_site=site,mole_bond_atom=[j])    
                if not os.path.exists(obpath + "/" +formula ):
                    os.mkdir(obpath + "/" + formula)                 
                #if compare_graph(compare_structure,tmp):
                structure.append(tmp)
                bridge_num+=1
    
    for j in double_site: #双氮位点1
            flag = get_smx_bond_order(mol,mole_bond_atom=j)
            if flag == 1:
                tmp=double_ad_structure(slab,mol,ad_site=site,mole_bond_atom=j)    
                #if compare_graph(compare_structure,tmp):
                structure.append(tmp)
                double_num+=1

        
    for j in double_cross_site:
            flag = get_smx_bond_order(mol,mole_bond_atom=j)
            if flag == 1:
                tmp=double_cross_ad_structure(slab,mol,ad_site=site,mole_bond_atom=j,type="normal")       
                #if compare_graph(compare_structure,tmp):
                structure.append(tmp)
                double_cross_num+=1    
   
    for i in structure:
        count+=1
        total_count+=1
        name=str(count)+"-"+formula
        write(obpath+"/"+formula+"/"+"{}.vasp".format(name),i)



for file in os.listdir(obpath):
    for subfile in os.listdir(obpath+"/"+file):
        tmpf = open(obpath + "/" + file+"/"+subfile, "r")
        elementlist = tmpf.readlines()[0]
        fp = open(obpath + "/" + file+"/"+subfile, "r")
        s = fp.read()
        fp.close()
        tmp = s.split('\n')
        tmp.insert(5, elementlist.strip())
        # print(tmp)
        s = "\n".join(tmp)
        fp = open(obpath + "/" + file+"/"+subfile, "w")
        fp.write(s)
        fp.close()
  
