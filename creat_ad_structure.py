import numpy as np
import os
import math
import networkx as nx
from ase.io import read
from ase.io import write
from ase import neighborlist
import networkx.algorithms.isomorphism as iso
from ase.data import covalent_radii
import shutil
from numpy.linalg import norm
from itertools import combinations
from mole_operation import dou_ad_entension,sin_ad_extension,sin_neg_extension,mole_extension, get_base_position
from default_func import get_neighbor_list,periodic_deal,get_radii
from ase import Atom,Atoms


def get_single_bond_atom(mol):   #判断单位点吸附的原子
    element=nx.get_node_attributes(mol,"symbol")
    C_N_list = []
    O_list = []              
    mole_bond_atom=[]      

    for key,value in element.items():
        if value == "C" or value == "N":
            C_N_list.append(key)
        if value == "O":
            O_list.append(key)
    for index in C_N_list:
        if mol.degree(index) <= 3:
            mole_bond_atom.append(index)
    for index in O_list:
        if mol.degree(index) <= 1:
            mole_bond_atom.append(index)  
    return mole_bond_atom  

def get_double_bond_atom(mol):    #判断双位点吸附的原子
    element=nx.get_node_attributes(mol,"symbol")
    C_N_list = []
    O_list = []
    mole_bond_atom=[]
    for key,value in element.items():
        if value == "C" or value == "N":
            C_N_list.append(key)
        if value == "O":
            O_list.append(key)    
    for index1 in C_N_list:
        if mol.degree(index1) <= 3:
            for index2 in C_N_list:
                if index2 != index1 and mol.degree(index2) <= 3:
                    #if ([index2,index1] not in mole_bond_atom):
                    if ((index1,index2) in mol.edges) or ((index2,index1) in mol.edges): 
                        mole_bond_atom.append([index1,index2])
            for index2 in O_list:
                if mol.degree(index2) <= 1:
                    if ((index1,index2) in mol.edges) or ((index2,index1) in mol.edges): 
                        mole_bond_atom.append([index1,index2])
                        mole_bond_atom.append([index2,index1])
    #for pair in mole_bond_atom:
    #    print(pair,":",(element[pair[0]],element[pair[1]]))                
    return mole_bond_atom

def get_cross_bond_atom(mol):     #判断双位点（中间隔着一个原子）吸附的原子
    element=nx.get_node_attributes(mol,"symbol")
    C_N_list = []
    O_list = []
    mole_bond_atom=[]
    for key,value in element.items():
        if value == "C" or value == "N":
            C_N_list.append(key)
        if value == "O":
            O_list.append(key) 
    for index1 in C_N_list:
        if mol.degree(index1) <= 3:
            index1_neighbor = mol[index1]
            for key, _ in index1_neighbor.items():
                if mol.degree(key) > 1:
                    index1_neighbor_2 = mol[key]
                    for index2, _ in index1_neighbor_2.items():
                        if index2 != index1:
                            if element[index2]=="C" or element[index2]=="N":
                                if mol.degree(index2) <= 3 and ([index2,index1] not in mole_bond_atom):
                                    mole_bond_atom.append([index1,index2])
                                    mole_bond_atom.append([index2,index1])

                            if element[index2]=="O":
                                if mol.degree(index2) <= 1:   
                                    mole_bond_atom.append([index1,index2]) 
                                    mole_bond_atom.append([index2,index1]) 
    for index1 in O_list:                                
        if mol.degree(index1) == 1:
            index1_neighbor = mol[index1]
            for key, _ in index1_neighbor.items():
                if mol.degree(key) > 1:
                    index1_neighbor_2 = mol[key]
                    for index2, _ in index1_neighbor_2.items():
                        if index2 != index1:
                            if element[index2]=="O":
                                if mol.degree(index2) <= 1 and ([index2,index1] not in mole_bond_atom):   
                                    mole_bond_atom.append([index1,index2])
                                    mole_bond_atom.append([index2,index1])  
    #for pair in mole_bond_atom:
    #    print(pair,":",(element[pair[0]],element[pair[1]]))                                 
    return mole_bond_atom


def single_ad_structure(slab,mol,ad_site,mole_bond_atom=[],type="auto"):  #如果和吸附原子一个原子相连，该原子就要竖直朝上，如果两个，就是要z方向分布
    element=nx.get_node_attributes(mol,"symbol")
    slab=slab.copy()  #必须copy
    new_mol = mol.copy()
    positions = {}
    for i in list(new_mol.nodes):   #初始化坐标信息
        positions[i]=(0,0,0)
    


    #print(element)
    if type=="auto":
        if len(ad_site)==1:
            R=get_radii(slab[ad_site[0]].symbol)+get_radii(element[mole_bond_atom[0]])
            #positions[mole_bond_atom[0]]=tuple(slab[ad_site[0]].position+np.array([0,0,R]))
            pos_ad_site_neigh=[slab[i].position for i in ad_site ]
    
            positions[mole_bond_atom[0]]=tuple(get_base_position(pos_ad_site_neigh,[R*0.9]))
            #positions[mole_bond_atom[0]]=tuple(get_base_position(pos_ad_site_neigh,[R]))
        elif len(ad_site)==2:
            #site1_coor=slab[ad_site[0]].position
            #print(site1)
            pos_ad_site_neigh=[slab[i].position for i in ad_site ]
            a=np.linalg.norm(slab[ad_site[0]].position-slab[ad_site[1]].position)
            R1=get_radii(slab[ad_site[0]].symbol)
            R2=get_radii(slab[ad_site[1]].symbol)
            R3=get_radii(element[mole_bond_atom[0]])
            R=[(R2+R3)*0.9,(R1+R3)*0.9]

            positions[mole_bond_atom[0]]= tuple(get_base_position(pos_ad_site_neigh,R))   #计算inter atom的位置,001应该确实不行？？？？？
    else:
        R=-1  # 人为指定
        positions[mole_bond_atom[0]]=tuple(slab[ad_site[0]].position+np.array([0,0,R]))

    total_branches = list(nx.bfs_successors(new_mol, mole_bond_atom[0]))

    if len(total_branches[0][1]) != 0:  #中心原子所连的原子不等于零

        sin_ad_extension(positions, element, total_branches[0])   #延展单位点吸附的函数，双位点，就只能延伸两个方向。
        center_atom = total_branches[0][0]
        graph_copy=new_mol.copy()
        #print("reduced_mole_graph_copy:",reduced_mole_graph_copy)
        
        if len(total_branches[1:]) > 0:  # 第一层以外的处理方式
            if len(total_branches[1][1]) == 2:
                #print("nodes1:",reduced_mole_graph_copy.nodes,"edges1:",reduced_mole_graph_copy.edges)
                graph_copy.remove_edge(*[total_branches[1][0],total_branches[0][0]])
                #print("nodes2:",reduced_mole_graph_copy.nodes,"edges2:",reduced_mole_graph_copy.edges)
                bran_tmp0=list(nx.bfs_successors(graph_copy, total_branches[1][1][0]))
                #print("bran_tmp0:",bran_tmp0)
                bran_tmp1 = list(nx.bfs_successors(graph_copy, total_branches[1][1][1]))
                #print("bran_tmp1:",bran_tmp1)
                if len(bran_tmp0[0][1])>len(bran_tmp1[0][1]):  #延展方向的改变
                    total_branches[1][1][1]=bran_tmp0[0][0]
                    total_branches[1][1][0] = bran_tmp1[0][0]
                else:
                    total_branches[1][1][1] = bran_tmp1[0][0]
                    total_branches[1][1][0] = bran_tmp0[0][0]
        
        for branch in total_branches[1:]: #  非中央位点的处理方式
            mole_extension(positions,element, branch, center_atom)
        #write("tmp2.vasp", atoms)

    atom_list=[]
    for i in list(new_mol.nodes):
        a = Atom(element[i],positions[i])
        atom_list.append(a)
    atoms = Atoms(atom_list)
    slab += atoms           #atoms这个类是否必须？调用ase.Atoms这个类可以把原子符号和对应的原子坐标写进去
    #total_structure.append(slab)
    #print("slab2:",slab)
    return slab


def double_ad_structure(slab,mol,ad_site,mole_bond_atom=[],type="auto"):
    element=nx.get_node_attributes(mol,"symbol")
    new_mol = mol.copy()  #对图进行复制，否则报错，直接用XX.copy()也行
    slab = slab.copy()   
    positions = {}    
    for i in list(new_mol.nodes):   #初始化坐标信息
        positions[i]=(0,0,0)
    r_bond=[]
    r_site=[]
    
    for i in mole_bond_atom:
        r_bond.append(get_radii(element[i]))

    for i in ad_site:
        r_site.append(get_radii(slab[i].symbol))

    #print("rs",r_site)
    
    #for i in [[0,1],[1,0]]:
    if type=="auto":
            #new_site_a = slab[ad_site[0]].position + np.array([0, 0, r_site[0]+r_bond[0]])
            #new_site_b = slab[ad_site[1]].position + np.array([0, 0, r_site[1]+r_bond[1]])    #与N原子相连的吸附原子的坐标
            new_site_a=get_base_position([slab[ad_site[0]].position],[(r_site[0]+r_bond[0])*0.9])
            #new_site_a=get_base_position([slab[ad_site[0]].position],[(r_site[0]+r_bond[0])])
            new_site_b=get_base_position([slab[ad_site[1]].position],[(r_site[1]+r_bond[1])*0.9])
            #new_site_b=get_base_position([slab[ad_site[1]].position],[(r_site[1]+r_bond[1])])
    else:
            #print("222")
            H=1.5  # 手动调整
            new_site_a = slab[ad_site[0]].position + np.array([0, 0, H])
            new_site_b = slab[ad_site[1]].position + np.array([0, 0, H])    #与N原子相连的吸附原子的坐标
        #print("2222")
    vec_site = new_site_b - new_site_a  # 一个位点指向另一个位点的矢量
    len_vec = np.linalg.norm(vec_site)  #范数
    uvec0 = vec_site / len_vec  #单位方向矢量
    d = np.sum(r_bond)
    dn = (d - len_vec) / 2  # 这个操作整体上是个什么意义？相比于原本的那两个矢量，做了微小的移动  ,调整位置
    base_position0 = new_site_a - uvec0 * dn    
    base_position1 = new_site_b + uvec0 * dn
        
    positions[mole_bond_atom[0]] = tuple(base_position0)
    positions[mole_bond_atom[1]] = tuple(base_position1)

    #print("base,position",base_position0,base_position1)   #这个地方对两个基本位置进行了一种处理

    if tuple(mole_bond_atom) in new_mol.edges:
            new_mol.remove_edge(*mole_bond_atom)

    uvec1  = np.array([[0,0,1]])  #这个不应该指定的
    uvec2 = np.cross(uvec1, uvec0)         #叉乘
    uvec2 = uvec2/np.linalg.norm(uvec2)
    uvec1 = -np.cross(uvec2, uvec0)
    #print("uc1", uvec1)
    

    branches0 = list(nx.bfs_successors(new_mol, mole_bond_atom[0]))
    #("brance0",branches0,branches0[1:])
    if len(branches0[0][1]) != 0:
            uvec = [-uvec0, uvec1[0], uvec2[0]]  #可能是这个地方出了问题。
            dou_ad_entension(positions,element, uvec, branches0[0])
            #print("first branch",atoms.positions)
            #write("tmp1.vasp",atoms)
            root = branches0[0][0]

            graph_copy = new_mol.copy()
            
            if len(branches0[1:])>0:
                if len(branches0[1][1]) == 2:
                    graph_copy.remove_edge(*[branches0[1][0], branches0[0][0]])
                    bran_tmp0 = list(nx.bfs_successors(graph_copy, branches0[1][1][0]))
                    bran_tmp1 = list(nx.bfs_successors(graph_copy, branches0[1][1][1]))
                    if len(bran_tmp0[0][1]) > len(bran_tmp1[0][1]):
                        branches0[1][1][1] = bran_tmp0[0][0]
                        branches0[1][1][0] = bran_tmp1[0][0]
                    else:
                        branches0[1][1][1] = bran_tmp1[0][0]
                        branches0[1][1][0] = bran_tmp0[0][0]
            
            for branch in branches0[1:]:
                mole_extension(positions,element, branch, root)


    branches1 = list(nx.bfs_successors(new_mol, mole_bond_atom[1]))
        #print("brance1", branches1,branches1[1:])
    if len(branches1[0][1]) != 0:
            uvec = [uvec0, uvec1[0], uvec2[0]]
            dou_ad_entension(positions,element, uvec, branches1[0])
            #print("second branch", atoms.positions)
            #write("tmp2.vasp", atoms)
            root = branches1[0][0]

            graph_copy = new_mol.copy()
            
            if len(branches1[1:]) > 0:
                if len(branches1[1][1]) == 2:
                    #reduced_mole_graph_copy.remove_edge(*[branches1[1][0], branches0[0][0]])
                    graph_copy.remove_edge(*[branches1[1][0], branches1[0][0]])
                    bran_tmp0 = list(nx.bfs_successors(graph_copy, branches1[1][1][0]))
                    bran_tmp1 = list(nx.bfs_successors(graph_copy, branches1[1][1][1]))
                    if len(bran_tmp0[0][1]) > len(bran_tmp1[0][1]):
                        branches1[1][1][1] = bran_tmp0[0][0]
                        branches1[1][1][0] = bran_tmp1[0][0]
                    else:
                        branches1[1][1][1] = bran_tmp1[0][0]
                        branches1[1][1][0] = bran_tmp0[0][0]
    
            for branch in branches1[1:]:
                mole_extension(positions,element, branch, root)

        # print("afteratom",atoms.positions)
        
    atom_list=[]
    for i in list(new_mol.nodes):
            a = Atom(element[i],positions[i])    #将处理后的“分子图”的元素以及坐标信息写到列表
            atom_list.append(a)
    atoms = Atoms(atom_list)        #将元素、原子坐标转化为ase包的分子图对象
        #print("after",atoms.positions)
    slab += atoms
        #total_structure.append(slab)
    return slab

def double_cross_ad_structure(slab,mol,ad_site,mole_bond_atom=[],type="normal"):
    element=nx.get_node_attributes(mol,"symbol")
    slab = slab.copy()
    new_mol = mol.copy()
    positions = {}
    for i in list(new_mol.nodes):   #初始化坐标信息
        positions[i]=(0,0,0)    
    vec_site_0=slab[ad_site[0]].position-slab[ad_site[1]].position

    r_bond=[]
    r_site=[]
    
    for i in mole_bond_atom:
        r_bond.append(get_radii(element[i]))

    for i in ad_site:
        r_site.append(get_radii(slab[i].symbol))

    inter_atom_index = 0
    connec_list_of_bond = []
    for i in mole_bond_atom:
        tmp = []
        for j in new_mol.edges:
            if i in j:
                tmp.append(j[0])
                tmp.append(j[1])
        connec_list_of_bond.append(tmp)
    for i in connec_list_of_bond[0]:
        for j in connec_list_of_bond[1]:
            if i == j:
                inter_atom_index = i
    
    #for i in [[0,1],[1,0]]:
    if type=="normal":
            #new_site_a = slab[ad_site[0]].position + np.array([0, 0, r_site[0]+r_bond[0]])
            #new_site_b = slab[ad_site[1]].position + np.array([0, 0, r_site[1]+r_bond[1]])    #与N原子相连的吸附原子的坐标
            new_site_a=get_base_position([slab[ad_site[0]].position],[(r_site[0]+r_bond[0])*0.9])
            #new_site_a=get_base_position([slab[ad_site[0]].position],[(r_site[0]+r_bond[0])])
            new_site_b=get_base_position([slab[ad_site[1]].position],[(r_site[1]+r_bond[1])*0.9])
            #new_site_b=get_base_position([slab[ad_site[1]].position],[(r_site[1]+r_bond[1])])
    else:
            #print("222")
            H=1.5  # 手动调整
            new_site_a = slab[ad_site[0]].position + np.array([0, 0, H])
            new_site_b = slab[ad_site[1]].position + np.array([0, 0, H])    #与N原子相连的吸附原子的坐标  
    vec_site = new_site_b- new_site_a  # 一个位点指向另一个位点的矢量
    len_vec = np.linalg.norm(vec_site)  # 范数
    uvec0 = vec_site / len_vec  # 单位方向矢量


    d = 2.1  # 可调整
    dn = (d - len_vec) / 2  
    base_position0 = new_site_a - uvec0 * dn
    base_position1 = new_site_b + uvec0 * dn

    positions[mole_bond_atom[0]] = tuple(base_position0)
    positions[mole_bond_atom[1]] = tuple(base_position1)
        #write("tmp.vasp",atoms)
    
    uvec1 = np.array([[0, 0, 1]])  
    uvec2 = np.cross(uvec1, uvec0)  
    uvec2 = uvec2/np.linalg.norm(uvec2)  #新增
    b=r_bond[0]+get_radii(element[inter_atom_index])
    c=r_bond[1]+get_radii(element[inter_atom_index])
    positions[inter_atom_index]= get_base_position([base_position1,base_position0],[b,c])  #计算inter atom的位置
    new_mol.remove_edge(*[mole_bond_atom[0],inter_atom_index])
    new_mol.remove_edge(*[mole_bond_atom[1], inter_atom_index])

        #write("tmp0.vasp", atoms)
    uvec1=np.cross(uvec0,uvec2)
    uvec1 = uvec1/np.linalg.norm(uvec1)  #新增

    branches0 = list(nx.bfs_successors(new_mol, mole_bond_atom[0]))
        #print("brance0", branches0, branches0[1:])
    if len(branches0[0][1]) != 0:
            uvec = [-uvec0, uvec1[0], uvec2[0]]  # 可能是这个地方出了问题。
            dou_ad_entension(positions,element, uvec, branches0[0])
            root = branches0[0][0]
            #print("before:",positions)
            for branch in branches0[1:]:
                mole_extension(positions,element, branch, root)
            #print("after:",positions)        

    branches1 = list(nx.bfs_successors(new_mol, mole_bond_atom[1]))
        #print("brance1", branches1, branches1[1:])
    if len(branches1[0][1]) != 0:
            uvec = [uvec0, uvec1[0], uvec2[0]]
            dou_ad_entension(positions,element, uvec, branches1[0])
            root = branches1[0][0]
            #print("before:",positions)
            for branch in branches1[1:]:
                mole_extension(positions,element, branch, root)
            #print("after:",positions)    

    branches2=list(nx.bfs_successors(new_mol, inter_atom_index))
    if len(branches2[0][1]) != 0:
            sin_ad_extension(positions,element, branches2[0],cross="y",vec=uvec2[0])
            root = branches2[0][0]
            #print("before:",positions)
            for branch in branches2[1:]:
                mole_extension(positions,element, branch, root)
            #print("after:",positions)

    
    atom_list=[]
    for i in list(new_mol.nodes):
            #print(positions)
            a = Atom(element[i],list(positions[i]))
            atom_list.append(a)
    atoms = Atoms(atom_list)      
        #print("after", atoms.positions)
    slab += atoms
        #total_structure.append(slab)
    return slab


def single_ad_neg_structure(slab, mol, ad_site, mole_bond_atom=[],neigh_atom=[]): #专门处理基底加氢的反
    element=nx.get_node_attributes(mol,"symbol")
    total_structure = []   
    slab = slab.copy()
    new_mol = mol.copy()
    neigh_position=[slab.positions[i] for i in neigh_atom]
    vec0 = neigh_position[1] - neigh_position[0]  # 一个位点指向另一个位点的矢量
    vec1=[0,0,1]
    vec=np.cross(vec0,vec1)
    len_vec = np.linalg.norm(vec)  # 范数
    # print("n222", n)
    vec = vec / len_vec  # 单位方向矢量
    positions = {}
    for i in list(new_mol.nodes):   #初始化坐标信息
        positions[i]=(0,0,0)
    for index,atom in enumerate(slab):
        if atom.index==ad_site[0]:        
            slab[index].position=slab[index].position+np.array([0,0,0.4])

    positions[mole_bond_atom[0]] = slab[ad_site[0]].position
    #slab.pop(ad_site[0])
    
    branches0 = list(nx.bfs_successors(new_mol, mole_bond_atom[0]))

    if len(branches0[0][1]) != 0:
        sin_neg_extension(positions,element, branches0[0],vec)
        root = branches0[0][0]


    
       
   
    atom_list=[]
    for i in list(new_mol.nodes):
            if element[i]!="N":
            #print(positions)
                a = Atom(element[i],list(positions[i]))
                atom_list.append(a)
    atoms = Atoms(atom_list)      
        #print("after", atoms.positions)
    slab += atoms
        #total_structure.append(slab)
    return slab

