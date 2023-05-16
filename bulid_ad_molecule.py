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
from geometry_operation import dou_ad_entension,sin_ad_extension,sin_neg_extension,mole_extension, get_base_position 
from basic_func import get_neighbor_list,get_radii
from ase import Atom,Atoms


class ad_modes:
    def __init__(self):
        pass
    def single_bond_index(mol):   #判断单位点吸附的原子
        element=nx.get_node_attributes(mol,"symbol")
        C_N_list = []
        O_list = []              
        mole_bond_index=[]      

        for key,value in element.items():
            if value == "C" or value == "N":
                C_N_list.append(key)
            if value == "O":
                O_list.append(key)
        for index in C_N_list:
            if mol.degree(index) <= 3:
                mole_bond_index.append(index)
        for index in O_list:
            if mol.degree(index) <= 1:
                mole_bond_index.append(index)  
        return mole_bond_index  

    def double_bond_index(mol):    #判断双位点吸附的原子
        element=nx.get_node_attributes(mol,"symbol")
        C_N_list = []
        O_list = []
        mole_bond_index=[]
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
                            mole_bond_index.append([index1,index2])
                for index2 in O_list:
                    if mol.degree(index2) <= 1:
                        if ((index1,index2) in mol.edges) or ((index2,index1) in mol.edges): 
                            mole_bond_index.append([index1,index2])
                            mole_bond_index.append([index2,index1])
        #for pair in mole_bond_atom:
        #    print(pair,":",(element[pair[0]],element[pair[1]]))                
        return mole_bond_index

    def cross_bond_index(mol):     #判断双位点（中间隔着一个原子）吸附的原子
        element=nx.get_node_attributes(mol,"symbol")
        C_N_list = []
        O_list = []
        mole_bond_index=[]
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
                                    if mol.degree(index2) <= 3 and ([index2,index1] not in mole_bond_index):
                                        mole_bond_index.append([index1,index2])
                                        mole_bond_index.append([index2,index1])

                                if element[index2]=="O":
                                    if mol.degree(index2) <= 1:   
                                        mole_bond_index.append([index1,index2]) 
                                        mole_bond_index.append([index2,index1]) 
        for index1 in O_list:                                
            if mol.degree(index1) == 1:
                index1_neighbor = mol[index1]
                for key, _ in index1_neighbor.items():
                    if mol.degree(key) > 1:
                        index1_neighbor_2 = mol[key]
                        for index2, _ in index1_neighbor_2.items():
                            if index2 != index1:
                                if element[index2]=="O":
                                    if mol.degree(index2) <= 1 and ([index2,index1] not in mole_bond_index):   
                                        mole_bond_index.append([index1,index2])
                                        mole_bond_index.append([index2,index1])  
        #for pair in mole_bond_atom:
        #    print(pair,":",(element[pair[0]],element[pair[1]]))                                 
        return mole_bond_index


class geometry_rules:  
    def __init__(self):
        pass
    def judge_bond_order(mol,mole_bond_atom=[],bond_to_basis=1):   #判断键级，flag=0说明没有多余的电子跟基底成键
        element=nx.get_node_attributes(mol,"symbol")
        flag = 1
        C_list = []
        O_list = []
        N_list = []
        new_mol = mol.copy()
        for key,value in element.items():
            if value == "C":
                C_list.append(key)
            if value == "O":
                O_list.append(key)  
            if value == "N": 
                N_list.append(key)      
        if bond_to_basis==2 and O_list:   
            for i in O_list:
                if i in mole_bond_atom: 
                    flag = 0
        if O_list:
            for i in O_list:
                tmp_neigh = new_mol[i]
                tmp = list(tmp_neigh.keys())
                for k in tmp:            
                    if "H" in element[k]:
                        new_mol.add_node(i, htag="yes")
                    elif i in mole_bond_atom:
                        new_mol.add_node(i, htag="yes")
        o_deal_dict = nx.get_node_attributes(new_mol, "htag")
        #print(o_deal_dict)

        N_elec_list = []
        C_elec_list = [] 
        if N_list:
            for i in N_list:
                N_elec = 5
                tmp_neigh = new_mol[i]
                tmp = list(tmp_neigh.keys())
                for k in tmp:
                    if k in [key for key in o_deal_dict.keys()]:
                        N_elec = N_elec - 1
                    elif "O" in element[k]:
                        N_elec = N_elec - 2
                    elif "H" in element[k]:
                        N_elec -= 1
                if i in  mole_bond_atom:
                    N_elec -= bond_to_basis
                if N_elec < 0:
                    flag = 0    
                N_elec_list.append(N_elec)
        if C_list:
            for i in C_list:
                C_elec = 4
                tmp_neigh = new_mol[i]
                tmp = list(tmp_neigh.keys())
                for k in tmp:
                    if k in [key for key in o_deal_dict.keys()]:
                        C_elec = C_elec - 1
                    elif "O" in element[k]:
                        C_elec = C_elec - 2  # 好像还没考虑C和界面成键的情况啊
                    elif "H" in element[k]:
                        C_elec = C_elec - 1
                if i in  mole_bond_atom:
                    C_elec -= bond_to_basis   
                C_elec_list.append(C_elec)     
        if C_elec_list:
            if C_elec_list[0] == 1:
                if len(N_elec_list) > 1:
                    print("the bond order of  C > 4,elimination")
                    flag=0
                if len(N_elec_list) == 1:
                    #c_elec = 4
                    n_elec = N_elec_list[0] - 1
                    if n_elec<0:
                        print("C and N can not satisfy the bond order together,elimination")
                        flag=0
            if C_elec_list[0] == 2:
                if len(N_elec_list) == 2:
                    n1_elec = N_elec_list[0] - 1
                    n2_elec = N_elec_list[1] - 1
                    if n1_elec<0 or n2_elec<0:
                        print("the bond order of one N would exceed 5,elimination")
                        flag=0
                if len(N_elec_list) == 1:
                    tmp_ele = []
                    for i in [1, 2]:
                        n_elec = N_elec_list[0] - i
                        tmp_ele.append(n_elec)
                    # if 2 in tmp_ele or 0 in tmp_ele or 4 in tmp_ele:
                    if tmp_ele[0]<0 and tmp_ele[1]<0:
                        print("the bond order of N would exceed 5,elimination")
                        flag=0
            if C_elec_list[0] == 3:
                if len(N_elec_list) == 2:
                    tmp_ele = []
                    for i in [[1, 1], [1, 2], [2, 1]]:  # 成键数目，其实判断条件只需要满足两个单键即可
                        n1_elec = N_elec_list[0] - i[0]
                        n2_elec = N_elec_list[1] - i[1]
                        tmp_ele.append([n1_elec, n2_elec])
                    tmp_flag=0
                    for i in tmp_ele:
                        if i[0]>0 and i[1]>0:
                            tmp_flag=1
                    if tmp_flag==0:
                        flag=0
                        print("the two N can not satisfy the bond order together,elimination")
                if len(N_elec_list) == 1:
                    tmp_ele = []
                    for i in [1, 2]:
                        n_elec = N_elec_list[0] - i
                        tmp_ele.append(n_elec)
                    # if 2 in tmp_ele or 0 in tmp_ele or 4 in tmp_ele:
                    if tmp_ele[0]<0 and tmp_ele[1]<0:
                        print("the bond order of N would exceed 5,elimination")
                        flag=0
            if C_elec_list[0] < 0:              
                flag=0
                print("the bond order of C would denifinetly exceed 4,elimanation")
            if C_elec_list[0] == 0:
                if len(N_elec_list)!=0:
                    flag=0
                    print("the bond order of C would exceed 4,elimanation")                    
        return flag

    def single_ad_structure(slab,mol,ad_site,mole_bond_index=[],type="auto"):  #如果和吸附原子一个原子相连，该原子就要竖直朝上，如果两个，就是要z方向分布
        element=nx.get_node_attributes(mol,"symbol")
        slab=slab.copy()  #必须copy
        new_mol = mol.copy()
        positions = {}
        for i in list(new_mol.nodes):   #初始化坐标信息
            positions[i]=(0,0,0)
        
        #print(element)
        if type=="auto":
            if len(ad_site)==1:
                R=get_radii(slab[ad_site[0]].symbol)+get_radii(element[mole_bond_index[0]])
                #positions[mole_bond_atom[0]]=tuple(slab[ad_site[0]].position+np.array([0,0,R]))
                pos_ad_site_neigh=[slab[i].position for i in ad_site ]
        
                positions[mole_bond_index[0]]=tuple(get_base_position(pos_ad_site_neigh,[R*0.9]))
                #positions[mole_bond_atom[0]]=tuple(get_base_position(pos_ad_site_neigh,[R]))
            elif len(ad_site)==2:
                #site1_coor=slab[ad_site[0]].position
                #print(site1)
                pos_ad_site_neigh=[slab[i].position for i in ad_site ]
                a=np.linalg.norm(slab[ad_site[0]].position-slab[ad_site[1]].position)
                R1=get_radii(slab[ad_site[0]].symbol)
                R2=get_radii(slab[ad_site[1]].symbol)
                R3=get_radii(element[mole_bond_index[0]])
                R=[(R2+R3)*0.9,(R1+R3)*0.9]

                positions[mole_bond_index[0]]= tuple(get_base_position(pos_ad_site_neigh,R))   #计算inter atom的位置,001应该确实不行？？？？？
        else:
            R=-1  # 人为指定
            positions[mole_bond_index[0]]=tuple(slab[ad_site[0]].position+np.array([0,0,R]))

        total_branches = list(nx.bfs_successors(new_mol, mole_bond_index[0]))

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


    def double_ad_structure(slab,mol,ad_site,mole_bond_index=[],type="auto"):  #双位点搭吸附构型的几何规则
        element=nx.get_node_attributes(mol,"symbol")
        new_mol = mol.copy()  #对图进行复制，否则报错，直接用XX.copy()也行
        slab = slab.copy()   
        positions = {}    
        for i in list(new_mol.nodes):   #初始化坐标信息
            positions[i]=(0,0,0)
        r_bond=[]
        r_site=[]
        
        for i in mole_bond_index:
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
            
        positions[mole_bond_index[0]] = tuple(base_position0)
        positions[mole_bond_index[1]] = tuple(base_position1)

        #print("base,position",base_position0,base_position1)   #这个地方对两个基本位置进行了一种处理

        if tuple(mole_bond_index) in new_mol.edges:
                new_mol.remove_edge(*mole_bond_index)

        uvec1  = np.array([[0,0,1]])  #这个不应该指定的
        uvec2 = np.cross(uvec1, uvec0)         #叉乘
        uvec2 = uvec2/np.linalg.norm(uvec2)
        uvec1 = -np.cross(uvec2, uvec0)
        #print("uc1", uvec1)
        

        branches0 = list(nx.bfs_successors(new_mol, mole_bond_index[0]))
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


        branches1 = list(nx.bfs_successors(new_mol, mole_bond_index[1]))
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

    def double_cross_ad_structure(slab,mol,ad_site,mole_bond_index=[],type="normal"):  #arc位点搭吸附构型的几何规则
        element=nx.get_node_attributes(mol,"symbol")
        slab = slab.copy()
        new_mol = mol.copy()
        positions = {}
        for i in list(new_mol.nodes):   #初始化坐标信息
            positions[i]=(0,0,0)    
        vec_site_0=slab[ad_site[0]].position-slab[ad_site[1]].position

        r_bond=[]
        r_site=[]
        
        for i in mole_bond_index:
            r_bond.append(get_radii(element[i]))

        for i in ad_site:
            r_site.append(get_radii(slab[i].symbol))

        inter_atom_index = 0
        connec_list_of_bond = []
        for i in mole_bond_index:
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

        positions[mole_bond_index[0]] = tuple(base_position0)
        positions[mole_bond_index[1]] = tuple(base_position1)
            #write("tmp.vasp",atoms)
        
        uvec1 = np.array([[0, 0, 1]])  
        uvec2 = np.cross(uvec1, uvec0)  
        uvec2 = uvec2/np.linalg.norm(uvec2)  #新增
        b=r_bond[0]+get_radii(element[inter_atom_index])
        c=r_bond[1]+get_radii(element[inter_atom_index])
        positions[inter_atom_index]= get_base_position([base_position1,base_position0],[b,c])  #计算inter atom的位置
        new_mol.remove_edge(*[mole_bond_index[0],inter_atom_index])
        new_mol.remove_edge(*[mole_bond_index[1], inter_atom_index])

            #write("tmp0.vasp", atoms)
        uvec1=np.cross(uvec0,uvec2)
        uvec1 = uvec1/np.linalg.norm(uvec1)  #新增

        branches0 = list(nx.bfs_successors(new_mol, mole_bond_index[0]))
            #print("brance0", branches0, branches0[1:])
        if len(branches0[0][1]) != 0:
                uvec = [-uvec0, uvec1[0], uvec2[0]]  # 可能是这个地方出了问题。
                dou_ad_entension(positions,element, uvec, branches0[0])
                root = branches0[0][0]
                #print("before:",positions)
                for branch in branches0[1:]:
                    mole_extension(positions,element, branch, root)
                #print("after:",positions)        

        branches1 = list(nx.bfs_successors(new_mol, mole_bond_index[1]))
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

