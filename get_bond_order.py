
import numpy as np
import os
import math
import networkx as nx
from ase.io import read
from ase.io import write
from ase import neighborlist
import networkx.algorithms.isomorphism as iso
from ase.data import covalent_radii
#from numba import jit
import shutil
#import scipy
from numpy.linalg import norm
from itertools import combinations
from default_func import get_neighbor_list,periodic_deal
import copy



def get_smx_bond_order(mol,mole_bond_atom=[],bond_to_basis=1):   #判断键级，跟之前没啥区别
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


