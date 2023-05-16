#import math
import os
from ase.io import read
from basic_func import get_neighbor_list,atom_to_graph

def judge_smx_connect(path): #判断界面上的分子是否分解，也可通过图连通性判断分解
    count=0
    for file_1 in os.listdir(path):  
        try:      
            atoms=read(path+"/"+file_1+"/CONTCAR")
        except IndexError:
            os.system("cp -r {}".format(file_1+" "+index_error))
            continue
        ad_list=[]
        flag=0
        for index,atom in enumerate(atoms):
            if index>=49:
                if atom.symbol =="C" :
                    flag=1
        if flag==1:
            for index,atom in enumerate(atoms):              
                if index>=49:
                    if atom.symbol =="C" :
                        ad_list.append(index)
                if index>=49:
                    if atom.symbol =="O" or atom.symbol =="N" :
                        ad_list.append(index)
            #print(file_1,ad_list)
            neighbor_list = get_neighbor_list(atoms,0.001,1.25)
            #full,cpmpare_graph,reduce_full,ad_chem1=atom_to_graph(atoms,neighbor_list,grid=(0,0,0),adsorbate_atoms=ad_list)  
            full,ad_chem1=atom_to_graph(atoms,neighbor_list,adsorbate_atoms=ad_list)
            if len(ad_chem1)!=1:
                print("{} is not connected ".format(file_1))
                os.chdir(path)
                os.system("cp -r {}".format(file_1+" "+decomposition))
        
        if flag==0:
            for index,atom in enumerate(atoms):
                if index>=49:
                    if atom.symbol =="O" or atom.symbol =="N" :
                        ad_list.append(index)
            neighbor_list = get_neighbor_list(atoms,0.001,1.25)
            full,ad_chem1=atom_to_graph(atoms,neighbor_list,adsorbate_atoms=ad_list)  
            #full,cpmpare_graph,reduce_full,ad_chem1=atom_to_graph(atoms,neighbor_list,grid=(0,0,0),adsorbate_atoms=ad_list)
            if len(ad_chem1)!=1:
                print("{} is not connected ".format(file_1))
                os.chdir(path)
                os.system("cp -r {}".format(file_1+" "+decomposition))



path = "/home/all"
index_error = "/home/index_error"
decomposition = "/home/decomposition"

judge_smx_connect(path)

