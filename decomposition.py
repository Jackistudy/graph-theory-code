#import math
import os
from ase.io import read
from struc_to_graph import atom_to_graph
from default_func import get_neighbor_list

def judge_smx_connect(path):
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
            full,cpmpare_graph,reduce_full,ad_chem1=atom_to_graph(atoms,neighbor_list,grid=(0,0,0),adsorbate_atoms=ad_list)  #全图的生成在边界上会出现重复的问题，这种重复具有不同oxoyoz
            #full,ad_chem1=atom_to_graph(atoms,neighbor_list,adsorbate_atoms=ad_list)
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
            #full,ad_chem1=atom_to_graph(atoms,neighbor_list,adsorbate_atoms=ad_list)  
            full,cpmpare_graph,reduce_full,ad_chem1=atom_to_graph(atoms,neighbor_list,grid=(0,0,0),adsorbate_atoms=ad_list)
            if len(ad_chem1)!=1:
                print("{} is not connected ".format(file_1))
                os.chdir(path)
                os.system("cp -r {}".format(file_1+" "+decomposition))



path = "/home/all"
index_error = "/home/index_error"
decomposition = "/home/decomposition"

judge_smx_connect(path)

