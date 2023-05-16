import networkx as nx
from create_reaction_network import add_h,drop_oh,judge_coupling_3,get_reaction_network2
from ase.io import read
from basic_func import get_neighbor_list,get_thermal_correction_from_vaspkit_output,atom_to_graph
import os
import numpy as np


least_energy_file = '/homw/least_energy'

reactnet_0210='/home/reactnet.txt'

def generate_ad_graph(file_path):  #返回子图，即吸附界面上的小分子的图
    atoms = read(file_path)
    ad_list = []
    for index,atom in enumerate(atoms):
        if index > 48:   #序号大于48的原子为界面上的原子
            ad_list.append(index)    
    neighbor_list = get_neighbor_list(atoms,0.001,1.25)  #成键阈值为1.25
    full,ad_graph=atom_to_graph(atoms,neighbor_list,ad_atoms=ad_list)
    return ad_graph


def create_total_graph(path):  #返回文件名、分子图和自由能
    total_graph=[]
    for file in os.listdir(path):
        file_path = path+'/'+file+'/'+'CONTCAR'
        ad_graph = generate_ad_graph(file_path)
        energy=0
        if os.path.exists(path+"/"+file+"/scf/output.relax"):    #文件名后缀写错了，正常应该是.scf
            for line in open(path+"/"+file+"/scf/output.relax"):
                    if "1 F" in line:
                        energy=float(line.split()[4])
        else:  
            for line in open(path+"/"+file+"/scf/output.scf"):   
                    if "1 F" in line:
                        energy=float(line.split()[4])                          
        thermal_correction=get_thermal_correction_from_vaspkit_output(path+"/"+file) #零点能和熵校正
        free_energy=energy+thermal_correction        
        total_graph.append([file.split('-')[1],ad_graph,free_energy])
    return total_graph    

def molecular_formula(before_formula):  #把smiles变成最简分子式
    C_num=0
    H_num=0
    N_num=0
    O_num=0
    length=len(before_formula)
    count=1
    for for_str in before_formula:
        if for_str.isdigit():
            count+=1
            continue
        if count<length:
            if before_formula[count].isdigit():
                if for_str=='C':
                    C_num+=int(before_formula[count])
                elif for_str=='H':
                    H_num+=int(before_formula[count])
                elif for_str=='N':
                    N_num+=int(before_formula[count])
                elif for_str=='O':
                    O_num+=int(before_formula[count])
            else:
                if for_str=='C':
                    C_num+=1
                elif for_str=='H':
                    H_num+=1
                elif for_str=='N':
                    N_num+=1
                elif for_str=='O':
                    O_num+=1  
        else:
            if for_str.isdigit():
                continue
            else:                     
                if for_str=='C':
                    C_num+=1
                elif for_str=='H':
                    H_num+=1
                elif for_str=='N':
                    N_num+=1
                elif for_str=='O':
                    O_num+=1  
        count+=1

    if C_num>1:
        C_num_str="C"+str(C_num)
    elif C_num==1:
        C_num_str="C"
    elif C_num==0:
        C_num_str=""

    if H_num>1:
        H_num_str="H"+str(H_num)
    elif H_num==1:
        H_num_str="H"
    elif H_num==0:
        H_num_str=""

    if N_num>1:
        N_num_str="N"+str(N_num)
    elif N_num==1:
        N_num_str="N"
    elif N_num==0:
        N_num_str=""

    if O_num>1:
        O_num_str="O"+str(O_num)
    elif O_num==1:
        O_num_str="O"
    elif O_num==0:
        O_num_str=""

    after_formula=C_num_str+H_num_str+N_num_str+O_num_str
    return after_formula

# create reaction path

total_graph=create_total_graph(least_energy_file)
reaction_network,file_list=get_reaction_network2(total_graph)
with open(reactnet_0210,'w') as f:
    for each in file_list:
        f.write(each)
        f.write('\n')


# use reaction path
cou_mod_path='/home/reactnet.txt' #该文件存放反应网络包含的所有化学反应的反应物、产物、能量变化、反应类型信息
path_energy_dict={}
reaction_network=nx.DiGraph()
pH_mod=0.4012  #ph=6.8对应的自由能校正
with open(cou_mod_path,'r') as f:
    lines = f.readlines()
for line in lines:
    line=line.strip('\n')
    line=line.split(',')

    if line[3][:3]=='add':
        tmpedge=[(line[0],line[1],{"deltE":float(line[2]),"type":"add_H"})]
        pair=line[0]+'_to_'+line[1]
        path_energy_dict[pair]=float(line[2])+pH_mod
        reaction_network.add_edges_from(tmpedge)
    elif line[3][:4]=='drop':
        tmpedge=[(line[0],line[1],{"deltE":float(line[2]),"type":"drop_H2O"})]
        pair=line[0]+'_to_'+line[1]
        path_energy_dict[pair]=float(line[2])+pH_mod
        reaction_network.add_edges_from(tmpedge)
    else:
        tmpedge=[(line[0],line[1],{"deltE":float(line[2]),"type":"couple"})]
        pair=line[0]+'_to_'+line[1]
        path_energy_dict[pair]=float(line[2])
        reaction_network.add_edges_from(tmpedge)        

molecule_path='/home/molecule.txt'
#for file in os.listdir(molecule_path):
with open(molecule_path,'r') as f:
    lines=f.readlines()
for line in lines:
    file=line.strip('\n')
    CO2_X_path=nx.all_simple_edge_paths(reaction_network,'CO2',file)
    max_energy=-100000
    abs_energy_sum=0    
    total=0
    for i in CO2_X_path: #通过两次列表排序找出最优路径
        count=1
        length=len(i)
        path_str=''
        max_energy_tmp=-100000
        abs_energy_sum_tmp=0   
        eng_list=[]
        total+=1
        for j in i:
            pair=j[0]+'_to_'+j[1]
            energy=round(path_energy_dict[pair],4)
            
            if energy>max_energy_tmp:
                max_energy_tmp=energy  #决速步能量
            #abs_energy=abs(energy)
            if energy>0:
                abs_energy=energy
            abs_energy_sum_tmp+=abs_energy
            #打印反应路径        
            if count==1:
                path_str+=j[0]+' -> '+j[1]+' -> '
                count+=1
                eng_list.append(energy)
            elif count == length:
                path_str+=j[1]
                eng_list.append(energy)
                #print(path_str)
                if total==1:
                    reaction_path_list=[path_str,max_energy_tmp,abs_energy_sum_tmp,eng_list]
                else:
                    if max_energy_tmp<reaction_path_list[1]:  #通过不断排序找出决速步能量最低的反应路径
                        reaction_path_list=[path_str,max_energy_tmp,abs_energy_sum_tmp,eng_list]
                    elif max_energy_tmp==reaction_path_list[1]:
                        if abs_energy_sum_tmp<reaction_path_list[2]:  #若决速步能量相同，找出绝对值为正的反应步骤，其累加的自由能最小
                            reaction_path_list=[path_str,max_energy_tmp,abs_energy_sum_tmp,eng_list]
                #reaction_path_list.append(tmp)
            else:
                path_str+=j[1]+' -> '
                count+=1  
                eng_list.append(energy)
    print(file) #产物名称
    print(str(total))  #某产物的反应路径数目
    print(reaction_path_list[0])  #最优路径
    print(reaction_path_list[1])  #决速步
    print(reaction_path_list[3])  #最优路径每一步的能量

