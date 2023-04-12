import os
import ase
import numpy as np
import networkx as nx
from ase.io import read
from ase.io import write
from ase import neighborlist
import networkx.algorithms.isomorphism as iso
from ase.data import covalent_radii
from ase import Atoms
#from numba import jit
import shutil
import scipy
from ase.neighborlist import NeighborList, natural_cutoffs
from numpy.linalg import norm
from itertools import combinations
#from default_func import graphstruc
import itertools
from default_func import get_radicial,get_radii,get_atomic_number
import math




def _get_basis_vectors(coordinates):  #获得一个正交的坐标系，basis1是键的方向，
    if len(coordinates) == 3:  #会有等于3的时候？
        c0, c1, c2 = coordinates
    else:
        c0, c1 = coordinates
        c2 = np.array([0, 1, 0])

    basis1 = c0 - c1   
    basis2 = np.cross(basis1, c0 - c2)
    basis3 = np.cross(basis1, basis2)
    #print("3 basis",basis1,basis2,basis3)  这个矢量构建的根据

    basis1 /= np.linalg.norm(basis1)
    basis2 /= np.linalg.norm(basis2)
    basis3 /= np.linalg.norm(basis3)

    basis_vectors = np.vstack([basis1, basis2, basis3])
    #print("basis_vec", basis_vectors)

    return basis_vectors


def get_exten_position(coordinates, distance, angle=109.47, dihedral=0):  
    #print("coordinate",coordinates)
    #print(dihedral)
    v0, v1 = coordinates
    v2 = np.array([0, 1, 0])
    basis1 = v0 - v1   
    basis2 = np.cross(basis1, v0 - v2)
    basis3 = np.cross(basis1, basis2)
    #print("3 basis",basis1,basis2,basis3)  这个矢量构建的根据，构建的是一个正交的矩阵
    basis1 /= np.linalg.norm(basis1)
    basis2 /= np.linalg.norm(basis2)
    basis3 /= np.linalg.norm(basis3)
    basis = np.vstack([basis1, basis2, basis3])
    inter_ang = np.deg2rad(angle)  
    rot_ang = np.deg2rad(dihedral)
    #print("diii",tor)

    basis[1] *= -np.sin(rot_ang)   #这边交换cos和sin可能可以解决延展的问题
    basis[2] *= np.cos(rot_ang)
    #basis[1] *=  np.cos(rot_ang)  
    #basis[2] *= -np.sin(rot_ang)

    vector = basis[1] + basis[2]
    vector /= np.linalg.norm(vector)  #单位向量

    vector *= distance * np.sin(inter_ang)
    basis[0] *= distance * np.cos(inter_ang)

    exten_atom_posi = coordinates[0] + vector - basis[0]
    
    return exten_atom_posi

def mole_extension(positions,element,branch,center_atom=None):  #目前最多只能考虑4个键的？  
    center, ex_atom = branch
    tmp_num = len(ex_atom)
    #print(tmp_num)
    if tmp_num == 0:
        return
    tmp_branch = [center] + ex_atom
    rlist = np.array([get_radii(element[i]) for i in tmp_branch])
    d=rlist[0]+rlist[1:]  #  中心原子和其他原子相加
    #print("d:", d,rlist)
    center_pos = np.array(positions[center])   
    root_center_pos = np.array(positions[center_atom])
    inter_angle = 120
    dihedral = np.array([])
    if tmp_num == 1:   #一个延伸，和原来那个键成180度
            inter_angle = 180
            dihedral = np.array([0, 0, 0])  #实际上只用的到第一个
    elif tmp_num == 2:
        inter_angle=120
        #dihedral[1] = 180  #两个延伸，延伸的两个原子180，但是和中间那个是120
        dihedral = np.array([0, 180, 0])
    elif tmp_num==3:  #三个，四面体
        inter_angle = 109.47 
        dihedral = np.array([0, 120, -120])

    for k in range(tmp_num):
        #print(k)
        c = get_exten_position([center_pos, root_center_pos],distance=d[k],angle=inter_angle,dihedral=dihedral[k])  #获得的是延伸分子的位置，中心原子是已经知道了的
        #print("before",nodes[k],"position",atoms[nodes[k]].position)
        positions[ex_atom[k]]=tuple(c)  
        #print("after",nodes[k],"position",atoms[nodes[k]].position)
    # print("root",root)
    return center





def get_base_position(ad_site_neigh_coor, r, zvector=None):#  理解起来应该是返回桥位或者3个位点的中心原子的坐标，以此作为吸附的位点。
    
    if zvector is None:
        zvector = np.array([0,0,1])

    if len(ad_site_neigh_coor) == 1:

        return ad_site_neigh_coor[0] + r[0] * zvector

    
    vec_site = ad_site_neigh_coor[1]-ad_site_neigh_coor[0]  # 一个位点指向另一个位点的矢量
    len_vec = np.linalg.norm(vec_site)  # 范数
    uvec1 = vec_site / len_vec  # 单位方向矢量

    if len(ad_site_neigh_coor)==2:
        
        a=np.linalg.norm(ad_site_neigh_coor[0]-ad_site_neigh_coor[1])        
        c=r[0]
        b=r[1]
        #print(a,b,c)
        #print("aaaaa", (b*b-a*a-c*c)/(-2*a*c))
        angle_1=math.degrees(math.acos((b*b-a*a-c*c)/(-2*a*c)))   #余弦公式求角度
        #angle_2=math.degrees(math.atan((H2-H1)/))     #2.007是Cu-N键的长度
        delt_h=abs(ad_site_neigh_coor[1][2]-ad_site_neigh_coor[0][2])
        a0=np.array([ad_site_neigh_coor[0][0],ad_site_neigh_coor[0][1],0])
        b0=np.array([ad_site_neigh_coor[1][0],ad_site_neigh_coor[1][1],0])
        vec_a0b0=(a0-b0)/np.linalg.norm(a0-b0)
        
        delt_ad_site=np.linalg.norm(np.array(a0-b0))
        
        angle_2=math.degrees(math.atan(delt_h/delt_ad_site))        
        #print(delt_ad_site,delt_h,angle_2,angle_1)
        dx=c*np.cos((angle_1+angle_2)/180*np.pi)
        dz=c*np.sin((angle_1+angle_2)/180*np.pi)

        base_position= tuple(ad_site_neigh_coor[1]+vec_a0b0*dx+np.array([0,0,1])*dz)  

        return base_position
    
    if len(ad_site_neigh_coor)==3:

        vec2 = ad_site_neigh_coor[2] - ad_site_neigh_coor[0]
        i = np.dot(uvec1, vec2)
        vec2 = vec2 - i * uvec1

        uvec2 = vec2 / np.linalg.norm(vec2)
        uvec3 = np.cross(uvec1, uvec2)
        j = np.dot(uvec2, vec2)

        x = (r[0]**2 - r[1]**2 + d**2) / (2 * d)
        y = (r[0]**2 - r[2]**2 - 2 * i * x + i**2 + j**2) / (2 * j)
        z = np.sqrt(r[0]**2 - x**2 - y**2)
        if np.isnan(z):
            z = 0.01
        base_position = ad_site_neigh_coor[0] + x * uvec1 + y * uvec2 + z * uvec3
        return base_position

    

def sin_ad_extension(positions,element,branch,cross=None,vec=None):   #处理和界面接触的那个原子或者cross吸附中间的那个原子


    center, ex_atom = branch
    tmp_num = len(ex_atom)
    #print(tmp_num)
    if tmp_num == 0:
        return
    tmp_branch = [center] + ex_atom
    rlist = np.array([get_radii(element[i]) for i in tmp_branch])
    d=rlist[0]+rlist[1:]  #  中心原子和其他原子相加
    center_pos = np.array(positions[center])   
    

    
    if cross=="y":  #处理cross吸附，中间那个原子的延展方式。
        if len(ex_atom)==1:
            coord0 = center_pos + np.array([0, 0, d[0]])
            positions[ex_atom[0]] = tuple(coord0)
        if len(ex_atom)==2:
            v1=vec   #这个vec是和底下两个双原子连键方式垂直的地方
            v2=-vec
            angle=(60/180)* np.pi
            coord0 = center_pos + np.array([d[0] * v1[0] * np.sin(angle), d[0] * v1[1] * np.sin(angle), d[0] * np.cos(angle)])
            positions[ex_atom[0]] = tuple(coord0)    

            coord1 = center_pos + np.array(
                [d[1] * v2[0] * np.sin(angle), d[1] * v2[1] * np.sin(angle), d[1] * np.cos(angle)])
            positions[ex_atom[1]] = tuple(coord1)    

    else:  #如果不是single时的处理方式，之后要把每个代码的思想都记下来
        if len(ex_atom) == 1:
            coord0 = center_pos + np.array([0,0,d[0]])
            positions[ex_atom[0]] = tuple(coord0) 
        elif len(ex_atom) == 2:
            coord0 = center_pos + np.array([d[0]*np.cos(1/6.*np.pi),0,d[0]*np.cos(1/3.*np.pi)])
            positions[ex_atom[0]] = tuple(coord0) 
            coord1 = center_pos + np.array([-d[1] * np.cos(1/6. * np.pi), 0, d[1] * np.cos(1 / 3. * np.pi)])
            positions[ex_atom[1]] = tuple(coord1)
        elif len(ex_atom) == 3:
            v_y_0=np.array([0,1,0])
            v_y_1=np.array([np.cos(1/6*np.pi),-np.sin(1/6*np.pi),0])
            v_y_2=np.array([-np.cos(1 / 6 * np.pi), -np.sin(1 / 6 * np.pi), 0])
            angle=(70.53/180)* np.pi
            #print(np.cos(angle),np.sin(angle))

            coord0 = center_pos + np.array([d[0]*v_y_0[0]*np.sin(angle),d[0]*v_y_0[1]*np.sin(angle) ,d[0]*np.cos(angle)])
            positions[ex_atom[0]] = tuple(coord0)

            coord1 = center_pos + np.array([d[1] * v_y_1[0] * np.sin(angle), d[1] * v_y_1[1] * np.sin(angle ),d[1] * np.cos(angle )])
            positions[ex_atom[1]] = tuple(coord1)

            coord2 = center_pos + np.array([d[2] * v_y_2[0] * np.sin(angle), d[2] * v_y_2[1] * np.sin(angle ),d[2] * np.cos(angle )])
            positions[ex_atom[2]] = tuple(coord2)

        else:
            raise ValueError('exceed the limitation of bond')

def dou_ad_entension(positions,element, uvec, branch):  #处理和界面吸附的双原子的延展方式
   
    center, ex_atom = branch
    tmp_num = len(ex_atom)
    #print(uvec)
    if tmp_num == 0:
        return
    tmp_branch = [center] + ex_atom
    rlist = np.array([get_radii(element[i]) for i in tmp_branch])
    d=rlist[0]+rlist[1:]  #  中心原子和其他原子相加
    center_pos = np.array(positions[center])   
    #print(atoms[nodes[0]].position)
    # Single additional atom
    if len(ex_atom) == 1:   #和键的方向延伸的方向成60角，然后在键方向和发下方向的平面上。
        #print("nodes==1", uvec[0],uvec[1] )

        coord0 = center_pos + d[0] * uvec[0] * np.cos(1 / 6. * np.pi) + d[0] * uvec[1] * np.sin(1 / 6. * np.pi)
        #print(center_pos)
        positions[ex_atom[0]] = coord0

    # Two branch system
    elif len(ex_atom) == 2:
        #print("ttotototototottoto")
        coord0 = center_pos + d[0] * uvec[1] * np.cos(1 / 3. * np.pi) + np.sin(1 / 3. * np.pi) * d[0] * uvec[0] * np.cos(1 / 3. * np.pi) + np.sin(1 / 3. * np.pi) * d[0] * uvec[2] * np.sin(1 / 3. * np.pi)
        positions[ex_atom[0]] = coord0

        coord1 = center_pos + d[1] * uvec[1] * np.cos(1 / 3. * np.pi) + np.sin(1 / 3. * np.pi)* d[1] * uvec[0] * np.cos(1 / 3. * np.pi) + np.sin(1 / 3. * np.pi) * d[1] * -uvec[2] * np.sin(1 / 3. * np.pi)
        positions[ex_atom[1]] = coord1

    else:
        #print("element:",element,"branch:",branch)
        raise ValueError('exceed the limitation of bond.')


def sin_neg_extension(positions,element,  branch ,vec=None):  #处理界面负向延展时的方式
    center, ex_atom = branch
    tmp_num = len(ex_atom)
    #print(tmp_num)
    if tmp_num == 0:
        return
    tmp_branch = [center] + ex_atom
    rlist = np.array([get_radii(element[i]) for i in tmp_branch])
    d=rlist[0]+rlist[1:]  #  中心原子和其他原子相加
    center_pos = np.array(positions[center])     
    if len(ex_atom) == 1:
            coord0 = center_pos - np.array([0, 0, d[0]])
            positions[ex_atom[0]] = coord0
    if len(ex_atom) == 2:
            v1 = vec
            v2 = -vec
            angle = (50 / 180) * np.pi
            coord0 = center_pos - np.array([d[0] * v1[0] * np.sin(angle), d[0] * v1[1] * np.sin(angle), d[0] * np.cos(angle)])
            positions[ex_atom[0]] = coord0
            coord1 = center_pos - np.array([d[1] * v2[0] * np.sin(angle), d[1] * v2[1] * np.sin(angle), d[1] * np.cos(angle)])
            positions[ex_atom[1]] = coord1



#print(molecule("CO2"))