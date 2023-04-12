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


em = iso.numerical_edge_match('bonds', 1)
nm = iso.numerical_node_match('number', 1)


radicial={1:1,6:4,7:4,8:2}  # the number of bonds that can be formed, this can be change in specific reaction
#print(radicial)

class graphstruc(Atoms):   #  还需不需要这个类，
 

    def __init__(self,
                 symbols=None,
                 positions=None,
                 numbers=None,
                 tags=None,
                 momenta=None,
                 masses=None,
                 magmoms=None,
                 charges=None,
                 scaled_positions=None,
                 cell=None,
                 pbc=None,
                 celldisp=None,
                 constraint=None,
                 calculator=None,
                 info=None,
                 edges=None):
        super().__init__(
            symbols, positions, numbers, tags, momenta,
            masses, magmoms, charges, scaled_positions, cell,
            pbc, celldisp, constraint, calculator, info)

        if isinstance(edges, np.ndarray):
            if self.pbc.any():
                self._graph = MultiGraph(edges)
            else:
                self._graph = Graph(edges)
        else:
            if self.pbc.any():
                self._graph = MultiGraph()
            else:
                self._graph = Graph()

        nodes = [[i, {'number': n}]
                 for i, n in enumerate(self.arrays['numbers'])]
        self._graph.add_nodes_from(nodes)

        if isinstance(edges, list):
            self._graph.add_edges_from(edges)

    @property
    def graph(self):
        return self._graph

    @property
    def nodes(self):
        return self._graph.nodes

    @property
    def edges(self):
        return self._graph.edges

    @property
    def adj(self):
        return self._graph.adj

    @property
    def degree(self):
        degree = self._graph.degree
        return np.array([_[1] for _ in degree])

    @property
    def connectivity(self):
        connectivity = nx.to_numpy_matrix(self._graph).astype(int)
        return connectivity
    def is_isomorph(self, other):
        """Check if isomorphic by bond count and atomic number."""
        isomorphic = nx.is_isomorphic(
            self._graph, other._graph, edge_match=em, node_match=nm)

        return isomorphic




    def copy(self):  #没有这个copy函数就会报错！！！！！！！！
        """Return a copy."""
        atoms = self.__class__(cell=self.cell, pbc=self.pbc, info=self.info)

        atoms.arrays = {}
        for name, a in self.arrays.items():
            atoms.arrays[name] = a.copy()
        atoms.constraints = copy.deepcopy(self.constraints)
        if hasattr(self, '_graph'):
            atoms._graph = self._graph.copy()

        return atoms


def get_neighbor_list(atoms,skin,cutoff=1):
    radii = neighborlist.natural_cutoffs(atoms,mult=cutoff)  #这个radii可以认为设置把，根据需要，未必采用ase的方法  获得近邻原子
    #print(radii)
    neighbor_list = neighborlist.NeighborList(radii,skin,bothways=True)
    neighbor_list.update(atoms)
    return neighbor_list


def get_atomic_number(formula,return_count=False,return_sym=False):  #获得字符串中的元素符号及对应个数，还需不需要？
    """Return the atomic numbers associated with a chemical formula.

        Parameters
        ----------
        formula : string
            A chemical formula to parse into atomic numbers.
        return_count : bool
            Return the count of each element in the formula.

        Returns
        -------
        numbers : ndarray (n,)
            Element numbers in associated species.
        counts : ndarray (n,)
            Count of each element in a species.
        """
    parse = re.findall('[A-Z][a-z]?|[0-9]+', formula)

    values = {}
    for i, e in enumerate(parse):
        if e.isdigit():
            values[parse[i - 1]] += int(e) - 1
        else:
            if e not in values:
                values[e] = 1
            else:
                values[e] += 1
    #print(type(sym))
    numbers = np.array([sym.index(k) for k in values.keys()])
    symbol=np.array([k for k in values.keys()])
    srt = np.argsort(numbers)
    numbers = numbers[srt]

    if return_count:
        counts = np.array([v for v in values.values()])[srt]
        #print(numbers,counts)
        return numbers, counts
    if return_sym:
        counts = np.array([v for v in values.values()])[srt]
        return symbol,counts
    return numbers

def get_radicial(ele_array):
    radicial_array=np.zeros(len(ele_array))
    for i in range(len(ele_array)):
        bond_number=radicial[ele_array[i]]
        radicial_array[i]=bond_number
    #print(radicial_array)
    return radicial_array

def get_radii(element):
    r_dict = {'C':0.76,'N':0.71,'O':0.66,'H':0.31,"Cu":1.32,"Ni":1.24} 
    return r_dict[element]

def periodic_deal(path):  
    vector=getlatvec(path+"/CONTCAR")
    re_vector=np.linalg.inv(vector)
    #print(vector)
    struc=read(path+"/CONTCAR")
    ad_list=[]
    for index,atom in enumerate(struc):
        if atom.symbol=="N" or atom.symbol=="O" or atom.symbol=="H" or atom.symbol=="C":
            #print("before",atom.position)
            position=np.dot(atom.position,re_vector)
            print(position)        
            if position[0]>0.85:
                position[0]-=1
            if position[1]>0.85:
                position[1]-=1
            re_position=np.dot(position,vector)
            atom.position=re_position
            #print("chuli",atom.position)
    return struc


def get_tri_inter_central(position):   #获得三角形的外接圆的圆心/外心
    a=position
    #hcp=np.array(tri[0][0]+)
    x1 = a[0][0]
    x2 = a[1][0]
    x3 = a[2][0]
    y1 = a[0][1]
    y2 = a[1][1]
    y3 = a[2][1]
    z1 = a[0][2]
    z2 = a[1][2]
    z3 = a[2][2]
    a1 = (y1 * z2 - y2 * z1 - y1 * z3 + y3 * z1 + y2 * z3 - y3 * z2)
    b1 = -(x1 * z2 - x2 * z1 - x1 * z3 + x3 * z1 + x2 * z3 - x3 * z2)
    c1 = (x1 * y2 - x2 * y1 - x1 * y3 + x3 * y1 + x2 * y3 - x3 * y2)
    d1 = -(x1 * y2 * z3 - x1 * y3 * z2 - x2 * y1 * z3 + x2 * y3 * z1 + x3 * y1 * z2 - x3 * y2 * z1)
    a2 = 2 * (x2 - x1)
    b2 = 2 * (y2 - y1)
    c2 = 2 * (z2 - z1)
    d2 = x1 * x1 + y1 * y1 + z1 * z1 - x2 * x2 - y2 * y2 - z2 * z2
    a3 = 2 * (x3 - x1)
    b3 = 2 * (y3 - y1)
    c3 = 2 * (z3 - z1)
    d3 = x1 * x1 + y1 * y1 + z1 * z1 - x3 * x3 - y3 * y3 - z3 * z3
    x = -(b1 * c2 * d3 - b1 * c3 * d2 - b2 * c1 * d3 + b2 * c3 * d1 + b3 * c1 * d2 - b3 * c2 * d1) / (
                a1 * b2 * c3 - a1 * b3 * c2 - a2 * b1 * c3 + a2 * b3 * c1 + a3 * b1 * c2 - a3 * b2 * c1)
    y = (a1 * c2 * d3 - a1 * c3 * d2 - a2 * c1 * d3 + a2 * c3 * d1 + a3 * c1 * d2 - a3 * c2 * d1) / (
                a1 * b2 * c3 - a1 * b3 * c2 - a2 * b1 * c3 + a2 * b3 * c1 + a3 * b1 * c2 - a3 * b2 * c1)
    z = -(a1 * b2 * d3 - a1 * b3 * d2 - a2 * b1 * d3 + a2 * b3 * d1 + a3 * b1 * d2 - a3 * b2 * d1) / (
                a1 * b2 * c3 - a1 * b3 * c2 - a2 * b1 * c3 + a2 * b3 * c1 + a3 * b1 * c2 - a3 * b2 * c1)
    center=np.array([x,y,z])

    return center


def judge_struc_proper(atoms):
    length=[]
    flag=1
    for i in atoms.get_all_distances():
        for j in i:
            if j!=0:
                length.append(j)
    for k in length:
        if k<0.9:
            flag=0
    return flag
    #print(adsorb_index)



def get_thermal_correction_from_vaspkit_output(path):
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


def smile_to_graph(smiles):
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
