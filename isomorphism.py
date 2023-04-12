import os
from default_func import get_neighbor_list
from struc_to_graph import atom_to_graph
import networkx.algorithms.isomorphism as iso
from ase.io import read

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
  