import networkx as nx
import numpy as np
import h5py
from scipy.spatial.distance import cdist
db = h5py.File('/home/jota/DataSets/ani1x-release.h5')

formula = 'C2H5N1O2'
coords = db[formula+'/coordinates'][50:100] 
species = db[formula+'/atomic_numbers']

max_valence = {'1':1,
               '6':4,
               '7':3,
               '8':2}

def check_max_valence(species, n_neighbors): 
    ok = True
    i=0
    while (ok == True) and i <= len(species) - 1:
        if max_valence[str(species[i])] < n_neighbors[i]:
            ok = False
        i+=1
    return ok

def get_neighbors(coord, species, cut_Heavy, cut_H, check_valence):
    n_atoms = len(species)
    bin_dist_matrix = np.zeros((n_atoms, n_atoms))
    dist_matrix = cdist(coord, coord)
    for i in range(n_atoms):
        for j in range(n_atoms):
            cut = cut_Heavy
            if (species[i] == 1) or (species[j]==1):
                cut = cut_H
            if dist_matrix[i][j] <= cut:
                bin_dist_matrix[i,j]=1
    bin_dist_matrix = bin_dist_matrix - np.identity(n_atoms)
    n_neighbors = np.sum(bin_dist_matrix, axis=1)
    if (check_valence):
        shorter_cut_Heavy = cut_Heavy
        shorter_cut_H = cut_H
        while (check_max_valence(species, n_neighbors) == False): #algun elemento excede la valencia
            shorter_cut_Heavy = shorter_cut_Heavy*0.95
            shorter_cut_H = shorter_cut_H*0.95
            bin_dist_matrix, n_neighbors = get_neighbors(coord, species, shorter_cut_Heavy, shorter_cut_H, False)               
    return bin_dist_matrix, n_neighbors

def find_aromatics(bin_dist, n_at):
    G = nx.Graph()
    for atom in range(n_at):
        G.add_node(atom)
    for u in G.nodes:
        for v in G.nodes:
            if u < v:
                if bin_dist[u][v] == 1:
                    G.add_edge(u, v)
    cycles = nx.cycle_basis(G)
    
    is_aromatic = np.full(n_at, False)
    bound_to_aromatic = np.full(n_at, False)
    for cycle in cycles:
        for i in cycle:
            is_aromatic[i] = True #todos los atomos del ciclo son aromaticos 
            for j, bond in enumerate(bin_dist[i]):
                if bond == 1:
                    bound_to_aromatic[j] = True
    return is_aromatic, bound_to_aromatic

def bound_to_carbonyl(i, bin_dist, atom_types):
    b2CO = False
    for j, bond1 in enumerate(bin_dist[i]):
        if (bond1 == 1) and (atom_types[j] == 'Car'):
            for k, bond2 in enumerate(bin_dist[j]):
                if bond2 == 1:
                    if (atom_types[k] == 'Osp2c'):
                        b2CO = True
                        return b2CO
                    else:
                        b2CO = None
    return b2CO

def bound_to_H(i, bin_dist):
    bounded = False
    for j, bond in enumerate(bin_dist[i]):
        if (bond == 1) and (species[j] == 1):
            bounded = True
            return bounded
    return bounded

                    
def asign_atom_types(species, coord):
    n_at = len(species)
    atom_types = [None]*n_at
    bin_dist, n_neighbors = get_neighbors(coord, species, 3, 2, True) 
    bound_to_het = np.full(n_at, False)
    is_aromatic, bound_to_aromatic = find_aromatics(bin_dist, n_at) 
    for i in range(n_at):
        if (species[i] == 7) or (species[i] == 8):
            for j, bond in enumerate(bin_dist[i]):
                if bond==1:
                    bound_to_het[j] = True
    for i in range(n_at):
        if species[i] == 6: #Carbon atom
            if n_neighbors[i] == 2:
                atom_types[i] = 'Csp'
            elif n_neighbors[i] == 3:
                atom_types[i] = 'Car'
            elif n_neighbors[i] == 4:
                if bound_to_het[i]: 
                    atom_types[i] = 'Cpolar'
                else:
                    atom_types[i] = 'Cnonpol'
            else:
                atom_types[i] = 'Cvac'
        elif species[i] == 1: # Hidrogen atom
            if bound_to_aromatic[i]: 
                atom_types[i] = 'Har'
            elif bound_to_het[i]: 
                atom_types[i] = 'Hpolar'
            else:
                atom_types[i] = 'Hnonpol'
        elif species[i] == 8: #Oxygen atom
            j = i 
            while j < n_at:
                if n_neighbors[j] == 1 :
                    atom_types[j] = 'Osp2c'
                j = j + 1
            if atom_types[i] == None: #No lo asigne en el loop anterior
                if bound_to_het[i]:
                    atom_types[i] = 'Opolar'
                elif is_aromatic[i]:
                    atom_types[i] = 'Oar'
                b2CO = bound_to_carbonyl(i, bin_dist, atom_types)
                if b2CO == True:
                    if bound_to_H(i, bin_dist) == True:
                        atom_types[i] = 'Ohacid'
                    else:
                        atom_types[i] = 'Oar'
                elif b2CO == False:
                    atom_types[i] = 'Ononpol'
        elif species[i] == 7: #Nitrogen atom
            if n_neighbors[i] == 1:
                atom_types[i] = 'Nsp'
            elif bound_to_het[i]: 
                atom_types[i] = 'Npolar'
            elif is_aromatic[i]:
                atom_types[i] = 'Nar'
            elif bound_to_carbonyl(i, bin_dist, atom_types):
                atom_types[i] = 'Nsp2'
            else:
                atom_types[i] = 'Nnonpol'
    return atom_types
