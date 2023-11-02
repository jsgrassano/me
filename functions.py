import numpy as np
import parmed as pmd
from parmed.amber import Rst7
import networkx as nx
import h5py
from scipy.spatial.distance import cdist
#########################################################################################################################
#DICCIONARIOS
##########################################################################################################################
periodic_table = {'H':1, 'HA':1, 'HB2':1,'HB3':1, 'He':2,'Li':3,
                   'Be':4,'B':5,'C':6,'CA':6,'CB':6,'N':7,'O':8,
                   'F':9,'Ne':10,'Na':11,'Mg':12,'Al':13,'Si':14,
                   'P':15,'S':16,'SG':16,'Cl':17,'Ar':18}

pol_table_1 = {
        '1':4.50711, # H 
        '6':11.3,    # C
        '7':7.4,     # N 
        '8':5.3,     # O
#https://doi.org/10.1080/00268976.2018.1535143
pol_table_litman = {
"HW":2.89031169275157;
"Har":2.91359344698923;
"Hnonpol":3.24149502333945;
"Hpolar":3.08601339503921;
"Car":13.93200418585700;
"Cch4":9.05106876411203;
"Cnonpol":9.54882592137880;
"Cpolar":10.92987258941900;
"Csp":14.86219462183380;
"Nnh3":7.79932018627293;
"Nar":11.48458569038600;
"Nnonpol":7.99448202179564;
"Npolar":9.74284054002603;
"Nsp":7.81450393903663;
"Nsp2":8.39006939379913;
"OW":6.58873644925985;
"Oar":5.61636892227275;
"Ononpol":5.48126726434866;
"Opolar":5.78689933664546;
"Ohacid":5.61684130569207;
"Osp2c":6.16642567239225;
"O−sulf":5.71253268977616;
"O=sulf":5.79546972153875}
# https://doi.org/10.1021/acs.jcim.1c01307

dicZ=[("H", 1),("He", 2), ("Li",3), ("Be",4), ("B",5), ("C",6), ("N",7), ("O",8), ("F",9),("Ne",10)]

def append_water_charges(charges, wat):
    charges_h20 = np.array([-0.834,0.417,0.417]*wat)
    total_charges = np.concatenate((charges, charges_h20))
    return total_charges

##########################################################################################################################
# FUNCIONES PARA COMPUTAR ENERGIAS 
##########################################################################################################################
def in_cutoff(r,nqm,j,cutoff):
    esta=False
    i=0
    while (not esta and i<nqm):
        distij = np.linalg.norm(r[j] - r[i]) 
        esta = esta or distij < cutoff
        i=i+1
    return esta

def atom_index_in_cutoff(r,atoms_qm,cutoff):
    lista = np.zeros(len(r), dtype=bool)
    for j in range(atoms_qm,len(r)):
        if (in_cutoff(r,atoms_qm,j,cutoff)):
            lista[i] = True
    return np.array(lista)

def E_elec(q, V_QM):
    nqm = len(V_QM) #numeros de atomos qm
    E_elec = 0.0
    for i in range(nqm):
        E_elec = E_elec + q[i]*V_QM[i]
    E_elec = E_elec*0.529177249*627.5094740631
    return E_elec 

def V_elec(r, sqm, q, in_cutoff):
    nqm = len(sqm)
    ntot = len(r)
    V_QM = np.empty(nqm)
    for i in range(nqm):
        V_i = 0.0
        for j in range(nqm, ntot):
            if in_cutoff[j]:
                distij = np.linalg.norm(r[j]-r[i])
                V_i = V_i + q[j]/distij
        V_QM[i] = V_i
    return V_QM

def E_pol(nqm, p, E_field_QM):
    E_pol = 0.0
    for i in range(nqm):
        alpha=p[i]
        E_field_ri = E_field_QM[i]
        E_field_2 = E_field_ri[0]**2 + E_field_ri[1]**2 + E_field_ri[2]**2
        E_pol = E_pol + alpha*E_field_2
    E_pol = -0.5*0.529177249*627.5094740631*E_pol
    return E_pol

def E_field(r, nqm, q, in_cutoff):
    ntot = len(r)
    E_field_QM = np.empty((nqm,3))
    for i in range(nqm):
        E_field_ri = np.zeros(3)
        for j in range(nqm, ntot):
            if in_cutoff[j]:
                distij = np.linalg.norm(r[j] - r[i])
                rij = (r[j]-r[i])/distij
                E_field_ri = E_field_ri + (q[j]/distij**2) * rij
        E_field_QM[i] = E_field_ri
    return E_field_QM

def get_pol_1(species, pol_table_1):
    p = []
    for z in sqm:
        alpha = 0.14818471*pol_table_1.items(str(z))
        p.append(alpha)
    p = np.array(p)
    return p

def get_pol_lit(atom_types, pol_table_litman):
    p = []
    for atom in atom_types:
        alpha = 0.14818471*pol_table_litman(atom)
        p.append(alpha)
    p = np.array(p)
    return p

##############################################################################################################
# FUNCIONES PARA ASIGNAR ATOM TYPES ##########################################################################
##############################################################################################################

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

###############################################################################
# FUNCIONES PARA CARGAR INFO DE ARCHIVOS ######################################
###############################################################################

def load_MK(nombre_archivo, nqm):
    inicio = 7 + nqm
    fin = 8 + 2*nqm
    q_ = []
    with open(nombre_archivo, 'r') as f:
        for renglon, line in enumerate(f):
            if renglon > inicio and renglon < fin:
                q = line.split()[-1]
                q_.append(float(q))
    q_ = np.array(q_)
    return q_ 

def load_q_dftb(nombre_archivo, nqm):
    q_mul = []
    q_cm3 = []
    inicio=0
    fin=0
    with open(nombre_archivo, 'r') as SE:
        for renglon, line in enumerate(SE):
            if 'Atom' in line:
                inicio = renglon
                fin = inicio + nqm
            if renglon > inicio and renglon <= fin: 
                q_m = line.split()[2]
                q_c = line.split()[3]
                q_mul.append(float(q_m))
                q_cm3.append(float(q_c))
    return q_mul, q_cm3

def load_q_am1(nombre_archivo, nqm):
    q_mul = []
    inicio=0
    fin=0
    with open(nombre_archivo, 'r') as SE:
        for renglon, line in enumerate(SE):
            if 'Atom' in line:
                inicio = renglon
                fin = inicio + nqm
            if renglon > inicio and renglon <= fin: 
                q_m = line.split()[2]
                q_mul.append(float(q_m))
    return q_mul

def sort_r_and_q(formula, indice, nqm, path_top):
    q_mul_sol = load_MK(f'../4-outputs-PBE/{formula}/{formula}_{indice}_sol_MK', nqm)
    q_mul_vac = load_MK(f'../4-outputs-PBE/{formula}/{formula}_{indice}_vac_MK', nqm)
    topfile = pmd.load_file(path_top)
    listZ_top=[]
    listS_top=[]
    for i in range(0,nqm):
        listZ_top.append(topfile.atoms[i].atomic_number)
        listS_top.append(dicZ[listZ_top[i]-1][0])
    listS_ord=sorted(listS_top)
    listZ_ord=[]
    for i in range(0,len(listS_ord)):
        for j in range(0,len(dicZ)):
            if (dicZ[j][0]==listS_ord[i]):
                listZ_ord.append(dicZ[j][1])
    a=''.join(listS_ord)
    q_mul_sol_ord = []
    q_mul_vac_ord = []
    rst7file=Rst7(f'../3-outputs-EQ+SE/{formula}/{formula}_{indice}_sol_term.rst7')
    coordenadas=rst7file.coordinates[:1][0]
    coordQM=coordenadas[0:nqm]
    coordQM_ord=[]
        
    for i in range(0,len(listZ_ord)):
        j=0
        encontre=False
        while (j<len(listZ_top) and not encontre):
            if (listZ_top[j]==listZ_ord[i]):
                q_mul_sol_ord.append(q_mul_sol[j])
                q_mul_vac_ord.append(q_mul_vac[j])
                q_mul_sol = np.delete(q_mul_sol, j, axis=0)
                q_mul_vac = np.delete(q_mul_vac, j, axis=0)

                coordQM_ord.append(coordQM[j])
                coordQM = np.delete(coordQM, j,axis=0)

                listZ_top = np.delete(listZ_top, j,axis=0)
                encontre = True
            j=j+1
    # cargar en un tensor las coordenadas ordenadas
    coordQM_ord = np.array(coordQM_ord)
    coordMM = coordenadas[nqm:]
    r = np.concatenate((coordQM_ord, coordMM))
    q_mul_sol_ord = np.array(q_mul_sol_ord)
    q_mul_vac_ord = np.array(q_mul_vac_ord)
    return r, q_mul_sol_ord, q_mul_vac_ord

def sort_q_array(q_not_ord, nqm, path_top):
    topfile = pmd.load_file(path_top)
    listZ_top = [topfile.atoms[i].atomic_number for i in range(nqm)]
    listS_top = [dicZ[z-1][0] for z in listZ_top]
    listS_ord = sorted(listS_top)
    listZ_ord = [dicZ[j][1] for i in listS_ord for j in range(len(dicZ)) if dicZ[j][0]==i]
    a=''.join(listS_ord)
    q_ord = np.empty(q_not_ord.shape)
    for i in range(0,len(listZ_ord)):
        j=0
        encontre=False
        while (j<len(listZ_top) and not encontre):
            if (listZ_top[j]==listZ_ord[i]):
                q_ord[:,i] = q_not_ord[:,j]
                q_not_ord = np.delete(q_not_ord, j,axis=1)
                listZ_top = np.delete(listZ_top, j,axis=0)
                encontre = True
            j=j+1
    return q_ord

def sort_q_dftb(formula, indice, nqm, path_top):
        q_mul_sol, q_cm3_sol = load_q_dftb(f'/media/jota/cachirulo/1-ME/3-outputs-EQ+SE/{formula}/charges_{formula}_{indice}_sp_sol_DFTB.out', nqm)
        q_mul_vac, q_cm3_vac = load_q_dftb(f'/media/jota/cachirulo/1-ME/3-outputs-EQ+SE/{formula}/charges_{formula}_{indice}_sp_vac_DFTB.out', nqm) 
        topfile = pmd.load_file(path_top)
        listZ_top=[]
        listS_top=[]
        for i in range(0,nqm):
            listZ_top.append(topfile.atoms[i].atomic_number)
            listS_top.append(dicZ[listZ_top[i]-1][0])
        listS_ord=sorted(listS_top)
        listZ_ord=[]
        for i in range(0,len(listS_ord)):
            for j in range(0,len(dicZ)):
                if (dicZ[j][0]==listS_ord[i]):
                    listZ_ord.append(dicZ[j][1])
        a=''.join(listS_ord)
        q_mul_sol_ord = []
        q_cm3_sol_ord = []
        q_mul_vac_ord = []
        q_cm3_vac_ord = []
        for i in range(0,len(listZ_ord)):
            j=0
            encontre=False
            while (j<len(listZ_top) and not encontre):
                if (listZ_top[j]==listZ_ord[i]):
                    q_mul_sol_ord.append(q_mul_sol[j])
                    q_cm3_sol_ord.append(q_cm3_sol[j])
                    q_mul_vac_ord.append(q_mul_vac[j])
                    q_cm3_vac_ord.append(q_cm3_vac[j])
                    q_mul_sol = np.delete(q_mul_sol, j, axis=0)
                    q_cm3_sol = np.delete(q_cm3_sol, j, axis=0)
                    q_mul_vac = np.delete(q_mul_vac, j, axis=0)
                    q_cm3_vac = np.delete(q_cm3_vac, j, axis=0)
                    listZ_top = np.delete(listZ_top, j,axis=0)
                    encontre = True
                j=j+1
        q_mul_sol = np.array(q_mul_sol_ord)
        q_cm3_sol = np.array(q_cm3_sol_ord)
        q_mul_vac = np.array(q_mul_vac_ord)
        q_cm3_vac = np.array(q_cm3_vac_ord)
        return q_mul_sol, q_cm3_sol, q_mul_vac, q_cm3_vac

def sort_q_am1(formula, indice, nqm, path_top):
        q_mul_sol = load_q_am1(f'/media/jota/cachirulo/1-ME/3-outputs-EQ+SE/{formula}/charges_{formula}_{indice}_sp_sol_AM1.out', nqm)
        q_mul_vac = load_q_am1(f'/media/jota/cachirulo/1-ME/3-outputs-EQ+SE/{formula}/charges_{formula}_{indice}_sp_vac_AM1.out', nqm) 
        topfile = pmd.load_file(path_top)
        listZ_top=[]
        listS_top=[]
        for i in range(0,nqm):
            listZ_top.append(topfile.atoms[i].atomic_number)
            listS_top.append(dicZ[listZ_top[i]-1][0])
        listS_ord=sorted(listS_top)
        listZ_ord=[]
        for i in range(0,len(listS_ord)):
            for j in range(0,len(dicZ)):
                if (dicZ[j][0]==listS_ord[i]):
                    listZ_ord.append(dicZ[j][1])
        a=''.join(listS_ord)
        q_mul_sol_ord = []
        q_mul_vac_ord = []
        for i in range(0,len(listZ_ord)):
            j=0
            encontre=False
            while (j<len(listZ_top) and not encontre):
                if (listZ_top[j]==listZ_ord[i]):
                    q_mul_sol_ord.append(q_mul_sol[j])
                    q_mul_vac_ord.append(q_mul_vac[j])
                    q_mul_sol = np.delete(q_mul_sol, j, axis=0)
                    q_mul_vac = np.delete(q_mul_vac, j, axis=0)
                    listZ_top = np.delete(listZ_top, j,axis=0)
                    encontre = True
                j=j+1
        q_mul_sol = np.array(q_mul_sol_ord)
        q_mul_vac = np.array(q_mul_vac_ord)
        return q_mul_sol, q_mul_vac

def e_dft(path_to_dft):
    nuc = 0.0
    elec = 0.0
    with open(path_to_dft, 'r') as f:
        for line in f:
            if 'Total energy = ' in line: #ENERGIA TOTAL (QM + QMMM sin LJ)
                linea = line.split(' ')
                tot = float(linea[-1])*627.509391
            elif 'QM-MM nuc.' in line:
                linea = line.split(' ')
                nuc = float(linea[-1])
            elif 'QM-MM elec.' in line:
                linea = line.split(' ')
                elec = float(linea[-1])
    QMMM = (nuc + elec)*627.509391
    QM = tot - QMMM
    return QM, QMMM
