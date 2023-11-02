import numpy as np
import os
import h5py
import math
import sys
import parmed as pmd
from parmed.amber import Rst7

from functions import *

database = h5py.File('/home/jota/DataSets/ani1x-release.h5')
dicZ=[("H", 1),("He", 2), ("Li",3), ("Be",4), ("B",5), ("C",6), ("N",7), ("O",8), ("F",9),("Ne",10)]

formulas = np.array(['C2H5N1O2'], dtype=str) 
cantidades = np.array([13717]) 
cutoff=15.0

for i in range(len(formulas)):
    formula = formulas[i]
    cantidad = cantidades[i]
    sqm = database[formula+'/atomic_numbers']
    nqm = sqm.shape[0]
   # if os.path.exists(f'energia_{formula}.csv'):
   #     os.remove(f'energia_{formula}.csv')
    energia = np.zeros((cantidad+1,16)) 
    lista_indices_file=f'../2-Estructuras/indices/indices_{formula}'
    lista_indices = np.loadtxt(lista_indices_file, dtype=int)
    for indice in lista_indices[10000:]:
        formula_indice = f'{formula}_{indice}'
        path = f'../4-outputs-PBE/{formula}/{formula}_{indice}'
        if os.path.exists(path+'_sol_MK'):
            path_top = f'/media/jota/cachirulo/1-ME/2-Estructuras/inputs_dft/{formula}/{formula}_{indice}_vac.prmtop'
            
            r, q_mul_sol, q_mul_vac = sort_r_and_q(formula, indice, nqm, path_top)

            q_mul_sol_dftb, q_cm3_sol_dftb, q_mul_vac_dftb, q_cm3_vac_dftb = sort_q_dftb(formula, indice, nqm, path_top)
            q_mul_sol_am1, q_mul_vac_am1 = sort_q_am1(formula, indice, nqm, path_top)
            q_hir = np.array(database[formula+'/wb97x_dz.hirshfeld_charges'][indice])
            q_cm5 = np.array(database[formula+'/wb97x_dz.cm5_charges'][indice])
            q_mbis= np.array(database[formula+'/wb97x_tz.mbis_charges'][indice])
         
            wat = int((r.shape[0] - nqm)/3)
            
            q_hir     = append_water_charges(q_hir,     wat)
            q_cm5     = append_water_charges(q_cm5,     wat)
            q_mbis    = append_water_charges(q_mbis,    wat)
            
            q_mul_vac = append_water_charges(q_mul_vac, wat)
            q_mul_vac_dftb = append_water_charges(q_mul_vac_dftb, wat)
            q_mul_vac_am1  = append_water_charges(q_mul_vac_am1, wat)
            
            q_mul_sol = append_water_charges(q_mul_sol, wat)
            q_mul_sol_dftb = append_water_charges(q_mul_sol_dftb, wat)
            q_mul_sol_am1  = append_water_charges(q_mul_sol_am1, wat)
            
            q_cm3_vac = append_water_charges(q_cm3_vac_dftb, wat)
            q_cm3_sol = append_water_charges(q_cm3_sol_dftb, wat)

            indice_a = atom_index_in_cutoff(r, nqm, cutoff)
        
            V_QM = V_elec(r, sqm, q_hir, indice_a)

            E_hir = E_elec(q_hir, V_QM)
            E_cm5 = E_elec(q_cm5, V_QM)
            E_mbis = E_elec(q_mbis, V_QM)
            E_mul_sol = E_elec(q_mul_sol, V_QM)
            E_mul_vac = E_elec(q_mul_vac, V_QM)
            
            E_mul_sol_dftb  = E_elec(q_mul_sol_dftb, V_QM)
            E_mul_sol_am1   = E_elec(q_mul_sol_am1, V_QM)
            E_mul_vac_dftb  = E_elec(q_mul_vac_dftb, V_QM)
            E_mul_vac_am1   = E_elec(q_mul_vac_am1, V_QM)
            E_cm3_sol =  E_elec(q_cm3_sol, V_QM)
            E_cm3_vac =  E_elec(q_cm3_vac, V_QM)

# Polarization energy        
            pol_1 = get_pol_1(sqm, pol_table_1)
            pol_2 = get_pol_lit(atom_types, pol_table_litman)

            Elec_field_QM = E_field(r, nqm, q_hir, indice_a)
            
            E_p_1 = E_pol(nqm, pol_1, Elec_fiel_QM)
            E_p_lit = E_pol(nqm, pol_2, Elec_field_QM)
# Electrostatic embedding        
            E_QM_sol, E_dft_sol = e_dft(path+'_sp_sol_DFT.out')
            E_QM_vac, E_dft_vac = e_dft(path+'_sp_vac_DFT.out')
            
            energia[indice] = np.array([indice, E_dft_sol, E_hir,E_cm5,E_mbis,E_mul_sol,E_mul_vac, E_mul_sol_dftb, E_mul_sol_am1, E_mul_vac_dftb, E_mul_vac_am1, E_cm3_sol, E_cm3_vac, E_p, E_QM_sol, E_QM_vac])
        else:
            print(f'no entre en {formula} {indice}')
    f1 = open(f'7-energies/energies_{formula}.csv', 'a')
    np.savetxt(f1, energia, delimiter=';')
    f1.close()
    
    f2 = open(f'electrostatics_{formula}.csv', 'a')
    np.savetxt(f2, , delimiter=';')
    f2.close()
    
