import numpy as np
import h5py
import re

from functions import *

formula = 'C4H8N2O3'
indices = np.loadtxt(f'../2-Estructuras/indices/indices_{formula}', dtype=int)
n_str = int(indices[-1])
n_atoms = sum([int(x) for x in re.findall(r'\d+', formula)])

AM1_mul_vac_charges_notsorted  = np.empty((n_str, n_atoms))
#AM1_mul_sol_charges_notsorted  = np.empty((n_str, n_atoms))
#DFTB_mul_vac_charges_notsorted = np.empty((n_str, n_atoms))
#DFTB_mul_sol_charges_notsorted = np.empty((n_str, n_atoms))
#DFTB_cm3_vac_charges_notsorted = np.empty((n_str, n_atoms))
#DFTB_cm3_sol_charges_notsorted = np.empty((n_str, n_atoms))
#PBE_mul_vac_charges_notsorted = np.empty((n_str, n_atoms))
#PBE_mul_sol_charges_notsorted = np.empty((n_str, n_atoms))

for indice in indices:
    try:
        AM1_mul_vac_charges_notsorted[indice] = load_q_am1(f'../3-outputs-EQ+SE/{formula}/charges_{formula}_{indice}_sp_vac_AM1.out', n_atoms)
#        AM1_mul_sol_charges_notsorted[indice] = load_q_am1(f'../3-outputs-EQ+SE/{formula}/charges_{formula}_{indice}_sp_sol_AM1.out', n_atoms)
#        DFTB_mul_vac_charges_notsorted[indice], DFTB_cm3_vac_charges_notsorted[indice] = load_q_dftb(f'../3-outputs-EQ+SE/{formula}/charges_{formula}_{indice}_sp_vac_DFTB.out', n_atoms)
#        DFTB_mul_sol_charges_notsorted[indice], DFTB_cm3_sol_charges_notsorted[indice] = load_q_dftb(f'../3-outputs-EQ+SE/{formula}/charges_{formula}_{indice}_sp_sol_DFTB.out', n_atoms)
#        PBE_mul_vac_charges_notsorted[indice] = load_MK(f'../4-outputs-PBE/{formula}/{formula}_{indice}_vac_MK', n_atoms) 
#        PBE_mul_sol_charges_notsorted[indice] = load_MK(f'../4-outputs-PBE/{formula}/{formula}_{indice}_sol_MK', n_atoms)
    except:
        pass

path_top = f'/media/jota/cachirulo/1-ME/2-Estructuras/inputs_dft/{formula}/{formula}_0_vac.prmtop'
AM1_mul_vac_charges_sorted = sort_q_array(AM1_mul_vac_charges_notsorted, n_atoms, path_top)

