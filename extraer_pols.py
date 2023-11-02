import numpy as np
import h5py

def extract_molcas_info(filename):
    #RECIBE EL PATH CON EL NOMBRE_INDICE DONDE ESTÁN LOS OUTPUTS DE MOLCAS ej '/C3H4N2O3/C3H4N2O3_0'
    path_to_log = filename+'.log'
    path_to_scf = filename+'.scf.molden'
    #las cargas de mulliken las saco de este archivo porque es mas facil :p

# GENERO LAS VARIABLES DONDE VOY A GUARDAR LA INFO
    molec_dip_magnitude = []
    loprop_molecular_dipole_vector = [0,0,0]
    loprop_molecular_polarizability_tensor = []
    atoms_ = []
    bonds_ = []
    mulliken_charges = []
    loprop_charges = []
    loprop_dipoles_magnitudes = []
    loprop_dipoles_vectors = []
    loprop_atomic_polarizability = []
    loprop_bond_polarizability = []

# ABRO EL ARCHIVO .LOG
    with open(path_to_log, 'r') as f:
# DEFINO ESTOS CONDICIONALES PARA ARRANCAR A LEER
        molecular_prop = False
        atomic_domain = True
        mol_pol_tensor= False
        mulliken = False
        linea_ref_1=0
        linea_ref_2=0
        linea_ref_3=0
        atoms_dipoles = False

        for l_n, line in enumerate(f): #ITERO SOBRE LAS LINEAS DEL ARCHIVO
            line.replace(","," ")

#  ENERGIA
            if 'Total KS-DFT energy' in line:
                energy = float(line.split()[-1])

#  PROPIEDADES MOLECULARES
            if 'Molecular properties' in line:
                molecular_prop = True
                linea_ref_1 = l_n #GUARDO LA LINEA COMO REFERENCIA

##  DIPOLO MOLECULAR (VECTOR Y MAGNITUD)
            if l_n == linea_ref_1 + 7 and molecular_prop == True: # 7 LINEAS MAS ABAJO ESTAN LOS VALORES
                loprop_molecular_dipole_vector[0] = float(line.split()[-7])
                loprop_molecular_dipole_vector[1] = float(line.split()[-5])
                loprop_molecular_dipole_vector[2] = float(line.split()[-3])

                molec_dip_magnitude = float(line.split()[-1])

# PROPIEDADES LOCALIZADAS
            if 'ATOMIC DOMAIN' in line:
                atoms_.append(line.split()[-1])
                atomic_domain = True
                bond_domain = False
            if 'BOND DOMAIN' in line:
                bonds_.append([line.split()[-2], line.split()[-1][1:]])
                bond_domain = True
                atomic_domain = False

# ATOMIC DIPOLE VECTOR
            if 'Electronic Dipole' in line:
                atoms_dipoles = True
                linea_ref_2 = l_n
            if l_n == linea_ref_2 + 3 and atoms_dipoles == True and atomic_domain== True:
                vector = [float(i) for i in line.split()]
                loprop_dipoles_vectors.append(vector)

# DIPOLE MAGNITUDE
            if 'Dipole magnitude' in line and atoms_dipoles == True and atomic_domain== True:
                loprop_dipoles_magnitudes.append(float(line.split()[-1]))

# CARGA LOPROP POR ATOMO
            if 'Total charge' in line:
                charge = float(line.split()[-1])
                if atomic_domain == True:
                    loprop_charges.append(charge)

# POLARIZABILIDAD POR DOMINIO
            if 'Isotropic' in line:
                polarizability = float(line.split()[-1])
                if atomic_domain == True:
                    loprop_atomic_polarizability.append(polarizability)
                elif atomic_domain == False:
                    loprop_bond_polarizability.append(polarizability)

            if 'Molecular Polarizability Tensor' in line:
                mol_pol_tensor = True
                linea_ref_3=l_n
            if linea_ref_3+2<l_n and l_n<linea_ref_3 + 6 and mol_pol_tensor == True:
                loprop_molecular_polarizability_tensor= loprop_molecular_polarizability_tensor +[float(i) for i in line.split()]
    f.close()

# CARGAS MULLIKEN
    with open(path_to_scf, 'r') as s:
        mulliken = False
        for l_n, line in enumerate(s): #ITERO SOBRE LAS LINEAS DEL ARCHIVO
            line.replace(","," ")
            if '[GTO]' in line:
                mulliken = False
            if mulliken == True:
                mulliken_charges.append(float(line.split()[0]))
            if '(Mulliken)' in line:
                mulliken = True
    s.close()

# SUMO A LA POLARIZABILIDAD ATOMICA 0.5 DE LA POLARIZABILIDAD DE LOS ENLACES

    loprop_total_atomic_polarizability=[0]*len(atoms_)
    for i, atom in enumerate(atoms_):
        loprop_total_atomic_polarizability[i] = loprop_atomic_polarizability[i]
        for j, bond in enumerate(bonds_):
            if atom in bond:
                loprop_total_atomic_polarizability[i] = loprop_total_atomic_polarizability[i] + loprop_bond_polarizability[j]*0.5

    loprop_molecular_dipole_vector = np.array(loprop_molecular_dipole_vector,  dtype=object)
    loprop_molecular_polarizability_tensor = np.array(loprop_molecular_polarizability_tensor,  dtype=object)
    mulliken_charges = np.array(mulliken_charges,  dtype=object)
    loprop_charges = np.array(loprop_charges,  dtype=object)
    loprop_dipoles_magnitudes = np.array(loprop_dipoles_magnitudes,  dtype=object)
    loprop_dipoles_vectors = np.array(loprop_dipoles_vectors,  dtype=object)
    loprop_atomic_polarizability = np.array(loprop_atomic_polarizability, dtype=object)
    loprop_bond_polarizability = np.array(loprop_bond_polarizability, dtype=object)
    loprop_total_atomic_polarizability = np.array(loprop_total_atomic_polarizability, dtype=object)

    info = [atoms_, bonds_, molec_dip_magnitude, loprop_molecular_dipole_vector,
            mulliken_charges, loprop_charges,  loprop_dipoles_magnitudes,
            loprop_dipoles_vectors, loprop_atomic_polarizability, loprop_bond_polarizability,
            loprop_total_atomic_polarizability]
    return info

def append_to_hdf5(formula): #recibe str tipo 'C2H5N1O2'
    for filename in lista_archivos: # LISTA ARCHIVOS DEBERIA SER TIPO ['C2H5N1O2_1', 'C2H5N1O2_4'] una lista de los archivos en ese directorio (sin la extensión)
    info = extract_molcas_info(filename)
    
    atoms_ = info[0]
    bonds_.append(info[1])
    molec_dip.magnitude.append(infor[2])
    loprop_molecular_dipole_vector.append(info[3])
    mulliken_charges.append(info[4])
    loprop_charges.append(info[5])
    loprop_dipoles_magnitudes.append(info[6])
    loprop_dipoles_vectors.append(info[7])
    loprop_atomic_polarizability.append(info[8])
    loprop_bond_polarizability.append(info[9])
    loprop_total_atomic_polarizability.append(info[10])
    
    with h5py.File('/home/santi/Descargas/ani1x_copia_prueba.h5', 'a') as db:
        db.create_dataset(formula+"/atomic_char", data = atoms_)
        db.create_dataset(formula+"/bonds", data=bonds_)
        db.create_dataset(formula+"/pbe0.loprop_molecular_dipole_magnitude", data = molec_dip.magnitude)
        db.create_dataset(formula+"/pbe0.loprop_molecular_dipole_vector" , data = bonds_)
        db.create_dataset(formula+"/pbe0.loprop_molecular_polarizability_tensor" , data = loprop_molecular_dipole_vector)
        db.create_dataset(formula+"/pbe0.mulliken_charges", data = mulliken_charges)
        db.create_dataset(formula+"/pbe0.loprop_charges", data = loprop_charges)
        db.create_dataset(formula+"/pbe0.loprop_dipoles_magnitudes", data = loprop_dipoles_magnitudes)
        db.create_dataset(formula+"/pbe0.loprop_dipoles_vectors", data = loprop_dipoles_vectors)
        db.create_dataset(formula+"/pbe0.loprop_atomic_polarizability", data =  loprop_atomic_polarizability)
        db.create_dataset(formula+"/pbe0.loprop_bond_polarizability", data = loprop_bond_polarizability)
        db.create_dataset(formula+"/pbe0.loprop_total_atomic_polarizability", data = loprop_total_atomic_polarizability)
    db.close()
    return
                          
