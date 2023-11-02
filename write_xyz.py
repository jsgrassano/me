#---------------------------------------------------------------------
import os
import h5py
import numpy as np
import subprocess
#---------------------------------------------------------------------

def ejecutar(comando):
    comando=comando.split()
    subprocess.call(comando)
  
#CARGA LA BASE DE DATOS
db = h5py.File('/home/jota/DataSets/ani1x-release.h5', 'r')
#periodic_table = {'H':1, 'He':2, 'Li':3, 'Be':4, 'B':5, 'C':6,  'N':7,  'O':8,  'F':9,
#                 'Ne':10,'Na':11,'Mg':12,'Al':13,'Si':14,'P':15,'S':16,'Cl':17,'Ar':18}

periodic_table = {'1':'H', '6':'C', '7':'N', '8':'O'}

#os.makedirs(path+str(formula), exist_ok=True)  #GENERO UNA CARPETA PARA CADA FORMULA EN EL ARCHIVO formulas.txt 

def writeXYZ(formula, index, database):
    atoms =  database[formula+'/atomic_numbers'] #EXTRAIGO LOS NUMEROS ATÓMICOS
    coord = np.array(database[formula+'/coordinates'][index]) #GUARDO LA I-ESIMA ESTRUCTURA
    with open(formula+"_"+str(index)+'.xyz', "w+") as XYZ: #ABRO UN ARCHIVO VACIO
        XYZ.write(str(len(atoms))+'\n'+'\n')
        for i, atom in enumerate(atoms):
            el = periodic_table[str(atom)]
            #list(periodic_table.keys())[list(periodic_table.values()).index(atom)] #EXTRAIGO EL SIMBOLO DEL ATOMO
            line = f'{el} {coord[i][0]} {coord[i][1]} {coord[i][2]} \n'
            #atom_charac+" "+str(float(coord[i][0]))+" "+str(float(coord[i][1]))+" "+str(float(structure[i][2]))+"\n"
            XYZ.write(line)
    XYZ.close()
    return

writeXYZ('H2O1',0,db)
