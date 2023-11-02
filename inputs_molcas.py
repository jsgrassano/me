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
periodic_table = {'H':1, 'He':2, 'Li':3, 'Be':4, 'B':5, 'C':6,  'N':7,  'O':8,  'F':9,
                 'Ne':10,'Na':11,'Mg':12,'Al':13,'Si':14,'P':15,'S':16,'Cl':17,'Ar':18}

formulas = []
indices = []
with open("lista_estructuras_corridas.txt", "r") as f:
    for line in f:
        formulas.append(line.split()[0])
        indices.append(line.split()[1])
f.close()

path = '/home/santi/0-Pol/molecules/'

#GENERO LOS ARCHIVOS .XYZ DE TODAS LAS ESTRUCTURAS CON CADA FORMULA MOLECULAR
for line, formula in enumerate(formulas[0:20]): #in formulas  #ITERO SOBRE LAS FORMULAS MOLECULARES
    index = int(indices[line])
    os.makedirs(path+str(formula), exist_ok=True)  #GENERO UNA CARPETA PARA CADA FORMULA EN EL ARCHIVO formulas.txt 
#    index_list = # len(database[formula+'/coordinates']) #ITERO SOBRE TODAS LAS ESTRUCTURAS DE CADA FORMULA
    database = h5py.File('/home/santi/Descargas/ani1x-release.h5', 'r')
    atoms =  database[formula+'/atomic_numbers'] #EXTRAIGO LOS NUMEROS ATÓMICOS
    structure = np.array(database[formula+'/coordinates'][index]) #GUARDO LA I-ESIMA ESTRUCTURA
    
    with open(path+formula+"/"+formula+"_"+str(index)+'.xyz', "w") as XYZ: #ABRO UN ARCHIVO VACIO
        XYZ.write(str(len(atoms))+'\n'+'\n')
        for i, atom in enumerate(atoms):
            atom_charac = list(periodic_table.keys())[list(periodic_table.values()).index(atom)] #EXTRAIGO EL SIMBOLO DEL ATOMO
            a_xyz = atom_charac+" "+str(float(structure[i][0]))+" "+str(float(structure[i][1]))+" "+str(float(structure[i][2]))+"\n"
            XYZ.write(a_xyz)
    XYZ.close()
    database.close()

# generar un input de molcas
    with open(path+formula+"/"+formula+"_"+str(index)+'.input', "w") as MOLCAS_INPUT: #ABRO UN ARCHIVO VACIO
        MOLCAS_INPUT.write('&GATEWAY\n')
        MOLCAS_INPUT.write('coord='+path+formula+"/"+formula+"_"+str(index)+'.xyz'+'\n')
        MOLCAS_INPUT.write('basis=ANO-L\ngroup=NOSYM\n&SEWARD\n&SFC\nKSDFT=PBE0\nLOPROP')
    MOLCAS_INPUT.close()        
        # correr molcas
     #   ejecutar('export MOLCAS=/home/santi/programas/molcas/build')
     #   ejecutar('export PATH=$PATH:$MOLCAS')
     #   ejecutar('pymolcas '+path+formula+"/"+formula+"_"+str(index)+'.input -f')
        
#escribo script de corrida ?

    with open(path+formula+"/run_"+formula+"_"+str(index)+'.sh', "w") as sh:
        sh.write('#!/bin/sh\n\n')
        sh.write('#/opt/intel/oneapi/compiler/2021.2.0/linux/bin/intel64/libcilkrts.so.5\n')
        sh.write('#/opt/intel/oneapi/compiler/2022.1.0/linux/bin/intel64/libcilkrts.so.5\n')
        sh.write('#/opt/intel/oneapi/compiler/2022.1.0/linux/compiler/lib/intel64_lin/libcilkrts.so.5\n')
        sh.write('export LIBEXTRA=/snap/gnome-3-34-1804/72/usr/lib/x86_64-linux-gnu\n') # estas cosas preguntarle a joni si sirven para cecar
        sh.write('export MOLCAS=/home/jota/molcas/build\n')
        sh.write('export PATH=$PATH:$MOLCAS\n')
        sh.write('export PATH=$PATH:$LIBEXTRA\n')
        sh.write('export MOLCAS_MEM=6Gb\n')
        sh.write('export MOLCAS_DISK=50Gb\n')
        sh.write('export MOLCAS_NPROCS=2\n\n')
        sh.write('pymolcas C3H7N1O3_0.input -f')
    sh.close()

# ejecuto el script q cree recien :ssssss

    ejecutar('\.'+path+formula+"/run_"+formula+"_"+str(index)+'.sh'
