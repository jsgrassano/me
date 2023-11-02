#---------------------------------------------------------------------
import torch
import numpy as np
import h5py
import random
import math
import parmed as pmd
from parmed.amber import Rst7

from functions import *

#---------------------------------------------------------------------
#CARGA LA BASE DE DATOS
database = h5py.File('C:/Users/juans/OneDrive/Escritorio/Labo/ani1x-release.h5', 'r')

# FORMULAS DE LOS COMPUESTOS QUE SE QUIEREN ANALIZAR
#formulas ya corridas "C2H5N1O2", "C3H7N1O2", "C5H11N1O2", "C6H13N1O2", "C9H11N1O2", "C9H11N1O3", "C11H12N2O2", "C3H7N1O3", "C4H9N1O3", "C6H9N3O2",
# "C6H14N2O2" , "C6H14N4O2" , "C4H8N2O3"
# "C6H13N1O2" estaba dos veces, deje solo una
#[ "C5H10N2O3"],  "C5H9N1O2"]
formulas = ['C4H7N1O4','C5H9N1O4']

#CANTIDAD DE ESTRUCTURAS A EXTRAER DE CADA FORMULA
N = 1

#CANTIDAD DE AGUAS EN EL SISTEMA MM
wat = 4000

#---------------------------------------------------------------------
random.seed(40118)
#---------------------------------------------------------------------
dicZ=[("H",1),("He",2),("Li",3),("Be",4),("B",5),("C",6),("N",7),("O",8),("F",9),("Ne",10)]


for formula in formulas:
    index_list = random.sample(range(0,len(database[formula+'/coordinates'])),N)
    for index in index_list:
        if not math.isnan(database[formula+'/wb97x_tz.mbis_charges'][index][0]):
            structure = torch.tensor(database[formula+'/coordinates'][index])
            nombretop=formula+".prmtop"
            topfile = pmd.load_file(nombretop)
            NQM=len(topfile.residues[0])
            listZ_top=[]
            listS_top=[]
            for i in range(0,NQM):
                listZ_top.append(topfile.atoms[i].atomic_number)
                listS_top.append(dicZ[listZ_top[i]-1][0])

            listS_rst7=sorted(listS_top)
            listZ_rst7=[]
            for i in range(0,len(listS_rst7)):
                for j in range(0,len(dicZ)):
                    if (dicZ[j][0]==listS_rst7[i]):
                        listZ_rst7.append(dicZ[j][1])

            a=''.join(listS_top)
            Z_in_top=[]
            for i in range(0,len(a)):
                for j in range(0,len(dicZ)):
                    if (dicZ[j][0]==a[i]):
                        Z_in_top.append(dicZ[j][1])

            a=sorted(a)
            Z_in_rst7=[]
            for i in range(0,len(a)):
                for j in range(0,len(dicZ)):
                    if (dicZ[j][0]==a[i]):
                        Z_in_rst7.append(dicZ[j][1])

            coordQM = structure.numpy()

            ordenada=[]
            
            #EXTRAIGO CARGAS DE LA BASE DE DATOS, GENERO UNA LISTA VACIA PARA ORDENARLAS
            cargas_mbis=database[formula+'/wb97x_tz.mbis_charges'][index]
            cargas_ordenadas = []
            
            for i in range(0,len(Z_in_top)):
                j=0
                encontre=False
                while (j<len(Z_in_rst7) and not encontre):
                    if (Z_in_rst7[j]==Z_in_top[i]):
                        cargas_ordenadas.append(cargas_mbis[j])
                        cargas_mbis = np.delete(cargas_mbis, j, axis=0)
                        ordenada.append(coordQM[j])
                        coordQM = np.delete(coordQM, j,axis=0)
                        Z_in_rst7 = np.delete(Z_in_rst7, j,axis=0)
                        encontre = True
                    j=j+1
                    
            ordenada=np.array(ordenada)
            coordQM = ordenada

        #se guarda .rst7 del soluto

            rst7 = Rst7(natom=len(coordQM))
            rst7.coordinates = coordQM
            fname=str(formula)+'_'+str(index)+"_vac.rst7"
            rst7.write(fname)

        #SOLVATAR
            tensor_coordQM = torch.Tensor(coordQM)
            solvated = solvate(tensor_coordQM, wat)

        #GENERO EL ARCHIVO .rst7
            rst7 = Rst7(natom = len(solvated))
            rst7.coordinates = solvated.numpy()
            fname=str(formula)+'_'+str(index)+'_sol.rst7'
            rst7.write(fname)

        #
        for i, carga in enumerate(cargas_ordenadas):
            topfile.atoms[i].charge = carga
        topfile.save(str(formula)+'_'+str(index)+'_ord.prmtop', format='amber')
            
            
            

        




