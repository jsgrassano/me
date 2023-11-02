#!/bin/sh 
export MOLCAS=/home/santi/programas/molcas/build
export PATH=$PATH:$MOLCAS
pymolcas /home/santi/0-Pol/molecules/C3H7N1O3/C3H7N1O3_0.input -f
