import os
import sys

import unittest

import ase
import ase.build
import ase.io
import ase.optimize
from ase.io import write

import atomistica

def molecule(mol):  #load molecule
    if os.path.exists('{0}.xyz'.format(mol)):
        a = ase.io.read('{0}.xyz'.format(mol))
    else:
        a = ase.build.molecule(mol)
    return a

fileName = sys.argv[1]
fileOut = fileName + '-done.xyz'
a = molecule(fileName)
write(fileOut, a) # escreve antes e depois sobrescreve se der certo

try:
    a.center(vacuum=5.0)
    a.set_calculator(atomistica.Rebo2())
    a.rattle(0.05)
    ase.optimize.QuasiNewton(a, logfile='QuasiNewton.log') \
        .run(fmax=0.001)
    #ase.optimize.FIRE(a, logfile='FIRE.log').run(fmax=0.001)
    e = a.get_potential_energy()
    write(fileOut, a)
    file = open(fileOut,'a+')
    file.write(str(e))
    file.close()
except:
    pass




