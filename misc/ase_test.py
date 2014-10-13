from __future__ import print_function
import os, time, itertools, math, sys
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from ase import Atom, Atoms
from ase.io import *
from ase.visualize import view
from ase.calculators.emt import EMT
#from ase.constraints import FixAtoms
#from ase.optimize import QuasiNewton
#from ase.lattice.surface import fcc111, add_adsorbate
from ase.lattice.cubic import FaceCenteredCubic as fcc

#LS = 3
#d = 1.1
#flake = fcc111('Au', size=(LS,LS,LS), vacuum=10.0)
#molecule = Atoms('2N', positions=[(0.,0.,0.),(0.,0.,d)])

#class Flake(object):
    #def __init__(self):
        #self.size = (LATTICE_SIZE, LATTICE_SIZE, LATTICE_SIZE)
        #self.packing = 'fcc111'
        #self.atom = 'Au'

    #def lattice(self):
        #lattice = self.packing(self.atom, size=self.size, vacuum=10.0)

# ::: main :::
LS = 2
LATTICE_SIZE = (LS,LS,LS)
SPACE = 0
GOLD = 'Au'

#flake = fcc111(GOLD, size=LATTICE_SIZE, vacuum=SPACE)
#for i in range(10):
    #add_adsorbate(flake, GOLD, height=2, position='hcp')

#view(flake)
#print(flake.get_cell())

flake = fcc(directions=[[1,-1,0], [1,1,-2], [1,1,1]], size=LATTICE_SIZE, 
        symbol=GOLD, pbc=(1,1,0))

with open('gold_cell.txt', 'w') as f:
    f.write(flake.get_cell())
    
view(flake)
