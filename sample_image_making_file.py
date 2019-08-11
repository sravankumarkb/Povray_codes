#!/usr/bin/env python

from ase.io import read, write
from ase.atoms import Atoms
import os, glob, subprocess
from ase.io import pov
import numpy as np

contcars = ['CONTCAR_IS','CONTCAR_TS','CONTCAR_FS']
x = -8
y = 25
box = 12
res = 80     # Resolution
top='-90z'
side='90x,0y,-180z'

#rad = {}
#rad['C'] = 0.6
#rad['Pd'] = 1.1
#rad['Zr'] = 1.2
#rad['H'] = 0.3
#rad['O'] = 0.6

for cont in contcars:
  name= cont+'_OOH_formation_from_bridge_OH_side'
  atoms = read(cont,format='vasp')
  atoms = atoms*(1,1,1)
  colors=[None]*len(atoms)


  O2=[75,76]
  cus_Oxygen=[73,74]
  cus_Hydrogen=[0,1]
  bridge_Oxygen=[72]
  bridge_Hydrogen=[2]
  H2=[3,4]
  donot_show=[1,3,4,74,0,73]

  rad={}
  rad['Ti']=0.2
  rad['O']=0.1
  rad['H']=0.28*0.7
  rad['Au']=0.2
  radii=np.zeros(len(atoms))
  for i in range(0,int(len(atoms))):
	if (atoms[i].symbol=='H'):
		colors[i]=(1,1,1)       # Hydrogen
		radii[i]=rad['H']
	if (atoms[i].symbol=='O'):
		colors[i]=(1,0,0)       #Oxygen
		radii[i]=rad['O']
	if (atoms[i].symbol=='Au'):
		colors[i]=(1.00, 0.82, 0.14)       #Gold
		radii[i]=rad['Au']
	if (atoms[i].symbol=='Ti'):
		colors[i]=(0.75, 0.76, 0.78)       #Titanium
		radii[i]=rad['Ti']

  for i in range(0,int(len(atoms))):
	if i in cus_Oxygen:
		colors[i]=(0,0,1,0,0.7)	# blue
		radii[i]=0.59*0.7
	if i in cus_Hydrogen:
		colors[i]=(0,0,1,0,0.7)	# blue
		radii[i]=0.28*0.7
	if i in bridge_Oxygen:
		colors[i]=(1,0.5,0,0,0.7)     #orange
		radii[i]=0.59*0.7
	if i in bridge_Hydrogen:
		colors[i]=(1,0.5,0,0,0.7)     #orange
		radii[i]=0.28*0.7
	if i in H2:
		colors[i]=(0,1,1)	#cyan
		radii[i]=0.28*0.7
	if i in O2:
		colors[i]=(1,0,0)	#red
		radii[i]=0.59*0.7
	if i in donot_show:
		colors[i]=(1,1,1,0,1)	#transperent

  bond_pairs=pov.get_bondpairs(atoms,radius=0.9)

  bonds_to_be_removed=[]
  for i in range(0,len(bond_pairs)-1):
        if (bond_pairs[i][0]== 73 and bond_pairs[i][1]==95):bonds_to_be_removed.append(i)
        if (bond_pairs[i][0]== 74 and bond_pairs[i][1]==98):bonds_to_be_removed.append(i)

  bond_pairs_updated=[]

  for i in range(0,len(bond_pairs)-1):
        if(i not in bonds_to_be_removed): bond_pairs_updated.append(bond_pairs[i])
  bond_pairs=bond_pairs_updated

  bond_pairs.append((11, 62, [0, 0, 0]))
  bond_pairs.append((14, 52, [0, 0, 0]))
  bond_pairs.append((17, 42, [0, 0, 0]))

  #h_h = atoms.get_distance(9,10)
  #radii=0.7
  write(name+'.pov',atoms*(1,1,1),format='pov',rotation=side, radii=radii, bondatoms=bond_pairs, bbox=(x-box/2.,y-box/2.,x+box/2.,y+box/2.),run_povray=False,canvas_width=box*res,colors=colors,transparent=True,display=False,pause=False)
  fo=open(name+'.pov',"r")
  line=fo.readlines()
  line.insert(10,'shadowless')
  fo.close()
  f = open(name+'.pov', "w")
  line = "".join(line)
  f.write(line)
  f.close()
  subprocess.call(["povray",name+".ini"])
