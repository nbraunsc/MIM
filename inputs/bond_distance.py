import os
import numpy as np
import sys
from sys import argv
import xml.etree.ElementTree as ET

opt_coords = sys.argv[1]
cml = sys.argv[2]

with open(opt_coords) as f:
    lines = f.readlines()
lines_new = lines[2:]

tree = ET.parse(cml)
root = tree.getroot()

#parse the cml for C-C bonds
rows = []
for node in root:
    rows.append(str(node.tag))
for i in range(0, len(rows)):
    if rows[i].endswith('atomArray'):
        atomArray = root[i]
    if rows[i].endswith('bondArray'):
        bondArray = root[i]

natoms = len(atomArray)
carbons = []
for atomi in range(0, natoms):
    atom = atomArray[atomi].attrib['elementType']
    if atom == 'C':
        atom_num = atomArray[atomi].attrib['id'].replace('a', "")
        atom_num = int(atom_num)-1
        carbons.append(atom_num)
    else:
        continue
carbon_bonds = []
for bondi in bondArray:
    a12 = bondi.attrib['atomRefs2'].replace("a", "").split()
    a12_new = []
    for atom in a12:
        atomi = int(atom)-1
        a12_new.append(atomi)
    if all(item in carbons for item in a12_new):
        carbon_bonds.append(a12_new)
    else:
        continue

distances = []
for bond in carbon_bonds:
    single = []
    for atom in bond:
        x = np.fromstring(lines_new[atom].replace("C", "").replace("\n", ""),  dtype=float, sep=' ')
        single.append(x)
    distance = np.sqrt((single[1][0]-single[0][0])**2+ (single[1][1]-single[0][1])**2 + (single[1][2]-single[0][2])**2)
    distances.append(distance)

print(distances)


    
