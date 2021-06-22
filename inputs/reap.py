import pickle
import os
import glob
import dill
import pickle
import numpy as np
from mendeleev import element
import os
import sys

folder = sys.argv[1]

np.set_printoptions(precision=8, linewidth = 200)

os.chdir(folder)
#levels = os.listdir()

levels = []
for filename in os.listdir():
    if os.path.isdir(os.path.join(filename)):
        levels.append(filename)

e = 0
g = 0
h = 0
apt = 0
aptgrad = 0
step = 0.001

for level in levels:
    os.chdir(level)
    frags = glob.glob('*.dill')
    for i in frags:
        #undill and run e, g, hess, apt etc
        infile = open(i, 'rb')
        new_class_loop = dill.load(infile)
        outfile = open(i, "wb")
        e += new_class_loop.energy
        g += new_class_loop.grad
        h += new_class_loop.hessian
        apt += new_class_loop.apt 
        #aptgrad += new_class_loop.aptgrad
        dill.dump(new_class_loop, outfile)
        infile.close()
        outfile.close()
    print("############## Energy for level(", level, ") is :", e, "################")
    os.chdir('../')

###############         ONLY FOR DOING WITHOUT FINITE GRAD METHOD ####################
#aptgrad = apt

##Opening just frag1 for IR mass-weighting Hessian, pickling/repickiling issued solved with following
os.chdir('frag1')
infile = open('fragment0.dill', 'rb')
new_class = dill.load(infile)
#

#mass-weight H, compute freq and normal modes
freq, modes_unweight, M, evalues = new_class.mw_hessian(h)
pq = np.dot(apt.T, modes_unweight)   #shape 3x3N    unit: D/A/np.sqrt(amu)
pq_pq = np.dot(pq.T, pq)    #shape 3Nx3N
intense = np.diagonal(pq_pq)
intense_kmmol = intense*42.256078   #D**2/A**2/amu -> km/mol

outfile = open('fragment0.dill', 'wb')
dill.dump(new_class, outfile)
infile.close()
outfile.close()

print(os.getcwd())
os.chdir('../')
print(os.getcwd())
np.save('energy.npy', e)
np.save('gradient.npy', g)
np.save('hessian.npy', h)
np.save('apt.npy', apt)
np.save('freq.npy', freq)
np.save('intensity.npy', intense_kmmol)
np.save('normalmodes.npy', modes_unweight)

print("\nMIM Energy (Hartree)\n", e)
print("\nMIM Gradient:\n", g)
print("\nMIM Hessian shape:\n", h.shape, "\n")
print("MIM frequencies:\n", freq)
print("\nMIM intensities:\n", intense_kmmol)
exit()

outfile = open('fragment0.dill', 'wb')
dill.dump(new_class, outfile)
infile.close()
outfile.close()

shape = int(g.shape[0])
diff = np.divide(apt, aptgrad)

print("\n apt from atomic coord disp [units=D/A/np.sqrt(amu)]:\n", apt)
print("\n apt from atomic coord disp [units=D/A/np.sqrt(amu)]:\n", apt.T)
print("\n apt from grad disp [units=D/A/np.sqrt(amu)]:\n", aptgrad)
print("\n apt from grad disp [units=D/A/np.sqrt(amu)]:\n", aptgrad.T)
#print("subtracting:\n", np.subtract(apt, aptgrad))
#print("\n apt/aptgrad:\n", np.divide(apt, aptgrad))
#print("\n apt/aptgrad:\n", np.divide(apt, aptgrad).T)

print("\nAPT from atomic coords in electric dipole x direction (units Debye**2/A**2 amu):\n", apt[:,0].reshape((shape,3)))
print("\nAPT grad component for x direction (units Debye**2/A**2 amu):\n", aptgrad[:,0].reshape((shape,3)))
print("\nAPT/APTGRAD:\n", diff[:,0].reshape((shape,3)))

print("\nAPT from atomic coords in electric dipole y direction (units Debye**2/A**2 amu):\n", apt[:,1].reshape((shape,3)))
print("\nAPT grad component for y direction (units Debye**2/A**2 amu):\n", aptgrad[:,1].reshape((shape,3)))
print("\nAPT/APTGRAD:\n", diff[:,1].reshape((shape,3)))

print("\nAPT from atomic coords in electric dipole z direction (units Debye**2/A**2 amu):\n", apt[:,2].reshape((shape,3)))
print("\nAPT grad component for z direction (units Debye**2/A**2 amu):\n", aptgrad[:,2].reshape((shape,3)))
print("\nAPT/APTGRAD:\n", diff[:,2].reshape((shape,3)))

#"""Start of apt from gradient derivates w.r.t applied electric field """
pqgrad = np.dot(aptgrad.T, modes_unweight)   #shape 3x3N, unit D/A/np.sqrt(amu)
pq_pqgrad = np.dot(pqgrad.T, pqgrad)    #shape 3Nx3N
intense_kmmolgrad = np.diagonal(pq_pqgrad)
intensegrad = np.diagonal(pq_pqgrad)
intense_kmmolgrad = intensegrad*42.256078  #D**2/A**2/amu -> km/mol

# """Start of apt from dipole moment derivatives w.r.t atomic coordinates """
pq = np.dot(apt.T, modes_unweight)   #shape 3x3N    unit: D/A/np.sqrt(amu)
pq_pq = np.dot(pq.T, pq)    #shape 3Nx3N
intense = np.diagonal(pq_pq)
intense_kmmol = intense*42.256078   #D**2/A**2/amu -> km/mol

print(np.dot(pq, pq.T))
print(np.dot(pqgrad, pqgrad.T))

print("\n apt from atomic coord disp in norm coords basis [units=D/A/np.sqrt(amu)]:\n", pq.T)
print("\n apt from grad disp in norm coords basis [units=D/A/np.sqrt(amu)]:\n", pqgrad.T)
print("\n pq/pqgrad:\n", np.divide(pq.T, pqgrad.T))

print("\nintensity in in km/mol from gradient: \n", intense_kmmolgrad)
print("\nintensity in kmmol: \n", intense_kmmol)
print("\ndividing apt int/grad int:\n", np.divide(intense_kmmol, intense_kmmolgrad)[6:])

os.chdir('../../')
np.save('energy.npy', e)
np.save('gradient.npy', g)
np.save('hessian.npy', h)
np.save('apt.npy', apt)
np.save('freq.npy', freq)
np.save('intensity.npy', intense_kmmol)
print("MIM Energy:", e, "Hartree")
print("MIM Gradient:\n", g)
print("MIM Hessian shape:\n", h.shape, "\n")
print("MIM mass-weighted APT's shape:", apt.shape)



#for i in range(0, modes_unweight.shape[0]):
#    modes_trans = np.dot(M, modes_unweight).T
#    modes_trans_nm = modes_unweight.T
#    modej = modes_trans[i].reshape((g.shape[0],3))
#    modej_nm = modes_trans_nm[i].reshape((g.shape[0],3))
#    print("Mode # (mass-weighted):", i)
#    print("\n", modej, "\n")
    #print("Mode # (NOT mass-weighted):", i)
    #print("\n", modej_nm, "\n")


################### All Normal mode testing for IR intensites ##########################

#modes1 = np.dot(M, modes_unweight).T   #mass-weighted normal modes
#modes = modes1
#mole_coords = new_class.molecule.atomtable
#labels = []
#
#for i in range(0, len(mole_coords)):
#    label = mole_coords[i][0]
#    mole_coords[i].remove(label)
#    labels.append(label)
#
#atoms = np.array(labels)
#coords = np.array(mole_coords)
#coords = np.array(mole_coords)
#
#x = np.indices(modes.shape)[1]
#apt_norm = np.zeros((modes.shape[0], 3))
#step_mode = 0.001
#
#def normal_disp(modes_list):
#    for i in modes_list:
#        print("mode:", i)
#        square = modes[i].reshape((g.shape[0], 3))
#        pos_list = []
#        neg_list = []
#        no = []
#        pos = np.add(coords, step_mode*square) #getting new global coords for displacement along normal mode
#        neg = np.subtract(coords, step_mode*square) #getting new global coords for displacement along normal mode
#        for row in range(0, new_class.molecule.natoms):
#            pos_list.append(list(pos[row]))
#            neg_list.append(list(neg[row]))
#            no.append(list(coords[row]))
#            pos_list[row].insert(0, labels[row])
#            neg_list[row].insert(0, labels[row])
#            no[row].insert(0, labels[row])
#        new_class.molecule.atomtable = pos_list #setting positive disp global coords 
#        inputxyz = new_class.build_xyz()    #getting pyscf xyz input format
#        dip, nuc = new_class.qc_class.get_dipole(inputxyz)   #getting dipole with new coords
#        new_class.molecule.atomtable = neg_list  #setting negative disp  global coords 
#        inputxyz = new_class.build_xyz()    #getting pyscf xyz input format
#        dip2, nuc2 = new_class.qc_class.get_dipole(inputxyz)   #getting dipole with new coords
#        apt_comp = (dip-dip2)/(2*step_mode)    #3x1 vector
#        apt_norm[i] = apt_comp
#        new_class.molecule.atomtable = no
#
#print("normal modes:\n")
#normal_disp(list(mode for mode in range(len(freq))))   #telling which modes to compute
#print(apt_norm)
#norm_freq = np.dot(apt_norm, apt_norm.T)
#test_intense = np.diagonal(norm_freq)
#print("intensity in unknown units: \n", test_intense)
#test_intense_kmmol = test_intense*42.256078
#print("intensity in kmmol: \n", test_intense_kmmol)
#################################### End of Normal mode testing code #############################

print("\nFrequencies [cm-1] and IR intensities [km/mol] from APT w.r.t atomic coords:")
for i in range(0, len(freq)):
    print("Freq:", freq[i], "int :", intense_kmmol[i])

print("\nFrequencies [cm-1] and IR intensities [km/mol] from gradient deriv w.r.t applied field:")
for i in range(0, len(freq)):
    print("Freq:", freq[i], "int :", intense_kmmolgrad[i])

#print("\nFrequencies [cm-1] and IR intensities [km/mol] from dipole moment deriv w.r.t normal coordiant displacements:")
#for i in range(0, len(freq)):
#    print("Freq:", freq[i], "int :", test_intense_kmmol[i])
