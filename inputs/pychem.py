from pyqchem import Structure
from pyqchem import QchemInput
from pyqchem import get_output_from_qchem
from pyqchem.parsers.basic import basic_parser_qchem
import numpy as np

#Defining molecule
molecule = Structure(coordinates=[[0.0, 0.0, 0.0],
                                  [0.0, 0.0, 0.9],
                                  [0.0, 0.0, -0.9]],
                     symbols=['O', 'H', 'H'],
                     charge=0,
                     multiplicity=1)

#Preparing QCHEM input
qc_input1 = QchemInput(molecule,
                      jobtype='force',
                      exchange='hf',
                      basis='6-31G',
                      extra_rem_keywords={'multipole_field': 'x 0.001'},
                      extra_sections=['$multipole_field \n x 0.001 \n $end'])
                        

string = '$multipole_field \n   x 0.001 \n$end'
print(string)
exit()
print(type(qc_input1))
inp = qc_input1.get_txt()
qc_input_test = 
print(inp)
exit()

qc_input2 = QchemInput(molecule,
                      jobtype='freq',
                      exchange='hf',
                      basis='6-31G')

def get_gradient(contents_splitlines, natoms):
    grad = []
    for i in range(0, len(contents_splitlines)):
        if 'Gradient of SCF Energy' in contents_splitlines[i]:
            grad.append(i)
        if 'Max gradient component' in contents_splitlines[i]:
            grad.append(i)
            break
    gradient = contents_splitlines[grad[0]:grad[1]]
    
    N = natoms
    gradient_scf = np.zeros(shape=(N, 3))
    ncols_block = 3
    
    contents_iter = iter(gradient)
    for line in contents_iter:
        if 'Gradient of SCF Energy' in line:
            # skip the line containing the header
            line = next(contents_iter)
            
            # keep track of how many "columns" there are left; better than searching
            # for the gradient line
            ncols_remaining = N
            while ncols_remaining > 0:
                # this must be the column header
                if len(line.split()) == min(ncols_block, ncols_remaining):
                    ncols_remaining -= ncols_block
                    sline = line.split()
                    colmin, colmax = int(sline[0]) - 1, int(sline[-1])
                    line = next(contents_iter)
                # iterate over the rows
                for row in range(N):
                    sline = line.split()
                    rowidx = int(sline[0]) - 1
                    gradient_scf[rowidx, colmin:colmax] = list(map(float, sline[1:]))
                    line = next(contents_iter, 'break')
    return gradient_scf


def get_hessian(contents_splitlines2, natoms):
    hess = []
    for i in range(0, len(contents_splitlines2)):
        if 'Hessian of the SCF Energy' in contents_splitlines2[i]:
            hess.append(i)
        if '**********************************************************************' in contents_splitlines2[i]:
            hess.append(i)
            break
    hessian = contents_splitlines2[hess[0]:hess[1]]
    
    N = natoms
    hessian_scf = np.zeros(shape=(3*N, 3*N))
    ncols_block = 6
    
    contents_iter = iter(hessian)
    for line in contents_iter:
        if 'Hessian of the SCF Energy' in line:
            # skip the line containing the header
            line = next(contents_iter)
            # keep track of how many "columns" there are left; better than searching
            # for the gradient line
            ncols_remaining = 3*N
            while ncols_remaining > 0:
                # this must be the column header
                if len(line.split()) == min(ncols_block, ncols_remaining):
                    ncols_remaining -= ncols_block
                    sline = line.split()
                    colmin, colmax = int(sline[0]) - 1, int(sline[-1])
                    line = next(contents_iter)
                # iterate over the rows
                for row in range(3*N):
                    sline = line.split()
                    rowidx = int(sline[0]) - 1
                    hessian_scf[rowidx, colmin:colmax] = list(map(float, sline[1:]))
                    line = next(contents_iter, 'break')
    return hessian_scf


#Running calcs
#Getting gradient
output_grad = get_output_from_qchem(qc_input1, processors=4)
contents_splitlines = output_grad.splitlines()
gradient = get_gradient(contents_splitlines, 3)
print(gradient)

#Getting hessian
output_hess = get_output_from_qchem(qc_input2, processors=4)
contents_splitlines2 = output_hess.splitlines()
hessian = get_hessian(contents_splitlines2, 3)
print(hessian)

#Getting Energy and Dipole moment
parsed_data = basic_parser_qchem(output_grad)
print(parsed_data)
e = parsed_data['scf_energy']
print(e)
multipole = parsed_data['multipole']
dipole = multipole['dipole_moment']
print(dipole)




