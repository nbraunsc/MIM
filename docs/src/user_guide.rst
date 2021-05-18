User Guide
==========

Energy
------
The energy is calculated for all MIM calculations.

Gradient
--------
The first derivative of the energy is computed either analytically
(if available for specific method in quantum chemistry backend) or 
numerically.

Geometry Optmizations
---------------------
You are able to run geometery optimizations within the MIM framework. 
This is set in the input file that is edited before running calculation.

Hessian
-------
The second derivative of the energy is computed either analytically
(if available for specific method in quantum chemistry backend) or 
numerically.

Atomic Polar Tensors
--------------------
The atomic polar tensors can be computed using the derivative
of the dipole moment with respect to the mass-weighted atomic
coords or using the finite field method.

Infrared Intensities
---------------------
You can compute IR intensities from the APTs or tensors computed
using the finite field method.


Quantum Chemistry Backends
--------------------------

.. autosummary::
   :toctree: autosummary
   
    mim.Pyscf
    mim.Psi4

