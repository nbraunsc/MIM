.. mim documentation master file, created by
   sphinx-quickstart on Thu Mar 15 13:55:56 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


Welcome to MIM's documentation!
=========================================================
.. image:: fragment4.png
    :scale: 40 %
    :align: right
    :alt: alternate text
This is the documentation for the molecular fragmentation method called Molecules-in-Molecules (MIM). MIM is a linear scaling multi-theory method that can have an arbitrary number of levels, MIMn. Where the general idea is that each lower level of theory is attempting to correct the fragmentation error associated with the higher layer. MIM calculations are largely based on the principle of inclusion/exclusion to generate highly overlapping fragments verses disjoint fragments to better predict macromolecules properties.

.. toctree::
   :maxdepth: 2

   install
   quickstart
   example
   user_guide
   api

Indices and tables
==================
* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

References
==========
- Mayhall, N. J.; Raghavachari, K., Molecules-in-Molecules: An Extrapolated Fragment-Based Approach for Accurate Calculations on Large Molecules and Materials. J Chem Theory Comput 2011, 7 (5), 1336-43
- Mayhall, N. J.; Raghavachari, K., Many-Overlapping-Body (MOB) Expansion: A Generalized Many Body Expansion for Nondisjoint Monomers in Molecular Fragmentation Calculations of Covalent Molecules. J Chem Theory Comput 2012, 8 (8), 2669-75.
- Jovan Jose, K. V.; Raghavachari, K., Molecules-in-molecules fragment-based method for the evaluation of Raman spectra of large molecules. Molecular Physics 2015, 113 (19-20), 3057-3066.

