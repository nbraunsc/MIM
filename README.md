Molecules-in-Molecules (MIM)
==============================
[//]: # (Badges)
[![codecov](https://codecov.io/gh/REPLACE_WITH_OWNER_ACCOUNT/MIM/branch/master/graph/badge.svg)](https://codecov.io/gh/REPLACE_WITH_OWNER_ACCOUNT/MIM/branch/master)
[![docs](https://readthedocs.org/projects/pip/badge/)](https://readthedocs.org/projects/pip/badge/?version=latest&style=plastic)

[Documentation](https://nbraunsc.github.io/MIM/)

Molecular Fragmentation Code

### Installation

1. Download

    ```
    git clone https://github.com/nbraunsc/MIM.git
    cd MIM/
    ```

2. Create a conda environment and install required modules

    ```
    conda create -n mim_env pip python=3.7
    conda activate mim_env
    cd docs/
    pip install -r requirements.txt
    ```

3. Install MIM package

    ```
    cd ../
    pip install -e .
    ```


### Copyright

Copyright (c) 2021, Nicole Braunscheidel


#### Acknowledgements
 
Project based on the 
[Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms) version 1.5.
