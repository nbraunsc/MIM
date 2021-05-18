Installation Instructions
=========================

1. Download::

    git clone https://github.com/nbraunsc/MIM.git 
    cd MIM/

2. Create a conda environment and install required modules::

    conda create -n mim_env pip python=3.7 
    conda activate mim_env 
    cd docs/ 
    pip install -r requirements.txt 

3. Install MIM package::

    cd ../ 
    pip install -e .

