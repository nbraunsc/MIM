Quick Start
===========

Set up
------

1. Make changes to input file: 

    It is prefered that you edit the input file already in the `inputs/`,
but you may have your own. ::

   $ cd inputs/
   $ vi input_file.py

2. Obtain molecular coordinates in Chemical Markup Language (.cml):

    You will next have to put your molecular coordinates file that is a `*.cml` file in `inputs/`.
We will use `molecule.cml` for reference.

3. Create a directory to store results:

    Finally, make a directory within `inputs/` where you would like to store all generated files.
We will use `to_run/` for reference.

Run 
---
At this point you have everything required to run your MIM calculation! ::
    
    $ python sow.py input_file.py molecule.cml to_run/

Collect Results
---------------
Once all of your calculations have finished, you may collect all your results
by running the following command::
    
    $ python reap to_run/
