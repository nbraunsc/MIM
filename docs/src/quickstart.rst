Quick Start
===========

Set up
------

Change into the inputs directory ::
   
    $ cd inputs/

Edit the input file with desired parameters for calculation ::

   $ vi input_file.py

You will need to obtain molecular coordinates in Chemical Markup Language (.cml) for the 
molecule of interest (See `MIM/inputs/coords/water.cml for example` .cml file).

Create a directory to store results ::

    $ mkdir to_run/

Run 
---

At this point you have everything required to run your MIM calculation! Run the
following command ::
    
    $ python sow.py input_file.py <path_to_coordinates> to_run/

Collect Results
---------------

Once all of your calculations have finished, you may collect all your results
by running the following command ::
    
    $ python reap.py to_run/
