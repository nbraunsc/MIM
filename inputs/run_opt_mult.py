import pickle
import dill
import numpy as np
import os
import sys

directory = sys.argv[1] #should be the directory of mim level
batch = sys.argv[2].split("_")
folder = sys.argv[3]
batch_list = []
for i in batch:
    if not i:
        continue
    else:
        entry = "fragment" + i + ".dill"
        batch_list.append(entry)

#print("Batch list:", batch_list)
os.chdir(folder)
os.chdir(directory)
for i in batch_list:
    #undill and run e, g, hess, apt etc
    infile = open(i, 'rb')
    new_class = dill.load(infile)
    outfile = open(i, "wb")
    e, g, h = new_class.qc_backend()
    new_class.hess_apt(h)  #this should only be done at optimized geometry
    
    #try:
    #    e, g, h = new_class.qc_backend()
    #    new_class.hess_apt(h)  #this should only be done at optimized geometry
    #
    #except:
    #    print("Job died:", i)
    #    name = i.replace("fragment", "").replace(".dill", "_")
    #    filename = i.replace(".dill", ".py")
    #    submit_line ='sbatch %s -J %s -o "%s" --export=LEVEL="%s",BATCH="%s",FOLDER="%s"  slurm_hess_apt.sh'%(name+"redo", directory+"/"+name+"redo.log", directory, name, folder)     ##For TinkerCliffs/Huckleberry
    #    lines = """import os
#impo#rt sys
#impo#rt dill
#
#old_#frag = open('{name}', 'rb')
#frag#ment_obj = dill.load(old_frag)
#
##mak#e changes to fragment_obj
## e.#g. fragment_obj.qc_class.spin = 0
#
#new_#frag = open('{name}', 'wb')
#dill#.dump(fragment_obj, new_frag)
#old_#frag.close()
#new_#frag.close()
#
#cmd #= '{submit_line}'
#os.c#hdir('../../')
#os.s#ystem(cmd)
#"""
    #    context = {
    #    "name":i,
    #    "submit_line":submit_line
    #    }

    #    with open(filename, "w") as myfile:
    #        myfile.write(lines.format(**context))
    #    
    #    myfile.close()
    #
    #finally:
    #    ##redill with updated fragment e, g, hess, apt, etc
    #    dill.dump(new_class, outfile)
    #    infile.close()
    #    outfile.close()
    

    print("##################################this is from run_opt", e, g.shape)
    infile.close()
    
    ##redill with updated fragment e, g, hess, apt, etc
    dill.dump(new_class, outfile)
    outfile.close()
