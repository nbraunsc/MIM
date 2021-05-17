import pickle
import dill
import numpy as np
import os
import sys

directory = sys.argv[1] #should be the directory of mim level
batch = sys.argv[2].split("_")
print("\n" + directory + "_" + sys.argv[2])
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
status = 1
for i in batch_list:
    #undill and run e, g
    infile = open(i, 'rb')
    new_class = dill.load(infile)
    outfile = open(i, "wb")
    try:
        new_class.qc_backend()
    
    except:
        print("Job died:", i)
        status = -1
        name = i.replace("fragment", "").replace(".dill", "_")
        filename = i.replace(".dill", ".py")
        #submit_line = 'python run.py %s %s %s'%(directory, name, folder) 

        #submit_line ='sbatch -e %s -J %s -o "%s" --export=LEVEL="%s",BATCH="%s",FOLDER="%s"  slurm_pbs.sh'%(directory+"/"+name+"redo.error", name+"redo", directory+"/"+name+"redo.log", directory, name, folder)     ##For TinkerCliffs/Huckleberry
        submit_line ='sbatch %s -J %s -o "%s" --export=LEVEL="%s",BATCH="%s",FOLDER="%s"  slurm_pbs.sh'%(name+"redo", directory+"/"+name+"redo.log", directory, name, folder)     ##For TinkerCliffs/Huckleberry
        lines = """import os
import sys
import dill

old_frag = open('{name}', 'rb')
fragment_obj = dill.load(old_frag)

#make changes to fragment_obj
# e.g. fragment_obj.qc_class.spin = 0

new_frag = open('{name}', 'wb')
dill.dump(fragment_obj, new_frag)
old_frag.close()
new_frag.close()

cmd = '{submit_line}'
os.chdir('../../')
os.system(cmd)
"""
        context = {
        "name":i,
        "submit_line":submit_line
        }

        with open(filename, "w") as myfile:
            myfile.write(lines.format(**context))
        
        myfile.close()
    
    finally:
        ##redill with updated fragment e, g, hess, apt, etc
        dill.dump(new_class, outfile)
        infile.close()
        outfile.close()
        
        #update status of calculation
        status_name = i.replace(".dill", ".status")
        out_stat = open(status_name, "wb")
        dill.dump(status, out_stat)
        out_stat.close()


