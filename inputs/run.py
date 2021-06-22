import numpy as np
import pickle
import dill
import os
import sys

directory = sys.argv[1] #should be the directory of mim level
single_frag = sys.argv[2]
folder = sys.argv[3]
tmpdir = sys.argv[4]
frag_name = "fragment" + single_frag + ".dill"
print("fragment name:", frag_name)

os.chdir(folder)
os.chdir(directory)
status = 1

#undill and run e, g
infile = open(frag_name, 'rb')
new_class = dill.load(infile)
try:
    new_class.qc_backend()

except:
    print("Job died:", frag_name)
    status = -1
    name = "redo_frag" + single_frag + ".py"
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
    "name":frag_name,
    "submit_line":submit_line
    }

    with open(filename, "w") as myfile:
        myfile.write(lines.format(**context))
    
    myfile.close()

finally:
    ##redill with updated fragment e, g, hess, apt, etc
    infile.close()
    outfile = open(frag_name, "wb")
    dill.dump(new_class, outfile)
    outfile.close()
    
    #update status of calculation
    status_name = frag_name.replace(".dill", ".status")
    out_stat = open(status_name, "wb")
    dill.dump(status, out_stat)
    out_stat.close()

    #put status file in $TMPDIR from Slurm
    print("Temperary directoy:", tmpdir)
    os.chdir(tmpdir)
    out_stat = open(status_name, "w")
    out_stat.write("status is done")
    out_stat.close()
    os.chdir('../')


