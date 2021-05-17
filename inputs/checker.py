import dill
import os
import glob
import sys

folder = sys.argv[1]
os.path.abspath(os.curdir)
os.chdir(folder)
status_list = []
status = 0

#def opt_fnc():
for i in os.listdir():
    os.chdir(i)
    stat_files = glob.glob('*.status')
    frag_files = glob.glob('*.dill')
    if len(stat_files) != len(frag_files):
        print(status)
        exit()

    else: 
        for j in stat_files:
            infile = open(j, 'rb')
            var = dill.load(infile)
            if var == -1:
                infile.close()
                print(status)
                exit()
            status_list.append(var)
            infile.close()
        os.chdir('../')

#if len(status_list) != 0 and -1 not in status_list:
status = 1
print(str(status))
exit()


