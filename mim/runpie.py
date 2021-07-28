import numpy as np
from copy import deepcopy

def connected(fragment, prim_dict):
    out = []
    for prim in fragment:
        out.extend(prim_dict[prim])
    final = list(set(out))
    #print(final)
    return final

def pie(f_curr, old_coeff, prim_dict, derv_list, signlist, fraglist, level_dict, level, seen):
    """
    Parameters
    ----------
    index : int
    f_curr : current fragment (or derivative)
    """
    #Make sure g is a sorted list
    g = connected(f_curr, prim_dict)

    small = list(set(g)-set(level_dict[level]))
    small.sort()
    
    #print("\ndecreased list by:", len(g)-len(small), "\n")
    print("\n ############## \n current frag:", f_curr)
    #print("connected:", g)
    print("smaller list:", small)
    checked = seen
    for fk in small:
        print("frag index:", fk)
        print("\nchecking fragment:", fraglist[fk], "with", f_curr)
        print("Level:", level)
        print("Original coeff:", old_coeff)
        print("level dict entry:", level_dict[level])
        clone = deepcopy(checked)
        print("checked", checked, type(checked))
        print("clone:", clone, type(clone))
        if fk <= max(clone):
            print("CONTINUE \n")
            continue #or break?
        
        checked.append(fk)
        f_new = f_curr.intersection(fraglist[fk])
        if len(f_new) == 0 and not f_new:
            level = level-1
            print("end of recursion\n")
            return "end of recursion"
            #continue
        else:
            print("new intersection:", f_new)
            derv_list.append(f_new)
            coeff = old_coeff*-1
            print("coeff =", coeff)
            signlist.append(coeff)
            new_entry = level_dict[level]
            level_new = level+1
            print("new level", level_new)
            new_list = new_entry + list([fk])
            level_dict[level_new] = new_list
            print("dictionary:", level_dict)
            pie(f_new, coeff, prim_dict, derv_list, signlist, fraglist, level_dict, level_new, checked)


def start_pie(fraglist, prim_dict):
    derv_list = []
    signlist = []
    for fi in range(0, len(fraglist)):
        derv_list.append(fraglist[fi])
        signlist.append(1)
        level_dict = {}
        level_dict[0] = [fi] 
        seen = [0]
        pie(fraglist[fi], 1, prim_dict, derv_list, signlist, fraglist, level_dict, 0, seen)
    return derv_list, signlist


