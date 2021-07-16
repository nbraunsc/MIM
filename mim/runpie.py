import numpy as np

def connected(fragment, prim_dict):
    out = []
    for prim in fragment:
        out.extend(prim_dict[prim])
    final = list(set(out))
    print(final)
    return final

def pie(f_curr, old_coeff, prim_dict, derv_list, signlist, fraglist, level_dict, level):
    """
    Parameters
    ----------
    index : int
    f_curr : current fragment (or derivative)
    """
    #Make sure g is a sorted list
    g = connected(f_curr, prim_dict)
    small = list(set(g)-set(level_dict[level]))
    print("connected:", g)
    print("previous checked:", level_dict[level])
    print("smaller list:", small)
    
    print("\ndecreased list by:", len(g)-len(small), "\n")
    print("\n ############## \n current frag:", f_curr)
    print("connected:", g)
    for fk in small:
        #print("\nchecking fragment:", fk)
        #print("Level:", level)
        #print("Original coeff:", old_coeff)
        #print("level dict entry:", level_dict[level])
        if fk <= max(level_dict[level]):
        #    print("CONTINUE \n")
            continue #or break?
        f_new = f_curr.intersection(fraglist[fk])
        if len(f_new) == 0:
            level_new = level-1
            print("end of recursion\n")
            return
        else:
        #    print("new intersection:", f_new)
            derv_list.append(f_new)
            coeff = old_coeff*-1
        #    print("coeff =", coeff)
            signlist.append(coeff)
            new_entry = level_dict[level]
            level_new = level+1
        #    print("new level", level_new)
            new_list = new_entry + list([fk])
            level_dict[level_new] = new_list
        #    print("dictionary:", level_dict)
            pie(f_new, coeff, prim_dict, derv_list, signlist, fraglist, level_dict, level_new)


def start_pie(fraglist, prim_dict):
    derv_list = []
    signlist = []
    for fi in range(0, len(fraglist)):
        derv_list.append(fraglist[fi])
        signlist.append(1)
        level_dict = {}
        level_dict[0] = [fi] 
        pie(fraglist[fi], 1, prim_dict, derv_list, signlist, fraglist, level_dict, 0)
    return derv_list, signlist


