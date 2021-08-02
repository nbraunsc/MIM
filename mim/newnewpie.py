import numpy as np
import time

# Using a Python dictonary to build an adjacency list
def connected(fragment, prim_dict, max_ind):
    #out = []
    #for prim in fragment:
    #    out.extend(prim_dict[prim])
    #final = list(set(out))
    #final.sort()
    
    #new = [x for x in final if x > max_ind]
    #return new 
    
    #return final
    out = {}
    for prim in fragment:
        for f2 in prim_dict[prim]:
            if f2 > max_ind:
                out[f2]=1
            else:
                continue
    keylist = list(out.keys())
    #keylist.sort()
    return keylist
    

def recurse(f_curr, indices, old_coeff, fraglist, dervlist, signlist, prim_dict, starttime):
    """
    f_curr : current fragment with list of primitives within fragment
    indices : indices of fragments from which f_curr was made from
    old_coeff : previous coeff
    fraglist : original frag list
    dervlist : new derivative list
    signlist : coeff list
    """
    max_ind = max(indices)
    adj_list = connected(f_curr, prim_dict, max_ind)
    
    if not adj_list:
        return

    #print("Adj list:", adj_list)
    for fk in adj_list:
        #max_ind = max(indices)
        #if max_ind >= fk: 
        #    continue
        #f_new = sorted(set(f_curr).intersection(fraglist[fk]))
        f_new = f_curr.intersection(fraglist[fk])
        #if len(f_new) == 0:
        #    raise ValueError("Len of intersection = 0")
        #    exit()
        #    return
        #else:
       #print("\nNew intersection found:", f_new)
       new_coeff = old_coeff*-1
       end = time.time()
       if end-starttime > 30:
           exit()
       #print("New coeff:", new_coeff)
       indexfk = indices + [fk]
       #print("index list:", indexfk)
       dervlist.append(f_new)
       signlist.append(new_coeff)
       recurse(f_new, indexfk, new_coeff, fraglist, dervlist, signlist, prim_dict, starttime)


def start_pie(fraglist, prim_dict):
    dervlist = []
    signlist = []
    start = time.time()
    for fi in range(0, len(fraglist)):
        signlist.append(1)
        #dervlist.append(sorted(fraglist[fi]))
        dervlist.append(fraglist[fi])
        location = [fi]
        recurse(fraglist[fi], location, 1, fraglist, dervlist, signlist, prim_dict, start)
    end = time.time()
    total_time = end-start
    return dervlist, signlist, total_time
