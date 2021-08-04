import numpy as np
import time

# Using a Python dictonary to build an adjacency list
def connected(fragment, prim_dict, max_ind):
    """ Find fragments that are connected to current fragment.

    Parameters
    ----------
    fragment : list
        Current fragment
    prim_dict : dict 
        Where prims know which fragments they are in
    max_ind : int
        max value of indices used to make fragment 

    Returns
    -------
    final : numpy array
        Array of connected fragment indices
    """
    'Numpy array method'
    #arr = np.zeros((len(prim_dict)))
    #for prim in fragment:
    #    for f2 in prim_dict[prim]:
    #        arr[f2]=1

    #pull out nonzero entries
    #x = np.array(np.nonzero(arr))
    #pull out fragments that are higher than max index
    #final = x[x > max_ind]
    #return final
    #return x
    
    'List method'
    #out = []
    #for prim in fragment:
    #    out.extend(prim_dict[prim])
    #final = list(set(out))
    #final.sort()
    
    #new = [x for x in final if x > max_ind]
    #return new 
    
    'Dict method'
    out = {}
    for prim in fragment:
        for f2 in prim_dict[prim]:
    ##        out[f2]=1
            if f2 > max_ind:
                out[f2]=1
            else:
                continue
    keylist = list(out.keys())
    ##keylist.sort()
    return keylist

def recurse(f_curr, indices, old_coeff, fraglist, dervlist, signlist, prim_dict, starttime, derv_dict):
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
    
    for fk in adj_list:
        #if max_ind >= fk:
        #    continue
        #end = time.time()
        #if end-starttime > 30:
        #    exit()
        f_new = sorted(set(f_curr).intersection(set(fraglist[fk])))
        try:
            derv_dict[tuple(f_new)]+=old_coeff*-1
        except:
            derv_dict[tuple(f_new)] = old_coeff*-1
        #f_new = f_curr.intersection(fraglist[fk])
        new_coeff = old_coeff*-1
        indexfk = indices + [fk]
        #dervlist.append(f_new)
        #signlist.append(new_coeff)
        recurse(f_new, indexfk, new_coeff, fraglist, dervlist, signlist, prim_dict, starttime, derv_dict)
    
    #if not adj_list:
    #    return

    #print("Adj list:", adj_list)
        #max_ind = max(indices)
        #if max_ind >= fk: 
        #    continue
        #f_new = sorted(set(f_curr).intersection(fraglist[fk]))
        #if len(f_new) == 0:
        #    raise ValueError("Len of intersection = 0")
        #    exit()
        #    return
        #else:
        #print("\nNew intersection found:", f_new)
        #end = time.time()
        #if end-starttime > 30:
        #    exit()
        #print("New coeff:", new_coeff)
        #print("index list:", indexfk)

def start_pie(fraglist, prim_dict):
    dervlist = []
    signlist = []
    derv_dict = {}
    start = time.time()
    for fi in range(0, len(fraglist)):
        derv_dict[tuple((sorted(fraglist[fi])))] = 1
        signlist.append(1)
        dervlist.append(fraglist[fi])
        location = [fi]
        recurse(fraglist[fi], location, 1, fraglist, dervlist, signlist, prim_dict, start, derv_dict)
    end = time.time()
    total_time = end-start
    return dervlist, signlist, total_time, derv_dict
