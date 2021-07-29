import numpy as np

# Using a Python dictonary to build an adjacency list
def connected(fragment, prim_dict):
    out = []
    for prim in fragment:
        out.extend(prim_dict[prim])
    final = list(set(out))
    final.sort()
    return final


def recurse(f_curr, indexes, old_coeff, fraglist, dervlist, signlist, prim_dict):
    """
    f_curr : current fragment with list of primitives within fragment
    indexes : indexes of fragments from which f_curr was made from
    old_coeff : previous coeff
    fraglist : original frag list
    dervlist : new derivative list
    signlist : coeff list
    """
    
    adj_list = connected(f_curr, prim_dict)
    print("Adj list:", adj_list)
    for fk in adj_list:
        if max(indexes) >= fk:
            continue
        f_new = sorted(set(f_curr).intersection(fraglist[fk]))
        if len(f_new) == 0:
            return
        else:
            print("\nNew intersection found:", f_new)
            new_coeff = old_coeff*-1
            print("New coeff:", new_coeff)
            indexfk = indexes + [fk]
            print("index list:", indexfk)
            dervlist.append(f_new)
            signlist.append(new_coeff)
            recurse(f_new, indexfk, new_coeff, fraglist, dervlist, signlist, prim_dict)


def start_pie(fraglist, prim_dict):
    dervlist = []
    signlist = []
    for fi in range(0, len(fraglist)):
        signlist.append(1)
        dervlist.append(sorted(fraglist[fi]))
        location = [fi]
        recurse(fraglist[fi], location, 1, fraglist, dervlist, signlist, prim_dict)
    return dervlist, signlist
