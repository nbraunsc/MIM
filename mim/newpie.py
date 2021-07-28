import numpy as np

# Using a Python dictonary to build an adjacency list
def connected(fragment, prim_dict, visited):
    out = []
    for prim in fragment:
        out.extend(prim_dict[prim])
    final = list(set(out))
    final.sort()
    return final

def dfs(visited, graph, f_curr, index, oldcoeff, signlist, dervlist, fraglist, static_frag, count, range_list):
    #Find list of attached fragments according to the prims in current frag/derivative
    adj_list = connected(f_curr, graph, visited)
    
    #when on last fragment, should not have anymore intersections so end recursive function
    if len(visited) == len(fraglist):
        return
     
    for neighbour in adj_list:
        # When back to top layer reset visited list
        #print("count:", count)
        if f_curr == static_frag:
            print("Back to top")
            count = count+1
            if count == 0:
                visited = set(range_list)
            else:
                #visited now everything higher than already checked within this node
                visited = set(adj_list[:count+1])
                print(visited)

        print("\n Current node:", index, "\nCurrent frag:", f_curr)
        #print("adjacency list:", adj_list)
        print("neighbour:", neighbour)
        print("visted:", visited)
        
        if neighbour <= max(visited):
            print("neighbour < max visited")
            count = count-1
            continue
        else:
            f_new = sorted(set(f_curr).intersection(fraglist[neighbour]))
            visited.add(neighbour)
            if len(f_new) == 0:
                return
            else:
                print("intersection found:", f_new)
                coeff = oldcoeff*-1
                print("coeff =", coeff)
                signlist.append(coeff)
                dervlist.append(f_new)
                print("\n Starting new recursive func")
                dfs(visited, graph, f_new, index, coeff, signlist, dervlist, fraglist, static_frag, count, range_list)


# Driver code
def start_pie(fraglist, graph):
    dervlist = []
    signlist = []
    for fi in range(0, len(fraglist)):
        signlist.append(1)
        dervlist.append(sorted(fraglist[fi]))
        visited = set()
        if fi == 0:
            range_list = [0]
        else:
            range_list = list(range(0, fi+1))
        #print("\nRange list:", range_list)
        for value in range_list:
            visited.add(value)

        dfs(visited, graph, sorted(fraglist[fi]), fi, 1, signlist, dervlist, fraglist, sorted(fraglist[fi]), -1, range_list)
    return dervlist, signlist

