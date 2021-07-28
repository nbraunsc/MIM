def recurse(f_old, start, derivs, fraglist, signlist, sign, att_list):    #1st depth and on, finding intersection adding new
    """ 
    Recursive part of the pie
    
    Where funciton  builds layers and layers of fragments from interestions

    Parameters
    ----------
    f_old : list
        A list of primiatives for a specific fragment
    start : int
        The depth of intersections.  Always increasing.
    derivs : list
        List of derivates that keeps getting appended tol.l
    fraglist : list of lists
        List of fragment that gets generated by functions in the Fragmentation class. Does not change in this function.
    signlist : list
        List of coefficents
    sign : int
        New sign that is determined based on PIE
    """
    for f1 in att_list:
        if f1 > start:
            df_new = fraglist[f1].intersection(f_old)
            print("intersetction between", f_old, " and", fraglist[f1], " = ", df_new)


    exit()
    print("start of recurse with frag:", f_old)
    for fj in range(start+1, len(fraglist)):
        print("\nchecking", f_old, "with", fraglist[fj])
        print("checked list:", checked_list)
        
        df_new = fraglist[fj].intersection(f_old)
        
        #if fj in checked_list and len(df_new) > 0:
        #    print("fj in checked list")
        #    continue
        
        if len(df_new) > 0:
            print("interstion:", df_new)
            print("temp list:", temp_list)
            #if df_new not in temp_list:
            #    print("intersection not in temp_list")
            #    temp_list.append(df_new)
            #    level = level+1
            #    checked_list = []
            #    print("New level so checked list = 0", checked_list)
            
            print("intersection found:", df_new)
            
            #add fj to checked list
            checked_list.append(fj)

            #add new deriv to list
            derivs.append(df_new)

            #change sign and add to coeff list
            df_newcoeff = sign * -1
            signlist.append(df_newcoeff)
            print("coeff:", df_newcoeff)

            #check for additional overlaps
            recurse(df_new, fj, derivs, fraglist, signlist, df_newcoeff, att_dict)

def find_inter(frag_list, sign, att_dict, full_derv, full_sign):
    temp_derv = []
    temp_sign = []
    print("\n", frag_list)
    
    for fj in range(0, len(frag_list)):
        print("fragment:", frag_list[fj])
        att_list = []
        for prim in frag_list[fj]:
            att_list.extend(att_dict[prim])
        att_list = list(set(att_list))
        print("attached list:", att_list)
        
        ### NEED TO FIX HERE WITH FRAG_LIST INTERSEÇTION and attached list calls
        for f1 in att_list:
            if f1 > fj:
                df_new = frag_list[f1].intersection(frag_list[fj])
                if len(df_new) > 0:
                    temp_derv.append(df_new)
                    temp_sign.append(sign*-1)
    print("temp derv", temp_derv)
    print("temp sign", temp_sign)
    full_derv.extend(temp_derv)
    full_sign.extend(temp_sign)

    if len(temp_derv) != 0:
        print("Going onto next level")
        new_sign = sign*-1
        print("new sign:", new_sign)
        print(len(temp_derv))
        find_inter(temp_derv, new_sign, att_dict, full_derv, full_sign)
    
    print("Done recursing full derv:", full_derv)
    return full_derv, full_sign

def runpie(fraglist, att_dict):
    """ 
    Runs the principle of inculsion-exculsion
    
    Parameters
    ----------
    fraglist : list of lists
        List of fragments that gets generated by functions in the Fragmentation class.

    Returns
    -------
    derivs : list
        List of derivatives created by prinicple on inclusion-exculusion
    signlist : list
        List of coefficents (1 or -1)
    """
    
    derivs = []
    signlist = []
    for fi in range(0, len(fraglist)):  #0th depth, just initial frags
        dfi = fraglist[fi]
        dfi_coeff = 1
        derivs.append(dfi)
        signlist.append(dfi_coeff)

    derivs, signlist = find_inter(fraglist, 1, att_dict, derivs, signlist)
    print(derivs)
    print(signlist)
    exit()
    print(new_derv)
    print(new_sign)
    signlist.extend(new_sign_list)
    derivs.extend(new_derv)
    print(derivs)
    print(signlist)
    exit()
        
    temp_frag, temp_sign = recurse(fraglist[fi], fi, derivs, fraglist, signlist, dfi_coeff, att_list)
    derivs.extend(temp_frag)
    signlist.extend(temp_sign)
    return derivs, signlist