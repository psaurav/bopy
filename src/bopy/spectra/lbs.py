import numpy as np

def read_optical_properties_fdpm(filename, prop='MUA', nwl=6):
    f = open(filename, 'r')
    insection = False
    count = 0
    for l in f:
        a = [v for v in l.split()]
        if prop == a[0] or insection:
            if not insection:
                insection = True
                continue
            count = count + 1 
            if count > nwl:
                break
            d = np.array([float(v) for v in a])
            try:
                data = np.vstack((data, d))
            except NameError:
                data = d
    return data

def read_chromophores_fdpm(filename, prop='(SSFDPM)'):
    f = open(filename,'r')
    insection = False
    datadic = {}
    for l in f:
        if prop in l or insection:          # enter only if in SSFDPM section
            if not insection:                # if in first line of SSFDPM 
                insection = True
                continue
            s = l.strip()
            if not s:                       # if empty line in SSFDPM end reached
                break
            a = [v for v in s.split()]
            if 'hromophore' in a[0]:        # the header line in section
                colnames = a[:]             # keep header names
                continue
            n = a.pop(0)                    # n is the row name
            d = np.array([float(v) for v in a])       # rest the floats
            datadic[n] = d
    colnames.pop(0)
    return (datadic, colnames)

def read_mua(filename):
    d = np.loadtxt(filename, comments='%')
    return d

def read_musp(filename):
    d = np.loadtxt(filename, skiprows=5)
    return d

