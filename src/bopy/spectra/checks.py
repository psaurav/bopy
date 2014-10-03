import os
import numpy as np
import bopy as bp


def read_namedrowstable(filename, skip=0):
    f = open(filename, 'r')
    rname = []
    c = 0 
    for l in f:
        c = c + 1
        if c <= skip:
            continue
        a = l.split()
        n = a.pop(0)
        val = np.array([float(v) for v in a])
        rname.append(n)
        try:
            data = np.vstack((data, val))
        except NameError:
            data = val
    return (rname, data)

def get_namedrow((names, data), name, cindex=None):
    if cindex is None:
        return data[names.index(name)]
    else:
        return data[names.index(name), cindex]

def get_hb(data):
    return get_namedow(data, 'Hb')

def get_hbo2(data):
    return get_namedow(data, 'HbO2')

def show_names((names, data)):
    print names

def show_data((names, data)):
    print data

def read_hbleff2008():
    """
    Read the values given in Leff 2008
    """
    fp = os.path.join(bp.spectra._DATA_DIR, 'hbleff2008.txt')
    f = open(fp, 'r')
    study = []
    for l in f:
        a = l.split()
        n = a.pop(0)
        val = np.array([float(v) for v in a])
        study.append(n)
        try:
            data = np.vstack((data, val))
        except NameError:
            data = val
    return (study, data)

def ssfdpm_data(filename):
    rname, data = read_namedrowstable(filename, skip=2) 
    print rname
    print data
    return (rname, data)


def show_thc_leff2008(hbo2=None, hb=None):
    study, data = read_hbleff2008()
    y_pos = np.arange(len(study))
    y_pos = y_pos[::-1]

    import matplotlib.pyplot as plt
    plt.barh(y_pos, data[:,0], xerr=data[:,1], align='center', alpha=0.0)
    plt.yticks(y_pos, study)
    plt.xlabel('THC [uM]')
    try:
        thb = hbo2 + hb
        plt.axvline(thb, color='r')
        print "THb:", thb
    except:
        pass
    plt.show()

def show_sto2_leff2008(hbo2=None, hb=None):
    study, data = read_hbleff2008()
    y_pos = np.arange(len(study))
    y_pos = y_pos[::-1]

    import matplotlib.pyplot as plt
    plt.barh(y_pos, data[:,2], xerr=data[:,3], align='center', alpha=0.0)
    plt.yticks(y_pos, study)
    plt.xlabel('StO2 [%]')
    try:
        sto2 = 100.*hbo2/(hbo2 + hb)
        plt.axvline(sto2, color='r')
        print "StO2:", sto2
    except:
        pass
    plt.show()


