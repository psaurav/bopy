import numpy as np

def read_fem(filename):
    """
    """
    f = open(filename, 'r')
    s = f.readline()
    s = s.replace('[', '')
    s = s.replace(']', '')
    return np.array([float(v) for v in s.split()])

def write_fem(filename, d):
    """
    """
    f = open(filename, 'w')
    print >>f, '[',
    f.softspace=False
    print >>f, ' '.join(str(i) for i in d),
    f.softspace=False
    print >>f, ']',
    f.close()
