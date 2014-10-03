import os
import sys
import numpy as np

def read_nim(filename, image=0):
    """
    Read a NIM file

    Args:
        filename - NIM file name
        image - the n-th image [default 0]
    Returns:
        np.array (float, dim (n))
    """

    f = open(filename, 'r')
    l = f.readline().rstrip()   # NIM
    if l != 'NIM':
        raise Exception('Not a NIM file')
    l = f.readline()                # Mesh =
    l = f.readline()                # SolutionType = 
    l = f.readline()                # ImageSize =
    n = int(l.strip().split('=')[1].strip())
    l = f.readline()                # EndHeader
    for i in range(image+1):
        try:
            l = f.readline()
            l = f.readline().rstrip()
        except:
            raise Exception('Image {0} not found'.format(image))
    return np.array([float(x) for x in l.split()])

def write_nim_header(f, meshname, imagesize, soltype='N/A'):
    """
    Write NIM header

    Args:
        f - file type
        meshname - name of the mesh file
        imagesize - size of images (number of nodes in meah)
        soltype - solution type [default: 'N/A']

    """
    print >>f, 'NIM'
    print >>f, 'Mesh = {0}'.format(meshname)
    print >>f, 'SolutionType = {0}'.format(soltype)
    print >>f, 'ImageSize = {0}'.format(imagesize)
    print >>f, 'EndHeader'

def write_nim_one(filename, vect, meshname=''):
    """
    Write a NIM file
    
    Args:
        filename - NIM file name
        vect - np.array (dim - n)
    """
    n = vect.shape
    if len(n) > 1:
        raise Exception('Cannot handle more than 1 image while writing to NIM')
    n = n[0]
    f = open(filename, 'w')
    write_nim_header(f, meshname, n) 
    print >>f, 'Image 0'
    print >>f, ' '.join(map(str, vect))
    f.close()
    
    
def write_nim(filename, vect, meshname=''):
    write_nim_one(filename, vect, meshname)
