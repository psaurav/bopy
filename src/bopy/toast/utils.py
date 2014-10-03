import os
import sys
import re
import numpy as np

def read_mesh_nodes(filename):
    """
    Read the mesh nodes

    Args:
        filename
    Returns:
        np.array - size (n, 3)
    """
    f = open(filename, 'r')
    
    # header
    l = f.readline().rstrip().split()    # header line
    if l[0] != 'MeshData' or l[1] != '5.0':
        print >>sys.stderr, l[0], l[1]
        raise Exception('{0} does not have recognizable header'.
                         format(filename))
    l = f.readline()            # empty line

    # Nodelist 
    l = f.readline().rstrip().split()    # NodeList
    if l[0] == 'NodeList':
        nn = int(l[1])
    else:
        raise Exception('{0} could not find NodeList'.
                         format(filename))

    # read the nodes
    nodes = []
    for i in range(nn):
        l = f.readline().rstrip()
        node = [float(x) for x in l.split('[')[1].split(']')[0].split()]
        nodes.append(node)
    f.close()
    return np.array(nodes)

def read_mesh_bbox(filename=None, nodes=None):
    """
    Get the bounding box of the mesh.  
    Essentially the min and max of node positions.

    Args:
        filename - mesh file name [default: None]
        nodes - np.array (n, 3) [default: None]
    Returns:
        np.array - (2, 3) 
    """
    if filename != None:
        nodes = read_mesh_nodes(filename)
    return np.array([[np.min(nodes[:,0]), np.min(nodes[:,1]), np.min(nodes[:,2])],
                     [np.max(nodes[:,0]), np.max(nodes[:,1]), np.max(nodes[:,2])]])

def node_in_breast(nodes, bs):
    """
    Returns an bool array 
    
    Args:
        nodes - mesh nodes np.array (nn, 3)
        bs - breast and chest shape 
    """
    from scipy.interpolate import griddata
    zn = griddata(bs[:,0:2], bs[:,2], nodes[:,0:2], method='cubic')
    return self.nodes[:,2] < zn

def mshbr2inbreast(mshfilename, brfilename, zchestwall=-np.inf):
    """
    Creates a file 'inbreast.npy' which holds truthvalues for 
    whether the nodes in a mesh file are in the breast.
    
    Args:
        mshfilename - name of the mesh file
        brfilename - name of the breast shape file
        zchestwall - the z position of the chestwall (default=-np.inf)
    """
    import bopy.io as bio
    from bopy.breastshape.breast import BreastResult

    nodes = bio.read_mesh_nodes(mshfilename)
    br = BreastResult(brfilename)
    inbreast = []
    for p in nodes:
        ib = br.isInsideBreast(p[0], p[1], p[2], zchestwall=zchestwall)
        inbreast.append(ib)
    inbreast = np.array(inbreast)
    np.save('inbreast.npy', inbreast)

def inbreast2nim(ibfilename, nimfilename, inval, outval, mshfilename=''):
    """
    Create NIM file from given inbreast.npy and in-breast and out-breast values

    Args:
        ibfilename - name of inbreast.npy file
        nimfilename - name of the output NIM file
        inval - the value in the breast
        outval - the value outside the breast
        mshfilename - name of the mesh file (default: '') 
    """
    inbreast = np.load(ibfilename)
    from bopy.io import write_nim_one
    
    x = np.ones(inbreast.shape)*outval
    x[inbreast==True] = inval
    write_nim_one(nimfilename, x, meshname=mshfilename)
    
class Mesh:
    def __init__(self, filename=None):
        if filename is not None:
            self.read_file(filename)

    def read_file(self, filename):
        f = open(filename, 'r')
        l = f.readline().rstrip().split()    # header line
        if l[0] != 'MeshData' or l[1] != '5.0':
            print >>sys.stderr, l[0], l[1]
            raise Exception('{0} does not have recognizable header'.
                             format(filename))
        l = f.readline()            # empty line
        l = f.readline().rstrip().split()    # NodeList
        if l[0] == 'NodeList':
            self.nn = int(l[1])
        else:
            raise Exception('{0} could not find NodeList'.
                             format(filename))
        nodes = []
        for i in range(self.nn):
            l = f.readline().rstrip()
            node = [float(x) for x in l.split('[')[1].split(']')[0].split()]
            nodes.append(node)
        self.nodes = np.array(nodes)
        l = f.readline()            # emply line
        l = f.readline().rstrip().split()    # ElementList
        if l[0] == 'ElementList':
            self.ne = int(l[1])
        else:
            raise Exception('{0} could not find ElementList'.
                             format(filename))
        print >>sys.stderr, "Did not element list..."
    
    def get_nodes(self):
        return self.nodes

    def get_bbox(self):
        n = self.nodes
        return np.array([[np.min(n[:,0]), np.min(n[:,1]), np.min(n[:,2])],
                         [np.max(n[:,0]), np.max(n[:,1]), np.max(n[:,2])]])

    def get_bbox2(self):
        n = self.nodes
        return np.array([[np.min(n[:,0]), np.max(n[:,0])],
                         [np.min(n[:,1]), np.max(n[:,1])],
                         [np.min(n[:,2]), np.max(n[:,2])]])

    def node_in_breast(self, p):
        xy = p[:,1:2]
        from scipy.interpolate import griddata
        zn = griddata(p[:,0:2], p[:,2], self.nodes[:,0:2], method='cubic')
        return self.nodes[:,2] < zn


class QM:
    def __init__(self, filename):
        self.filename = filename
        self.read_file()
        #self.print_stats()

    def print_stats(self):
        print >>sys.stderr, "QM filename    :", self.filename
        print >>sys.stderr, "QM nSources    :", self.nsrc
        print >>sys.stderr, "QM nDetectors  :", self.ndet
        print >>sys.stderr, "QM nLinks      :", 
        for s in range(self.nsrc):
            print >>sys.stderr, len(self.ll[s]),

    def write_file(self, filename):
        f = open(filename, 'w')
        print >>f, 'QM file %dD'%(self.dim)
        print >>f, 'Dimension %d'%(self.dim)
        print >>f, ''
        print >>f, 'SourceList %d fixed'%(self.nsrc)
        for s in range(0, self.nsrc):
            for d in range(self.dim):
                print >>f, self.src[s, d],
            f.softspace=False
            print >>f, ''
        print >>f, ""
        print >>f, 'MeasurementList %d'%(self.ndet)
        for s in range(0, self.ndet):
            for d in range(self.dim):
                print >>f, self.det[s, d],
            f.softspace=False
            print >>f, ''
        print >>f, ""
        print >>f, 'LinkList'
        for s in range(0, self.nsrc):
            print >>f, '%d:'%(len(self.ll[s])),
            for d in range(len(self.ll[s])):
                print >>f, self.ll[s][d],
            f.softspace=False
            print >>f, ''
        f.close()

    def read_file(self):
        f = open(self.filename, 'r')
        #
        # read the first 3 lines and throw them away
        #
        f.readline()         # QM file 3D
        self.dim = int(f.readline().split()[1])         # Dimension 3
        f.readline()         # empty line
        #
        # read the sources
        words = [word for word in f.readline().split()] # SourceList header
        self.nsrc = int(words[1])
        self.src = None
        for s in range(0, self.nsrc):
            pos = np.array([float(val) for val in f.readline().split()])
            try:
                self.src = np.vstack((self.src, pos))
            except:
                self.src = pos

        f.readline()         # empty line
        #
        # read measurements
        words = [word for word in f.readline().split()] # MeasurementList header
        self.ndet = int(words[1])
        self.det = None
        for d in range(0, self.ndet):
            pos = np.array([float(val) for val in f.readline().split()])
            try:
                self.det = np.vstack((self.det, pos))
            except:
                self.det = pos

        f.readline()         # empty line
        f.readline()         # LinkList 

        self.ll = []
        for s in range(0, self.nsrc):
            line = f.readline()
            links = np.array([val for val in line.split(':').pop().strip().split()]).astype(int)
            self.ll.append(links)
        f.close()

    def get_detbbox(self):
        bbox = []
        for i in range(0, self.dim):
            b = np.array([np.min(self.det[:,i]), np.max(self.det[:,i])])
            try:
                bbox = np.vstack((bbox, b))
            except:
                bbox = b
        return bbox

    def get_detpos(self, sindex=None):
        if sindex is None:
            return self.det

        pos = []
        for i in self.ll[sindex-1]: # the sindex is in gen3 notation
            try:
                pos = np.vstack((pos, self.det[i,:]))
            except:
                pos = self.det[i,:]
        return pos

    def get_srcpos(self, sindex=None):
        if sindex is None:
            return self.src

        return self.src[sindex - 1,:] # the sindex is in gen3 notation
    
    def get_Yrho(self, sindex):
        s = self.get_srcpos(sindex)
        d = self.get_detpos(sindex)
        s = s[0:3:2]
        d = d[:,0:3:2]
        d = d - s
        d = d[:,0]**2 + d[:,1]**2
        d = np.sqrt(d)
        return d

