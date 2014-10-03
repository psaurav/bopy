import os
import sys
import re
import numpy as np

class struct():
    def __init__(self, *args, **kwargs):
        for k, v in kwargs.items():
            setattr(self, k, v)

class ParamList:

    def __init__(self, filename=None):
        #self._fp  = re.compile('^\s+')
        #self._bp  = re.compile('\s+$')
        self._pl = {}
        if not filename is None: 
            self.read_from_file(filename)

    def set_val(self, key, val):
        #key = str(key)
        #val = str(val)
        #key = key.strip()
        #val = val.strip()
        #key = self._fp.sub('', key)
        #key = self._bp.sub('', key)
        #val = self._fp.sub('', val)
        #val = self._bp.sub('', val)
        self._pl[key] = val

    def get_val(self, key):
        if key in self._pl:
            return self._pl[key]
        else:
            None

    def write_to_file(self, filename):
        f = open(filename, 'w')
        for key, val in sorted(self._pl.items()):
            line = '%s=%s'%(key, val)
            f.write(line + '\n')
        f.close()

    def read_from_file(self, filename):
        if not os.path.exists(filename):
            sys.exit("ERROR: "+__file__+" File does not exist: "+filename)
        f = file(filename, mode='rt')
        for line in f:
            if line.find('=') > -1:
                key, val = line.split('=')
                self.set_val(key.strip(), val.strip()) 
        f.close()

class QM:
    def __init__(self, filename):
        self.filename = filename
        self.read_file()
        self.print_stats()

    def print_stats(self):
        print >>sys.stderr, "QM filename    :", self.filename
        print >>sys.stderr, "QM nSources    :", self.nsrc
        print >>sys.stderr, "QM nDetectors  :", self.ndet
        #print >>sys.stderr, "QM nLinks      :", self.nll
        #print >>sys.stderr, "QM begin_index :", self.begin_index

    def read_file(self):
        f = open(self.filename, 'r')
        #
        # read the first 3 lines and throw them away
        #
        f.readline()         # QM file 3D
        self.dim = float(f.readline().split()[1])         # Dimension 3
        f.readline()         # empty line
        #
        # read the sources
        words = [word for word in f.readline().split()] # SourceList header
        self.nsrc = int(words[1])
        self.src = []
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
        self.det = []
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
            words = [val for val in line.split()]
            words.pop(0)
            words = np.array([int(val) for val in words])
            self.ll.append(words)
        f.close()

        print >>sys.stderr, "Diagnostics..."
        print >>sys.stderr, self.src.shape
        print >>sys.stderr, self.det.shape

    def get_detbbox(self):
        bbox = []
        for i in range(0, self.dim):
            b = np.array([np.min(self.det[:,i]), np.max(self.det[:,i])])
            try:
                bbox = np.vstack((bbox, b))
            except:
                bbox = b
        return bbox

    def get_detpos(self, sindex):
        pos = []
        for i in self.ll[sindex-1]: # the sindex is in gen3 notation
            try:
                pos = np.vstack((pos, self.det[i,:]))
            except:
                pos = self.det[i,:]
        return pos

    def get_srcpos(self, sindex):
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
