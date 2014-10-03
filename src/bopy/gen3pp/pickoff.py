"""
This module implements the pickoff methods to correct 
the data.   
"""
import os
import sys
import numpy as np
from scipy.ndimage.filters import gaussian_filter
import h5py
from bopy.utils import ParamList, struct

def get_pickoff_params(cfgfile):
    pl = ParamList(cfgfile)
    p = struct()
    p.root = pl.get_val('ROOT')
    p.wls = [int(v) for v in pl.get_val('WAVELENGTHS').strip().split()]
    p.sigma = float(pl.get_val('SIGMA'))
    p.dkinds = [v for v in pl.get_val('DKINDS').strip().split()]
    p.pickoffpixels = np.array([int(s) for s in pl.get_val('PICKOFF_PIXELS').strip().split()])
    return p
    

def read_cfg(fname):
    pl = ParamList(fname)
    p = struct()
    p.root = pl.get_val('ROOT')
    p.ddirs = [v for v in pl.get_val('PICKOFF_DDIRS').strip().split()]
    ref = pl.get_val('REF_NAME')
    trg = pl.get_val('TRG_NAME')
    p.exs = [v for v in pl.get_val('PICKOFF_EXS').strip().split()]
    p.wls = [int(v) for v in pl.get_val('WAVELENGTHS').strip().split()]

def read_smoothed(fname, sigma):
    d = np.load(fname)
    return gaussian_filter(d, sigma=sigma)

def get(d, pixels, ckind='median'):
    if pixels.size == 2:
        r, c = pixels
        return d[r,c]
    else:
        rmin, rmax, cmin, cmax = pixels
        if ckind is 'median':
            return np.median(d[rmin:rmax, cmin:cmax].flatten())
        if ckind is 'mean':
            return np.mean(d[rmin:rmax, cmin:cmax].flatten())
    return None
     
def get_line_filenames(ddir, rootex, wl, dkind='phi'):
    fpaths = []
    for i in np.arange(209)+1:
        fpath = '{0}_wl{1}_s{2}_{3}.npy'.format(rootex, wl, i, dkind)
        fpath = os.path.join(ddir, fpath)
        fpaths.append(fpath)
    return fpaths 

def get_line(fpaths, sigma, pixels, ckind='median'):     
    pline = []
    #for i in np.arange(209)+1:
    for fname in fpaths:
        d = read_smoothed(fname, sigma)
        d =  get(d, pixels, ckind=ckind)
        pline.append(d)
    return np.array(pline)

def get_alllines(ddir, rootex, sigma, pixels, wls, ckind='median', dkind='phi'): 
    lines = {}
    for w in wls:
        fpaths = get_line_filenames(ddir, rootex, w, dkind=dkind)
        line = get_line(fpaths, sigma, pixels, ckind=ckind)
        lines[w] = line
    return lines

def get_allrootex(ddirs, root, exs, sigma, pixels, wls, ckind='median', dkind='phi'):
    all_lines = {}
    for ddir, ex in zip(ddirs, exs):
        rootex = '{0}_{1}'.format(root, ex) 
        lines = get_alllines(ddir, rootex, sigma, pixels, wls, ckind=ckind, dkind=dkind)
        all_lines[ex] = lines
    return all_lines

def get_alldkinds(ddirs, root, exs, sigma, pixels, wls, dkinds, ckind='median'):
    all_dkinds = {}
    for dkind in dkinds:
        all_rootex = get_allrootex(ddirs, root, exs, sigma, pixels, wls, ckind=ckind, dkind=dkind)
        all_dkinds[dkind] = all_rootex 
    return all_dkinds

def save_alldkinds(fname, d):
    # save to hdf5 files
    f = h5py.File(fname, 'w')
    for dk in d.keys():
        for ex in d[dk].keys():
            for w in d[dk][ex].keys():
                gname = '/{0}/{1}/{2}'.format(dk, ex, w)
                #dset = f.create_dataset(gname, shape=d[dk][ex][w].shape, 
                #                        dtype=d[dk][ex][w].dtype, data=d[dk][ex][w], 
                #                        chunks=True, compression='gzip',compression_opts=9)
                f[gname] = d[dk][ex][w]
    f.close()

def load_alldkinds(fname):
    f = h5py.File(fname, 'r')
    d = {}
    for dk in f.keys():
        d[str(dk)] = {}
        for ex in f[dk].keys():
            d[str(dk)][str(ex)] = {}
            for w in f[dk][ex].keys():
                #print dk, ex, w
                #print np.array(f[dk][ex][w])
                d[str(dk)][str(ex)][int(w)] = np.array(f[dk][ex][w])
    f.close()
    return d

def display_alllines(ddir, rootex, sigma, pixels, wls, ckind='median', dkind='phi'):
    all_lines =  get_alllines(ddir, rootex, sigma, pixels, wls, ckind=ckind, dkind=dkind) 
    import matplotlib.pyplot as plt
    for key in sorted(all_lines.keys()):
        line = all_lines[key]
        if dkind is 'phi':
            line = line - line[0]
            plt.ylabel(r'$\phi - \phi_0$ (rad)')
        else:
            line = line/line[0]
            plt.ylabel(r'$A/A_0$')
        plt.plot(line, label='{0}'.format(key))
        plt.legend()
        plt.show()

class RowMax:
    def __init__(self, ifile, ofile):
        self.pl = ParamList(ifile)
        self.ofile = ofile
        return

    def read_data(self):
        #SourcePlate/20130711_G3040_TwoTarget_IndiaNigro_SourcePlate_picture.fits
        ddir = self.pl.get_val('DATA_DIRECTORY')
        sname = self.pl.get_val('STUDY_NAME')
        root = self.pl.get_val('ROOT')
        fname = "{0}_SourcePlate_picture.fits".format(root)
        fname = os.path.join(ddir, sname, 'SourcePlate', fname)
        import bopy.io as bio
        data = bio.read_frame(fname)
        nd = np.shape(data)
        if len(nd) == 3:
            data = data[0]
        # now reorient
        from bopy.utils import data_reorient_udrot
        data = data_reorient_udrot(data, self.pl)
        return data

    def set_rowmax(self, rowmax):
        self.pl.set_val('GEN3FIT_PICKOFF_MASK_ROWS_UPTO', rowmax)
        return

    def write_cfg(self):
        self.pl.write_to_file(self.ofile)
        return
        
class ColMin:
    def __init__(self, ifile, ofile):
        self.pl = ParamList(ifile)
        self.ofile = ofile
        return

    def read_data(self):
        #SourcePlate/20130711_G3040_TwoTarget_IndiaNigro_SourcePlate_picture.fits
        ddir = self.pl.get_val('DATA_DIRECTORY')
        sname = self.pl.get_val('STUDY_NAME')
        root = self.pl.get_val('ROOT')
        fname = "{0}_SourcePlate_picture.fits".format(root)
        fname = os.path.join(ddir, sname, 'SourcePlate', fname)
        import bopy.io as bio
        data = bio.read_frame(fname)
        nd = np.shape(data)
        if len(nd) == 3:
            data = data[0]
        # now reorient
        from bopy.utils import data_reorient_udrot
        data = data_reorient_udrot(data, self.pl)
        return data

    def set_colmin(self, colmin):
        self.pl.set_val('GEN3FIT_PICKOFF_MASK_COLS_FROM', colmin)
        return

    def write_cfg(self):
        self.pl.write_to_file(self.ofile)
        return
        

class CCDRegion:
    def __init__(self, ifile, ofile, imgfname):
        self.pl = ParamList(ifile)
        self.ofile = ofile
        self.imgfname = imgfname
        return

    def read_data(self):
        #SourcePlate/20130711_G3040_TwoTarget_IndiaNigro_SourcePlate_picture.fits
        import bopy.io as bio
        data = bio.read_frame(self.imgfname)
        nd = np.shape(data)
        if len(nd) == 3:
            data = data[0]
        self.img_max_index = np.array(data.shape) - 1
        return data

    def set(self, x0, z0):
        s = '{0} {1} {2} {3}'.format(z0, self.img_max_index[0], 0, x0)
        self.pl.set_val('GEN3FIT_PICKOFF', s)
        return

    def write_cfg(self):
        self.pl.write_to_file(self.ofile)
        return
        
