import os
import sys
import numpy as np

class Gen3Info:
    def __init__(self, filename):

        from bopy.utils import ParamList
        self.prm = ParamList(filename)
        self.wls = np.sort(np.array([int(w) for w in self.prm.get_val('WAVELENGTHS').split()]))

    def read_cfg(self, filename):
        """
        Read more cfg files
        """
        sefl.prm.read_from_file(filename)

    def get_wavelengths(self):
        return self.wls

    def get_root(self):
        return self.prm.get_val('ROOT')

    def get_sources(self):
        return np.array([int(v) for v in self.prm.get_val('SOURCES').split()])

    def get_slabthickness(self):
        return float(self.prm.get_val('BREAST_TANK_THICKNESS_MM'))

    def get_refmua(self):
        tmp = []
        for w in self.wls:
            mua = float(self.prm.get_val('CALCULATED_MUA_MM_{0}'.format(w)))
            mua = tmp.append(mua)
        return np.array(tmp)

    def get_refmusp(self):
        tmp = []
        for w in self.wls:
            musp = float(self.prm.get_val('CALCULATED_MUSP_MM_{0}'.format(w)))
            musp = tmp.append(musp)
        return np.array(tmp)

    def get_power_b(self):
        return float(self.prm.get_val('CALCULATED_SCATTERING_POWER_B'))
    
    def get_prefactor_a(self):
        return float(self.prm.get_val('CALCULATED_SCATTERING_PREFACTOR_A'))

    def get_breast_tank_thickness(self):
        return float(self.prm.get_val('BREAST_TANK_THICKNESS_MM'))

    def get_list_data_with_dark(self, file_format=None):
        ddir = self.prm.get_val('DATA_DIRECTORY')
        if ddir.endswith('/'):
            ddir = ddir[:-1]
        sname = self.prm.get_val('STUDY_NAME')
        mname = self.prm.get_val('MEASUREMENT_NAME')
        if mname is None:
            mname = self.prm.get_val('GEN3FIT_MEASUREMENT')
        srcs = self.prm.get_val('SOURCES')
        print ddir, sname, mname
        if srcs == '0':
            srcs = np.arange(209).astype(int) + 1
        else:
            srcs = np.sort(np.array([int(v) for v in srcs.split()]))
        
        if file_format is None:
            file_format = '{0}_{1}_wl{2}_s{3}.fits'
        f = '{0}_{1}_dark.fits'.format(sname, mname)
        f = os.path.join(ddir, sname, mname, f)
        l = [f, ] 
        for w in self.wls:
            for s in srcs:
                f = file_format.format(sname, mname, w, s)
                f = os.path.join(ddir, sname, mname, f)
                l.append(f)
        return l
