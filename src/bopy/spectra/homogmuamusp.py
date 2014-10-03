import os
import sys
import bopy as bp
import numpy as np
import scipy
 
class HomogMuaMusp:

    def __init__(self, filename):
        self.file_dirname = os.path.dirname(filename)
        self.p = bp.utils.ParamList(filename)
        self._get_values()

    def _get_values(self):
        self.wl_nm = np.array([np.double(v) for v in self.p.get_val("WAVELENGTHS").split()])
        print "Wavelengths:", self.wl_nm

        #waterfile = p.get_val("WATER_SPECT_FILE")
        self.mua_water_cm = bp.spectra.Spectras().muaH2O_cm(self.wl_nm)
        print "mua_water_cm =", self.mua_water_cm

        spect_filename = self.p.get_val("INK_SPECTROPHOTOMETER_FILE")
        full_spect_filename = os.path.join(self.file_dirname, spect_filename)
        if not os.path.exists(full_spect_filename):
            full_spect_filename = os.path.join(self.p.get_val('DATA_DIRECTORY'),
                                               self.p.get_val('STUDY_NAME'),
                                               spect_filename)

        self.mua_spect_cm = np.log(10.0)*bp.spectra.Spectras().absorbance_spectfile(full_spect_filename, self.wl_nm)
        print "mua_spect_cm =", self.mua_spect_cm

        spect_dilution = float(self.p.get_val("INK_SPECT_DILUTION_PERCENTAGE"))/100.0
        self.mua_stock_cm = self.mua_spect_cm/spect_dilution
        print "mua_stock_cm =", self.mua_stock_cm

        #
        # volumes
        #
        self.vol_water_mL = float(self.p.get_val("HOMOG_WATER_VOLUME_ORIGINAL_ML"))
        self.vol_water_removed_mL = float(self.p.get_val("HOMOG_WATER_VOLUME_REMOVED_ML"))
        self.vol_water_mL = self.vol_water_mL - self.vol_water_removed_mL
        self.vol_intralipid_added_mL = float(self.p.get_val("HOMOG_INTRALIP_VOLUME_ADDED_ML"))
        self.vol_water_and_intralipid_mL = self.vol_water_mL + self.vol_intralipid_added_mL
        self.vol_stock_mL = float(self.p.get_val("HOMOG_INK_VOLUME_ADDED_ML"))
        print "final water+intralip volume (mL):", self.vol_water_and_intralipid_mL
        print "ink stock volume (mL):", self.vol_stock_mL

    def get_mua(self):
        mua_bath_cm = (self.mua_water_cm*self.vol_water_and_intralipid_mL + self.mua_stock_cm*self.vol_stock_mL 
                      )/(self.vol_water_and_intralipid_mL + self.vol_stock_mL)
        mua_bath_mm = mua_bath_cm/10.
        print "mua_bath_mm = ", mua_bath_mm
        for w, m in zip(self.wl_nm, mua_bath_mm):
            self.p.set_val('CALCULATED_MUA_MM_%d'%(int(w)), m)

    def get_musp(self):
        print "==Scattering=="
        intralipid_pct_vanstaveren = 10.0
        vanstaveren_unit_conversion = 1000.0*10.0 # from mL^-1 L mm^-1 to cm^-1
        intralipid_pct = 20.0
        print >>sys.stderr, 'Intralipid percentage:', intralipid_pct
        intralipid_factor = intralipid_pct/intralipid_pct_vanstaveren
        wl_um = self.wl_nm*1.0e-3
        g = 1.1-0.58*wl_um
        print >>sys.stderr, "g:", g

        musp_intralipid_10pct = (1.0-g)*0.016*np.power(wl_um, -2.4) * vanstaveren_unit_conversion

        musp_intralipid = musp_intralipid_10pct*intralipid_factor*self.vol_intralipid_added_mL/self.vol_water_and_intralipid_mL
        musp_intralipid_mm = musp_intralipid/10.0
        print "musp intralipid:", musp_intralipid_mm
        for w, m in zip(self.wl_nm, musp_intralipid_mm):
            self.p.set_val('CALCULATED_MUSP_MM_%d'%(int(w)), m)

        fit_func = lambda p, x: p[0]*(x**(-p[1]))
        err_func = lambda p, x, y: (y - fit_func(p,x))
        p = np.array([1000.0, 1.0])
        p1, success = scipy.optimize.leastsq(err_func, p[:], args=(self.wl_nm, musp_intralipid_mm))
        print "fitting success:", success
        #print "p1:", p1
        print "musp_check:", fit_func(p1,self.wl_nm)

        print "SCATTERING_PREFACTOR_A = HOMOG %g" % (p1[0])
        print "SCATTERING_POWER_B = HOMOG %g" %(p1[1])

        self.p.set_val('CALCULATED_SCATTERING_PREFACTOR_A', p1[0])
        self.p.set_val('CALCULATED_SCATTERING_POWER_B', p1[1])

    def output_values(self):
        filename = self.p.get_val('STUDY_NAME')+'-homog.cfg'
        self.p.write_to_file(filename)
