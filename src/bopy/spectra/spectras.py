import os
import sys
import numpy as np
import scipy.interpolate
from scipy import stats
import bopy as bp

def fit_scattering_Ab(wls, musp):
    """
    Formula:
        \mu_s^\prime = A \lambda^{-b}
    """
    fit_func = lambda p, x: p[0]*(x**(-p[1]))
    err_func = lambda p, x, y: (y - fit_func(p,x))
    #div_func = lambda p, x: [x**(-p[1]), -p[0]*x**(-p[1]-1.0)]
    p = np.array([1000.0, 1.0])
    from scipy.optimize import leastsq
    #p1, success = leastsq(err_func, p[:], Dfun=div_func, args=(wls, musp))
    p1, success = leastsq(err_func, p[:], args=(wls, musp))
    return p1[0], p1[1], success

def musp2Ab(wls, musp):
    x = np.log(wls)
    y = np.log(musp)
    slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
    A = np.exp(intercept)
    b = -slope
    return A, b

def A2al0(A, b, l0):
    return A*np.power(l0, -b)

def al02A(al0, b, l0):
    return al0*np.power(l0, b)

def get_spectra_filename(spectra):
    """
    """
    if spectra == 'sophie36zijlstra':
        return os.path.join(bp.spectra._DATA_DIR, 
               'chromophores_vanveen_water36c_zijlstra_sophie.txt')
    elif spectra == 'zijlstra2000':
        return os.path.join(bp.spectra._DATA_DIR, 
               'zijlstra2000.dat')
    elif spectra == 'taroni2007':
        return os.path.join(bp.spectra._DATA_DIR,
                'collagen-absorption-taroni2007.txt')
    elif spectra == 'segelstein81':
        return os.path.join(bp.spectra._DATA_DIR,
                'segelstein81.dat')
    elif spectra == 'vanveen':
        return os.path.join(bp.spectra._DATA_DIR,
                'fat.txt')
    elif spectra == 'prahl':
        return os.path.join(bp.spectra._DATA_DIR,
                'hemoglobin.dat')

def lbs_spectra(filename, colnum, wl_nm):
    d = bp.io.read_float_table(filename, comments='%', colstart=0, colend=colnum+1, colskip=colnum)
    return bp.utils.interpolate1d(d[:,0], d[:,1], wl_nm)

def get_spect_wavelengths_nm(spectra):
    """
    """
    filename = get_spectra_filename(spectra)
    if spectra == 'sophie36zijlstra':
        d = np.loadtxt(filename, comments='%')
        return d[:,0]
        
def get_spect_muaHBO2_mm(cHBO2, spectra='sophie36zijlstra'):
    """
    """
    filename = get_spectra_filename(spectra)
    if spectra == 'sophie36zijlstra':
        d = np.loadtxt(filename, comments='%')
        return cHBO2*d[:,3]/1000.

def get_spect_muaHB_mm(cHBO2, spectra='sophie36zijlstra'):
    """
    """
    filename = get_spectra_filename(spectra)
    if spectra == 'sophie36zijlstra':
        d = np.loadtxt(filename, comments='%')
        return cHBO2*d[:,4]/1000.

def get_spect_muaH2O_mm(pH2O, spectra='sophie36zijlstra'):
    """
    """
    filename = get_spectra_filename(spectra)
    if spectra == 'sophie36zijlstra':
        d = np.loadtxt(filename, comments='%')
        return pH2O*d[:,2]/100.

def get_spect_muaFAT_mm(pFAT, spectra='sophie36zijlstra'):
    """
    """
    filename = get_spectra_filename(spectra)
    if spectra == 'sophie36zijlstra':
        d = np.loadtxt(filename, comments='%')
        return pFAT*d[:,1]/100.

def get_eHBO2_mm_uM(wl_nm, spectra='zijlstra2000'):
    filename = get_spectra_filename(spectra) 
    if spectra == 'zijlstra2000':
        d = np.loadtxt(filename, comments='%')
        w = d[:,0]
        hbo2 = d[:,1]/1000.
    elif spectra == 'sophie36zijlstra':
        d = np.loadtxt(filename, comments='%')
        w = d[:,0]
        hbo2 = d[:,3]/1000.
    elif spectra == 'prahl':
        d = np.loadtxt(filename, comments='#')
        w = d[:,0]
        hbo2 = d[:,1]*np.log(10.)/10000000.0
    else:
        sys.exit('In module {0}, spectra \'{1}\' not found in get_eHBO2_mm_uM().  Exiting'.format(__name__, spectra))
    f = scipy.interpolate.interp1d(w, hbo2)
    return f(wl_nm)
    
def get_eHB_mm_uM(wl_nm, spectra='zijlstra2000'):
    filename = get_spectra_filename(spectra) 
    if spectra == 'zijlstra2000':
        d = np.loadtxt(filename, comments='%')
        w = d[:,0]
        hb = d[:,2]/1000.
    elif spectra == 'sophie36zijlstra':
        d = np.loadtxt(filename, comments='%')
        w = d[:,0]
        hb = d[:,4]/1000.
    elif spectra == 'prahl':
        d = np.loadtxt(filename, comments='#')
        w = d[:,0]
        hb = d[:,2]*np.log(10.)/10000000.0
    else:
        sys.exit('In module {0}, spectra \'{1}\' not found in get_eHB_mm_uM().  Exiting'.format(__name__, spectra))
    f = scipy.interpolate.interp1d(w, hb)
    return f(wl_nm)

def get_mua_fh2o_mm(wl_nm, fh2o, spectra='zijlstra2000'):
    filename = get_spectra_filename(spectra) 
    if spectra == 'zijlstra2000':
        d = np.loadtxt(filename, comments='%')
        w = d[:,0]
        h2o = d[:,3]
    elif spectra == 'sophie36zijlstra':
        d = np.loadtxt(filename, comments='%')
        w = d[:,0]
        h2o = d[:,2]
    elif spectra == 'segelstein81':
        d = np.loadtxt(filename, comments='#')
        w = d[:,0]
        h2o = d[:,1]/10.
    else:
        sys.exit('In module {0}, spectra \'{1}\' not found in get_mua_fh2o_mm().  Exiting'.format(__name__, spectra))
    f = scipy.interpolate.interp1d(w, h2o)
    return f(wl_nm)*fh2o
    
def get_mua_ffat_mm(wl_nm, ffat, spectra='zijlstra2000'):
    filename = get_spectra_filename(spectra) 
    if spectra == 'zijlstra2000':
        d = np.loadtxt(filename, comments='%')
        w = d[:,0]
        fat = d[:,4]
    elif spectra == 'sophie36zijlstra':
        d = np.loadtxt(filename, comments='%')
        w = d[:,0]
        fat = d[:,1]
    elif spectra == 'vanveen':
        d = np.loadtxt(filename, comments='#')
        w = d[:,0]
        fat = d[:,1]/1000.
    else:
        sys.exit('In module {0}, spectra \'{1}\' not found in get_mua_ffat_mm().  Exiting'.format(__name__, spectra))
    f = scipy.interpolate.interp1d(w, fat)
    return f(wl_nm)*ffat

def get_mua_fcoll_mm(wl_nm, fcoll, spectra='taroni2007'):
    filename = get_spectra_filename(spectra)
    if spectra == 'taroni2007':
        d = np.loadtxt(filename, comments='%')
        w = d[:,0]
        coll = d[:,1]/10.
    f = scipy.interpolate.interp1d(w, coll)
    return f(wl_nm)*fcoll

def get_mua_ffat_fh2o_mm(wl_nm, ffat, fh2o, spectra='sophie36zijlstra'):
    mua_fat = get_mua_ffat_mm(wl_nm, ffat, spectra)
    mua_h2o = get_mua_fh2o_mm(wl_nm, fh2o, spectra)
    return mua_fat + mua_h2o

def get_emat_tissue_mm(wl_nm, blood_spec='prahl', h2o_spec='segelstein81', fat_spec='vanveen'):
    ehb = get_eHB_mm_uM(wl_nm, spectra=blood_spec)
    ehbo2 = get_eHBO2_mm_uM(wl_nm, spectra=blood_spec)
    mh2o = get_mua_fh2o_mm(wl_nm, 1.0, spectra=h2o_spec)
    mfat = get_mua_ffat_mm(wl_nm, 1.0, spectra=fat_spec)
    bl = np.ones(len(wl_nm))  
    print >>sys.stderr, 'wl (nm):', wl_nm
    print >>sys.stderr, 'eHb:', ehb
    print >>sys.stderr, 'eHbO2:', ehbo2
    print >>sys.stderr, 'mH2O:', mh2o
    print >>sys.stderr, 'mfat:', mfat
    print >>sys.stderr, 'basel:', bl
    return np.vstack((ehb, ehbo2, mh2o, mfat, bl)).T

def mua_tissue_mm(emat, conc):
    return emat.dot(conc)

def musp_Ab_mm(wl_nm, A, b):
    return A*np.power(wl_nm, -b)

def read_spect_mm(filename):
    x = np.loadtxt(filename, comments='>>', skiprows=20)
    w = x[:,0]
    v = np.log(10.0)*x[:,1]/10.0 
    return (w, v)

def read_stock_mm(filename, p):
    w, v = read_spect_mm(filename)
    v = 100.*v/p
    return (w, v)

def get_stock_conc_uMmL(wgt_gm, mol_wgt_gm, h2o_ml):
    return wgt_gm/mol_wgt_gm/h2o_ml*1.e6

def get_diluted_conc_uMmL(stock_conc_uMmL, stock_vol_mL, liquid_vol_mL):
    return stock_conc_uMmL*stock_vol_mL/(stock_vol_mL+liquid_vol_mL)

def get_extinction_spectra_uMmLmm(filename, wl_nm, p, wgt_gm, mol_wgt_gm, h2o_ml):
    (w, mua_stock_mm) = read_stock_mm(filename, p)
    mua_stock_mm = bp.utils.interpolate1d(w, mua_stock_mm, wl_nm)
    conc_stock_uMmL = get_stock_conc_uMmL(wgt_gm, mol_wgt_gm, h2o_ml)
    return mua_stock_mm/conc_stock_mm

def get_toast_lbs_breast_spectral_matlab(wl_nm, hbo2, hb, ph2o, pfat):
    ehbo2 = get_eHBO2_mm_uM(wl_nm, spectra='sophie36zijlstra')
    ehb = get_eHB_mm_uM(wl_nm, spectra='sophie36zijlstra')
    mua_h2o = get_mua_fh2o_mm(wl_nm, ph2o/100., spectra='sophie36zijlstra')
    mua_fat = get_mua_fh2o_mm(wl_nm, pfat/100., spectra='sophie36zijlstra')
    mua_bkg = mua_h2o + mua_fat
    
    print "Wavelengths"
    print wl_nm
    print "HbO2 extinction coefficients"
    print ehbo2
    print "Hb extinction coefficients"
    print ehb
    print "Background mua"
    print mua_bkg

def get_tissue_breakup(spectra, hbo2, hb, ph2o, pfat, baseline, font_size=None, 
                       g3_wls=None, mua_maxlim=None, w_meas=None, mua_meas=None):
    wls = get_spect_wavelengths_nm(spectra)
    mua_hbo2 = get_spect_muaHBO2_mm(hbo2, spectra)
    mua_hb = get_spect_muaHB_mm(hb, spectra)
    mua_h2o = get_spect_muaH2O_mm(ph2o, spectra)
    mua_fat = get_spect_muaFAT_mm(pfat, spectra)
    mua_base = np.ones(wls.shape[0])*baseline
    mua_tot = mua_hbo2 + mua_hb + mua_h2o + mua_fat + mua_base
    outputfile = 'tissue_spectra.png'
    import matplotlib 
    import matplotlib.pyplot as plt 
    if font_size is not None:
        matplotlib.rcParams.update({'font.size': font_size})
    plt.plot(wls, mua_hbo2, 'r', label=r'$\mathrm{HbO_2}$') 
    plt.plot(wls, mua_hb, 'b', label=r'$\mathrm{Hb}$') 
    plt.plot(wls, mua_h2o, 'c', label=r'$\mathrm{H_2O}$') 
    plt.plot(wls, mua_fat, 'y', label=r'$\mathrm{Fat}$') 
    plt.plot(wls, mua_tot, 'k', label='Total', linewidth=2) 
    if w_meas is not None:
        plt.plot(w_meas, mua_meas, 'k--', label='Measured', linewidth=2) 
    plt.xlabel('Wavelength [nm]')
    plt.ylabel(r'$\mu_a$ [mm$^{-1}$]')
    plt.legend(loc='upper left')
    if mua_maxlim is not None:
        plt.ylim((0, mua_maxlim))
    if g3_wls is not None:
        for w in g3_wls:
            plt.axvline(w, color='m', linewidth=0.5)
    plt.title('Tissue Absorption')
    plt.savefig(outputfile, bbox_inches='tight')
    print 'Saved to file:', outputfile


class Spectras:
    def __init__(self):
        pass

    def get_musp_vanstaveren_cm(self, ilipid_frac, wl_nm):
        ilipid_frac_vanstaveren = .1     # 10% intralipd in vanstaveren 
        vanstaveren_unit_conversion = 1000.0*10.0 # from mL^-1 L mm^-1 to cm^-1
        wl_um = wl_nm*1.0e-3
        g = 1.1-0.58*wl_um                  # eq.13
        mus = 0.016*np.power(wl_um, -2.4)   # eq.12
        return (1.0-g)*mus*vanstaveren_unit_conversion*ilipid_frac/ilipid_frac_vanstaveren

    def get_mua_stock_cm(self, filename, stock_frac, wl_nm):
        return self.get_absorbance_twocol_file(filename, wl_nm, comments='>>', skiprows=20)/stock_frac

    def get_mua_spect_and_h20_cm(self, filename, stock_frac, wl_nm, vol_h20, vol_stock):
        mua_stock = self.get_mua_stock_cm(filename, stock_frac, wl_nm)
        mua_h2o = self.muaH2O_mm(wl_nm)
        return (mua_stock*vol_stock + mua_h2o*vol_h20)/(vol_h2o + vol_stock)

    def get_absorbance_twocol_file(self, filename, wl_nm, comments='', skiprows=0):
        d = np.loadtxt(filename, comments=comments, skiprows=skiprows).astype(float)
        f_interpolate = scipy.interpolate.interp1d(d[:,0], d[:,1])
        d = np.array(f_interpolate(wl_nm))
        return d
    
    def get_absorbance_threecol_file(self, filename, wl_nm, comments='', skiprows=0):
        d = np.loadtxt(filename, comments=comments, skiprows=skiprows).astype(float)
        f1_interpolate = scipy.interpolate.interp1d(d[:,0], d[:,1])
        f2_interpolate = scipy.interpolate.interp1d(d[:,0], d[:,2])
        d1 = np.array(f1_interpolate(wl_nm))
        d2 = np.array(f2_interpolate(wl_nm))
        return np.vstack((d1, d2)).T 
        
    def absorbance_spectfile(self, filename, wl_nm):
        return self.get_absorbance_twocol_file(filename, wl_nm, comments='>>', skiprows=20)

    def muaH2O_cm(self, wl_nm):
        return self.get_absorbance_twocol_file(os.path.join(bp.spectra._DATA_DIR, 'segelstein81.dat'), 
                                               wl_nm, comments='#')
    
    def muaH2O_mm(self, wl_nm):
        return self.get_absorbance_twocol_file(os.path.join(bp.spectra._DATA_DIR, 'segelstein81.dat'), 
                                               wl_nm, comments='#')/10.0

    def muaLipid_cm(self, wl_nm):
        return self.get_absorbance_twocol_file(os.path.join(bp.spectra._DATA_DIR, 'lipid.dat'), 
                                               wl_nm, comments='#')
    
    def muaLipid_mm(self, wl_nm):
        return self.get_absorbance_twocol_file(os.path.join(bp.spectra._DATA_DIR, 'lipid.dat'), 
                                               wl_nm, comments='#')/10.0

    def extCoeffHB_cmpuM(self, wl_nm):
        cmpM2cmpuM  = 1000000.0        # cm^-1/M -> cm^-1/uM
        return self.get_absorbance_threecol_file(os.path.join(bp.spectra._DATA_DIR, 'hemoglobin.dat'), 
                                                 wl_nm, comments='#')*np.log(10.)/cmpM2cmpuM
    
    def extCoeffHB_mmpuM(self, wl_nm):
        cmpM2mmpuM  = 10000000.0       # cm^-1/M -> mm^-1/uM
        return self.get_absorbance_threecol_file(os.path.join(bp.spectra._DATA_DIR, 'hemoglobin.dat'), 
                                                 wl_nm, comments='#')*np.log(10.)/cmpM2mmpuM
    
    def extCoeffHB_cmpmM(self, wl_nm):
        cmpM2cmpmM  = 1000.0        # cm^-1/M -> cm^-1/mM
        return self.get_absorbance_threecol_file(os.path.join(bp.spectra._DATA_DIR, 'hemoglobin.dat'), 
                                                 wl_nm, comments='#')*np.log(10.)/cmpM2cmpmM
    
    def extCoeffHB_mmpmM(self, wl_nm):
        cmpM2mmpmM  = 10000.0       # cm^-1/M -> mm^-1/uM
        return self.get_absorbance_threecol_file(os.path.join(bp.spectra._DATA_DIR, 'hemoglobin.dat'), 
                                                 wl_nm, comments='#')*np.log(10.)/cmpM2mmpmM
