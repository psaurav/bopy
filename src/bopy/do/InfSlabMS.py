import os, sys
import numpy as np
import bopy as bp
import scipy.optimize as opt
import numpy.ma as ma
import matplotlib.pyplot as plt

import bopy.spectra.spectras as spectras

def get_Contini_A(n):
    n2 = n*n
    n3 = n2*n
    n4 = n2*n2
    if n < 1:                      # Eq (A2)
        A =  (3.084635 - 
              6.531194*n + 
              8.357854*n2 - 
              5.082751*n3 + 
              1.171382*n4)
    elif n > 1:                    # Eq (A3)
        n5    = n2*n3
        n6    = n3*n3
        n7    = n3*n4
        A =  (504.332889 - 
              2641.00214*n + 
              5923.699064*n2 - 
              7376.355814*n3 + 
              5507.53041*n4 - 
              2463.357945*n5 + 
              610.956547*n6 - 
              64.8047*n7)
    else:
        A = 1
    return A

def greens_func11(S, mua, musp, ls, ld, x, y, L, xs, ys, mx, wbyv): 
    Rho2 = (x-xs)**2 + (y-ys)**2 
    k = np.sqrt((-mua + 1j*wbyv)*3*musp)
    zs0 = 1/musp
    fluence = np.zeros_like(Rho2, dtype=complex) 
    #print fluence.shape
    if mx != 0: 
        if np.abs(mx) % 2 == 0: 
            x_odd = np.arange(-np.abs(mx)+1, np.abs(mx), 2)
            x_even = np.arange(-np.abs(mx), np.abs(mx)+1, 2) 
        else: 
            x_odd = np.arange(-np.abs(mx), np.abs(mx)+1, 2)
            x_even = np.arange(-np.abs(mx)+1, np.abs(mx), 2)
        for i in range(len(x_even)): 
            l = ls 
            pm = 2*L + 4*l
            zpm = x_even[i]*pm + zs0
            zmm = x_even[i]*pm - 2*l - zs0 
            rp = np.sqrt(Rho2 + (L-zpm)**2) 
            rm = np.sqrt(Rho2 + (L-zmm)**2)
            fluence += (np.exp(1j*k*rp)/rp - np.exp(1j*k*rm)/rm)
        for i in range(len(x_odd)): 
            l = ld 
            pm = 2*L + 4*l
            zpm = x_odd[i]*pm + zs0
            zmm = x_odd[i]*pm - 2*l - zs0 
            rp = np.sqrt(Rho2 + (L-zpm)**2) 
            rm = np.sqrt(Rho2 + (L-zmm)**2) 
            fluence += (np.exp(1j*k*rp)/rp - np.exp(1j*k*rm)/rm) 
    else:
        zp = zs0 
        zm = -2*ld - zs0
        rp = np.sqrt(Rho2 + (L-zp)**2) 
        rm = np.sqrt(Rho2 + (L-zm)**2) 
        fluence = (np.exp(1j*k*rp)/rp - np.exp(1j*k*rm)/rm)

    #print fluence.shape
    return (3*musp*S)/(4*np.pi)*fluence 
    
def greens_func(S, mua, musp, l, x, y, L, xs, ys, mx, wbyv):
    """This function calculates the Green's function for the slab"""
    Rho2 = (x - xs)**2 + (y - ys)**2
    k = np.sqrt((-mua + 1j*wbyv)*3*musp)
    fluence = np.zeros_like(Rho2, dtype = complex)
    #l = 2/(3*musp)*self.A
    zs0 = 1/musp
    pm = 2*L + 4*l
    #compute the fluence rate
    for i in range(-np.abs(mx),np.abs(mx)+1):
        zpm = i*pm + zs0
        zmm = i*pm - 2*l - zs0
        rp = np.sqrt(Rho2 + (L - zpm)**2)
        rm = np.sqrt(Rho2 + (L - zmm)**2)
        fluence += (np.exp(1j*k*rp)/rp - np.exp(1j*k*rm)/rm)
    return (3*musp*S)/(4*np.pi)*fluence

def trans_contini_slab(S, mua, musp, l, x, y, L, xs, ys, mx=7, wbyv=0):
    """ Transmittance, from Contini (1997) eq 46
    """
    l = np.abs(l)
    z0 = 1.0/musp
    if wbyv == 0:
        k = 3.0*mua*musp
    else:
        k = 3.0*musp*(mua - 1j*wbyv)
    rho2 = (x-xs)**2 + (y-ys)**2
    trans = np.zeros_like(rho2, dtype = complex)
    mx = np.abs(mx)
    for m in range(-mx, mx+1):
        z1m = L*(1.0 - 2.0*m) - 4.0*m      *l - z0
        z2m = L*(1.0 - 2.0*m) -(4.0*m -2.0)*l + z0
        r12 = rho2 + z1m**2
        r22 = rho2 + z2m**2
        k1sq = np.sqrt(k*r12)
        k2sq = np.sqrt(k*r22)
        trans += z1m*r12**(-1.5)*(1.0+k1sq)*np.exp(-k1sq) - z2m*r22**(-1.5)*(1.0+k2sq)*np.exp(-k2sq)
    return S*trans/(4.0*np.pi)

def read_results_file(filename, src_index):
    d = np.loadtxt(filename)[src_index-1,:]
    pos_x = d[8]
    pos_z = d[9]
    S = d[5]
    p = d[10]
    return (S, p, pos_x, pos_z)

def get_emat(wl_nm, blood_spec='prahl', h2o_spec='segelstein81', fat_spec='vanveen'):
    return spectras.get_emat_tissue_mm(wl_nm, blood_spec=blood_spec, h2o_spec=h2o_spec, fat_spec=fat_spec)

class InfSlabMS(object):
    """Optimization class for the slab solution experiment."""
    #def __init__(self, func, slab_d, n_in, n_out, wl_nm, seriessum=7, blood_spec='prahl', h2o_spec='segelstein81', fat_spec='vanveen', A=None):
    def __init__(self, func, slab_d, n_in, n_out, wl_nm, emat, seriessum=7, A=None):
        """Initialization method that reads the images of the specified source"""
        
        c_light = 299792458000 # mm/s
        n = n_in/n_out
        self.v = c_light/n_in

        if not A == None:
            self.A = A
        else:
            self.A = get_Contini_A(n)
        self.slab_d = slab_d
        self.seriessum = seriessum
        if func == 'GREENS_FUNC':
            self.func = greens_func
        if func == 'GREENS_FUNC_LSLD':
            self.func = greens_func11
        elif func == 'CONTINI':
            self.func = trans_contini_slab

        self.wls = wl_nm
        self.nwls = len(wl_nm)
        #self.emat = spectras.get_emat_tissue_mm(wl_nm, 
        #                               blood_spec=blood_spec, 
        #                               h2o_spec=h2o_spec, 
        #                               fat_spec=fat_spec)
        self.emat = emat
        

    def set_gen3data(self, Apaths, Ppaths, Xdet, Zdet, omega, mask=None):

        self.omega = omega
        self.wbyv = omega/self.v
        self.amshape = mask.shape
        self.mask = mask

        for Apath, Ppath in zip(Apaths, Ppaths):
            # read the data
            Amp = bp.io.read_frame(Apath)
            Pha = bp.io.read_frame(Ppath)

            # mask all data arrays
            Amp = ma.compressed(ma.masked_array(Amp, self.mask))
            Pha = ma.compressed(ma.masked_array(Pha, self.mask))

            try:
                Amps = np.vstack((Amps, Amp))
                Phas = np.vstack((Phas, Pha))
            except:
                Amps = Amp
                Phas = Pha
        self.Amps = Amps
        self.Phas = Phas

        # masked detector positions
        self.dIdx = np.arange(len(Xdet))
        self.dIdx = ma.compressed(ma.masked_array(self.dIdx.reshape(self.amshape), self.mask))
        self.Xdet = ma.compressed(ma.masked_array(Xdet.reshape(self.amshape), self.mask))
        self.Zdet = ma.compressed(ma.masked_array(Zdet.reshape(self.amshape), self.mask))


    def opt_hbhbo2A(self, hb, hbo2, fh2o, ffat, A, b, src, S=1., phi0=1.):

        self.src_pos = src
        self.fh2o = fh2o
        self.ffat = ffat
        self.b = b
        self.xsrc = self.src_pos[0]
        self.zsrc = self.src_pos[1]
        Ss = np.ones(len(self.wls))*S
        Phi0s = np.ones(len(self.wls))*phi0
        hbsq = np.sqrt(hb)
        hbo2sq = np.sqrt(hbo2)
        mua_bl = 1.
        mua_blsq = np.sqrt(mua_bl)
        self.npar = 4
        self.conc = np.array([hb, hbo2, fh2o, ffat, mua_bl])
        #               0     1       2         3
        xin = np.array([hbsq, hbo2sq, mua_blsq, A])
        xin = np.append(xin, Ss)
        xin = np.append(xin, Phi0s)
        out = opt.leastsq(self.err_hbhbo2A, xin, factor = 100, maxfev = 2000, xtol=1e-12, ftol=1e-12)
        xopt = out[0]
        succ = out[1]
        #print >>sys.stderr, 'succ:', succ
        #print >>sys.stderr, 'Cs:', out[0][:4]
        #print >>sys.stderr, 'Ss:', out[0][4:4+len(self.wls)]
        #print >>sys.stderr, 'Ps:', out[0][4+len(self.wls):]
        ov = {'hb' : xopt[0]*xopt[0], 'hbo2' : xopt[1]*xopt[1], 'mua_bl' : xopt[2]*xopt[2], 'A' : xopt[3], 'succ' : succ}
        print >>sys.stderr, ov
        return ov
        

    def err_hbhbo2A(self, x):
        
        hb = x[0]*x[0]
        hbo2 = x[1]*x[1]
        mua_bl = x[2]*x[2] 
        A = x[3]
        Ss = x[self.npar:self.npar+self.nwls] 
        Phi0s = x[self.npar+self.nwls:]
        self.conc[0] = hb
        self.conc[1] = hbo2
        self.conc[4] = mua_bl 
        muas = spectras.mua_tissue_mm(self.emat, self.conc)
        musps = spectras.musp_Ab_mm(self.wls, A, self.b)

        diff = np.array([])
        
        for mua, musp, S, Phi0, Amp, Pha in zip(muas, musps, Ss, Phi0s, self.Amps, self.Phas): 
            l = 2./(3.*np.abs(musp))*self.A
            g = self.func(S, mua, musp, l, self.Xdet, self.Zdet, self.slab_d, self.xsrc, self.zsrc, self.seriessum, self.wbyv)
            a = np.abs(g)
            p = np.angle(g)
            # phase unwrap
            p0 = p[0]
            p[p-p0<-np.pi] += 2.*np.pi
            p[p-p0>np.pi] -= 2.*np.pi

            diff_a = (np.log(a) - np.log(Amp))
            diff_p = p - (Pha - Phi0)

            diff = np.append(diff, diff_a.ravel()) 
            diff = np.append(diff, diff_p.ravel())
        return diff
        
    def opt_hbhbo2al0b(self, hb, hbo2, fh2o, ffat, al0, b, src, l0=785, S=1., phi0=1.):

        self.wlsbyl0 = np.wls/l0
        self.l0 = l0
        self.src_pos = src
        self.fh2o = fh2o
        self.ffat = ffat
        self.xsrc = self.src_pos[0]
        self.zsrc = self.src_pos[1]
        Ss = np.ones(len(self.wls))*S
        Phi0s = np.ones(len(self.wls))*phi0
        hbsq = np.sqrt(hb)
        hbo2sq = np.sqrt(hbo2)
        al0sq = np.sqrt(al0)
        bsq = np.sqrt(b)
        mua_bl = 1.
        mua_blsq = np.sqrt(mua_bl)
        self.conc = np.array([hb, hbo2, fh2o, ffat, mua_bl])
        #               0     1       2         3      4
        xin = np.array([hbsq, hbo2sq, mua_blsq, al0sq, bsq])
        self.npar = len(xin)
        xin = np.append(xin, Ss)
        xin = np.append(xin, Phi0s)
        out = opt.leastsq(self.err_hbhbo2alob, xin, factor = 100, maxfev = 2000, xtol=1e-12, ftol=1e-12)
        xopt = out[0]
        succ = out[1]
        ov = {'hb' : xopt[0]*xopt[0], 'hbo2' : xopt[1]*xopt[1], 'mua_bl' : xopt[2]*xopt[2], 'al0' : xopt[3]*xopt[3], 'b': xopt[4]*xopt[4], 'succ' : succ}
        print >>sys.stderr, ov
        return ov

    def err_hbhbo2al0b(self, x):
        
        hb = x[0]*x[0]
        hbo2 = x[1]*x[1]
        mua_bl = x[2]*x[2] 
        al0 = x[3]*x[3]
        b = x[4]*x[4]
        Ss = x[self.npar:self.npar+self.nwls] 
        Phi0s = x[self.npar+self.nwls:]
        self.conc[0] = hb
        self.conc[1] = hbo2
        self.conc[4] = mua_bl 
        muas = spectras.mua_tissue_mm(self.emat, self.conc)
        musps = al0*np.power(self.wlsbyl0, -b)

        diff = np.array([])
        
        for mua, musp, S, Phi0, Amp, Pha in zip(muas, musps, Ss, Phi0s, self.Amps, self.Phas): 
            l = 2./(3.*np.abs(musp))*self.A
            g = self.func(S, mua, musp, l, self.Xdet, self.Zdet, self.slab_d, self.xsrc, self.zsrc, self.seriessum, self.wbyv)
            a = np.abs(g)
            p = np.angle(g)
            # phase unwrap
            p0 = p[0]
            p[p-p0<-np.pi] += 2.*np.pi
            p[p-p0>np.pi] -= 2.*np.pi

            diff_a = (np.log(a) - np.log(Amp))
            diff_p = p - (Pha - Phi0)

            diff = np.append(diff, diff_a.ravel()) 
            diff = np.append(diff, diff_p.ravel())
        
    def opt_run(self, mua, musp, src, S=1., phi0=1.):
        """A function that performs an ordinary optimization run.
        The parameters being optimized are:
        S, mua, musp
        """
        phi0 = 1.
        self.src_pos = src
        #               0  1    2     3                4                5
        xin = np.array([S, mua, musp, self.src_pos[0], self.src_pos[1], phi0])
        try:
            out = opt.leastsq(self.obj_fun_w_position, xin,
                factor = 100, maxfev = 2000, xtol=1e-12, ftol=1e-12)
            xopt = out[0]
            err = out[1]
        except:
            xopt = [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan]
            err = 100

        self.src_pos = [xopt[3], xopt[4]]
        self.l = 2./(3.*np.abs(xopt[2]))*self.A
        return {'S' : xopt[0]*1e9, 'mua' : xopt[1], 'musp' : xopt[2], 'l' : self.l,
                'xsrc' : xopt[3], 'zsrc' : xopt[4], 'phi0' : xopt[5], 'ier' : err}

    def obj_fun_onlylsld(self, x):
        """Objective function for opt_run_normal"""
        ls = np.abs(x[0])
        ld = np.abs(x[1])
        f = self.func(self.S, self.mua, self.musp, ls, ld, self.Xdet, self.Zdet, self.slab_d, self.src_pos[0], self.src_pos[1], self.seriessum, self.wbyv)
        a = np.abs(f)
        p = np.angle(f)
        diff_a = (np.log(a) - np.log(self.Amp))
        diff_p = p - (self.Pha - self.phi0)
        return np.concatenate((diff_a.ravel(), diff_p.ravel()))

    def obj_fun_muamusp_withlsld(self, x):
        """Objective function for opt_run_normal"""
        S = 1e9*np.abs(x[0])
        mua = np.abs(x[1])
        musp = np.abs(x[2])
        f = self.func(S, mua, musp, self.ls, self.ld, self.Xdet, self.Zdet, self.slab_d, self.src_pos[0], self.src_pos[1], self.seriessum, self.wbyv)
        a = np.abs(f)
        p = np.angle(f)
        diff_a = (np.log(a) - np.log(self.Amp))
        diff_p = p - (self.Pha - x[3])
        return np.concatenate((diff_a.ravel(), diff_p.ravel()))

    def obj_fun_lsld(self, x):
        """Objective function for opt_run_normal"""
        S = 1e9*np.abs(x[0])
        ls = np.abs(x[1])
        ld = np.abs(x[2])
        f = self.func(S, self.mua, self.musp, ls, ld, self.Xdet, self.Zdet, self.slab_d, self.src_pos[0], self.src_pos[1], self.seriessum, self.wbyv)
        a = np.abs(f)
        p = np.angle(f)
        diff_a = (np.log(a) - np.log(self.Amp))
        diff_p = p - (self.Pha - x[3])
        return np.concatenate((diff_a.ravel(), diff_p.ravel()))

    def obj_fun_lSphi(self, x):
        """Objective function for opt_run_normal"""
        l = np.abs(x[0])
        S = 1e9*np.abs(x[1])
        phi0 = np.abs(x[2])
         
        f = self.func(S, self.mua, self.musp, l, self.Xdet, self.Zdet, self.slab_d, self.src_pos[0], self.src_pos[1], self.seriessum, self.wbyv)
        a = np.abs(f)
        p = np.angle(f)
        diff_a = (np.log(a) - np.log(self.Amp))
        diff_p = p - (self.Pha - phi0)
        return np.concatenate((diff_a.ravel(), diff_p.ravel()))

    def obj_fun_l(self, x):
        """Objective function for opt_run_normal"""
        S = 1e9*np.abs(x[0])
        l = np.abs(x[1])
        f = self.func(S, np.abs(x[1]), np.abs(x[2]), l, self.Xdet, self.Zdet, self.slab_d, x[2], x[3], self.seriessum, self.wbyv)
        a = np.abs(f)
        p = np.angle(f)
        diff_a = (np.log(a) - np.log(self.Amp))
        diff_p = p - (self.Pha - x[4])
        return np.concatenate((diff_a.ravel(), diff_p.ravel()))

    def obj_fun_w_position(self, x):
        """Objective function for opt_run_normal"""
        S = 1e9*np.abs(x[0])
        l = 2./(3.*np.abs(x[2]))*self.A
        f = self.func(S, np.abs(x[1]), np.abs(x[2]), l, self.Xdet, self.Zdet, self.slab_d, x[3], x[4], self.seriessum, self.wbyv)
        a = np.abs(f)
        p = np.angle(f)
        diff_a = (np.log(a) - np.log(self.Amp))
        diff_p = p - (self.Pha - x[5])
        return np.concatenate((diff_a.ravel(), diff_p.ravel()))

    def print_A_phi(self, S, mua, musp, xsrc, zsrc, phi0, root, w, s):
        l = self.l
        g = self.func(S, mua, musp, l, self.Xdet, self.Zdet, self.slab_d, xsrc, zsrc, self.seriessum, self.wbyv)
        r = np.sqrt(self.slab_d**2 + (self.Xdet - xsrc)**2 + (self.Zdet - zsrc)**2)
        alnrA = np.log(np.abs(g)*r)
        dlnrA = np.log(self.Amp*r)

        filename = "lnrA_%s_%d.png" % (w, s)
        plt.cla()
        plt.plot(r, dlnrA, 'b,', label='Data')
        plt.plot(r, alnrA, 'r.', label='Fit')
        plt.xlabel(r'$r$ (mm)')
        plt.ylabel(r'$\log(rA)$')
        plt.title("%s (w=%s, s=%d)"%(root, w, s))
        plt.legend()
        plt.savefig(filename)

        filename = "phi_%s_%d.png" % (w, s)
        plt.cla()
        plt.plot(r, self.Pha, 'b,', label='Data')
        plt.plot(r, np.angle(g)+phi0, 'r.', label='Fit')
        plt.xlabel('$r$ (mm)')
        plt.ylabel(r'$\phi$ (radian)')
        plt.title("%s (w=%s, s=%d)"%(root, w, s))
        plt.legend()
        plt.savefig(filename)

        filename = "dat_%s_%d.txt" % (w, s)
        r = self.dIdx.astype(int)/int(self.amshape[0])
        c = self.dIdx.astype(int)%int(self.amshape[0])
        d = np.vstack((self.dIdx, r, c, self.Xdet, self.Zdet, self.Amp, self.Pha)).T
        #np.savetxt(filename, d, fmt='%.0f %.0f %.0f %f %f %f %f', header="index row col xpos zpos Amp Pha")
        np.savetxt(filename, d, fmt='%.0f %.0f %.0f %f %f %f %f')
    
    def print_A_phi_lsld(self, S, mua, musp, xsrc, zsrc, phi0, root, w, s):
        #print S/1.e9
        ls = self.ls
        ld = self.ld
        g = self.func(S, mua, musp, ls, ld, self.Xdet, self.Zdet, self.slab_d, xsrc, zsrc, self.seriessum, self.wbyv)
        r = np.sqrt(self.slab_d**2 + (self.Xdet - xsrc)**2 + (self.Zdet - zsrc)**2)
        alnrA = np.log(np.abs(g)*r)
        dlnrA = np.log(self.Amp*r)

        filename = "lnrA_lsld_%s_%d.png" % (w, s)
        plt.cla()
        plt.plot(r, dlnrA, 'b,', label='Data')
        plt.plot(r, alnrA, 'r.', label='Fit')
        plt.xlabel(r'$r$ (mm)')
        plt.ylabel(r'$\log(rA)$')
        plt.title("%s (w=%s, s=%d)"%(root, w, s))
        plt.legend()
        plt.savefig(filename)

        filename = "phi_lsld_%s_%d.png" % (w, s)
        plt.cla()
        plt.plot(r, self.Pha, 'b,', label='Data')
        plt.plot(r, np.angle(g)+phi0, 'r.', label='Fit')
        plt.xlabel('$r$ (mm)')
        plt.ylabel(r'$\phi$ (radian)')
        plt.title("%s (w=%s, s=%d)"%(root, w, s))
        plt.legend()
        plt.savefig(filename)

        filename = "dat_lsld_%s_%d.txt" % (w, s)
        r = self.dIdx.astype(int)/int(self.amshape[0])
        c = self.dIdx.astype(int)%int(self.amshape[0])
        d = np.vstack((self.dIdx, r, c, self.Xdet, self.Zdet, self.Amp, self.Pha)).T
        #np.savetxt(filename, d, fmt='%.0f %.0f %.0f %f %f %f %f', header="index row col xpos zpos Amp Pha")
        np.savetxt(filename, d, fmt='%.0f %.0f %.0f %f %f %f %f')
