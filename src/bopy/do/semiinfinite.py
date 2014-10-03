import numpy as np
import bc

class SemiInfinite:
    def __init__(self, nin, nout, solntype='fluence'):
        print 'Setting up interface:'
        self._n = nin/nout
        self._A = bc.integrate_A(self._n)
        self._flc = bc.fluence_coeff(self._n)
        self._fxc = bc.flux_coeff(self._n)
        print '     A =', self._A
        print '  Reff =', (self._A - 1)/(self._A + 1)
        print '   flc =', self._flc 
        print '   fxc =', self._fxc 
        if solntype is 'fluence':
            print '  type = fluence'
            self.reflectance = self.fluence
        elif solntype is 'kienle':
            print '  type = kienle'
            self.reflectance = self.kienle
        elif solntype is 'flux':
            print '  type = flux'
            self.reflectnace = self.flux

    def kienle(self, mua, musp, rho, wbyv=0, S=1.):
        return self._flc * self.fluence(mua, musp, rho, wbyv=wbyv, S=S) + \
               self._fxc * self.flux(mua, musp, rho, wbyv=wbyv, S=S)
        

    def fluence(self, mua, musp, rho, wbyv=0, S=1.):
        """
        Returns the DOS Green's function for semi-infinte homogeneous medium 
        with extrapolated boundary conditions

        Keyword arguments:
        mua     -- Absorption coefficient
        musp    -- Scattering coefficient
        rho     -- Source-Detector distance
        wbyv    -- angular frequency of RF over speed of light in medium
        S       -- Multiplicative coefficient (default = 1)
        """
        l = 1./(musp+mua)                       # l_str in durduran et al
        D = l/3.                                # D = 1/(3*(mua + musp))
        ik = np.sqrt((mua - 1j*wbyv)/D)         # iK
        r = np.sqrt(np.power(-l, 2) + rho**2)    # isotropic source-detector distance
        zb = 2.*D*self._A                             #  
        rb = np.sqrt(np.power((2.*zb + l), 2) + rho**2)   # source image-detector position 
        return (0.25*S/np.pi/D)*(np.exp(-ik*r)/r - np.exp(-ik*rb)/rb)

    def flux(self, mua, musp, rho, wbyv=0, S=1.):
        """ 
        Returns the DOS Green's function derivative wrt z, times -D, 
        
            3D d\phi / dz; z=0

        for semi-infinite 
        homogeneous medium with extrapolated boundary conditions

        Keyword arguments:
        mua     -- Absorption coefficient
        musp    -- Scattering coefficient
        rho     -- Source-Detector distance
        wbyv    -- angular frequency of RF over speed of light in medium
        S       -- Multiplicative coefficient (default = 1)
        """
        l = 1./(musp+mua)                       # l_str in durduran et al
        D = l/3.                                # D = 1/(3*(mua + musp))
        ik = np.sqrt((mua - 1j*wbyv)/D)         # iK
        r = np.sqrt(np.power(-l, 2) + rho**2)    # isotropic source-detector distance
        zb = 2.*D*self._A                             #  
        rb = np.sqrt(np.power((2.*zb + l), 2) + rho**2)   # source image-detector position 
        return (0.25*3.*S/np.pi)*(l*(ik+1./r)*np.exp(-ik*r)/np.power(r,2) \
                           +(l+2.*zb)*(ik+1./rb)*np.exp(-ik*rb)/np.power(rb,2))

