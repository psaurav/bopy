#!/usr/bin/python

import numpy as np
from numpy import inf as Inf
from scipy.integrate import quad

def n2v_cm_sec(n):
    from scipy.constants import c as c
    return c*100./n

def n2v_mm_sec(n):
    from scipy.constants import c as c
    return c*1000./n

def wbyv_mm(nin, fMHz):
    return fMHz*np.pi*1e6/n2v_mm_sec(nin)


def R(n, theta_i):
    """
    Returns Fresnel's Reflection Coefficient 

    Arguments:
    ----------
    n = n_in/n_out
    theta_i = angle of incidence 
    """
    costheta_i   = np.cos(theta_i)
    ntsintheta_i = n*np.sin(theta_i)
    r            = 1
    if ntsintheta_i < 1:
        theta_t    = np.arcsin(ntsintheta_i)
        costheta_t = np.cos(theta_t)
        r1 = (n*costheta_i - costheta_t)/(n*costheta_i + costheta_t)
        r1 = r1**2
        r2 = (costheta_i - n*costheta_t)/(costheta_i + n*costheta_t)
        r2 = r2**2
        r  = (r1+r2)/2.0
    return r

def integrand_rj(theta_i, n):
    return R(n, theta_i)*np.cos(theta_i)*np.cos(theta_i)*np.sin(theta_i)

def integrand_rphi(theta_i, n):
    return R(n, theta_i)*np.cos(theta_i)*np.sin(theta_i)

def integrate_A(n):
    R_J   = 3*quad(integrand_rj,   0, np.pi/2, args=(n))[0]
    R_Phi = 2*quad(integrand_rphi, 0, np.pi/2, args=(n))[0]
    return (1+R_J)/(1-R_Phi)

def integrand_fluencecoeff(theta_i, n):
    return (1.-R(n, theta_i))*np.cos(theta_i)*np.sin(theta_i)

def integrand_fluxcoeff(theta_i, n):
    return (1.-R(n, theta_i))*np.cos(theta_i)*np.cos(theta_i)*np.sin(theta_i)

def fluence_coeff(n):
    return 0.5*quad(integrand_fluencecoeff, 0, np.pi/2., args=(n))[0]

def flux_coeff(n):
    return 3.*0.5*quad(integrand_fluxcoeff, 0, np.pi/2., args=(n))[0]

def fluflx_coeff(n):
    return fluence_coeff(n), flux_coeff(n)

def cfl(n):
    return 0.5*quad(integrand_fluencecoeff, 0, np.pi/2., args=(n))[0]

def cfx(n):
    return 0.5*quad(integrand_fluxcoeff, 0, np.pi/2., args=(n))[0]


def keijzer_A(n):
    Ro = (n-1)/(n+1)
    Ro = Ro*Ro

    costheta_c = 0
    if n > 1:
        theta_c = np.arcsin(1/n)
        costheta_c = np.cos(theta_c)

    costheta_c2 = costheta_c*costheta_c
    costheta_c3 = np.fabs(costheta_c)*costheta_c2
    return (2/(1-Ro) -1 + costheta_c3)/(1-costheta_c2)

def contini_A(n):
    n2 = n*n
    n3 = n*n2
    n4 = n2*n2
    if n < 1:                      # Eq (A2)
        return 3.084635 - 6.531194*n + 8.357854*n2 - 5.082751*n3 + 1.171382*n4
    elif n > 1:                    # Eq (A3)
        n5    = n2*n3
        n6    = n3*n3
        n7    = n3*n4
        return 504.332889 - 2641.00214*n + 5923.699064*n2 - 7376.355814*n3 + 5507.53041*n4 - 2463.357945*n5 + 610.956547*n6 - 64.8047*n7
    else:
        return 1

def egan_Reff(n):
    return -1.440/n/n + 0.710/n + 0.668 + 0.0636*n

def egan_A(n):
    Reff = egan_Reff(n)
    return (1+Reff_egan)/(1-Reff_egan)

def A(n, type='contini'):
    if type is 'contini':
        return contini_A(n)
    elif type is 'egan':
        return egan_A(n)
    elif type is 'keitjzer':
        return keitjzer_A(n)
    elif type is 'exact':
        return integrate_A(n)
    else:
        raise NameError('Type: {0} no defined in A()'.format(type))
    

if __name__ == '__main__':
    min = 0.5
    dn  = 0.01 

    for i in xrange(0,150):
        n = min + i*dn
        A_contini = contini_A(n)
        A_keijzer = keijzer_A(n)
        A_integrt  = integrate_A(n)
        Reff_contini = (A_contini-1)/(A_contini+1)
        Reff_keijzer = (A_keijzer-1)/(A_keijzer+1)
        Reff_integrt = (A_integrt-1)/(A_integrt+1)
        Reff_egan    = egan_Reff(n)
        A_egan       = (1+Reff_egan)/(1-Reff_egan)
        print n, A_contini, A_keijzer, A_integrt, A_egan, Reff_contini, Reff_keijzer, Reff_integrt, Reff_egan
