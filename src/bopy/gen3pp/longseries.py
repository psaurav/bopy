import os
import sys
import re
import numpy as np
from scipy.ndimage.filters import gaussian_filter
import matplotlib.pyplot as plt
from bopy.utils import ParamList, struct

def get_pickoff_value(x, r):
    rm, rx, cm, cx = r
    return np.mean(x[rm:rx, cm:cx].flatten())

def get_pickoff_normalized(pickoff_ddir, fname, pickoff_pixels, norm_index, nts, w, s, sigma):
    a = []
    phi = []
    for i in np.arange(nts)+1:
        # amplitude
        fa = os.path.join(pickoff_ddir, fname.format(i, w, s, 'A'))
        aa = gaussian_filter(np.load(fa), sigma)
        aa = get_pickoff_value(aa, pickoff_pixels)
        a.append(aa)
        #phase
        fp = os.path.join(pickoff_ddir, fname.format(i, w, s, 'phi'))
        pp = gaussian_filter(np.load(fp), sigma)
        pp = get_pickoff_value(pp, pickoff_pixels)
        phi.append(pp)

    # pickoff normalization
    a0 = a[norm_index-1]
    p0 = phi[norm_index-1]

    a = np.array(a)/a0
    phi = np.array(phi) - p0
    return a, phi

def get_src_normalized(ddir, fname, px, norm_index, nts, w, s, sigma):
    a = []
    phi = []
    for i in np.arange(nts)+1:
        fa = os.path.join(ddir, fname.format(i, w, s, 'A'))
        aa = gaussian_filter(np.load(fa), sigma)
        a.append(aa[px])
        fp = os.path.join(ddir, fname.format(i, w, s, 'phi'))
        pp = gaussian_filter(np.load(fp), sigma)
        phi.append(pp[px])

    a0 = a[norm_index-1]
    p0 = phi[norm_index-1]

    a = np.array(a)/a0
    phi = np.array(phi) - p0
    return a, phi

def get_pixels(ddir, fname, truncate, sigma, srcs):
    pixels = {}
    for s in srcs:
        d = np.load(os.path.join(ddir, fname.format(truncate, 785, s, 'A')))
        d = gaussian_filter(d, sigma)
        px = d.sum(1)[:100].argmax(), d.sum(0)[35:].argmax()+35
        pixels[s] = px
    return pixels

def ls_corr(a, phi, ac, phic):
    acorr = a/ac
    pcorr = phi - phic
    return acorr, pcorr

def get_lonseries_images(wls, srcs, nts, norm_index, sigma, ddir, pixels, pickoff_ddir, pickoff_pixels, fname, ofname_a, ofname_p, maxradians, maxpercent, calibname, correction=None):   
    maxp = 1.0 + maxpercent/100.
    minp = 1.0 - maxpercent/100.
    if not correction is None:
        srcs.remove(calibname)
        print >>sys.stderr, srcs

    # metrics file
    if correction is 'calib':
        fcsv = open('calib.csv', 'w')
    elif correction is 'pickoff':
        fcsv = open('pickoff.csv', 'w')

    figa = plt.figure()
    figp = plt.figure()
    for w in wls:
        figa.clf()
        figp.clf()
        axa = figa.add_subplot(111)
        axp = figp.add_subplot(111)
        print >>sys.stderr, "w =", w

        a_calib, phi_calib = get_src_normalized(ddir, fname, pixels[calibname], norm_index, nts, w, calibname, sigma)
        for s in srcs:
            a, phi = get_src_normalized(ddir, fname, pixels[s], norm_index, nts, w, s, sigma)

            # "area" under the curve (AUC)
            auc_a = np.sum(np.abs(a-1.))
            auc_phi = np.sum(np.abs(phi))

            a_pickoff, phi_pickoff = get_pickoff_normalized(pickoff_ddir, fname, pickoff_pixels, norm_index, nts, w, s, sigma)

            # correction
            if correction is 'calib':
                a, phi = ls_corr(a, phi, a_calib, phi_calib)
            elif correction is 'pickoff':
                a, phi = ls_corr(a, phi, a_pickoff, phi_pickoff)
            
            # corrected AUC
            auc_a_corr = np.sum(np.abs(a-1.))
            auc_phi_corr = np.sum(np.abs(phi))

            # the ratio of AUC
            r_auc_a = auc_a_corr/auc_a
            r_auc_phi = auc_phi_corr/auc_phi

            # draw the plots
            if correction is None:
                base_line, = axa.plot(a, label=r'{0}'.format(s))
                axp.plot(phi, color=base_line.get_color(), label=r'{0}'.format(s))
                axa.plot(a_pickoff, '--', color=base_line.get_color(), label=r'{0} '.format(s))
                axp.plot(phi_pickoff, '--', color=base_line.get_color(), label=r'{0} '.format(s))
            else:
                base_line, = axa.plot(a, label=r'{0} {1:1.3f}'.format(s, r_auc_a))
                axp.plot(phi, color=base_line.get_color(), label=r'{0} {1:1.3f}'.format(s, r_auc_phi))
                print >>fcsv, ','.join(map(str, (correction, w, s, r_auc_a, r_auc_phi)))

        axa.legend()
        axa.set_title(r'Amplitude, ($\lambda$={0}nm)'.format(w))
        axa.set_xlabel('Time index')
        axa.set_ylabel(r'A$_i$/A$_0$')
        axa.set_ylim((minp, maxp))
        axp.legend()
        axp.set_title(r'Phase, ($\lambda$={0}nm)'.format(w))
        axp.set_xlabel(r'$i$ [Time index]')
        axp.set_ylabel(r'$\phi - \phi_0$ (rad)')
        axp.set_ylim((-maxradians, maxradians))
        axp2 = axp.twinx()
        x1, x2 = axp.get_xlim()
        axp2.set_xlim(x1, x2)
        y1 = -maxradians*180./np.pi
        y2 = maxradians*180./np.pi
        axp2.set_ylim(y1, y2)
        axp2.set_ylabel(r'$\phi - \phi_0$ (degree)')

        figa.savefig(ofname_a.format(w))
        figp.savefig(ofname_p.format(w))
        
        #plt.close(figa)
        #plt.close(figp)
        
    if not correction is None:
        fcsv.close()


def driver(plfname, correction=None):
    pl = ParamList(plfname)
    nts = int(pl.get_val('LONGSERIES_NTS'))
    sigma = float(pl.get_val('LONGSERIES_SIGMA'))
    ddir = pl.get_val('LONGSERIES_DDIR')
    wls = [int(w) for w in pl.get_val('WAVELENGTHS').strip().split()]
    srcs = [s for s in pl.get_val('LONGSERIES_SRCS').strip().split()]
    norm_index = int(pl.get_val('LONGSERIES_NORM_INDEX'))
    truncate = int(pl.get_val('LONGSERIES_TRUNCATE'))
    fname = pl.get_val('LONGSERIES_FNAME')
    pickoff_ddir = pl.get_val('LONGSERIES_PICKOFF_DDIR')
    pickoff_pixels = [int(p) for p in pl.get_val('LONGSERIES_PICKOFF_PIXELS').strip().split()]
    print >>sys.stderr, pickoff_pixels
    maxradians = float(pl.get_val('LONGSERIES_MAXRADIANS'))
    maxpercent = float(pl.get_val('LONGSERIES_MAXPERCENT'))
    calibname = pl.get_val('LONGSERIES_CALIBNAME')
    
    pixels = get_pixels(ddir, fname, truncate, sigma, srcs)
    if correction is None:
        ofname_a = 'ls-amplitude-{0}.png'
        ofname_p = 'ls-phase-{0}.png'
    elif correction is 'calib':
        ofname_a = 'ls-amplitude-calibcorr-{0}.png'
        ofname_p = 'ls-phase-calibcorr-{0}.png'
    elif correction is 'pickoff':
        ofname_a = 'ls-amplitude-pickoff-{0}.png'
        ofname_p = 'ls-phase-pickoff-{0}.png'

    get_lonseries_images(wls, srcs, nts, norm_index, sigma, ddir, pixels, pickoff_ddir, pickoff_pixels, fname, ofname_a, ofname_p, maxradians, maxpercent, calibname, correction=correction)   

def print_defaults(plfname=None, nts=None, fname=None):
    if plfname is None:
        print 'LONGSERIES_NTS='
        print 'LONGSERIES_SIGMA=8.0'
        print 'LONGSERIES_DDIR=../pp-masked'
        print 'WAVELENGTHS=660 690 785 808 830'
        print 'LONGSERIES_SRCS=s101 s122 s177 s178 s199 calib0'
        print 'LONGSERIES_NORM_INDEX=1'
        print 'LONGSERIES_TRUNCATE=1'
        print 'LONGSERIES_FNAME=20140902_G3093_LS{0}_wl{1}_{2}_{3}.npy'
        print 'LONGSERIES_PICKOFF_DDIR=../pp-pickoff'
        print 'LONGSERIES_PICKOFF_PIXELS=18 28 3 12'
        print 'LONGSERIES_MAXRADIANS=0.1'
        print 'LONGSERIES_MAXPERCENT=15.0'
        print 'LONGSERIES_CALIBNAME=calib0'
    else:
        f = open(plfname, 'w')
        print >>f, 'LONGSERIES_NTS='
        print >>f, 'LONGSERIES_SIGMA=8.0'
        print >>f, 'LONGSERIES_DDIR=../pp-masked'
        print >>f, 'WAVELENGTHS=660 690 785 808 830'
        print >>f, 'LONGSERIES_SRCS=s101 s122 s177 s178 s199 calib0'
        print >>f, 'LONGSERIES_NORM_INDEX=1'
        print >>f, 'LONGSERIES_TRUNCATE=1'
        print >>f, 'LONGSERIES_FNAME=20140902_G3093_LS{0}_wl{1}_{2}_{3}.npy'
        print >>f, 'LONGSERIES_PICKOFF_DDIR=../pp-pickoff'
        print >>f, 'LONGSERIES_PICKOFF_PIXELS=18 28 3 12'
        print >>f, 'LONGSERIES_MAXRADIANS=0.1'
        print >>f, 'LONGSERIES_MAXPERCENT=15.0'
        print >>f, 'LONGSERIES_CALIBNAME=calib0'
        f.close()
    

def get_fitsfiles(ddir, wls):
    import fnmatch

    # get all fits files
    matches = []
    for root, dirnames, filenames in os.walk(ddir):
        for filename in fnmatch.filter(filenames, '*.fits'):
            matches.append(os.path.join(root, filename))

    # get the darks and the calib fits files
    darks = [s for s in matches if 'dark' in os.path.basename(s)]
    calibs = [s for s in matches if 'calib' in os.path.basename(s)]

    # remove the dark and the calibs from the master list
    for s in darks:
        matches.remove(s)
    for s in calibs:
        matches.remove(s)

    # pick the first dark from the list
    darkfits = darks[0]

    # collect the general and calib files according to their 
    # wavelengths and insert the dark at the begging
    srcfits = {} 
    calibfits = {}
    for w in wls:
        wl = 'wl{0}'.format(w)
        srcfits[w] = [s for s in matches if wl in os.path.basename(s)]
        calibfits[w] = [s for s in calibs if wl in os.path.basename(s)]

    # return the src and the calib fits files
    return darkfits, srcfits, calibfits 

def print_fitsdefaults():
    print 'GEN3FIT_SAMPLING_INTERVAL_SEC=0.1'
    print 'GEN3FIT_BEAT_FREQ_HZ=1.0'
    print 'GEN3FIT_N_DATAPOINTS=17'
    print 'GEN3FIT_ABS_ERROR=0.0'
    print 'GEN3FIT_REL_ERROR=1.0e-5'
    print 'GEN3FIT_MAX_ITERATIONS=200'
    print 'GEN3FIT_OUTPUTS=A phi DC'
    print 'GEN3FIT_DROP_FIRST_FRAMES=0'
    print 'GEN3FIT_CALL=aphidc_flat_andro_udrot2'
    print 'GEN3FIT_PICKOFF_MASK_ROWS_UPTO=110'
    print 'CORRECTION_FLIPUD=-1'
    print 'CORRECTION_ROTATION=0'
    print 'WAVELENGTHS=660 690 785 808 830'
    print 'G3ID=G3093'
    print 'ROOT=20140902_G3093_LS'

def driver_fits(ddir, plfname):
    import bopy.utils as butils
    p = butils.struct()
    pl = ParamList(plfname)
    p.rebin = int(pl.get_val('REBIN'))
    p.delt = float(pl.get_val('GEN3FIT_SAMPLING_INTERVAL_SEC'))
    p.freq = float(pl.get_val('GEN3FIT_BEAT_FREQ_HZ'))
    p.nOfData = int (pl.get_val('GEN3FIT_N_DATAPOINTS'))
    p.absErr = float(pl.get_val('GEN3FIT_ABS_ERROR'))
    p.relErr = float(pl.get_val('GEN3FIT_REL_ERROR'))
    p.maxIter = int (pl.get_val('GEN3FIT_MAX_ITERATIONS'))
    p.outputs = [v for v in pl.get_val('GEN3FIT_OUTPUTS').split()]
    wavelengths = [int(v) for v in pl.get_val('WAVELENGTHS').split()]
    p.drop_first_frames = int(pl.get_val('GEN3FIT_DROP_FIRST_FRAMES')) 

def display_rootandnts(darkfit, srcfit):
    root = os.path.basename(darkfit)
    root = re.sub('_dark.*', '', root)
    root = re.sub('\d+$', '', root)
    print 'LONGSERIES_FNAME='+root+'{0}_wl{1}_{2}_{3}.npy'
    nts = map(os.path.basename, srcfit)
    pattern = r'({0})(\d+)'.format(root)
    p = re.compile(pattern)
    nts = np.max(np.array([int(p.match(s).group(2)) for s in nts])) - 1
    print ''.join(('LONGSERIES_NTS=', str(nts)))

def parse_cfgfile_fit(cfgfile):
    pl = ParamList(cfgfile)
    p = struct()

    p.wls = [int(v) for v in pl.get_val('WAVELENGTHS').split()]
    p.rebin = int(pl.get_val('REBIN'))
    p.delt = float(pl.get_val('GEN3FIT_SAMPLING_INTERVAL_SEC'))
    p.freq = float(pl.get_val('GEN3FIT_BEAT_FREQ_HZ'))
    p.nOfData = int (pl.get_val('GEN3FIT_N_DATAPOINTS'))
    p.absErr = float(pl.get_val('GEN3FIT_ABS_ERROR'))
    p.relErr = float(pl.get_val('GEN3FIT_REL_ERROR'))
    p.maxIter = int (pl.get_val('GEN3FIT_MAX_ITERATIONS'))
    p.outputs = [v for v in pl.get_val('GEN3FIT_OUTPUTS').split()]
    p.drop_first_frames = int(pl.get_val('GEN3FIT_DROP_FIRST_FRAMES'))
    p.pickoff_mask_rows_upto = int(pl.get_val('GEN3FIT_PICKOFF_MASK_ROWS_UPTO'))
    p.pickoff_mask_cols_from = int(pl.get_val('GEN3FIT_PICKOFF_MASK_COLS_FROM'))
    p.pickoff = [int(v) for v in pl.get_val('GEN3FIT_PICKOFF').split()]

    return p

def fitit((darkfile, filelist, cfgfile, p)):
    #print darkfile, filelist
    import bopy.gen3fit as bg3f
    bg3f.aphidc_flat_andro_udrot2(filelist, darkfile, cfgfile, p)
    return 0

def lsfit(ddir, cfgfile, w, kind):
    p = parse_cfgfile_fit(cfgfile)
    wls = p.wls
    darkfits, srcfits, calibfits = get_fitsfiles(ddir, wls)
    w = int(w)
    import bopy.gen3fit as bg3f
    print kind
    if kind == 'src':
        print 'in src'
        p.pickoff_mask_cols_from = None
        bg3f.aphidc_flat_andro_udrot2(srcfits[w], darkfits, cfgfile, p)
    elif kind == 'calib':
        print 'in calib'
        p.pickoff_mask_rows_upto = None
        bg3f.aphidc_flat_andro_udrot2(calibfits[w], darkfits, cfgfile, p)
    elif kind == 'pickoff':
        print 'in pickoff'
        p.pickoff_mask_cols_from = None
        bg3f.aphidc_flat_andro_udrot2_pickoff(srcfits[w], darkfits, cfgfile, p)
        bg3f.aphidc_flat_andro_udrot2_pickoff(calibfits[w], darkfits, cfgfile, p)
    return 0
    
