#!/usr/bin/python

import os
import sys
from datetime import date

from traits.api import *
from traitsui.api import *
from traitsui.menu import OKButton, CancelButton

import bopy as bp

class G3Study(HasTraits):

    SOURCE_NCOLUMNS = CInt(19,
        desc='Number of Columns (Source)',
        label='Number of Columns (Source)')
        
    SOURCE_NROWS= CInt(11,
        desc='Number of Rows (Source)',
        label='Number of Rows (Source)')
    SOURCE_SPACING=CInt(8,
        desc='Spacing between sources (mm)',
        label='Spacing between sources (mm)') 
     
    DATA_DIRECTORY = Directory('',
        desc='Data directory',
        label='Data directory')

    STUDY_DIRECTORY = Directory('',
        desc='Study directory',
        label='Study directory')

    MEASUREMENT_REF = Directory('',
        desc='Measurement REF',
        label='Measurement REF')

    MEASUREMENT_TRG = Directory('',
        desc='Measurement TRG',
        label='Measurement TRG')

    #SOURCE_PLATE_IMAGE_FILE = File('',
    #    desc='Source Plate Image file',
    #    label='Source Plate Image file')

    INK_SPECTROPHOTOMETER_FILE = String('',
        desc='Ink Absorbance file',
        label='Ink Absorbance file')

    #SOURCE_PLATE_IMAGE_FILE = File('',
    #    desc='Source Plate file',
    #    label='Source Plate file')

    SOURCES = String('',
        desc='All source indexes',
        label='All source indexes')

    #BREAST_IMAGE_FILE = File('',
    #    desc='Breast image file',
    #    label='Breast image file')

    INIT_CONFIG_FILE = File('',
        desc='Init config file',
        label='Init config file')

    WAVELENGTHS = String('660 690 785 808 830',
        desc='Wavelengths (nm)',
        label='Wavelengths (nm)')

    INK_SPECT_DILUTION_PERCENTAGE = Float(0.0,
        desc='Ink spect dilution percentage',
        label='Ink spect dilution percentage')

    SAMPLING_RATE = Float(0.1,
        desc='Sampling rate (sec)',
        label='Sampling rate (sec)')

    BREAST_TANK_THICKNESS_MM = Float(60,
        desc='Breast tank thickness (mm)',
        label='Breast tank thickness (mm)')

    RF_FREQUENCY_MHZ = Float(70,
        desc='RF Frequency (MHz)',
        label='RF Frequency (MHz)')

    HOMOG_WATER_VOLUME_ORIGINAL_ML = Float(0.0,
        desc='Homog Water volume original (ml)',
        label='Homog Water volume original (ml)')

    HOMOG_WATER_VOLUME_REMOVED_ML = Float(0.0,
        desc='Homog Water volume removed (ml)',
        label='Homog Water volume removed (ml)')

    HOMOG_INTRALIP_VOLUME_ADDED_ML = Float(0.0,
        desc='Homog Intralipid volume added (ml)',
        label='Homog Intrlipid volume added (ml)')

    HOMOG_INK_VOLUME_ADDED_ML = Float(0.0,
        desc='Homog Ink volume added (ml)',
        label='Homog Ink volume added (ml)')

    def _DATA_DIRECTORY_changed(self):
        self.STUDY_DIRECTORY = self.DATA_DIRECTORY

    def _STUDY_DIRECTORY_changed(self):
        self.MEASUREMENT_REF = self.STUDY_DIRECTORY
        self.MEASUREMENT_TRG = self.STUDY_DIRECTORY
        self.INK_ABSORBANCE_FILE = self.STUDY_DIRECTORY
        #self.SOURCE_PLATE_IMAGE_FILE = self.STUDY_DIRECTORY
        #self.BREAST_IMAGE_FILE = self.STUDY_DIRECTORY
        self.INIT_CONFIG_FILE = self.STUDY_DIRECTORY

    def _MEASUREMENT_REF_changed(self):
        refname = os.path.basename(self.MEASUREMENT_REF)
        fname = os.path.join(self.MEASUREMENT_REF, '%s.ini'%(refname))
        if os.path.exists(fname):
            pl = bp.utils.ParamList(fname)
            self.SOURCES = pl.get_val('Source_SourceIndex_NONE_String').replace(',', ' ')
            try:
                self.WAVELENGTHS = pl.get_val('Switch_Wavelength_NONE_String').replace(',', ' ')
            except:
                self.WAVELENGTHS = pl.get_val('Switch6_Wavelength_NONE_String').replace(',', ' ')

            print >>sys.stderr, 'unsorted wavelengths:',  self.WAVELENGTHS
            self.WAVELENGTHS = ' '.join(sorted([v for v in self.WAVELENGTHS.strip().split()]))
            print >>sys.stderr, '  sorted wavelengths:',  self.WAVELENGTHS
        

    def _INIT_CONFIG_FILE_changed(self):
        fpath = self.INIT_CONFIG_FILE
        if os.path.isdir(fpath):
            return
        if os.path.exists(fpath):
            pl = bp.utils.ParamList(fpath)
            self.INK_SPECTROPHOTOMETER_FILE = pl.get_val('INK_SPECTROPHOTOMETER_FILE')
            self.INK_SPECT_DILUTION_PERCENTAGE = float(pl.get_val('INK_SPECT_DILUTION_PERCENTAGE'))
            self.HOMOG_WATER_VOLUME_ORIGINAL_ML = float(pl.get_val('HOMOG_WATER_VOLUME_ORIGINAL_ML'))
            self.HOMOG_WATER_VOLUME_REMOVED_ML = float(pl.get_val('HOMOG_WATER_VOLUME_REMOVED_ML'))
            self.HOMOG_INTRALIP_VOLUME_ADDED_ML = float(pl.get_val('HOMOG_INTRALIP_VOLUME_ADDED_ML'))
            self.HOMOG_INK_VOLUME_ADDED_ML = float(pl.get_val('HOMOG_INK_VOLUME_ADDED_ML'))
            self.BREAST_TANK_THICKNESS_MM = float(pl.get_val('BREAST_TANK_THICKNESS_MM'))
            self.RF_FREQUENCY_MHZ = float(pl.get_val('RF_FREQUENCY_MHZ'))

    def print_values(self):
        filename = os.path.basename(self.STUDY_DIRECTORY) + '-init.cfg'
        f = open(filename, 'w')
        #for attr, value in sorted(self.__dict__.iteritems()):
        #    print >>f, attr, '=', value
        
        
        print >>f, 'DATA_DIRECTORY =', self.DATA_DIRECTORY
        print >>f, 'STUDY_NAME =', os.path.basename(self.STUDY_DIRECTORY)
        print >>f, 'REF_NAME =', os.path.basename(self.MEASUREMENT_REF)
        print >>f, 'TRG_NAME =', os.path.basename(self.MEASUREMENT_TRG)
        print >>f, 'SOURCE_PLATE_NAME = SourcePlate'
        print >>f, 'SOURCE_PLATE_FILE =', os.path.basename(self.STUDY_DIRECTORY)+'_SourcePlate_picture.fits'
        print >>f, 'BREAST_IMAGE_NAME = BreastPicture'
        print >>f, 'BREAST_IMAGE_FILE =', os.path.basename(self.STUDY_DIRECTORY)+'_BreastPicture_picture.fits'
        print '*** DEBUG', os.path.basename(self.STUDY_DIRECTORY)
        #root, g3id, t = [v for v in os.path.basename(self.STUDY_DIRECTORY).split('_')]
        strs = [v for v in os.path.basename(self.STUDY_DIRECTORY).split('_')]
        root = os.path.basename(self.STUDY_DIRECTORY)
        date = strs[0]
        g3id = strs[1]
        t = '_'.join([str(v) for v in strs[2:]])
        print >>f, 'ROOT =', root
        print >>f, 'ROOT_DATE =', date
        print >>f, 'G3ID =', g3id
        print >>f, 'TYPE =', t
        print >>f, 'SOURCES =', self.SOURCES
        print >>f, 'WAVELENGTHS =', self.WAVELENGTHS
        print >>f, 'INK_SPECT_DILUTION_PERCENTAGE =', self.INK_SPECT_DILUTION_PERCENTAGE
        print >>f, 'INK_SPECTROPHOTOMETER_FILE =', self.INK_SPECTROPHOTOMETER_FILE
        print >>f, 'INK_SPECT_DILUTION_PERCENTAGE =', self.INK_SPECT_DILUTION_PERCENTAGE
        print >>f, 'HOMOG_WATER_VOLUME_ORIGINAL_ML =', self.HOMOG_WATER_VOLUME_ORIGINAL_ML
        print >>f, 'HOMOG_WATER_VOLUME_REMOVED_ML =', self.HOMOG_WATER_VOLUME_REMOVED_ML
        print >>f, 'HOMOG_INTRALIP_VOLUME_ADDED_ML =', self.HOMOG_INTRALIP_VOLUME_ADDED_ML
        print >>f, 'HOMOG_INK_VOLUME_ADDED_ML =', self.HOMOG_INK_VOLUME_ADDED_ML
        print >>f, 'SAMPLING_RATE =', self.SAMPLING_RATE
        print >>f, 'BREAST_TANK_THICKNESS_MM =', self.BREAST_TANK_THICKNESS_MM
        print >>f, 'RF_FREQUENCY_MHZ =', self.RF_FREQUENCY_MHZ
        print >>f, 'SOURCE_NCOLUMNS =', self.SOURCE_NCOLUMNS
        print >>f, 'SOURCE_NROWS =', self.SOURCE_NROWS
        print >>f, 'SOURCE_SPACING =', self.SOURCE_SPACING
        f.close()
        return filename

    def read_values(self, filename):
        print >>sys.stderr, "========================"
        print >>sys.stderr, sys.argv[0], ":  reading values"
        print >>sys.stderr, "========================"
        f = open(filename, 'r')
        for line in f:
            line = line.strip()
            print >>sys.stderr, line
            k, v = line.split('=')
            v = v.replace(' ', '')
            k = k.replace(' ', '')
            if v:
                print >>sys.stderr, k, v, "'self.'+k+'='+v"
                d = 'isinstance('+k+', Date)'
                isDate = False
                try:
                    isDate = eval(d)
                except:
                    print >>sys.stderr, "***eval failed to test"
                    print >>sys.stderr, "*** d =", d
                if isDate:
                    print >>sys.stderr, "*** isDate: "
                try:
                    exec 'self.'+k+'='+v
                except:
                    print >>sys.stderr, "***exec failed to assign value"
            else:
                print >>sys.stderr, k
        f.close()

study_info = Group(Item(name='DATA_DIRECTORY'),
    Item(name='STUDY_DIRECTORY'),
    Item(name='MEASUREMENT_REF'),
    Item(name='MEASUREMENT_TRG'),
    #Item(name='SOURCE_PLATE_IMAGE_FILE'),
    #Item(name='BREAST_IMAGE_FILE'),
    Item(name='INIT_CONFIG_FILE'),
    Item(name='SOURCES'),
    Item(name='SOURCE_NCOLUMNS'),
    Item(name='SOURCE_NROWS'),
    Item(name='SOURCE_SPACING'),
    Item(name='WAVELENGTHS'),
    Item(name='INK_SPECT_DILUTION_PERCENTAGE'),
    Item(name='INK_SPECTROPHOTOMETER_FILE'),
    Item(name='INK_SPECT_DILUTION_PERCENTAGE'),
    Item(name='HOMOG_WATER_VOLUME_ORIGINAL_ML'),
    Item(name='HOMOG_WATER_VOLUME_REMOVED_ML'),
    Item(name='HOMOG_INTRALIP_VOLUME_ADDED_ML'),
    Item(name='HOMOG_INK_VOLUME_ADDED_ML'),
    Item(name='SAMPLING_RATE'),
    Item(name='BREAST_TANK_THICKNESS_MM'),
    Item(name='RF_FREQUENCY_MHZ'),
    label = 'G3 Data',
    dock='tab'
    )

#ccd_info = Group(Item(name='SOURCE_PLATE_IMAGE_FILE'),
#    label = 'CCD',
#    dock = 'tab',
#    #buttons = [OKButton]
#    )

viewtabbed = View(study_info,
    #ccd_info, 
    #layout='tabbed',
    buttons = [OKButton, CancelButton],
    title = 'G3 Configuration',
    width=700,
    height=250,
    resizable=True,
    scrollable=True
    )

def run_g3study(filename=None):

    study = G3Study()
    if not filename == None:
        study.read_values(filename)
    study.configure_traits(view=viewtabbed)

    ofilename = study.print_values()
    homog = bp.spectra.HomogMuaMusp(ofilename)
    homog.get_mua()
    homog.get_musp()
    homog.output_values()

if __name__ == "__main__":
    
    try:
        filename = sys.argv[1]
        print sys.argv[0], filename
        run_g3study(filename)
    except:
        run_g3study()
