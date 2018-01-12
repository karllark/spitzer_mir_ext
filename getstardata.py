from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

#import string
import math

import numpy as np
from astropy.io import fits as pyfits

def read_line_val(line):
    eqpos = line.find('=')
    if eqpos >= 0:
        colpos = line.find(';')
        if colpos == 0: colpos = len(line)
        return line[eqpos+1:colpos]

# -----------------------------------------------------------
# the object to store the data for photometric band obsservations
class BandData(): 
    # -----------------------------------------------------------
    def __init__(self, type):
        self.type = type
        self.n_bands = 0
        self.bands = {}
        self.band_units = {}
        self.band_waves = {}
        self.band_fluxes = {}
    # -----------------------------------------------------------

    # -----------------------------------------------------------
    # extract the photometric band data from the DAT file
    def read_bands(self, lines):
        for line in lines:
            eqpos = line.find('=')
            pmpos = line.find('+/-')
            if (eqpos >= 0) & (pmpos >= 0):
                # check for reference or unit 
                colpos = max((line.find(';'), line.find('#'), line.find('mJy')))
                if colpos == -1: 
                    colpos = len(line)
                self.bands[line[0:eqpos].strip()] = (float(line[eqpos+1:pmpos]),
                                                     float(line[pmpos+3:colpos]))
                # units
                if line.find('mJy') >= 0:
                    self.band_units[line[0:eqpos].strip()] = 'mJy'
                else:
                    self.band_units[line[0:eqpos].strip()] = 'mag'

        self.n_bands = len(self.bands)
    # -----------------------------------------------------------

    # -----------------------------------------------------------
    # get the information for the supported photometric bands
    def get_poss_bands(self):
        _johnson_band_names = ['U','B','V','R','I',
                               'J','H','K','L','M']
        _johnson_band_waves = np.array([0.366,0.438,0.545,0.641,0.798,
                                        1.22,1.63,2.19,3.45,4.75])
        _johnson_band_zeromag_fluxes = np.array([417.5,632.0,363.1,
                                                 217.7,112.6,31.47,
                                                 11.38,3.961,0.699,
                                                 0.204])*1e-11
        
        _spitzer_band_names = ['IRAC1','IRAC2','IRAC3','IRAC4','IRS15','MIPS24']
        _spitzer_band_waves = np.array([3.52,4.45,5.66,7.67,15.4,23.36])
        _spitzer_band_zeromag_fluxes = np.array([0.650,0.266,0.102,
                                                 0.0304,1.941e-3,
                                                 3.83e-4])*1e-11

        _wfpc2_band_names = ['F170W','F255W','F336W','F439W','F555W','F814W']
        _wfpc2_band_waves = np.array([0.170,0.255,0.336,0.439,0.555,0.814])
        _wfpc2_photflam = np.array([1.551e-15,5.736e-16,5.613e-17,2.945e-17,
                                    3.483e-18,2.508e-18])
        _wfpc2_vegamag = np.array([16.335,17.019,19.429,20.884,22.545,21.639])
        _n_wfpc2_bands = len(_wfpc2_vegamag)
        _wfpc2_band_zeromag_fluxes = np.zeros(_n_wfpc2_bands)

        for i in range(_n_wfpc2_bands):
            _wfpc2_band_zeromag_fluxes[i] = _wfpc2_photflam[i] \
                                            *(10**(0.4*_wfpc2_vegamag[i]))

        _poss_band_names = np.concatenate([_johnson_band_names, 
                                           _spitzer_band_names, 
                                           _wfpc2_band_names])
        _poss_band_waves = np.concatenate([_johnson_band_waves, 
                                           _spitzer_band_waves, 
                                           _wfpc2_band_waves])
        _poss_band_zeromag_fluxes = np.concatenate([_johnson_band_zeromag_fluxes,
                                                    _spitzer_band_zeromag_fluxes,
                                                    _wfpc2_band_zeromag_fluxes])

        # zip everything together into a dictonary to pass back
        return dict(zip(_poss_band_names,
                        zip(_poss_band_zeromag_fluxes,_poss_band_waves)))
    # -----------------------------------------------------------

    # -----------------------------------------------------------
    # get the magnitude for a specific band
    #   supports deriving the magnitude from colors
    def get_band_mag(self, band_name):
        band_data = self.bands.get(band_name)
        if band_data != None:
            return self.bands.get(band_name) + (self.band_units.get(band_name),)
        else:  # check if the band can be made from colors
            _mag = 0.0
            _mag_unc = 0.0
            if band_name == 'U':
                if ((self.bands.get('V') != None) 
                    & (self.bands.get('(B-V)') != None) 
                    & (self.bands.get('(U-B)') != None)):
                    _mag = self.bands['V'][0] + self.bands['(B-V)'][0] \
                           + self.bands['(U-B)'][0]
                    _mag_unc = math.sqrt(self.bands['V'][1]**2 \
                                         + self.bands['(B-V)'][1]**2 \
                                         + self.bands['(U-B)'][1]**2)
            elif band_name == 'B':
                if ((self.bands.get('V') != None) 
                    & (self.bands.get('(B-V)') != None)):
                    _mag = self.bands['V'][0] + self.bands['(B-V)'][0]
                    _mag_unc = math.sqrt(self.bands['V'][1]**2 
                                         + self.bands['(B-V)'][1]**2)
            elif band_name == 'R':
                if ((self.bands.get('V') != None) \
                    & (self.bands.get('(V-R)') != None)):
                    _mag = self.bands['V'][0] - self.bands['(V-R)'][0]
                    _mag_unc = math.sqrt(self.bands['V'][1]**2 
                                         + self.bands['(V-R)'][1]**2)
            elif band_name == 'I':
                if ((self.bands.get('V') != None) 
                    & (self.bands.get('(V-I)') != None)):
                    _mag = self.bands['V'][0] - self.bands['(V-I)'][0]
                    _mag_unc = math.sqrt(self.bands['V'][1]**2 
                                         + self.bands['(V-I)'][1]**2)
                elif ((self.bands.get('V') != None)
                      & (self.bands.get('(V-R)') != None) 
                      & (self.bands.get('(R-I)') != None)):
                    _mag = self.bands['V'][0] - self.bands['(V-R)'][0] \
                      - self.bands['(R-I)'][0]
                    _mag_unc = math.sqrt(self.bands['V'][1]**2 
                                         + self.bands['(V-R)'][1]**2 
                                         + self.bands['(R-I)'][1]**2)
            if (_mag != 0.0) & (_mag_unc != 0.0):
                return (_mag, _mag_unc, 'mag')
    # -----------------------------------------------------------

    # -----------------------------------------------------------
    # get the fluxes in photometric bands using the photometric band data
    def get_band_fluxes(self):
        poss_bands = self.get_poss_bands()

        for pband_name in poss_bands.keys():
            _mag_vals = self.get_band_mag(pband_name)
            if _mag_vals != None:
                if _mag_vals[2] == 'mag':
                    _flux1 = poss_bands[pband_name][0] \
                             *(10**(-0.4*(_mag_vals[0] + _mag_vals[1])))
                    _flux2 = poss_bands[pband_name][0] \
                             *(10**(-0.4*(_mag_vals[0] - _mag_vals[1])))
                    self.band_fluxes[pband_name] = (0.5*(_flux1 + _flux2), 
                                                    0.5*(_flux2 - _flux1))
                    self.band_waves[pband_name] = poss_bands[pband_name][1]
                elif _mag_vals[2] == 'mJy':
                    light_speed = 2.998e-6
                    self.band_waves[pband_name] = poss_bands[pband_name][1]
                    mfac = (1e-3*1e-26*1e7
                            *(light_speed/
                              np.square(self.band_waves[pband_name]))
                            *(1e4)*1e8)
                    self.band_fluxes[pband_name] = (_mag_vals[0]*mfac,
                                                    _mag_vals[1]*mfac)

        # also store the band data in flat numpy vectors for 
        #   computational speed in fitting routines
        # mimics the format of the spectral data
        self.n_flat_bands = len(self.band_waves)
        self.flat_bands_waves = np.zeros(self.n_flat_bands)
        self.flat_bands_fluxes = np.zeros(self.n_flat_bands)
        self.flat_bands_uncs = np.zeros(self.n_flat_bands)

        for k, pband_name in enumerate(self.band_waves.keys()):
            self.flat_bands_waves[k] = self.band_waves[pband_name]
            self.flat_bands_fluxes[k] = self.band_fluxes[pband_name][0]
            self.flat_bands_uncs[k] = self.band_fluxes[pband_name][1]

        # useful for normalization 
        if np.sum(self.flat_bands_uncs) > 0:
            self.ave_band_fluxes = np.average(self.flat_bands_fluxes,
                                              weights=self.flat_bands_uncs**(-2))
        else:
            self.ave_band_fluxes = np.average(self.flat_bands_fluxes)

    # -----------------------------------------------------------
# -----------------------------------------------------------

# -----------------------------------------------------------
# the object to store the data for the spectroscopic observations
#   by this it means wavelength range (e.g. IUE, FUSE, NIR, etc.)
class SpecData(): 
    # -----------------------------------------------------------
    def __init__(self, type):
        self.type = type
        self.n_waves = 0
    # -----------------------------------------------------------

    # -----------------------------------------------------------
    # generic reading of spectra from FITS table
    #  assumes a homoginized format
    def read_spectra(self, line, path='./'):
        eqpos = line.find('=')
        self.file = line[eqpos+2:].rstrip()

        # open and read the spectrum
        datafile = pyfits.open(path + self.file)
        tdata = datafile[1].data  # data are in the 1st extension
        theader = datafile[1].header  # header

        self.wave_range = [theader['wmin'],theader['wmax']]
        self.waves = tdata['wavelength']
        self.flux = tdata['flux']
        self.uncs = tdata['sigma']
        self.npts = tdata['npts']
        self.n_waves = len(self.waves)

        # trim any data that is not finite
        indxs, = np.where(np.isfinite(self.flux) == False)
        if len(indxs) > 0:
            self.npts[indxs] = 0
    # -----------------------------------------------------------

    # -----------------------------------------------------------
    # specific details for IUE spectra
    def read_iue(self, line, path='./'):
        self.read_spectra(line, path)
        # trim the long wavelength data (bad)
        indxs, = np.where(self.waves > 3200.)
        if len(indxs) > 0:
            self.npts[indxs] = 0
        # convert wavelengths from Angstroms to microns (standardization)
        self.waves *= 1e-4

        # exclude regions
        # lya = [8.0, 8.55]
        #ex_regions = [[8.23-0.1,8.23+0.1],
        #              [6.4,6.6], 
        #              [7.1,7.3], 
        #              [7.45,7.55], 
        #              [8.7,10.0]]

        #x = 1.0/self.waves
        #for exreg in ex_regions:
        #    indxs, = np.where((x >= exreg[0]) & (x <= exreg[1]))
        #    if len(indxs) > 0:
        #        self.npts[indxs] = 0

    # -----------------------------------------------------------
    # -----------------------------------------------------------
    # specific details for IRS spectra
    def read_irs(self, line, path='./'):
        self.read_spectra(line, path)

        # convert units of spectra from Jy to ergs/(cm^2 s A)
        #   standardization
        light_speed = 2.998e-6
        mfac = 1e-26*1e7*(light_speed/np.square(self.waves))*(1e4)*1e8
        self.flux *= mfac
        self.uncs *= mfac

        # correct the spectra if correction factors are present
        

# -----------------------------------------------------------

# -----------------------------------------------------------
# the object to store the data for a single star/model
class StarData():
    # -----------------------------------------------------------
    def __init__(self, datfile, path='./'):
        self.file = datfile
        self.path = path
        self.sptype = ''
        self.data = {}
        self.corfac = {}

        # open and read all the lines in the file
        f = open(self.path + self.file, 'r')
        self.datfile_lines = list(f)

        # get the photometric band data
        self.data['BANDS'] = BandData('BANDS')
        self.data['BANDS'].read_bands(self.datfile_lines)

        # covert the photoemtric band data to fluxes in all possible bands
        self.data['BANDS'].get_band_fluxes()

        # extract and store the data in individual ObsData objects
        for line in self.datfile_lines:
            if line.find('sptype') == 0:
                self.sptype = read_line_val(line)
            elif line.find('corfac_irs_zerowave') == 0:
                self.corfac['IRS_zerowave'] = float(read_line_val(line))
            elif line.find('corfac_irs_slope') == 0:
                self.corfac['IRS_slope'] = float(read_line_val(line))
            elif line.find('corfac_irs') == 0:
                self.corfac['IRS'] = float(read_line_val(line))
            elif line.find('IUE') == 0 or line.find('STIS') == 0:
                self.data['STIS'] = SpecData('STIS')
                self.data['STIS'].read_iue(line, path=self.path)
            elif line.find('IRS') == 0 and line.find('IRS15') < 0:
                self.data['IRS'] = SpecData('IRS')
                self.data['IRS'].read_irs(line, path=self.path)

        # correct the spectra if corfacs defined
        if 'IRS' in self.corfac.keys():
            if (('IRS_zerowave' in self.corfac.keys()) 
                and ('IRS_slope' in self.corfac.keys())):
                mod_line = (self.corfac['IRS'] 
                            + (self.corfac['IRS_slope']
                               *(self.data['IRS'].waves 
                                 - self.corfac['IRS_zerowave']) ) )
                self.data['IRS'].flux *= mod_line
                self.data['IRS'].uncs *= mod_line
            else:
                self.data['IRS'].flux *= self.corfac['IRS']
                self.data['IRS'].uncs *= self.corfac['IRS']
                

    # -----------------------------------------------------------
# -----------------------------------------------------------


                
  
