#!/usr/bin/env python
'''############################################################################
## imaka_massdimm.py
##  
## Reads in an imaka fits frame and finds matching FWHM from MASS and DIMM 
## from the MKWC page. Fits file path/name can be given as the argument.
##
## df 2017-02-08. MASS/DIMM reading functions based on massdimm.py by jlu
## last modifield df 2017-02-15
## to run on onaga: cat imaka_massdimm.py | ssh imaka@onaga.ifa.hawaii.edu python -
'''############################################################################

import numpy as np
import os
from astropy.io import fits
from pandas import read_csv
from datetime import datetime
import socket
import pdb

def plotRo(dateSuffix):

    imgdir = '/Users/dorafohring/Desktop/plots'

    #dateSuffix = str(year) + str(month).zfill(2) + str(day).zfill(2)

    dimmfile = imgdir + 'results.TMTDIMM.T6-Hawaii.' + dateSuffix
    massfile = imgdir + 'results.TMTMASS.T6-Hawaii.' + dateSuffix
    proffile = imgdir + 'results.TMTMASS_profile.T6-Hawaii.' + dateSuffix

    dimm = DIMM(dimmfile)

    py.clf()
    py.plot(dimm.r0)

class Image(object):
    ''' Input a FITS image '''
    def __init__ (self, imgfile):
        self.infile = imgfile

        self.hdulist = fits.open(self.infile)

    def time(self):
        ''' Extracts time from FITS header'''

        datetime = self.hdulist[0].header['DATETIME']

        self.hour   = float(datetime[-6:-4])
        self.minute = float(datetime[-4:-2])
        self.second = float(datetime[-2:])

        # Convert from HST to UT
        self.hour += 10
        
        if self.hour >= 24:
            self.hour -= 24

        self.timeInHours = self.hour + (self.minute/60.0) + (self.second/3600.0)

        return np.array([self.timeInHours])


class Table(object):
    ''' Input a FITS table '''
    def __init__ (self, tablfile):
        self.infile = tablfile

        self.data = fits.open(self.infile)[1].data

    def time(self):
        ''' Extracts time from FITS table. Returns array of hh, mm, ss '''
        self.times   = self.data['TIME_UTC']

        self.times   = np.array([item for item in self.times.split(':')]).astype('float')
        self.time_24 = self.times[:,0] + self.times[:,1]/60. + self.times[:,2]/3600.

        return self.time_24

    def fwhm(self):
        self.fwhm  = self.data['FWHM']
        return self.fwhm

    def EE80(self):
        self.ee80  = self.data['EE80']
        return self.ee80
    
    def add_mass(self, masscolumn):
        ''' Adds new MASS column to table '''
        orig_cols  = self.data.columns
        new_cols   = fits.ColDefs(\
            [fits.Column(name='MASS', format='D', array=masscolumn)])
        self.data = fits.BinTableHDU.from_columns(orig_cols + new_cols)

        return

    def add_dimm(self, dimmcolumn):
        ''' Adds new DIMM column to table '''
        orig_cols  = self.data.columns
        new_cols   = fits.ColDefs(\
            [fits.Column(name='DIMM', format='D', array=dimmcolumn)])
        self.data = fits.BinTableHDU.from_columns(orig_cols + new_cols)

        return
    
    def add_profile(self, profcols):
        ''' Adds new profile columns to table '''
        orig_cols  = self.data.columns
        pdb.set_trace()
        new_cols   = fits.ColDefs(\
            [fits.Column(name='Cn2dh_05', format='D', array=profcols[0]), \
             fits.Column(name='Cn2dh_1' , format='D', array=profcols[1]), \
             fits.Column(name='Cn2dh_2' , format='D', array=profcols[2]), \
             fits.Column(name='Cn2dh_4' , format='D', array=profcols[3]), \
             fits.Column(name='Cn2dh_8' , format='D', array=profcols[4]), \
             fits.Column(name='Cn2dh_16', format='D', array=profcols[5])  ])
        self.data = fits.BinTableHDU.from_columns(orig_cols + new_cols)

        return
    

    def output(self):
        ''' Write modified table '''
        self.data.writeto(self.infile[:-5] + '_mdp' + '.fits')
        return

class DIMM(object):
    def __init__(self, dimmfile):
        table = read_csv(dimmfile, delim_whitespace=True, names= \
            ['year', 'month', 'day', 'hour', 'minute', 'second', 'seeing'])
        
        # Date and time are in HST
        self.year  = np.array(table['year'])
        self.month = np.array(table['month'])
        self.day   = np.array(table['day'])
    
        self.hour   = np.array(table['hour'])
        self.minute = np.array(table['minute'])
        self.second = np.array(table['second'])

        self.seeing = np.array(table['seeing'])

        # Convert from HST to UT
        self.hour += 10
        
        idx = np.where(self.hour >= 24)[0]
        self.day[idx] += 1
        self.hour[idx] -= 24
        

        self.r0 = 0.98 * 500e-7 * 206265.0 / self.seeing # in cm

        self.timeInHours = self.hour + (self.minute/60.0) + (self.second/3600.0)

    def indexTime(self, hhmmss):
        """
        Fetch the closest row of data for a specified time (UT).
        """
        closestVals = []
        for item in hhmmss:
            timeDiff = abs(self.timeInHours - item)
            closestIndex = timeDiff.argmin()

            ## Make sure we aren't off by more than an hour
            if timeDiff[closestIndex] > 1.0:
                print('Could not find DIMM data close to ', item)
                closestVal = -1.
            else:
                closestVal = self.seeing[closestIndex]
            closestVals.append(closestVal)

        return closestVals

class MASS(object):

    def __init__(self, massfile):
        self.file = massfile

        # redid this to not use asciidata module but pandas instead
        table = read_csv(massfile, delim_whitespace=True, names= \
            ['year', 'month', 'day', 'hour', 'minute', 'second', 'seeing'])
        
        # Date and time are in HST
        self.year  = np.array(table['year'])
        self.month = np.array(table['month'])
        self.day   = np.array(table['day'])
    
        self.hour   = np.array(table['hour'])
        self.minute = np.array(table['minute'])
        self.second = np.array(table['second'])

        self.free_seeing = np.array(table['seeing'])

        # Convert from HST to UT
        self.hour += 10
        
        idx = np.where(self.hour >= 24)[0]
        self.day[idx] += 1
        self.hour[idx] -= 24
        self.timeInHours = self.hour + (self.minute/60.0) + (self.second/3600.0)

    def indexTime(self, hhmmss):
        """
        Fetch the closest row of data for a specified time (UT).
        """

        closestVals = []
        for item in hhmmss:
            timeDiff = abs(self.timeInHours - item)
            closestIndex = timeDiff.argmin()

            ## Make sure we aren't off by more than an hour
            if timeDiff[closestIndex] > 1.0:
                print('Could not find MASS data close to ', item) 
                closestVal = -1.
            else:
                closestVal = self.free_seeing[closestIndex]
            closestVals.append(closestVal)

        return closestVals

class MASSPROF(object):

    def __init__(self, proffile):
        self.file = proffile

        # redid this to not use asciidata module but pandas instead
        table = read_csv(proffile, delim_whitespace=True, names= \
            ['year', 'month', 'day', 'hour', 'minute', 'second', \
            'cn2dh_05', 'cn2dh_1', 'cn2dh_2', 'cn2dh_4', 'cn2dh_8', \
            'cn2dh_16', 'seeing'])
        
        # Date and time are in HST
        self.year  = np.array(table['year'])
        self.month = np.array(table['month'])
        self.day   = np.array(table['day'])
    
        self.hour   = np.array(table['hour'])
        self.minute = np.array(table['minute'])
        self.second = np.array(table['second'])

        self.profs = np.array([table['cn2dh_05'], table['cn2dh_1'], table['cn2dh_2'], \
            table['cn2dh_4'], table['cn2dh_8'], table['cn2dh_16']]).T
        pdb.set_trace()
        # Convert from HST to UT
        self.hour += 10
        
        idx = np.where(self.hour >= 24)[0]
        self.day[idx] += 1
        self.hour[idx] -= 24
        self.timeInHours = self.hour + (self.minute/60.0) + (self.second/3600.0)

    def indexTime(self, hhmmss):
        """
        Fetch the closest row of data for a specified time (UT).
        """

        closestVals = []
        for item in hhmmss:
            timeDiff = abs(self.timeInHours - item)
            closestIndex = timeDiff.argmin()
            ## Make sure we aren't off by more than an hour
            if timeDiff[closestIndex] > 1.0:
                print('Could not find MASS data close to ', item) 
                closestVal = -1.
            else:
                closestVal = self.profs[closestIndex]
            closestVals.append(closestVal)

        closestVals = np.array(closestVals).T

        return closestVals

def fetch_data(utDate, saveTo):
    ''' Saves massdimm files to directory specified.'''
    import urllib

    print('Saving MASS/DIMM data to directory:')
    print(saveTo)

    urlRoot = 'http://mkwc.ifa.hawaii.edu/current/seeing/'
    
    # Save the MASS file
    massFile = utDate + '.mass.dat'
    url = urlRoot + 'mass/' + massFile
    urllib.urlretrieve(url, saveTo + massFile)

    # Save the DIMM file
    dimmFile = utDate + '.dimm.dat'
    url = urlRoot + 'dimm/' + dimmFile
    urllib.urlretrieve(url, saveTo + dimmFile)

    # Save the MASS profile
    massproFile = utDate + '.masspro.dat'
    url = urlRoot + 'masspro/' + massproFile
    urllib.urlretrieve(url, saveTo + massproFile)
    

if __name__ == "__main__":
    ## initialize ##
    realtime = False    #set by hand
    
    ## When running on onaga
    if socket.gethostname() == 'onaga.ifa.hawaii.edu':
        
        if realtime == True:
            ## get most recent mass/dimm files. 
            ## only need to do this when processing data in real time
            date  = datetime.utcnow().timetuple()
            td    = '%i%02d%02d' %(date.tm_year, date.tm_mon, date.tm_mday)

            massdimmdir = '/Volumes/DATA/imaka/' + td + '/mkwc/'

            try:
                fetch_data(td, massdimmdir)
            except:
                'Most recent MASS/DIMM data not updated.\n'

        #datadir     = '/Volumes/DATA/imaka/20170111/ao/'
        #filename    = 'aocb_20170111_205.fits'
        td          = '20170114'  # <-- set this by hand
        filename    = 'stats_ttf2.fits'  # <-- set this by hand
        ## The format bit is needed when the data directory is in HST on 201701 and before
        datadir     = '/Volumes/DATA4/imaka/{0}/fli/reduce/stats/'.format(int(float(td))-1)
        massdimmdir = '/Volumes/DATA4/imaka/{0}/mkwc/'.format(int(float(td))-1)

    else:
        ## When running on local files

        #datadir     = '/Users/dorafohring/Desktop/imaka/data/20170111/'
        #filename    = 'aocb_20170111_205.fits'

        td          = '20170112'  # <-- set this by hand
        filename    = 'stats_closed1.fits'
        datadir     = '/Users/dorafohring/Desktop/imaka/data/{0}/'.format(int(float(td))-1)
        massdimmdir = '/Users/dorafohring/Desktop/imaka/massdimm/'

    file = datadir + filename

    ## optional: can give filename as argument
    if len(os.sys.argv) == 2:
    	file = os.sys.argv[1]

    ## determine if data is image or table based on name
    if filename.startswith('stats'):
    	isfitsimg = False
    	datafile  = Table(file)	
    else:
    	isfitsimg = True
    	datafile  = Image(file)

    ## get mass/dimm times
    datatime = datafile.time()
    print(massdimmdir+ td +'.mass.dat')
    print(datatime)

    massvalue  = MASS(massdimmdir    +td+'.mass.dat')   .indexTime(datatime)
    dimmvalue  = DIMM(massdimmdir    +td+'.dimm.dat')   .indexTime(datatime)
    profvalues = MASSPROF(massdimmdir+td+'.masspro.dat').indexTime(datatime)

    ## display results
    if isfitsimg:
    	print('MASS:', massvalue, 'DIMM:', dimmvalue)
    else:
        if filename.endswith('md.fits'):
            pass
        else:
            datafile.add_dimm(dimmvalue)
            datafile.add_mass(massvalue)
            datafile.add_profile(profvalues)
            datafile.output()
            #output(td, datafile.r0(), datafile.EE80(), massvalue, dimmvalue)

    print('Done!')


