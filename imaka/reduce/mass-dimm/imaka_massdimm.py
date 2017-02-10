#!/usr/bin/env python
'''############################################################################
## imaka_massdimm.py
##  
## Reads in an imaka fits frame and finds matching FWHM from MASS and DIMM 
## from the MKWC page. Fits file path/name can be given as the argument.
##
## df 2017-02-08. MASS/DIMM reading functions based on massdimm.py by jlu
## last modifield df 2017-02-10
'''############################################################################

import numpy as np
import os
from astropy.io import fits
from datetime import datetime
import asciidata
#import pdb

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

        # tidy this v
		self.hour   = float(datetime[-6:-4])
		self.minute = float(datetime[-4:-2])
		self.second = float(datetime[-2:])

		return np.array([self.hour, self.minute, self.second])


class Table(object):
    ''' Input a FITS table '''
    def __init__ (self, tablfile):
		self.infile = tablfile

		self.data = fits.open(self.infile)[1].data

    def time(self):
        ''' Extracts time from FITS table. Returns array of hh, mm, ss '''
        self.times  = self.data['TIME']

        self.times = self.times.replace(' ',':')
        self.times = self.times.split(':')
        
        ## separate times from am/pm and convert to array
        hhmmss = np.array([item[0:3] for item in self.times]).astype('float')

        ampm   = np.array([item[3] for item in self.times])

        ## add 12 to hour if time is 'PM'
        self.time_24 = hhmmss[:,0] + (ampm[:]=='PM')*12 + hhmmss[:,1]/60. + hhmmss[:,2]/3600.

        return self.time_24

    def fwhm(self):
        self.fwhm  = self.data['FWHM']
        return self.fwhm

    def EE80(self):
        self.ee80  = self.data['EE80']
        return self.ee80
    
    def add_massdimm(self, masscolumn, dimmcolumn):
        ''' Adds new column to table '''
        orig_cols  = self.data.columns
        new_cols   = fits.ColDefs(\
            [fits.Column(name='MASS', format='D', array=masscolumn), \
             fits.Column(name='DIMM', format='D', array=dimmcolumn)])
        hdu = fits.BinTableHDU.from_columns(orig_cols + new_cols)
        hdu.writeto(self.infile+'.md')
        return

class DIMM(object):
    def __init__(self, dimmfile):
        self.file = dimmfile
        table = asciidata.open(dimmfile)

        # Date and time are in UT
        self.year = table[0].tonumpy()
        self.month = table[1].tonumpy()
        self.day = table[2].tonumpy()
    
        self.hour = table[3].tonumpy()
        self.minute = table[4].tonumpy()
        self.second = table[5].tonumpy()

        self.seeing = table[6].tonumpy()

        '''
        # No airmass in new file format
        self.airmass = np.zeros(len(self.hour))

        # Convert from HST to UT
        self.hour += 10
        
        idx = np.where(self.hour > 24)[0]
        self.day[idx] += 1
        self.hour[idx] -= 24
        '''

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
                print 'Could not find DIMM data close to ', item
                closestVal = -1.
            else:
                closestVal = self.seeing[closestIndex]
            closestVals.append(closestVal)

        return closestVals

class MASS(object):

    def __init__(self, massfile):
        self.file = massfile
        # can redo this to not use asciidata module but numpy array instead v
        table = asciidata.open(massfile)
        
        # Date and time are in UT
        self.year = table[0].tonumpy()
        self.month = table[1].tonumpy()
        self.day = table[2].tonumpy()

        self.hour = table[3].tonumpy()
        self.minute = table[4].tonumpy()
        self.second = table[5].tonumpy()

        self.free_seeing = table[6].tonumpy()

        '''
        ## NB all times in HST!!
        # Values Don't exist
        self.isoplanatic_angle = np.zeros(len(self.hour))
        self.tau0 = np.zeros(len(self.hour))

        # Convert from HST to UT
        self.hour += 10
        
        idx = np.where(self.hour > 24)[0]
        self.day[idx] += 1
        self.hour[idx] -= 24
        '''
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
                print 'Could not find MASS data close to ', item
                closestVal = -1.
            else:
                closestVal = self.free_seeing[closestIndex]
            closestVals.append(closestVal)

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
    datadir     = '/Users/dorafohring/Desktop/imaka/data/20170111/'
    filename    = 'aocb_20170111_205.fits'

    datadir     = '/Users/dorafohring/Desktop/imaka/data/20170113/'
    filename    = 'stats_closed1.fits'

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
    
    ## get most recent mass/dimm files. NB for the moment need to ensure massdimm dates match manually!

    #date  = datetime.utcnow().timetuple()
    #td    = '%i%02d%02d' %(date.tm_year, date.tm_mon, date.tm_mday)

    td    = '20170112'

    try:
    	fetch_data(td, massdimmdir)
    except:
    	'Most recent MASS/DIMM data not updated.\n'

    ## get mass/dimm times
    datatime = datafile.time()
    try:
        massvalue = MASS(massdimmdir+td+'.mass.dat').indexTime(datatime)
        dimmvalue = DIMM(massdimmdir+td+'.dimm.dat').indexTime(datatime)
    except:
        print('Error in opening MASS/DIMM files.')

    ## display results
    if isfitsimg:
    	print('MASS:', massvalue, 'DIMM:', dimmvalue)
    else:
        if filename.endswith('.md'):
            pass
        else:
            datafile.add_massdimm(massvalue, dimmvalue)
            #output(td, datafile.r0(), datafile.EE80(), massvalue, dimmvalue)


