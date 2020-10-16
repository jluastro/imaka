## Eden McEwen
## September 25, 2020
# This code takes in an imaka formated aocb file
# Will return auto correlations (acor) or cross correlations (ccor)

import math
import imageio
import time
import os
import sys

import logging
import logging.config

import pandas as pd
import numpy as np
from scipy.io import readsav

import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib as mpl
import matplotlib.animation as animation

from astropy.table import Table
from astropy.stats import sigma_clip
from astropy.io import fits

# Self-witten code
from pipeline.code.corr_code import *
from pipeline.code.corr_plots import *


mask_8_8_center = [[0,0,1,1,1,1,0,0],
           [0,1,1,1,1,1,1,0],
           [1,1,1,1,1,1,1,1],
           [1,1,1,0,0,1,1,1],
           [1,1,1,0,0,1,1,1],
           [1,1,1,1,1,1,1,1],
           [0,1,1,1,1,1,1,0],
           [0,0,1,1,1,1,0,0]]


############################################################## 
####################### Correlator Class #####################
##############################################################
        
class Correlator:
    def __init__(self, name, data_f, out_d, f_file=None, tmax=200, s_sub=False, tt_sub=False):
        # Filepath setups
        self.name = name
        self.out_dir = out_d
        self.data_file = data_f
        self.fits_file = None
        self.date = None
        self.path_valid = self.path_check(out_d, loud=True)
        # Variables in corr generation
        self.s_sub = s_sub
        self.tt_sub = tt_sub
        self.tmax = tmax
        self.n_wfs = self.n_wfs_gen()
        self.active_wfs = [True]*self.n_wfs
        # Target dependent vars
        self.target = None
        self.target_file = None
        self.wfs_ra = [None]*self.n_wfs
        self.wfs_dec = [None]*self.n_wfs
        self.wfs_mag = [None]*self.n_wfs
        # Data inputs
        self.x_slopes = None 
        self.y_slopes = None
        # Set container for autocorr
        self.acor = False
        self.acor_x = np.zeros(self.n_wfs)
        self.acor_y = np.zeros(self.n_wfs)     
        # Set container for crosscorr
        self.ccor = False
        self.ccor_x = np.zeros((self.n_wfs,self.n_wfs))
        self.ccor_y = np.zeros((self.n_wfs,self.n_wfs))
        # if there is a fits file, use it
        self.data_valid = self.slopes_gen()
        if self.fits_check(f_file, loud=True): self.fits_pull()

    ####################### Checking files #######################
    ### These files work with os.path functions to check for a files existence
    ### retunr based on if the request was successful or not
    
    def data_check(self, file, loud=True):
        """
        Tries to find file
            args: file (string)
            return: T/F 
        """
        #REMOVE
        logging.debug("datacheckfile: %s"% file)
        try:
            if os.path.isfile(file): 
                if loud: logging.debug("file found: %s" % file)
                return True
            if loud: logging.debug("data_check! file not found: %s" % file)
            return False
        except Exception as e: 
            logging.error("Error in data_check: %s", e)
            return False

    def path_check(self, path, loud=True):
        """
        Tries to find path, if not found, attempts to create 
            args: path (string)
            return: T/F (whether path is found/created)
        """
        if os.path.isdir(path):
            logging.debug("path found: %s" % path)
            return True
        else:
            try:
                os.makedirs(path)
                logging.debug("path created: %s" % path)
                return True
            except:
                logging.warning("WARNING! Out dir creation error ")
                return False
    
    def fits_check(self, f_file, loud=True):
        '''
        Tries to find either the fits file given, if one stored in class
        Or tries to find the one this function would have been named
            args: path (string)
            return: T/F (whether fits file is found)
        '''
        # Decide on f_file name
        if f_file:
            out_file = f_file
        elif self.fits_file:
            out_file = self.fits_file
        else:
            out_file = self.fits_name_gen()
        #check if we can find the given fits file
        out = os.path.isfile(out_file)
        # if found
        if out:
            self.fits_file = out_file
            if loud: logging.debug("fits found: %s", out_file)
        else:
            if loud: logging.debug("fits not found: %s", out_file)
        return out
    
    
    ####################### Naming files #######################
    ### This defines how we'll save things as we run thruough a pipeline
    ### Very important in terms of defining structure
    
    def fits_name_gen(self):
        """
        Gives back the fits file name associated with parameters
            args: none
            returns: string (fits file)
        """
        out_dir = self.out_dir + "fits/"
        out_file = self.name +"_tmax" + str(self.tmax)
        # generating file structure more finely
        if self.s_sub and self.tt_sub:
            #out_dir = out_dir + "s_tt_sub/"
            out_file = out_file + "_stt"
        elif self.tt_sub:
            #out_dir = out_dir + "tt_sub/"
            out_file = out_file + "_tt"
        elif self.s_sub:
            #out_dir = out_dir + "s_sub/"
            out_file = out_file + "_s"
        else:
            #out_dir = out_dir + "raw/"
            out_file = out_file + "_raw"
        out_file = out_file + ".fits"
        #checking to make sure this location exists
        path = self.path_check(out_dir, loud=False)
        if path: 
            fits_file = out_dir + out_file
            logging.debug("fits name: %s"% fits_file)
            return out_dir + out_file
        else:
            logging.warning("ERROR: fits path not created: %s"% out_dir)
            return out_file

    
    def plot_file_gen(self, plot_dir, plot_type, file_type, med_sub=False, avg_sub=False, avg_len=0):
        """
        Generates a file name for a plot, renumbers if name already taken
            args: plot_type (varies by plot func), file_type ("png" or "gif"), Graphing specifications
            retuns: out_file (descriptive string)
        """
        out_dir = self.out_dir + "plots/" + plot_dir + "/"
        out_file = self.name 
        # generating file structure more finely
        if self.s_sub and self.tt_sub:
            out_file = out_file + "_stt"
        elif self.tt_sub:
            out_file = out_file + "_tt"
        elif self.s_sub:
            out_file = out_file + "_s"
        else:
            out_file = out_file + "_raw"
        out_file = out_file + plot_type
        #check if path exists
        path = self.path_check(out_dir, loud=False)
        # Renumber if repeats 
        out_base = out_dir + out_file
        i = 0
        while os.path.exists(out_base + "_%s.%s" % (i, file_type)):
            i += 1
        out = out_base + "_%s.%s" % (i, file_type)
        return out
    
    
    ####################### Setting Data #######################
    ### Changing variables defined on input
    ### These functions will flush computations, be careful
    ### return True if a change was made, False if value was the same
            
    def set_tmax(self, tmax, loud=True):
        """
        Changes tmax if the requested state is new
        This flushes the old acor and ccor as they no longer reflect the tmax val
            args: tmax (int)
            return: T/F (whether or not tmax val changed)
        """
        if tmax == self.tmax:
            if loud: logging.info("Tmax remains: %s"%tmax)
            return False
        else:
            if loud: logging.info("Tmax now: %s" % tmax)
            self.tmax = tmax
            self.acor = False
            self.ccor = False
            self.fits_file = None
            return True
            
    def set_pre_subs(self, s_sub = -1, tt_sub = -1, loud=True):
        """
        Changes s_sub and or t_sub if the requested state is new
        This flushes the old acor and ccor as they no longer reflect the s_sub state
        Be careful!!!
            args: s_sub (Boolean, T/F) tt_sub (Boolean, T/F)
            return: T/F (whether or not s_sub val changed)
        """
        # Unused inputs
        if s_sub == -1: s_sub = self.s_sub
        if tt_sub == -1: tt_sub = self.tt_sub
        # Seeing who is changing
        change_s = s_sub != self.s_sub
        change_tt = tt_sub != self.tt_sub
        # Confirm changes:
        if change_s:
            if loud: logging.debug("Static Subtraction now: ", s_sub)
            self.s_sub = s_sub
        else:
            if loud: logging.debug("Static Subtraction remains: ", s_sub)
        if change_tt:
            if loud: logging.debug("Tip-tilt Subtraction now: ", tt_sub)
            self.tt_sub = tt_sub
        else:
            if loud: logging.debug("Tip-tilt Subtraction remains: ", tt_sub)
        # Making changes
        if change_tt or change_s:
            self.acor = False
            self.ccor = False
            self.fits_file = None
            if self.fits_check(None, loud=False): self.fits_pull()
            self.data_valid = self.slopes_gen()
            return True
        else: 
            return False
            
    def set_target(self, file):
        """
        Takes in target file and fills in wfs_ra, wfs_dec and wfs_mag
            args: file (in specified format)
            return: T/F (whether or not target was updated)
        """
        #REMOVE
        logging.debug("set_target_file: %s"% file)
        if self.data_check(file):
            self.target_file = file
            with open(file,'r') as i:
                entries = i.readlines()
                entries = [e.replace('\n', '') for e in entries]
            self.target = entries[0]      
            wfs_used = int(entries[1])
            active_wfs = [False]*self.n_wfs
            for w in range(wfs_used):
                wfs_i = entries[w+2].split()
                wfs_index = int(wfs_i[0])
                if wfs_index < self.n_wfs:
                    active_wfs[wfs_index] = True
                    self.wfs_ra[wfs_index] = float(wfs_i[1])
                    self.wfs_dec[wfs_index] = float(wfs_i[2])
                    self.wfs_mag[wfs_index] = float(wfs_i[3])
            self.active_wfs = active_wfs
            return True
        else:
            logging.warning("Target file not found")
            return False

    
    ####################### Generating computations #######################
    ## Generation functions work to define internal class variables
    ## Returns if the call was successful 
    ## Some are automatically calles, and some will need to be requested
    ## Usually required for other functions
    def n_wfs_gen(self):
        try:
            # check to see if data address works
            if self.data_check(self.data_file, loud=False):
                hdulist = fits.open(self.data_file)
                WFS_data = hdulist[3] #the slopes
                return hdulist[0].header['NWFS']
            else:
                logging.warning("No valid data")
                return 5
        except:
            logging.warning("Warning: Default number WFS")
            return 5
    
    def slopes_gen(self):
        """
        If data file exists, pulls x and y slopes from datafile
        If fits is formatted differently, you'll want to change this file
            args: N/A
            return: T/F (whether or not slopes were updated, and thus data valid)
        """
        try:
            # check to see if data address works
            if self.data_check(self.data_file, loud=False):
                hdulist = fits.open(self.data_file)
                WFS_data = hdulist[3] #the slopes
                WFS_list = np.arange(WFS_data.header['NAXIS2'])
                half_slopes = WFS_data.header['NAXIS1']//2
                WFS_shape = (WFS_data.header['NAXIS3'],8,8)
                x_wfs_slopes = np.array([np.array(WFS_data.data[:,wfs, :half_slopes]).reshape(WFS_shape) for wfs in WFS_list])
                y_wfs_slopes = np.array([np.array(WFS_data.data[:,wfs, half_slopes:]).reshape(WFS_shape) for wfs in WFS_list])
                self.n_wfs = hdulist[0].header['NWFS'] 
                self.date = hdulist[0].header['OBSDATE']
                self.x_slopes = x_wfs_slopes 
                self.y_slopes = y_wfs_slopes
                # extra steps if subtracting tiptilt / statics
                if self.tt_sub:
                    self.x_slopes, self.y_slopes = self.get_tiptilt_sub()
                if self.s_sub:
                    self.x_slopes, self.y_slopes = self.get_statics_sub()
                return True
            else:
                logging.warning("No valid data")
                return False
        except:
            return False
    
    def acor_gen(self, mask=mask_8_8_center):
        """
        If acor already calculated, does nothing
        If ccor already calculated, uses that to populate acor
        else, if data is valid, calculates acor
            args: N/A
            return: T/F (whether or not acor is generated)
        """
        if self.acor == True:
            logging.info("Auto corr already generated, tmax = %s" % self.tmax)
            return True
        elif self.ccor == True:
            return self.acor_from_ccor(self, self.ccor_x, self.ccor_y)
        elif self.data_valid:
            t_max = self.tmax
            logging.info("Generating auto corr tmax = %s" % t_max)
            self.acor_x = np.array([td_auto_corr(mat, t_max, mask) for mat in self.x_slopes])
            self.acor_y = np.array([td_auto_corr(mat, t_max, mask) for mat in self.y_slopes])
            self.acor = True
            return True
        else:
            logging.warning("Data invalid, auto corr not generated")
            return False

    def ccor_gen(self, mask=mask_8_8_center):
        """
        If ccor already calculated, does nothing
        else, if data is valid, calculates ccor, using acor vals if available
            args: N/A
            return: T/F (whether or not ccor is generated)
        """
        if self.ccor == True:
            logging.info("Cross corr already generated, tmax = %s" % self.tmax)
            return True
        elif self.data_valid:
            self.ccor = True
            tmax = self.tmax
            logging.info("Generating cross corr tmax = %s" % tmax)
            x_ccor = []
            y_ccor = []
            for i in range(self.n_wfs):
                x_corr_i = []
                y_corr_i = []
                for j in range(self.n_wfs):
                    if i == j and self.acor:
                        x_corr_i.append(self.acor_x[i])
                        y_corr_i.append(self.acor_y[i])
                    elif i > j:
                        x_corr_i.append(x_ccor[j][i])
                        y_corr_i.append(y_ccor[j][i])
                    else:
                        x_cc = td_cross_corr(self.x_slopes[i], self.x_slopes[j], tmax, mask)
                        y_cc = td_cross_corr(self.y_slopes[i], self.y_slopes[j], tmax, mask)
                        x_corr_i.append(x_cc)
                        y_corr_i.append(y_cc)
                    #print("finished cc %s and %s" % (i, j))
                x_ccor.append(x_corr_i)
                y_ccor.append(y_corr_i)
            self.ccor_x = np.array(x_ccor)
            self.ccor_y = np.array(y_ccor)
            if not self.acor:
                self.acor_from_ccor(self.ccor_x, self.ccor_y)
            return True
        else:
            logging.warning("Data invalid, cross corr not generated")
            return False
            
    def acor_from_ccor(self, ccorx, ccory):
        """
        Updates self.acor values from given ccor
            input: ccorx (wfs x wfs array), ccory (wfs x wfs array)
            return: T/F (whether or not function completes)
        """
        try:
            logging.info("Generating auto corr from ccor tmax = %s" % self.tmax)
            self.acor_x = np.array([ci[i] for i, ci in enumerate(ccorx)]) 
            self.acor_y = np.array([ci[i] for i, ci in enumerate(ccory)])
            self.acor = True
            return True
        except Exception as e: 
            logging.error("Error in acor_from_ccor: %s", e)
            return False
    
    
    ####################### Getting Data #######################
    ## Get functions will return you the requested data, used in graphing
    ## usually two returns, one for x slopes, one for y slopes
    ## med_sub = takes out image median
    ## avg_sub = takes out image average
    ## avg_len = total number of images averaged over, taken half from timesteps before t, half from after
        
    def get_statics_sub(self):
        """
        returns slopes with the average static slope subtracted
            input: nothing
            output: x slopes, y slopes
        """
        try:
            x_slopes = np.array([sub_data_avg(mat) for mat in self.x_slopes])
            y_slopes = np.array([sub_data_avg(mat) for mat in self.y_slopes])
            return x_slopes, y_slopes
        except Exception as e: 
            logging.warning("Error with get_statics_sub: %s", e)
        return False, False
    
    def get_tiptilt_sub(self):
        """
        returns slopes with the average tip tilt subtracted
            input: nothing
            output: x slopes, y slopes
        """
        try:
            x_slopes = np.array([sub_data_t_avg(mat) for mat in self.x_slopes])
            y_slopes = np.array([sub_data_t_avg(mat) for mat in self.y_slopes])
            return x_slopes, y_slopes
        except Exception as e: 
            logging.warning("Error with get_tiptilt_sub: %s", e)
        return False, False
    
    
    def data_get_ac(self, med_sub=False, avg_sub=False, avg_len=10):
        """
        Pulls and proccesses all auto corr data for graphing
            args: graphing inputs
            output: x slopes acor n_wfs array, y slopes acor n_wfs array
        """
        if self.data_valid and self.acor:
            data_x = np.array([np.array([mat / mat[7,7] for mat in wfs]) for wfs in self.acor_x])
            data_y = np.array([np.array([mat / mat[7,7] for mat in wfs]) for wfs in self.acor_y])
            if med_sub:
                xcor_med = np.array([running_med(wfs, avg_len) for wfs in data_x])
                ycor_med = np.array([running_med(wfs, avg_len) for wfs in data_y])
                data_x = data_x - xcor_med
                data_y = data_y - ycor_med
            elif avg_sub:
                xcor_avg = np.array([running_avg(wfs, avg_len) for wfs in data_x])
                ycor_avg = np.array([running_avg(wfs, avg_len) for wfs in data_y])
                data_x = data_x - xcor_avg
                data_y = data_y - ycor_avg
            return data_x, data_y 
        else:
            logging.warning("Error with data_get_ac: self.acor is %s, self.data_valid is %s"% (self.acor, self.data_valid))
            return None, None
    
    def data_get_cc(self, wfs_a, wfs_b, med_sub=False, avg_sub=False, avg_len=10):
        """
        Pulls and proccesses wfs_a and wfs_b cc data for graphing
            args: graphing inputs
            output: x slopes ccor wfsa/b, y slopes ccor wfsa/b
        """
        if self.data_valid and self.ccor:
            data_x = np.array([mat / mat[7,7] for mat in self.ccor_x[wfs_a][wfs_b]])
            data_y = np.array([mat / mat[7,7] for mat in self.ccor_y[wfs_a][wfs_b]])
            if med_sub:
                xcor_med = running_med(data_x, avg_len)
                ycor_med = running_med(data_y, avg_len)
                data_x = data_x - xcor_med
                data_y = data_y - ycor_med
            elif avg_sub:
                xcor_avg = running_avg(data_x, avg_len)
                ycor_avg = running_avg(data_y, avg_len)
                data_x = data_x - xcor_avg
                data_y = data_y - ycor_avg
            return data_x, data_y 
        else:
            logging.warning("Error with data_get_cc: self.acor is %s, self.data_valid is %s"% (self.ccor, self.data_valid))
            return None, None
    
    def data_get_cc_all(self, med_sub=False, avg_sub=False, avg_len=10):
        """
        Pulls and proccesses all cc data for graphing
            args: graphing inputs
            output: x slopes ccor n_wfs array, y slopes ccor n_wfs array
        """
        x_ccor = []
        y_ccor = []
        for i in range(self.n_wfs):
            x_corr_i = []
            y_corr_i = []
            for j in range(self.n_wfs):
                if i == j and self.acor:
                    x_acorr, y_acorr = self.data_get_cc(i,i,med_sub,avg_sub,avg_len)
                    x_corr_i.append(x_acorr)
                    y_corr_i.append(y_acorr)
                elif i > j:
                    x_corr_i.append(x_ccor[j][i])
                    y_corr_i.append(y_ccor[j][i])
                else:
                    x_cc, y_cc = self.data_get_cc(i,j,med_sub,avg_sub,avg_len)
                    x_corr_i.append(x_cc)
                    y_corr_i.append(y_cc)   
            x_ccor.append(x_corr_i)
            y_ccor.append(y_corr_i)
        return np.array(x_ccor), np.array(y_ccor)
    
    
    ####################### FITS files #######################
    
    def fits_pull(self, fits_file=None):
        """
        Takes a fits file, checks to see if format, pulls information into object.
            args: fits_file (string path to file)
            returns: Boolean (if pull was successful)
        """
        if fits_file and fits_file != self.fits_file:
            if self.data_check(fits_file):
                self.fits_file = fits_file
                logging.debug("fits_pull: fits file set")
                
        if self.fits_file:
            logging.info("---> Pulling data from fits file :  %s"% self.fits_file)
            hdulist = fits.open(self.fits_file)
            hdr = hdulist[0].header
            self.name = hdr.get('DATANAME')
            self.data_file = hdr.get('DATAFILE')
            self.out_dir = hdr.get('OUTPATH')
            self.target = hdr.get('TARGET')
            self.target_file = hdr.get('TFILE')
            self.n_wfs = hdr.get('NWFS')
            self.tmax = hdr.get('TMAX')
            self.s_sub = hdr.get('SSUB')
            self.tt_sub = hdr.get('TTSUB')
            #Retrieve info from all wfs
            self.ccor = hdr.get("CCOR") 
            self.acor = hdr.get("ACOR")

            x_cor = np.array([])
            y_cor = np.array([])
            for i in range(self.n_wfs):
                hdr_i = hdulist[i+1].header
                self.active_wfs[i] = hdr_i.get('WFSVALID')
                self.wfs_ra[i] = hdr_i.get('RA')
                self.wfs_dec[i] = hdr_i.get('DEC')
                self.wfs_mag[i] = hdr_i.get('MAG')
            if self.ccor or self.acor:
                x_cor = [hdulist[i+1].data[0] for i in range(self.n_wfs)]
                y_cor = [hdulist[i+1].data[1] for i in range(self.n_wfs)]
            if self.ccor:
                self.ccor_x = np.array(x_cor)
                self.ccor_y = np.array(y_cor)
                self.acor_from_ccor(self.ccor_x, self.ccor_y)
            elif self.acor:
                self.acor_x = np.array(x_cor)
                self.acor_y = np.array(y_cor)
            self.data_valid = self.slopes_gen()
            #REMOVE
            logging.debug("self.target_file: %s"% self.target_file)
            if self.target_file:
                self.set_target(self.target_file)
            return True
        else:
            return False
    
    
    def fits_gen(self):
        """
        Generates fits header given data file to pull from
           input: none
           returns: 1/-1 (if data pull was successful/unsuccessful)
        """
        if not self.data_valid:
            logging.warning("Not able to pull data.")
            return -1
        ### define file name
        self.fits_file = self.fits_name_gen()
        out_file = self.fits_file
        ### opening data fits file
        hdulist = fits.open(self.data_file)
        # general breakdown of what's in fits
        loop_state = hdulist[0]
        WFS_cam_raw = hdulist[1]
        WFS_camera_proc = hdulist[2]
        WFS_data = hdulist[3] #the slopes
        ##### set up primary
        hdr = fits.Header()
        hdr['DATANAME'] = (self.name, 'data file name')
        hdr['DATAFILE'] = (self.data_file, 'data file path')
        hdr['OUTPATH'] = (self.out_dir, 'out path')
        hdr['DATETIME'] = (loop_state.header['DATETIME'], loop_state.header.comments['DATETIME'])
        hdr['OBSDATE'] = (loop_state.header['OBSDATE'], loop_state.header.comments['OBSDATE'])
        hdr['TARGET'] = (self.target, "Target field name")
        hdr['TFILE'] = self.target_file
        hdr['NWFS'] = (self.n_wfs, "WFS used in correction")
        hdr['TMAX'] = (self.tmax, 'Max time taken in temporal correlations')
        hdr['SSUB'] = (self.s_sub, "If static slope subtracted before comp")
        hdr['TTSUB'] = (self.tt_sub, "If global tip/tilt subtracted before comp")
        hdr["CCOR"] = (False, "Contains cross correlations")
        hdr["ACOR"] = (False, "Contains only auto correlations")
        primary_hdu = fits.PrimaryHDU(header=hdr)
        hdl_list = [primary_hdu]
        ##### set up headers
        for i in range(self.n_wfs):
            wfs_h = WFS_cam_raw.header
            wfs_c = wfs_h.comments
            hdr = fits.Header()
            hdr['EXTNAME'] = ('WFS' + str(i), "Extension name")
            hdr['WFS'] = (i, "WFS name")
            hdr['TSAMPLE'] = (wfs_h['TSAMPLE'+ str(i)], wfs_c['TSAMPLE'+ str(i)])
            hdr['FSAMPLE'] = (wfs_h['FSAMPLE'+ str(i)], wfs_c['FSAMPLE'+ str(i)])
            hdr['TEXP'] = (wfs_h['TEXP'+ str(i)],  wfs_c['TEXP'+ str(i)])
            hdr['EMGAIN'] = (wfs_h['EMGAIN'+ str(i)], wfs_c['EMGAIN'+ str(i)])
            hdr['WFSVALID'] = (True if hdr['TSAMPLE']!=0 else False, "If WFS useable, based on TSAMPLE")
            hdr['RA'] = (self.wfs_ra[i], "WFS GS right ascension")
            hdr['DEC'] = (self.wfs_dec[i], "WFS GS declination")
            hdr['MAG'] = (self.wfs_mag[i], "WFS GS MAG")
            image_hdu = fits.ImageHDU(header=hdr)
            hdl_list.append(image_hdu)
        ### take list of headers, write to outfile
        hdul = fits.HDUList(hdl_list)
        hdul.writeto(out_file, overwrite=True)
        return 1                         
    
    def fits_write(self):
        """
        Writes data to self.fits_file (ccor, or acor if no ccor)
           input: none
           returns: string (outfile written to)
        """
        ### Set up fits_file if haven't already
        if not self.fits_file:
            self.fits_gen()
        out_file = self.fits_file      
        ### saving correlations
        f = fits.open(out_file, mode='update')
        if self.ccor:
            logging.info("---> saving ccor")
            #save for each wfs
            f[0].header["ccor"] = True
            for i in range(self.n_wfs):
                f[i+1].data = np.array([self.ccor_x[i], self.ccor_y[i]])
        elif self.acor:
            logging.info("---> saving acor")
            f[0].header["acor"] = True
            for i in range(self.n_wfs):
                f[i+1].data = np.array([self.acor_x[i], self.acor_y[i]])
        else:
            logging.warning("---> no cor found")
        f.flush()
        f.close()
        return out_file
    
    ###################### Plotting Helper Functions #############
    ### Consistently set up plot titles and file names
    ### Graphing Specifications 
    ### (optional but avaliable for all graphs in this section)
    ##### med_sub = takes out image median
    ##### avg_sub = takes out image average
    ##### avg_len = total number of images averaged over, taken half from timesteps before t, half from after
    
    
    def plot_title_gen(self, title_type, med_sub, avg_sub, avg_len):
        """
        Generates a Title for a plot
            args: title_type (varies by plot func), Graphing specifications
            retuns: title (nicely spaced string)
        """
        title = self.name + ", " + title_type + "\ntmax="+ str(self.tmax)
        if self.s_sub:
            title = title + ", s_sub"
        if self.tt_sub:
            title = title + ", tt_sub"
        if med_sub:
            title = title + ", med sub len "+ str(avg_len)
        elif avg_sub:
            title = title + ", avg sub len "+ str(avg_len) 
        return title
    
    def check_t_list(self, t_list, max_t, max_len=5):
        """
        Checks a t_list used for a graph, based on what it's going to be used on
            input: t_list (proposed times), max_t (max t length in data), max_len (how many t values desired)
            output: t_list_up
        """
        # Generating new
        if not t_list:
            if max_t <= max_len:
                t_list = np.arange(max_t)
            else:
                # want an even distribution of t's
                cuts = max_t//(max_len-1)
                t_list = [i*cuts for i in range(max_len)]
        # Checking existing
        else:
            if len(t_list) > max_len:
                t_list = tlist[:max_len]
            valid_values = sum([ x < max_t for x in t_list ])
            if valid_values > 2*max_len//3:
                #there are enough valid values to use
                t_list = [t if t < max_t else max_t - 1 for t in t_list]
            else:
                #should just regenerate the list
                t_list = self.check_t_list([], max_t, max_len)
        return t_list

    
    
    ####################### Plotting ACOR #######################
    ### Graphs of the auto correlations
    ### Either a graph, or an animation
    ### Either all WFS, or their average
    
    def acor_graph(self, t_list=[0,5,10,15,20], med_sub=False, avg_sub=False, avg_len=10):
        """
        Create a plot of WFS acor (x, y, and xy cor) t times in t_list 
            args: graph specifications
            outputs: png file location, or false
        """
        # making sure that there is acorr data to work with 
        if self.data_valid and not self.acor:
            logging.warning("Auto corr not available, generating")
            self.acor_gen()
        elif not self.data_valid:
            logging.error("Data not available")
            return None             
        # figure out what to graph (mean_sub or avg_sub)
        data_x, data_y = self.data_get_ac(med_sub, avg_sub, avg_len)
        # masking data, taking average over wfs axis
        x_out = np.average(data_x, axis = 0)
        y_out = np.average(data_y, axis = 0)
        # checking tmax
        t_list = self.check_t_list(t_list, x_out.shape[0])
        # generate names
        title = self.plot_title_gen(" Auto Corr, WFS avg", med_sub, avg_sub, avg_len)
        out_file = self.plot_file_gen("acor_png", "_acor_avg", "png", med_sub, avg_sub, avg_len)
        # graph
        mat_cube_1 = x_out
        mat_cube_2 = y_out
        mat_cube_3 = (x_out + y_out)/2
        label_1, label_2, label_3 = "Sx", "Sy", "Savg"
        try:
            logging.debug("Graphing: graph_3_rows_t_mat")
            graph_3_rows_t_mat(mat_cube_1, mat_cube_2, mat_cube_3, 
                           t_list, title, 
                           label_1, label_2, label_3).savefig(out_file)
            plt.close("all")
            return out_file
        except Exception as e:
            logging.error("Error in graph_3_rows_t_mat: %s"% e)
            return False
          
        
    def acor_animate(self, dt_max = 30, med_sub=False, avg_sub=False, avg_len=10):
        """
        Create an animation length dt_max, of the all WFS (x, y, xy avg)
            args: Graph specifications
            outputs: gif file location, or false
        """
        # making sure that there is acorr data to work with 
        if self.data_valid and not self.acor:
            logging.warning("Auto corr not available, generating...")
            self.acor_gen()
        elif not self.data_valid:
            logging.error("Data not available")
            return None       
        # retrieve data for graphing
        if dt_max > self.tmax: dt_max = self.tmax
        data_x, data_y = self.data_get_ac(med_sub, avg_sub, avg_len)
        # Plot title and file
        title = self.plot_title_gen(" Auto Corr, all WFS", med_sub, avg_sub, avg_len)
        out_file = self.plot_file_gen("acor_all_gif", "_acor_all", "gif", med_sub, avg_sub, avg_len)
        # graph
        label_1, label_2, label_3 = "Sx", "Sy", "Savg"
        try:
            logging.debug("Graphing: gif_3_rows_mat")
            gif_3_rows_mat(data_x, data_y, (data_x + data_y)/2, 
                  dt_max, title, out_file, label_1, label_2, label_3)
            plt.close("all")
            return out_file
        except Exception as e:
            logging.error("Error in gif_3_rows_mat: %s"% e)
            return None
        
        
    def acor_animate_avg(self, dt_max = 30, med_sub=False, avg_sub=False, avg_len=10):
        """
        Create an animation length dt_max, of the average WFS acor
            args: Graph specifications
            returns: gif file location, or false
        """
        # making sure that there is acorr data to work with 
        if self.data_valid and not self.acor:
            logging.warning("Auto corr not available, generating")
            self.acor_gen()
        elif not self.data_valid:
            logging.error("Data not available") 
            return None         
        # retrieve data for graphing
        if dt_max > self.tmax: dt_max = self.tmax
        data_x, data_y = self.data_get_ac(med_sub, avg_sub, avg_len)   
        # averaging over valid wfs
        x_out = np.average(data_x[self.active_wfs], axis=0)
        y_out = np.average(data_y[self.active_wfs], axis=0)     
        # Plot title and file
        title = self.plot_title_gen(" Auto Corr, WFS avg", med_sub, avg_sub, avg_len)
        out_file = self.plot_file_gen("acor_avg_gif", "_acor_avg", "gif", med_sub, avg_sub, avg_len)
        # graph
        label_1, label_2, label_3 = "Sx", "Sy", "Savg"
        try:
            logging.debug("Graphing: gif_3_mat")
            gif_3_mat(x_out, y_out, (x_out + y_out)/2, 
                  dt_max, title, out_file, True, label_1, label_2, label_3)
            plt.close("all")
            return out_file
        except Exception as e:
            logging.error("Error in gif_3_mat: %s"% e)
            return None
    
    ####################### Graphing CCOR #######################
    ### Graphs of the cross correlations
    ### Either a graph, or an animation
    ### Either all WFS, or pic a specific cross correlation
        
    def ccor_graph(self, wfs_a, wfs_b, t_list=[0,5,10,15,20], med_sub=False, avg_sub=False, avg_len=10):
        """
        Static plot of wfs_a and wfs_b, their xy-average slope cross corr, given timeslices in list t
            args: graphing preferences
            returns: string (file location)
        """
        # check to see if data is valid
        if self.data_valid and not self.ccor:
            logging.warning("Auto corr not available, generating")
            self.ccor_gen()
        elif not self.data_valid:
            logging.error("Data not available") 
            return None
        # ordering wfs
        if wfs_a > wfs_b: wfs_a, wfs_b = wfs_b, wfs_a
        # retrieve data for graphing
        data_cx, data_cy = self.data_get_cc(wfs_a, wfs_b, med_sub, avg_sub, avg_len)      
        # checking tmax
        t_list = self.check_t_list(t_list, data_cx.shape[0])
        # Plot title and file
        title = self.plot_title_gen("Cross Corr, WFS"  + str(wfs_a) + " WFS"+ str(wfs_b), med_sub, avg_sub, avg_len)
        out_file = self.plot_file_gen("ccor_png", "_ccor_wfs" + str(wfs_a) + str(wfs_b), "png", med_sub, avg_sub, avg_len)
        # graph
        label_1, label_2, label_3 = "Sx", "Sy", "Savg"
        try:
            logging.debug("Graphing: graph_3_rows_t_mat")
            graph_3_rows_t_mat(data_cx, data_cy, (data_cx + data_cy)/2, 
                           t_list, title, 
                           label_1, label_2, label_3).savefig(out_file)
            plt.close("all")
            return out_file
        except Exception as e:
            logging.error("Error in graph_3_rows_t_mat: %s"% e)
            return False
            
    
    def ccor_graph_all(self, t=0, med_sub=False, avg_sub=False, avg_len=10):
        """
        Graphs all wfs, their xy-average slope cross corr, given timeslice t
            args: graphing preferences
            returns: string (file location)
        """
        # check to see if data is valid
        if self.data_valid and not self.ccor:
            logging.warning("Cross corr not available, generating")
            self.ccor_gen()
        elif not self.data_valid:
            logging.error("Data not available") 
            return None
        # Retrieve data for graphing
        data_cx, data_cy = self.data_get_cc_all(med_sub, avg_sub, avg_len)
        # Plot title and file
        title = self.plot_title_gen(" Cross Corr, all WFS", med_sub, avg_sub, avg_len)
        out_file = self.plot_file_gen("ccor_all_png", "_ccor_all", "png", med_sub, avg_sub, avg_len)       
        # graph
        try:
            logging.debug("Graphing: graph_5_5_ccor")
            graph_5_5_ccor((data_cx+data_cy)/2, t, title).savefig(out_file)
            plt.close("all")
            return out_file
        except Exception as e:
            logging.error("Error in graph_5_5_ccor: %s"% e)
            return False
        
    def ccor_animate(self, wfs_a, wfs_b, dt_max = 30, med_sub=False, avg_sub=False, avg_len=10):
        """
        Static plot of wfs_a and wfs_b, their a, y, and xy-average slope cross corr, given timeslices in list t
            args: graphing preferences
            returns: string (file location)
        """
        # check to see if data is valid
        if self.data_valid and not self.ccor:
            logging.warning("Cross corr not available, generating")
            self.ccor_gen()
        elif not self.data_valid:
            logging.error("Data not available") 
            return None
        # ordering wfs
        if wfs_a > wfs_b: wfs_a, wfs_b = wfs_b, wfs_a
        if dt_max > self.tmax: dt_max = self.tmax
        # Retrieve data for graphing
        data_cx, data_cy = self.data_get_cc(wfs_a, wfs_b, med_sub, avg_sub, avg_len)
        # Plot title and file
        title = self.plot_title_gen(" Cross Corr, WFS"  + str(wfs_a) + " WFS"+ str(wfs_b), med_sub, avg_sub, avg_len)
        out_file = self.plot_file_gen("ccor_gif", "_ccor_wfs" + str(wfs_a) + str(wfs_b), "gif", med_sub, avg_sub, avg_len) 
        # graph
        label_1, label_2, label_3 = "Sx", "Sy", "Savg"
        try:
            logging.debug("Graphing:  gif_3_mat_vmin")
            gif_3_mat_vmin(data_cx, data_cy, (data_cx + data_cy)/2, 
                  dt_max, title, out_file, label_1, label_2, label_3)
            plt.close("all")
            return out_file
        except Exception as e:
            logging.error("Error in gif_3_mat_vmin: %s"% e)
            return False

        
    def cor_animate_all(self, dt_max = 30, med_sub=False, avg_sub=False, avg_len=10):
        # check to see if data is valid
        if self.data_valid and not self.ccor:
            logging.warning("Cross corr not available, generating")
            self.ccor_gen()
        elif not self.data_valid:
            logging.error("Data not available") 
            return None
        # Retrieve data for graphing
        if dt_max > self.tmax: dt_max = self.tmax
        data_cx, data_cy = self.data_get_cc_all(med_sub, avg_sub, avg_len)
        # Plot title and file
        title = self.plot_title_gen(" Cross Corr, all WFS", med_sub, avg_sub, avg_len)
        out_file = self.plot_file_gen("ccor_all_gif", "_ccor_all", "gif", med_sub, avg_sub, avg_len) 
        # animate all correlations 
        try:
            logging.debug("Graphing: cor_animate_all")
            gif_ccor_mat((data_cx+data_cy), dt_max, title, out_file)
            plt.close("all")
            return out_file
        except Exception as e:
            logging.error("Error in gif_ccor_mat: %s"% e)
            return False