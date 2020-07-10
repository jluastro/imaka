# Eden McEwen
# June 18th 2020
# Major edit July 10th
# removed target support
# deleted set_target
# A class for holding all the functions for processing data in the wind profiling pipeline

import math
import imageio
import time
import os
import sys

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
from corr_code import *
from graph_code import *


mask_8_8_center = [[0,0,1,1,1,1,0,0],
           [0,1,1,1,1,1,1,0],
           [1,1,1,1,1,1,1,1],
           [1,1,1,0,0,1,1,1],
           [1,1,1,0,0,1,1,1],
           [1,1,1,1,1,1,1,1],
           [0,1,1,1,1,1,1,0],
           [0,0,1,1,1,1,0,0]]

class DataPipe:
    def __init__(self, name, data_f, out_d, tmax=200, s_sub=False, tt_sub=False):
        self.name = name
        self.out_dir = out_d
        self.data_file = data_f
        self.fits_file = None
        self.path_valid = self.path_check(out_d)
        self.n_wfs = 5
        # TODO: modularize this?
        self.s_sub = s_sub
        self.tt_sub = tt_sub
        self.x_slopes = None 
        self.y_slopes = None
        self.data_valid = self.slopes_gen()
        self.tmax = tmax
        # Set container for autocorr
        self.acor = False
        self.acor_x = np.zeros(self.n_wfs)
        self.acor_y = np.zeros(self.n_wfs)     
        # Set container for crosscorr
        self.ccor = False
        self.ccor_x = np.zeros((self.n_wfs,self.n_wfs))
        self.ccor_y = np.zeros((self.n_wfs,self.n_wfs))
        # if there is a fits file, use it
        if self.fits_check(): self.fits_pull()

    ####################### Checking files #######################
    ### These files work with os.path functions to check for a files existence
    ### retunr based on if the request was successful or not
    
    def data_check(self, file):
        """
        Tries to find file
            args: file (string)
            return: T/F 
        """
        if os.path.isfile(file): 
            #print("file found: %s" % file)
            return True
        print("WARNING! file not found: %s" % file)
        return False 

    def path_check(self, path):
        """
        Tries to find path, if not found, attempts to create 
            args: path (string)
            return: T/F (whether path is found/created)
        """
        if os.path.isdir(path):
            print("path found: %s" % path)
            return True
        else:
            try:
                os.mkdir(path)
                print("path created: %s" % path)
                return True
            except:
                print("WARNING! Out dir creation error")
                return False
    
    def fits_check(self):
        """
        Tries to find either the fits file given, if one stored in class
        Or tries to find the one this function would have been named
            args: path (string)
            return: T/F (whether fits file is found)
        """
        if self.fits_file:
            out_file = self.fits_file
        else:
            out_file = self.fits_name_gen()
        out = os.path.isfile(out_file)
        if out:
            self.fits_file = out_file
            print("fits found: ", out_file)
        else:
            print("fits not found: ", out_file)
        return out
    
    ####################### Setting Data #######################
    ### Changing variables defined on input
    ### These functions will flush computations, be careful
    ### return True if a change was made, False if value was the same
            
    def set_tmax(self, tmax):
        """
        Changes tmax if the requested state is new
        This flushes the old acor and ccor as they no longer reflect the tmax val
            args: tmax (int)
            return: T/F (whether or not tmax val changed)
        """
        if tmax == self.tmax:
            print("Tmax remains: %s"%tmax)
            return False
        else:
            print("Tmax now: %s" % tmax)
            self.tmax = tmax
            self.acor = False
            self.ccor = False
            self.fits_file = None
            return True
            
    def set_ssub(self, s_sub):
        """
        Changes s_sub if the requested state is new
        This flushes the old acor and ccor as they no longer reflect the s_sub state
        Be careful!!!
            args: s_sub (Boolean, T/F)
            return: T/F (whether or not s_sub val changed)
        """
        if s_sub == self.s_sub:
            print("Static Subtraction remains: ", s_sub)
            return False
        else: 
            print("Static Subtraction now: ", s_sub)
            self.s_sub = s_sub
            self.acor = False
            self.ccor = False
            self.fits_file = None
            self.data_valid = self.slopes_gen()
            return True
        
    def set_ttsub(self, tt_sub):
        """
        Changes tt_sub if the requested state is new
        This flushes the old acor and ccor as they no longer reflect the s_sub state
        Be careful!!!
            args: tt_sub (Boolean, T/F)
            return: T/F (whether or not tt_sub val changed)
        """
        if tt_sub == self.tt_sub:
            print("Tip-tilt Subtraction remains: ", tt_sub)
            return False
        else: 
            print("Tip-tilt Subtraction now: ", tt_sub)
            self.tt_sub = tt_sub
            self.acor = False
            self.ccor = False
            self.fits_file = None
            self.data_valid = self.slopes_gen()
            return True
            
    def set_target(self, file):
        """
        Takes in target file and fills in wfs_ra, wfs_dec and wfs_mag
            args: file (in specified format)
            return: T/F (whether or not target was updated)
        """
        print("set_target depricated. Nothing updated")
        return False

    
    ####################### Generating computations #######################
    ## Generation functions work to define internal class variables
    ## Returns if the call was successful 
    ## Some are automatically calles, and some will need to be requested
    ## Usually required for other functions
    
    def slopes_gen(self):
        """
        If data file exists, pulls x and y slopes from datafile
        If fits is formatted differently, you'll want to change this file
            args: N/A
            return: T/F (whether or not slopes were updated, and thus data valid)
        """
        try:
            # check to see if data address works
            if self.data_check(self.data_file):
                hdulist = fits.open(self.data_file)
                WFS_data = hdulist[3] #the slopes
                WFS_list = np.arange(WFS_data.header['NAXIS2'])
                half_slopes = WFS_data.header['NAXIS1']//2
                WFS_shape = (WFS_data.header['NAXIS3'],8,8)
                x_wfs_slopes = np.array([np.array(WFS_data.data[:,wfs, :half_slopes]).reshape(WFS_shape) for wfs in WFS_list])
                y_wfs_slopes = np.array([np.array(WFS_data.data[:,wfs, half_slopes:]).reshape(WFS_shape) for wfs in WFS_list])
                self.n_wfs = hdulist[0].header['NWFS'] 
                self.x_slopes = x_wfs_slopes 
                self.y_slopes = y_wfs_slopes
                # extra steps if subtracting tiptilt / statics
                if self.tt_sub:
                    self.x_slopes, self.y_slopes = self.get_tiptilt_sub()
                if self.s_sub:
                    self.x_slopes, self.y_slopes = self.get_statics_sub()
                return True
            else: 
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
            print("Auto corr already generated, tmax = %s" % self.tmax)
            return True
        elif self.ccor == True:
            return self.acor_from_ccor(self, self.ccor_x, self.ccor_y)
        elif self.data_valid:
            t_max = self.tmax
            print("Generating auto corr tmax = %s" % t_max)
            self.acor_x = np.array([td_auto_corr(mat, t_max, mask) for mat in self.x_slopes])
            self.acor_y = np.array([td_auto_corr(mat, t_max, mask) for mat in self.y_slopes])
            self.acor = True
            return True
        else:
            print("Data invalid")
            return False

    def ccor_gen(self, mask=mask_8_8_center):
        """
        If ccor already calculated, does nothing
        else, if data is valid, calculates ccor, using acor vals if available
            args: N/A
            return: T/F (whether or not ccor is generated)
        """
        if self.ccor == True:
            print("Cross corr already generated, tmax = %s" % self.tmax)
            return True
        elif self.data_valid:
            self.ccor = True
            tmax = self.tmax
            print("Generating cross corr tmax = %s" % tmax)
            x_ccor = []
            y_ccor = []
            for i in range(len(self.active_wfs)):
                x_corr_i = []
                y_corr_i = []
                for j in range(len(self.active_wfs)):
                    if i == j and self.acor:
                        x_corr_i.append(self.acor_x[i])
                        y_corr_i.append(self.acor_y[i])
                    if i > j:
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
            print("Data invalid")
            return False
            
    def acor_from_ccor(self, ccorx, ccory):
        """
        Updates self.acor values from given ccor
            input: ccorx (wfs x wfs array), ccory (wfs x wfs array)
            return: T/F (whether or not function completes)
        """
        try:
            print("Generating auto corr from ccor tmax = %s" % self.tmax)
            self.acor_x = np.array([ci[i] for i, ci in enumerate(ccorx)]) 
            self.acor_y = np.array([ci[i] for i, ci in enumerate(ccory)])
            self.acor = True
            return True
        except: 
            print("Error in acor_from_ccor")
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
        except:
            print("Error with get_statics_sub")
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
        except:
            print("Error with get_tiptilt_sub")
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
            print("Error with data_get_ac: self.acor is %s, self.data_valid is %s"% (self.acor, self.data_valid))
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
            print("Error with data_get_cc: self.acor is %s, self.data_valid is %s"% (self.ccor, self.data_valid))
            return None, None
    
    def data_get_cc_all(self, med_sub=False, avg_sub=False, avg_len=10):
        """
        Pulls and proccesses all cc data for graphing
            args: graphing inputs
            output: x slopes ccor n_wfs array, y slopes ccor n_wfs array
        """
        x_ccor = []
        y_ccor = []
        for i in range(len(self.active_wfs)):
            x_corr_i = []
            y_corr_i = []
            for j in range(len(self.active_wfs)):
                if i == j and self.acor:
                    x_acorr, y_acorr = self.data_get_cc(i,i,med_sub,avg_sub,avg_len)
                    x_corr_i.append(x_acorr)
                    y_corr_i.append(y_acorr)
                if i > j:
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
    
    def fits_name_gen(self):
        """
        Gives back the fits file name associated with parameters
            args: none
            returns: string (fits file)
        """
        out_file = self.out_dir + self.name + "_tmax" + str(self.tmax) +"_corr"
        if self.s_sub:
            out_file = out_file + "_ssub"
        if self.tt_sub:
            out_file = out_file + "_ttsub"
        out_file = out_file +".fits"
        return out_file
    
    def fits_pull(self, fits_file=None):
        """
        Takes a fits file, checks to see if format, pulls information into object.
            args: fits_file (string path to file)
            returns: Boolean (if pull was successful)
        """
        if fits_file and fits_file != self.fits_file:
            if self.data_check(fits_file):
                self.fits_file = fits_file
        if self.fits_file:
            hdulist = fits.open(self.fits_file)
            hdr = hdulist[0].header
            self.name = hdr.get('DATANAME')
            self.data_file = hdr.get('DATAFILE')
            self.out_dir = hdr.get('OUTPATH')
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
            if self.ccor or self.acor:
                x_cor = [hdulist[i+1].data[0] for i in range(self.n_wfs)]
                y_cor = [hdulist[i+1].data[1] for i in range(self.n_wfs)]
                print(np.array(x_cor).shape)
                print(np.array(y_cor).shape)   
            if self.ccor:
                self.ccor_x = np.array(x_cor)
                self.ccor_y = np.array(y_cor)
                self.acor_from_ccor(self.ccor_x, self.ccor_y)
            elif self.acor:
                self.acor_x = np.array(x_cor)
                self.acor_y = np.array(y_cor)
            self.data_valid = self.slopes_gen()
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
            print("Not able to pull data.")
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
        hdr['NWFS'] = (self.n_wfs, "WFS used in correction")
        hdr['TMAX'] = (self.tmax, 'Max time taken in temporal correlations')
        hdr['SSUB'] = (self.s_sub, "If static slope subtracted before comp")
        hdr['TTSUB'] = (self.tt_sub, "If global tip/tilt subtracted before comp")
        hdr["CCOR"] = (False, "Contains cross correlations")
        hdr["CCOR"] = (False, "Contains only auto correlations")
        primary_hdu = fits.PrimaryHDU(header=hdr)
        hdl_list = [primary_hdu]
        ##### set up headers
        for i in range(len(self.active_wfs)):
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
            print("---> saving ccor")
            #save for each wfs
            f[0].header["ccor"] = True
            for i in range(len(self.active_wfs)):
                f[i+1].data = np.array([self.ccor_x[i], self.ccor_y[i]])
        elif self.acor:
            print("---> saving acor")
            f[0].header["acor"] = True
            for i in range(len(self.active_wfs)):
                f[i+1].data = np.array([self.acor_x[i], self.acor_y[i]])
        else:
            print("---> no cor found")
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
        title = self.name + ", tmax="+ str(self.tmax) + ", " + title_type 
        if self.s_sub:
            title = title + ", s_sub"
        if self.tt_sub:
            title = title + ", tt_sub"
        if med_sub:
            title = title + ", med sub, avg_len "+ str(avg_len)
        elif avg_sub:
            title = title + ", avg sub, avg_len "+ str(avg_len) 
        return title
    
    def plot_file_gen(self, plot_type, file_type, med_sub, avg_sub, avg_len):
        """
        Generates a file name for a plot, renumbers if name already taken
            args: plot_type (varies by plot func), file_type ("png" or "gif"), Graphing specifications
            retuns: out_file (descriptive string)
        """
        out_file_base = self.out_dir + self.name
        out_file_base = out_file_base + "_tmax" + str(self.tmax) + plot_type 
        if self.s_sub:
            out_file_base = out_file_base + "_ssub"
        if self.tt_sub:
            out_file_base = out_file_base + "_ttsub"
        # Specify for graphs specs
        if med_sub:
            out_file_base = out_file_base + "_medsub"+ str(avg_len)
        elif avg_sub:
            out_file_base = out_file_base + "_avgsub" + str(avg_len)
        # Renumber if repeats 
        i = 0
        while os.path.exists(out_file_base + "_%s.%s" % (i, file_type)):
            i += 1
        out_file = out_file_base + "_%s.%s" % (i, file_type)  
        return out_file
    
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
            print("Warning: Auto corr not available")
            self.acor_gen()
        elif not self.data_valid:
            print("Error: Data not available")
            return -1             
        # figure out what to graph (mean_sub or avg_sub)
        data_x, data_y = self.data_get_ac(med_sub, avg_sub, avg_len)
        # masking data, taking average over wfs axis
        x_out = np.average(data_x, axis = 0)
        y_out = np.average(data_y, axis = 0)
        # checking tmax
        t_list = self.check_t_list(t_list, x_out.shape[0])
        # generate names
        title = self.plot_title_gen(" Auto Corr, WFS avg", med_sub, avg_sub, avg_len)
        out_file = self.plot_file_gen("_acor_avg", "png", med_sub, avg_sub, avg_len)
        # graph
        mat_cube_1 = x_out
        mat_cube_2 = y_out
        mat_cube_3 = (x_out + y_out)/2
        label_1, label_2, label_3 = "Sx acor", "Sy acor", "Avg acor"
        try:
            print("Graphing: graph_3_rows_t_mat")
            graph_3_rows_t_mat(mat_cube_1, mat_cube_2, mat_cube_3, 
                           t_list, title, 
                           label_1, label_2, label_3).savefig(out_file)
            plt.close("all")
            return out_file
        except:
            print("Error in graphing: graph_3_rows_t_mat")
            return False
          
        
    def acor_animate(self, dt_max = 30, med_sub=False, avg_sub=False, avg_len=10):
        """
        Create an animation length dt_max, of the all WFS (x, y, xy avg)
            args: Graph specifications
            outputs: gif file location, or false
        """
        # making sure that there is acorr data to work with 
        if self.data_valid and not self.acor:
            print("Warning: Auto corr not available")
            self.acor_gen()
        elif not self.data_valid:
            print("Error: Data not available")
            return None       
        tmax = self.tmax   
        # retrieve data for graphing
        data_x, data_y = self.data_get_ac(med_sub, avg_sub, avg_len)
        # Plot title and file
        title = self.plot_title_gen(" Auto Corr, all WFS", med_sub, avg_sub, avg_len)
        out_file = self.plot_file_gen("_acor", "gif", med_sub, avg_sub, avg_len)
        # graph
        label_1, label_2, label_3 = "Sx acor", "Sy acor", "Avg acor"
        try:
            print("Graphing: gif_3_rows_mat")
            gif_3_rows_mat(data_x, data_y, (data_x + data_y)/2, 
                  tmax, title, out_file, label_1, label_2, label_3)
            plt.close("all")
            return out_file
        except:
            print("Error in graphing: gif_3_mat_vmin")
            return None
        
        
    def acor_animate_avg(self, dt_max = 100, med_sub=False, avg_sub=False, avg_len=10):
        """
        Create an animation length dt_max, of the average WFS acor
            args: Graph specifications
            returns: gif file location, or false
        """
        # making sure that there is acorr data to work with 
        if self.data_valid and not self.acor:
            print("Warning: Auto corr not available")
            self.acor_gen()
        elif not self.data_valid:
            print("Error: Data not available")
            return None         
        tmax = self.tmax
        # retrieve data for graphing
        data_x, data_y = self.data_get_ac(med_sub, avg_sub, avg_len)   
        # averaging over valid wfs
        x_out = np.average(data_x[self.active_wfs], axis=0)
        y_out = np.average(data_y[self.active_wfs], axis=0)     
        # Plot title and file
        title = self.plot_title_gen(" Auto Corr, WFS avg", med_sub, avg_sub, avg_len)
        out_file = self.plot_file_gen("_acor_avg", "gif", med_sub, avg_sub, avg_len)
        # graph
        label_1, label_2, label_3 = "Sx acor", "Sy acor", "Avg acor"
        try:
            print("Graphing: gif_3_mat_vmin")
            gif_3_mat_vmin(x_out, y_out, (x_out + y_out)/2, 
                  tmax, title, out_file, label_1, label_2, label_3)
            plt.close("all")
            return out_file
        except:
            print("Error in graphing: gif_3_mat_vmin")
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
            print("Warning: Cross corr not available")
            self.ccor_gen()
        elif not self.data_valid:
            print("Error: Data not available")
            return None
        # ordering wfs
        if wfs_a > wfs_b: wfs_a, wfs_b = wfs_b, wfs_a
        # retrieve data for graphing
        data_cx, data_cy = self.data_get_cc(wfs_a, wfs_b, med_sub, avg_sub, avg_len)      
        # checking tmax
        t_list = self.check_t_list(t_list, data_cx.shape[0])
        # Plot title and file
        title = self.plot_title_gen("Cross Corr, WFS"  + str(wfs_a) + " WFS"+ str(wfs_b), med_sub, avg_sub, avg_len)
        out_file = self.plot_file_gen("_ccor_wfs" + str(wfs_a) + str(wfs_b), "png", med_sub, avg_sub, avg_len)
        # graph
        label_1, label_2, label_3 = "Sx ccor", "Sy ccor", "Avg ccor"
        try:
            print("Graphing: graph_3_rows_t_mat")
            graph_3_rows_t_mat(data_cx, data_cy, (data_cx + data_cy)/2, 
                           t_list, title, 
                           label_1, label_2, label_3).savefig(out_file)
            plt.close("all")
            return out_file
        except:
            print("Error in graphing: graph_3_rows_t_mat")
            return False
            
    
    def ccor_graph_all(self, t=0, med_sub=False, avg_sub=False, avg_len=10):
        """
        Graphs all wfs, their xy-average slope cross corr, given timeslice t
            args: graphing preferences
            returns: string (file location)
        """
        # check to see if data is valid
        if self.data_valid and not self.ccor:
            print("Warning: Cross corr not available")
            self.ccor_gen()
        elif not self.data_valid:
            print("Error: Data not available")
            return None
        # Retrieve data for graphing
        data_cx, data_cy = self.data_get_cc_all(med_sub, avg_sub, avg_len)
        # Plot title and file
        title = self.plot_title_gen(" Cross Corr, all WFS,", med_sub, avg_sub, avg_len)
        out_file = self.plot_file_gen("_ccor", "png", med_sub, avg_sub, avg_len)       
        # graph
        try:
            print("Graphing: graph_5_5_ccor")
            graph_5_5_ccor((data_cx+data_cy)/2, t, title).savefig(out_file)
            plt.close("all")
            return out_file
        except:
            print("Error in graphing: graph_5_5_ccor")
            return False
        
    def ccor_animate(self, wfs_a, wfs_b, t_max, med_sub=False, avg_sub=False, avg_len=10):
        """
        Static plot of wfs_a and wfs_b, their a, y, and xy-average slope cross corr, given timeslices in list t
            args: graphing preferences
            returns: string (file location)
        """
        # check to see if data is valid
        if self.data_valid and not self.ccor:
            print("Warning: Cross corr not available")
            self.ccor_gen()
        elif not self.data_valid:
            print("Error: Data not available")
            return None
        # ordering wfs
        if wfs_a > wfs_b: wfs_a, wfs_b = wfs_b, wfs_a
        # Retrieve data for graphing
        data_cx, data_cy = self.data_get_cc(wfs_a, wfs_b, med_sub, avg_sub, avg_len)
        # Plot title and file
        title = self.plot_title_gen(" Cross Corr, WFS"  + str(wfs_a) + " WFS"+ str(wfs_b), med_sub, avg_sub, avg_len)
        out_file = self.plot_file_gen("_ccor_wfs" + str(wfs_a) + str(wfs_b), "gif", med_sub, avg_sub, avg_len) 
        # graph
        label_1, label_2, label_3 = "Sx ccor", "Sy ccor", "Avg ccor"
        try:
            print("Graphing:  gif_3_mat_vmin")
            gif_3_mat_vmin(data_cx, data_cy, (data_cx + data_cy)/2, 
                  t_max, title, out_file, label_1, label_2, label_3)
            plt.close("all")
            return out_file
        except:
            print("Error in graphing: gif_3_mat_vmin")
            return False

        
    def cor_animate_all(self, t_max, med_sub=False, avg_sub=False, avg_len=10):
        # check to see if data is valid
        if self.data_valid and not self.ccor:
            print("Warning: Cross corr not available")
            self.ccor_gen()
        elif not self.data_valid:
            print("Error: Data not available")
            return None 
        # Retrieve data for graphing
        data_cx, data_cy = self.data_get_cc_all(med_sub, avg_sub, avg_len)
        # Plot title and file
        title = self.plot_title_gen(" Cross Corr, all WFS", med_sub, avg_sub, avg_len)
        out_file = self.plot_file_gen("_ccor_wfsall", "gif", med_sub, avg_sub, avg_len) 
        # animate all correlations 
        try:
            print("Graphing: cor_animate_all")
            gif_ccor_mat((data_cx+data_cy), t_max, title, out_file)
            plt.close("all")
            return out_file
        except:
            print("Error in graphing: gif_3_mat_vmin")
            return False
        
    ####################### Graphing MISC #######################
    ### Functions that aren't stardard images reps of the correlation matricies
    
    
    def cor_graph_cpeak(self, wfs1, wfs2, avg_sub=False):
        # graph the center peak, choosing which cross corr, if want avg_mean
        #TODO
        t = np.arange(tmax)
        center_peak_x = x_slopes_td_a_corr[:, :, 7, 7]
        center_peak_y = y_slopes_td_a_corr[:, :, 7, 7]
        center_peak_avg = avg_slopes_td_a_corr[:, :, 7, 7]

        p_cmap = cm.get_cmap('Purples', 8)
        o_cmap = cm.get_cmap('Oranges', 8)
        g_cmap = cm.get_cmap('Greens', 8)

        plt.figure(figsize=(10,8))

        #for i in WFS_list: 
        #    c_float = (i+1)/len(WFS_list+1)
        #    plt.plot(t, center_peak_x[i][t], label = "Sx WFS" + str(i), linestyle=':', color = p_cmap(c_float))
        #    plt.plot(t, center_peak_y[i][t], label = "Sy WFS" + str(i), linestyle=':', color = o_cmap(c_float))
            #plt.plot(t, center_peak_avg[i][t], label = "WFS" + str(i), color = g_cmap(c_float))
    
        plt.title('WFS corr peak vs time, X and Y slopes')
        plt.ylabel('Intensity')
        plt.xlabel('time')
        plt.legend(loc = 'upper right')
        return 0
    
    def cor_graph_tslice(self, wfs1, wfs2, avg_sub=False):
        # graph the center peak, choosing which cross corr, if want avg_mean
        # TODO
        return 0
    
    def cor_animate_tslice(self, wfs1, wfs2, avg_sub=False):
        # graph the center peak, choosing which cross corr, if want avg_mean
        # TODO
        return 0