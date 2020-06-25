# Eden McEwen
# June 18th 2019
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
    def __init__(self, name, data_f, out_d, s_sub=False, tmax=200):
        self.name = name
        self.out_dir = out_d
        self.data_file = data_f
        self.fits_file = None
        self.path_valid = self.path_check(out_d)
        self.n_wfs = 5
        # TODO: modularize this?
        self.x_slopes = None 
        self.y_slopes = None
        self.data_valid = self.slopes_get()
        self.tmax = tmax
        # TODO: make a check on this
        self.target = None
        self.target_file = None
        self.active_wfs = [True]*self.n_wfs
        self.wfs_ra = [None]*self.n_wfs
        self.wfs_dec = [None]*self.n_wfs
        self.wfs_mag = [None]*self.n_wfs
        # Set container for autocorr
        self.acor = False
        self.acor_x = np.zeros(self.n_wfs)
        self.acor_y = np.zeros(self.n_wfs)     
        # Set container for crosscorr
        self.ccor = False
        self.ccor_x = np.zeros((self.n_wfs,self.n_wfs))
        self.ccor_y = np.zeros((self.n_wfs,self.n_wfs))

    ####################### Checking files #######################    
    
    def data_check(self, file):
        """
        Tries to find file
            args: file (string)
            return: T/F 
        """
        if os.path.isfile(file): 
            print("file found: %s" % file)
            return True
        print("WARNING! file not found: %s" % file)
        return False 

    def path_check(self, path):
        """
        Tries to find path, if not found, attempts to create 
            args: path (string)
            return: T/F (whether path is found/created)
        """
        #if path exists, return true, 
        #if path can be created, return true,
        #else, if error, false
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
        #returns true if there is a fits file in designated out file
        """
        Tries to find either the fits file given, if one
        Or tries to find the 
            args: path (string)
            return: T/F (whether path is found/created)
        """
        if self.fits_file:
            out_file = self.fits_file
        else:
            out_file = self.out_dir + self.name + "_corr_tmax" + str(self.tmax) +".fits"
        out = os.path.isfile(out_file)
        if out:
            self.fits_file = out_file
            print("Found fits file: ", out_file)
        return out
    
    ####################### Setting Data #######################
            
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
            args: s_sub (Boolean, T/F)
            return: T/F (whether or not s_sub val changed)
        """
        if s_sub == self.s_sub:
            print("Static Subtraction remains: ", s_sub)
        else: 
            print("Static Subtraction now: ", s_sub)
            self.s_sub = s_sub
            self.acor = False
            self.ccor = False
            self.fits_file = None
            self.data_valid = self.get_slopes()
            
    def set_target(self, file):
        """
        Takes in target file and fills in wfs_ra, wfs_dec and wfs_mag
            input: file (in specified format)
            return: T/F (whether or not target was updated)
        """
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
            print("Target file not found")
            return False
            
    ####################### Getting Data #######################
            
    def slopes_get(self):
        """
        If data is valid, pulls x and y slopes from datafile
        If fits is formatted differently, you'll want to change this file
            input: N/A
            return: T/F (whether or not slopes were updated)
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
                if self.s_sub:
                    
                return True
            else: 
                return False
        except:
            return False
        
        
    def data_get(self, med_sub=False, avg_sub=False, avg_len=10):
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
        else self.data_valid and not self.acor:
            print("Error with data_get: self.acor is %s, self.data_valid is %s"% (self.acor, self.data_valid))
            return None, None

    
    def data_get_cc(self, wfs_a, wfs_b, med_sub=False, avg_sub=False, avg_len=10):
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
    
    def data_get_cc_all(self, med_sub=False, avg_sub=False, avg_len=10):
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
        if not self.acor:
            
        return np.array(x_ccor), np.array(y_ccor)

    ####################### Generating correlations #######################
    
    def acor_gen(self, mask=mask_8_8_center):
        # create a list of autocorrelations for each WFS
        if self.acor == True:
            print("Auto corr already generated, tmax = %s" % self.tmax)
        elif self.ccor == True:
            self.acor_from_ccor(self, self.ccor_x, self.ccor_y)
        elif self.data_valid:
            t_max = self.tmax
            print("Generating auto corr tmax = %s" % t_max)
            self.acor = True
            self.acor_x = np.array([td_auto_corr(mat, t_max, mask) for mat in self.x_slopes])
            self.acor_y = np.array([td_auto_corr(mat, t_max, mask) for mat in self.y_slopes])
        else:
            print("Data invalid")

    def ccor_gen(self, mask=mask_8_8_center):
        if self.ccor == True:
            print("Cross corr already generated, tmax = %s" % self.tmax)
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
        else:
            print("Data invalid")
            
    def acor_from_ccor(self, ccorx, ccory):
        print("Generating auto corr from ccor tmax = %s" % tmax)
        self.acor_x = np.array([ci[i] for i, ci in enumerate(ccorx)]) 
        self.acor_y = np.array([ci[i] for i, ci in enumerate(ccory)])
        self.acor = True
        return self.acor_x, self.acor_y
    
    ####################### FITS files #######################
    
    def fits_name_gen(self):
        """
        Gives back the fits file name associated with parameters
            args: N/A
            returns: string (fits file)
        """
        out_file = self.out_dir + self.name + "_corr_tmax" + str(self.tmax)
        if self.s_sub:
            out_file = out_file + "_ssub"
        out_file = out_file +".fits"
        return out_file
    
    def fits_pull(self, fits_file):
        """
        Takes a fits file, checks to see if format, pulls information into object.
            args: fits_file (string path to file)
            returns: Boolean (if pull was successful)
        """
        if self.data_check(fits_file):
            self.fits_file = fits_file
            hdulist = fits.open(self.data_file)
            hdr = hdulist[0].header
            self.name = hdr['DATANAME']
            self.data_file = hdr['DATAFILE']
            self.out_dir = hdr['OUTPATH'] 
            self.target = hdr['TARGET']
            self.target_file = hdr['TFILE']
            self.n_wfs = hdr['NWFS']
            self.tmax = hdr['TMAX']
            #Retrieve info from all wfs
            self.ccor = hdr["CCOR"] 
            self.acor = hdr["ACOR"]
            x_cor = []
            y_cor = []
            for i in range(self.n_wfs):
                hdr_i = hdulist[i+1].header
                self.wfs_ra[i] = hdr_i['RA']
                self.wfs_dec[i] = hdr_i['DEC']
                self.wfs_mag[i] = hdr_i['MAG']
                self.active_wfs[i] = hdr_i['WFSVALID']
                if self.ccor or self.acor:
                    np.append(x_cor, hdulist[i+1].data[0])
                    np.append(y_cor, hdulist[i+1].data[1])
            if self.ccor:
                self.ccor_x = np.array(x_cor)
                self.ccor_y = np.array(y_cor)
                acor_from_ccor(self, self.ccor_x, self.ccor_y)
            elif self.acor:
                self.acor_x = np.array(x_cor)
                self.acor_y = np.array(y_cor)
            self.slopes_get()
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
            print("Not able to pull data.")
            return -1
        ### define file name
        self.fits_file = self.fits_name_gen()
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
        hdr['TFILE'] = (self.target_file, "Target field file")
        hdr['NWFS'] = (self.n_wfs, "WFS used in correction")
        hdr['TMAX'] = (self.tmax, 'Max time taken in temporal correlations')
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
            hdr['RA'] = (self.wfs_ra[i], "WFS GS right ascension")
            hdr['DEC'] = (self.wfs_dec[i], "WFS GS declination")
            hdr['MAG'] = (self.wfs_mag[i], "WFS GS MAG")
            #if self.target:
            #    hdr['WFSVALID'] = (self.active_wfs[i], "If WFS useable, based on TARGETFILE")
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
        # Set up if haven't already
        if not self.fits_file:
            self.fits_gen()
        out_file = self.fits_file      
        ##### saving correlations
        f = fits.open(out_file, mode='update')
        if self.ccor:
            print("---> saving ccor")
            #save for each wfs
            f[0].header["ccor"] = True
            for i in range(len(self.active_wfs)):
                f[i+1].data = np.array([self.ccor_x[i], self.ccor_y[i]])
        elif self.acor:
            print("  saving acor")
            f[0].header["acor"] = True
            for i in range(len(self.active_wfs)):
                f[i+1].data = np.array([self.acor_x[i], self.acor_y[i]])
        else:
            print("---> no cor found")
        f.flush()
        f.close()
        return out_file
    
    ####################### Graphing ACOR #######################
    def plot_title_gen(self, title_type, med_sub, avg_sub, avg_len):
        title = self.name + title_type +" tmax="+ str(self.tmax) 
        if med_sub:
            out_file_base = out_file_base + "_medsub"+ str(avg_len)
            title = title + ", med sub, avg_len "+ str(avg_len)
        elif avg_sub:
            title = title + ", avg sub, avg_len "+ str(avg_len) 
        return title
    
    def plot_file_gen(self, plot_type, file_type, med_sub, avg_sub, avg_len):
        """
        Graph naming scheme based on 
        """
        out_file_base = self.out_dir + self.name
        if self.s_sub:
            out_file_base = out_file_base + "_ssub"
        out_file_base = out_file_base + plot_type +"_tmax" + str(self.tmax)
        if med_sub:
            out_file_base = out_file_base + "_medsub"+ str(avg_len)
        elif avg_sub:
            out_file_base = out_file_base + "_avgsub" + str(avg_len)
        i = 0
        while os.path.exists(out_file_base + "_%s.%s" % (i, file_type)):
            i += 1
        out_file = out_file_base + "_%s.%s" % (i, file_type)  
        
        return out_file
    
    def acor_graph(self, t_list=[0,5,10,15,20], med_sub=False, avg_sub=False, avg_len=10):
        """
        Create a plot of WFS acor (x, y, and xy cor) t times in t_list 
        inputs: Graph specifications
        outputs: png file location, or false
        """
        # making sure that there is acorr data to work with 
        if self.data_valid and not self.acor:
            print("Auto corr not available")
            self.acor_gen()
        elif not self.data_valid:
            return -1            
        # base name of file
        out_file_base = self.out_dir + self.name + "_acor_avg_tmax" + str(self.tmax)
        title = self.name + " Auto Corr, WFS avg, tmax="+ str(self.tmax)       
        # figure out what to graph (mean_sub or avg_sub)
        data_x, data_y = self.data_get(med_sub, avg_sub, avg_len)
        if med_sub:
            print("Calculating median subtracted")
            out_file_base = out_file_base + "_medsub"+ str(avg_len)
            title = title + " median subtracted"
        elif avg_sub:
            print("Calculating average subtracted")
            out_file_base = out_file_base + "_avgsub" + str(avg_len)
            title = title + " average subtracted"        
        # masking data, taking average over wfs axis
        x_out = np.average(data_x, axis = 0)
        y_out = np.average(data_y, axis = 0)
        print(y_out.shape)
        # figure out what to name the file
        title, out_file_base = self.plot_name_gen(self, title, out_file_base, med_sub, avg_sub, avg_len)
        i = 0
        while os.path.exists(out_file_base + "_%s.png" % i):
            i += 1
        out_file = out_file_base + "_%s.png" % i
        # graph
        mat_cube_1 = x_out
        mat_cube_2 = y_out
        mat_cube_3 = (x_out + y_out)/2
        label_1 = "Sx acor"
        label_2 = "Sy acor"
        label_3 = "Avg acor"
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
        inputs: Graph specifications
        outputs: gif file location, or false
        """
        # create an animation max length dt_max, all wfs
        # create an image of average over valid wfs acor at times in t_list 
        # making sure that there is acorr data to work with 
        if self.data_valid and not self.acor:
            print("Auto corr not available")
            self.acor_gen()
        elif not self.data_valid:
            return -1        
        tmax = self.tmax
        # base name of file
        out_file_base = self.out_dir + self.name + "_acor_tmax" + str(tmax)
        title = self.name + " Auto Corr, all WFS, tmax="+ str(tmax)        
        # figure out what to graph (mean_sub or avg_sub)
        data_x, data_y = self.data_get(med_sub, avg_sub, avg_len)
        if med_sub:
            print("Calculating median subtracted")
            out_file_base = out_file_base + "_medsub"+ str(avg_len)
            title = title + " median subtracted"
        elif avg_sub:
            print("Calculating average subtracted")
            out_file_base = out_file_base + "_avgsub"+ str(avg_len)
            title = title + " average subtracted"        
        # figure out what to name the file
        i = 0
        while os.path.exists(out_file_base + "_%s.gif" % i):
            i += 1
        out_file = out_file_base + "_%s.gif" % i        
        # graph
        label_1 = "Sx acor"
        label_2 = "Sy acor"
        label_3 = "Avg acor"
        try:
            print("Graphing: gif_3_mat_vmin")
            gif_3_rows_mat(data_x, data_y, (data_x + data_y)/2, 
                  tmax, title, out_file, label_1, label_2, label_3)
            plt.close("all")
            return out_file
        except:
            print("Error in graphing: gif_3_mat_vmin")
            return False
        
        
    def acor_animate_avg(self, dt_max = 100, med_sub=False, avg_sub=False, avg_len=10):
        """
        Create an animation length dt_max, of the average WFS acor
        inputs: Graph specifications
        outputs: gif file location, or false
        """
        # making sure that there is acorr data to work with 
        if self.data_valid and not self.acor:
            print("Auto corr not available")
            self.acor_gen()
        elif not self.data_valid:
            return -1        
        tmax = self.tmax
        # base name of file
        out_file_base = self.out_dir + self.name + "_acor_avg_tmax" + str(tmax)
        title = self.name + " Auto Corr, WFS avg, tmax="+ str(tmax)
        # figure out what to graph (mean_sub or avg_sub)
        data_x, data_y = self.data_get(med_sub, avg_sub, avg_len)
        if med_sub:
            print("Calculating median subtracted")
            out_file_base = out_file_base + "_medsub"+ str(avg_len)
            title = title + " median subtracted"
        elif avg_sub:
            print("Calculating average subtracted")
            out_file_base = out_file_base + "_avgsub"+ str(avg_len)
            title = title + " average subtracted"      
        # masking data, taking average over wfs axis
        x_out = np.average(data_x[self.active_wfs], axis=0)
        y_out = np.average(data_y[self.active_wfs], axis=0)     
        # figure out what to name the file
        i = 0
        while os.path.exists(out_file_base + "_%s.gif" % i):
            i += 1
        out_file = out_file_base + "_%s.gif" % i  
        # graph
        label_1 = "Sx acor"
        label_2 = "Sy acor"
        label_3 = "Avg acor"
        try:
            print("Graphing:  gif_3_mat_vmin")
            gif_3_mat_vmin(x_out, y_out, (x_out + y_out)/2, 
                  tmax, title, out_file, label_1, label_2, label_3)
            return out_file
        except:
            print("Error in graphing: gif_3_mat_vmin")
            return False
    
    ####################### Graphing CCOR #######################
        
    def ccor_graph(self, wfs_a, wfs_b, t_list=[0,5,10,15,20], med_sub=False, avg_sub=False, avg_len=10):
        """
        Static plot of wfs_a and wfs_b, their xy-average slope cross corr, given timeslices in list t
            input: graphing preferences
            returns: string (file location)
        """
        # check to see if data is valid
        if self.data_valid and not self.ccor:
            print("cross corr not available")
            self.ccor_gen()
        elif not self.data_valid:
            return -1
        # ordering wfs
        if wfs_a > wfs_b: wfs_a, wfs_b = wfs_b, wfs_a
        # base name of file
        out_file_base = self.out_dir + self.name + "_ccor_tmax" + str(self.tmax) + "_wfs" + str(wfs_a) + str(wfs_b)
        title = self.name + ", Cross Corr, WFS"  + str(wfs_a) + " WFS"+ str(wfs_b)+ ", tmax="+ str(self.tmax)
        # figure out what to graph (mean_sub or avg_sub)
        data_cx, data_cy = self.data_get_cc(wfs_a, wfs_b, med_sub, avg_sub, avg_len)
        if med_sub:
            print("Calculating median subtracted")
            out_file_base = out_file_base + "_medsub"+ str(avg_len)
            title = title + " median sub, avg_len "+ str(avg_len)
        elif avg_sub:
            print("Calculating average subtracted")
            out_file_base = out_file_base + "_avgsub"+ str(avg_len)
            title = title + " average sub, avg_len "+ str(avg_len)
        # figure out what to name the file
        i = 0
        while os.path.exists(out_file_base + "_%s.png" % i):
            i += 1
        out_file = out_file_base + "_%s.png" % i
        # graph
        label_1 = "Sx ccor"
        label_2 = "Sy ccor"
        label_3 = "Avg ccor"
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
            input: graphing preferences
            returns: string (file location)
        """
        # check to see if data is valid
        if self.data_valid and not self.ccor:
            print("cross corr not available")
            self.ccor_gen()
        elif not self.data_valid:
            return -1
        # base name of file
        out_file_base = self.out_dir + self.name + "_ccor_tmax" + str(self.tmax) + "_allwfs"
        title = self.name + ", Cross Corr, all WFS, tmax="+ str(self.tmax)
        # Normalizing and taking average
        data_cx, data_cy = self.data_get_cc_all(med_sub, avg_sub, avg_len)
        if med_sub:
            print("Calculating median subtracted")
            out_file_base = out_file_base + "_medsub"+ str(avg_len)
            title = title + " median sub, avg_len "+ str(avg_len)
        elif avg_sub:
            print("Calculating average subtracted")
            out_file_base = out_file_base + "_avgsub"+ str(avg_len)
            title = title + " average sub, avg_len "+ str(avg_len)
        # figure out what to name the file
        i = 0
        while os.path.exists(out_file_base + "_%s.png" % i):
            i += 1
        out_file = out_file_base + "_%s.png" % i
        # graph
        try:
            print("Graphing: graph_5_5_ccor")
            graph_5_5_ccor((data_cx+data_cy)/2, t, title).savefig(out_file)
            plt.close("all")
            return out_file
        except:
            print("Error in graphing: graph_5_5_ccor")
            return False
        
    def ccor_animate(self, wfs_a, wfs_b, t_max):
        """
        Static plot of wfs_a and wfs_b, their a, y, and xy-average slope cross corr, given timeslices in list t
            input: graphing preferences
            returns: string (file location)
        """
        # check to see if data is valid
        if self.data_valid and not self.ccor:
            print("cross corr not available")
            self.ccor_gen()
        elif not self.data_valid:
            return -1
        # ordering wfs
        if wfs_a > wfs_b: wfs_a, wfs_b = wfs_b, wfs_a
        # base name of file
        out_file_base = self.out_dir + self.name + "_ccor_tmax" + str(self.tmax) + "_wfs" + str(wfs_a) + str(wfs_b)
        title = self.name + ", Cross Corr, WFS"  + str(wfs_a) + " WFS"+ str(wfs_b)+ ", tmax="+ str(self.tmax)
        # figure out what to graph (mean_sub or avg_sub)
        data_cx, data_cy = self.data_get_cc(wfs_a, wfs_b, med_sub, avg_sub, avg_len)
        if med_sub:
            print("Calculating median subtracted")
            out_file_base = out_file_base + "_medsub"+ str(avg_len)
            title = title + " median sub, avg_len "+ str(avg_len)
        elif avg_sub:
            print("Calculating average subtracted")
            out_file_base = out_file_base + "_avgsub"+ str(avg_len)
            title = title + " average sub, avg_len "+ str(avg_len)
        # figure out what to name the file
        i = 0
        while os.path.exists(out_file_base + "_%s.gif" % i):
            i += 1
        out_file = out_file_base + "_%s.gif" % i
        # graph
        label_1 = "Sx ccor"
        label_2 = "Sy ccor"
        label_3 = "Avg ccor"
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
            print("cross corr not available")
            self.ccor_gen()
        elif not self.data_valid:
            return -1
        # base name of file
        out_file_base = self.out_dir + self.name + "_ccor_tmax" + str(self.tmax) + "_allwfs"
        title = self.name + ", Cross Corr, all WFS, tmax="+ str(self.tmax)
        # Normalizing and taking average
        data_cx, data_cy = self.data_get_cc_all(med_sub, avg_sub, avg_len)
        if med_sub:
            print("Calculating median subtracted")
            out_file_base = out_file_base + "_medsub"+ str(avg_len)
            title = title + " median sub, avg_len "+ str(avg_len)
        elif avg_sub:
            print("Calculating average subtracted")
            out_file_base = out_file_base + "_avgsub"+ str(avg_len)
            title = title + " average sub, avg_len "+ str(avg_len)
        # figure out what to name the file
        i = 0
        while os.path.exists(out_file_base + "_%s.gif" % i):
            i += 1
        out_file = out_file_base + "_%s.gif" % i
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
        return 0
    
    def cor_animate_tslice(self, wfs1, wfs2, avg_sub=False):
        # graph the center peak, choosing which cross corr, if want avg_mean
        return 0