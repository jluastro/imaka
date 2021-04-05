### Def Estimate R
## Eden McEwen
## February 2, 2020

import math
import time
import os
import sys
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import matplotlib
import statistics as stat
from scipy import signal
from astropy.io import fits
import pandas as pd
from astropy.stats import sigma_clipped_stats
from celluloid import Camera

# Windplot for detections
from windrose import WindroseAxes

# Using Clustering
import hdbscan
import seaborn as sns
from numpy import unique
from numpy import where
from sklearn.datasets import make_classification
from sklearn.cluster import MeanShift
from sklearn.cluster import DBSCAN
from sklearn.cluster import AffinityPropagation

# Self-witten code
from pipeline.code.tpoint import TPoint
from pipeline.code.corr_code import *
from pipeline.code.Correlator import *
from pipeline.code.est_plots import *
from pipeline.code.cluster import *

#generating a mask for corr data
mask_data = mask_8_8_center
mask_cor = signal.correlate2d(mask_data, mask_data)
mask = ~np.array(mask_cor, dtype=bool)


#### Static code ####
#Cov map radii
rad_map = np.array([[np.sqrt(i**2 + j**2) for j in np.arange(-7, 8)] for i in np.arange(-7, 8)])
rad_map_boolean = np.zeros((8, 15, 15))
# boolean mask matrix
for i in np.arange(0, 8):
    rad_max = ma.masked_where(rad_map <=  i + 0.5, rad_map)
    rad_min = ma.masked_where(rad_map >  i - 0.5, rad_map)
    rad = rad_max.mask & rad_min.mask
    rad_map_boolean[i-1] = rad
rad_map_inv = np.ones(rad_map_boolean.shape) - rad_map_boolean


#### 

class Estimate_simple(object):
    
    def __init__(self, file):
        self.out_fits = file
        self.data = Correlator("", "", "", f_file = file)
        self.name = self.data.name
        self.date = self.data.date
        # these will change each run
        # auto correlation 
        self.a_dtct = None
        self.a_spds = None
        self.a_dirs = None
        self.a_xs = None
        self.a_ys = None
        self.a_cstr = None
        # cross correlation
        self.x_dtct = None
        self.x_spds = None
        self.x_dirs = None
        self.x_xs = None
        self.x_ys = None
        self.x_cstr = None
        
    def acor_map(self):
        data = self.data
        x_acor, y_acor = data.data_get_ac(avg_sub=False, avg_len=0)
        avg_wfs = np.average((x_acor + y_acor)/2, axis=0)
        return avg_wfs 
        
    def xcor_map(self):
        data = self.data
        x_cor, y_cor = data.data_get_cc_all(avg_sub=False, avg_len=0)
        avg_cor = (x_cor + y_cor)/2
        wfs_use = data.active_wfs
        avg_xcor = np.zeros_like(avg_cor[0][0])
        count = 0
        for i, row in enumerate(avg_cor):
            for j in range(len(row)):
                if i < j and wfs_use[i] and wfs_use[j]:
                    avg_xcor = avg_xcor + row[j]
                    count = count + 1
        avg_xcor = np.divide(avg_xcor, count)
        return avg_xcor
    
    def run(self, xcor=False, detect_clp = 4, clstr_method=detect_cluster_meanshift):
        #Cluster methods
        # detect_cluster_meanshift
        # detect_cluster_dbscan
        # detect_cluster_affprop
        
        # choosing between auto or crodd correlations 
        if not xcor:
            avg_awfs = self.acor_map()
        if xcor:    
            avg_awfs = self.xcor_map()
        #subtracting averages pixel value from frame
        avg_awfs_sub = np.subtract(avg_awfs, np.average(avg_awfs, axis = 0))
        # PART 1 : determining detection levels
        detect_lvl = detect_map(avg_awfs_sub)
        # PART 2 : Speed Map
        dtct, spds, dirs, xs, ys = speed_map(detect_lvl, detect_clp, 333)
        # PART 3: Clustering
        clusters, yhat = clstr_method(self.data, spds, xs, ys)
        # saving these parts based on a or x cor     
        if not xcor:
            self.a_dtct = dtct
            self.a_spds, self.a_dirs = spds, dirs
            self.a_xs, self.a_ys= xs, ys
            self.a_cstr = clusters
            self.a_yhat = yhat
        if xcor: 
            self.x_dtct = dtct
            self.x_spds, self.x_dirs = spds, dirs
            self.x_xs, self.x_ys = xs, ys
            self.x_cstr = clusters
            self.x_yhat = yhat
        return spds, dirs
    
    def return_table(self, sdv_cnd = 0, n_filter = 0, clstr_method = detect_cluster_meanshift):
        # Checking to see if system has been run before
        if self.a_dtct is None:
            self.run(clstr_method=clstr_method)
        if self.x_dtct is None:
            self.run(xcor = True, clstr_method=clstr_method)
        # Returning summary of cluster
        asum = self.summary_cluster(sdv_cnd = sdv_cnd, n_filter = n_filter)
        xsum = self.summary_cluster(xcor = True, sdv_cnd = sdv_cnd, n_filter = n_filter)
        # Classifying peaks
        aclass, xclass = self.classify_peaks(asum, xsum)
        # df for acor
        dfa = pd.DataFrame(asum, columns =['dir', 'dir_std', 'spd', 'spd_std', 'count'], dtype = float)
        dfa['class'] = aclass
        dfa['xcor'] = 0 
        # df for xcor
        dfx = pd.DataFrame(xsum, columns =['dir', 'dir_std', 'spd', 'spd_std', 'count'], dtype = float)
        dfx['class'] = xclass
        dfx['xcor'] = 1 
        df = pd.concat([dfa, dfx], ignore_index=True)
        df['name'] = self.name
        return df
    
    def summary_cluster(self, xcor=False, sdv_cnd = 0, n_filter = 0):
        # based off of cluster output
        # include counts per cluster     
        #choosing saved variables base on xcor/acor
        if not xcor and self.a_dtct is None:
            self.run()
        if xcor and self.x_dtct is None:
            self.run(xcor = True)
        dtct = self.x_dtct if xcor else self.a_dtct
        spds = self.x_spds if xcor else self.a_spds
        dirs = self.x_dirs if xcor else self.a_dirs
        xs = self.x_xs if xcor else self.a_xs
        ys = self.x_ys if xcor else self.a_ys
        clusters = self.x_cstr if xcor else self.a_cstr
        yhat = self.x_yhat if xcor else self.a_yhat        
        # for each cluster return summary of cluster
        summary = []
        for cluster in clusters:
            # get row indexes for samples with this cluster
            row_ix = where(yhat == cluster)
            ## Standardized return 5 element structure
            mean_x = np.mean(xs[row_ix])
            mean_y = np.mean(ys[row_ix])
            mean_dir = to_degrees(dir_c_proj(mean_x, mean_y, 0))
            std_dir = to_degrees(dir_std(dirs[row_ix]))
            summary.append([mean_dir,
                            std_dir,
                            np.mean(spds[row_ix]),
                            np.std(spds[row_ix]),
                            row_ix[0].shape[0]])
        # condensing peaks
        summary = self.condense_peaks(summary, sdv_cnd)
        # filtering peaks
        summary = self.filter_peaks(summary, n_filter)
        return summary
    
    def classify_peaks(self, acor_sum, xcor_sum):
        # standard deviation and means
        # mean within a standard deviation
        # setting up flag arrays
        acor_clas = np.array(["FL" for peak in acor_sum])
        xcor_clas = np.array(["NA" for peak in xcor_sum])
        #iteratively comparing between acor and xcor
        for a in range(len(acor_sum)):
            for x in range(len(xcor_sum)):
                if self.compare_peaks(acor_sum[a], xcor_sum[x]):
                    xcor_clas[x] = "GL"
                    acor_clas[a] = "GL"
        return acor_clas, xcor_clas
    
    def condense_peaks(self, summary, sdv_cnd):
        # takes in a summary based on func summary_cluster
        # looks at peaks within n standard deviations (sdv_cnd)
        if sdv_cnd == 0:
            return summary
        # compare all the peaks, and condense those that are closer
        # sort by count
        summary = sorted(summary, key=lambda x: x[4], reverse=True)
        sum_new = []
        for a1 in range(len(summary)):
            # seed with first array element
            if a1 == 0:
                sum_new = [summary[a1]]
            else:
                avged = 0
                for a2 in range(len(sum_new)):
                    if avged == 0 and self.compare_peaks(summary[a1], sum_new[a2], sd = sdv_cnd):
                        # average peaks
                        pk = self.avg_peaks(summary[a1], sum_new[a2])
                        # change sum_new
                        sum_new[a2] = pk
                        avged = 1
                        break
                if avged == 0 :
                    sum_new.append(summary[a1])
        return sum_new

    def avg_peaks(self, pk1, pk2):
        # combines two peaks
        dir_avg = np.average([pk1[0], pk2[0]])
        spds_avg = np.average([pk1[2], pk2[2]])
        dir_sdv = np.sqrt(pk1[1]**2 + pk2[1]**2)
        spds_sdv = np.sqrt(pk1[3]**2 + pk2[3]**2)
        counts = pk1[4] + pk2[4]
        return [dir_avg, dir_sdv, spds_avg, spds_sdv, counts]
    
    def filter_peaks(self, summary, n_filter):
        # only returns peaks that are greater or eaqual to filter number
        return [pk for pk in summary if pk[4] >= n_filter]
    
    def compare_peaks(self, peak1, peak2, sd=3):
        #returns if two peaks are similar
        spd_diff = np.abs(peak1[0] - peak2[0])
        dir_diff = np.abs(peak1[2] - peak2[2])
        #HARDCODED looking in range of 3 times stdev
        spd_stdev = sd*np.max([peak1[1], peak2[1]])
        dir_stdev = sd*np.max([peak1[3], peak2[3]])
        #boolean values if these match
        spd_q = spd_diff <= spd_stdev
        dir_q = dir_diff <= dir_stdev
        return spd_q & dir_q
    
    def plot_speeds(self, dirs, spds):
        plt.scatter(np.degrees(dirs), spds, alpha=0.5)
        plt.title(self.name)
        plt.xlim([-180, 180])
        plt.ylim([0, 30])
        plt.xlabel('Dirs (degrees)')
        plt.ylabel('speed (m/s?)')
        plt.show()
        
    def plot_windrose(self, dirs, spds):
        ax = WindroseAxes.from_ax()
        ax.contourf((np.degrees(dirs) + 360)%360, spds, bins=np.arange(0.1, 35, 5), cmap=cm.viridis)
        ax.set_legend()
        ax.set_title(self.name + " Detections speed  and direction")
        plt.show()
        
    def plot_spds_detect(self, acor=True, xcor=True):
        # pipeline iterating

        fig = plt.figure(figsize=(10,5))
        plt.xlim([0, 360])
        plt.ylim([0, 30])
        
        if acor:
            plt.scatter(np.degrees(self.a_spds), self.a_dirs, c="green",  alpha = .3, label="a")
        
        if xcor:
            plt.scatter(np.degrees(self.x_spds), self.x_dirs, c="red",  alpha = .3, label="x")
            
        plt.title( self.name + 'Detections')
        plt.ylabel('wind speed')
        plt.xlabel('wind dir')
        plt.legend()
        
        return fig


        
        
# static functions
 
##############
### Part 1 ###
##############
def detect_map(avg_wf):
    #detect_lvl = np.zeros_like(avg_wfs)
    t_mean = np.zeros_like(avg_wf)
    t_stdev = np.zeros_like(avg_wf)

    for t in range(avg_wf.shape[0]):
        # 1.1.  For each pixel in the cov map assign a radius
        t_slice = avg_wf[t]
        sigma = 3
        # 1.2. find the mean (m) and stddev (sd) of the pixels 
        # in a radial annuli one pixel wide (r from 1 to 7 +/- 0.5.  
        # I do a sigma clipping on this with a clip=3 si
        for r in np.arange(8):
            rad_vals = np.multiply(t_slice, rad_map_boolean[r])
            mean, median, std = sigma_clipped_stats(rad_vals, sigma=sigma, mask=rad_map_inv[r])
            # TRY REDOING MASKING => try to flatten?
            #print(t, r, mean, std)
            t_mean[t] = t_mean[t] + mean*rad_map_boolean[r]
            t_stdev[t] = t_stdev[t] + std*rad_map_boolean[r]
    # 1.3. For each pixel I assign a "detection level in sigma over the mean"
    # = (im[x,y] - m)/sd.   
    # This results in a detection map for each cov map time slice. 
    detect_map = np.divide(np.subtract(avg_wf, t_mean), t_stdev)
    return detect_map
    
    
##############
### Part 2 ###
##############

# For each detection greater than some value I calculate a speed (vspd), direction (vdir). 
# This results in a speed and a direction "detection" map for each time slice.

def dist_c(x, y, c):
    #BUG: hard-coded 2.2 meter telescope
    pix_to_m = 2.2/15.0
    return np.sqrt(np.square(x-c) + np.square(y-c))

# issues with wind projections

def dir_wind(x,y,c):
    # this projects from (0,7) as north, clockwise angle
    return 0

def dir_c(x,y,c):
    return np.arctan2(x-c, c-y) + np.pi

def dir_c_proj(xp,yp,c):
    return np.arctan2(xp, yp)

def proj_x(dir_rad, c):
    return c*np.sin(dir_rad)

def proj_y(dir_rad, c):
    return c*np.cos(dir_rad)

def dir_std(dirs):
    dirs_shift = (dirs + np.pi) % (2*np.pi)
    return np.min([np.std(dirs), np.std(dirs_shift)]) 

def to_degrees(rads):
    return (np.degrees(rads) + 360)%360

    
##############
### Part 3 ###
##############

# I also calculate the x,y location at the edge of the "circle" where the vdir projects to.
# Namely, xedgei = r0 * cos(thetai) 
# This makes the clustering analysis easier as it does away with the 2pi wrapping
    
def speed_map(detection_map, detect_val, hz):
    radius = 7
    
    #maps of the distance and direction for each pixel
    dist_map = np.array([[dist_c(x,y, radius) for y in range(15)] for x in range(15)])
    dir_map = np.array([[dir_c(x,y, radius) for y in range(15)] for x in range(15)])
    
    #storing arrays
    detect_list = np.array([])
    speed_list = np.array([])
    dir_list = np.array([])
    x_proj = np.array([])
    y_proj = np.array([])
    
    #values above detection value => make a map for detections
    i = 1
    for t_slice in detection_map:
        t_slice = ma.masked_invalid(t_slice)
        detect_t = ma.masked_where(t_slice <= detect_val, t_slice)
        #BUG: positive detections only
        
        t = i/hz
        speed_t = np.divide(dist_map, t*(8/2.2))
        
        #detections
        detect_det = t_slice[detect_t.mask == False]
        detect_list = np.append(detect_list, detect_det)
        # speeds 
        speed_det = speed_t[detect_t.mask == False]
        speed_list = np.append(speed_list, speed_det)
        # directions 
        dir_det = dir_map[detect_t.mask == False]
        dir_list = np.append(dir_list, dir_det)
        # x-y projections
        x_t = proj_x(dir_det, radius)
        y_t = proj_y(dir_det, radius)
        x_proj = np.append(x_proj, x_t)
        y_proj = np.append(y_proj, y_t)
        
        i+=1
        #if i>100: break
        
    return [detect_list,speed_list, dir_list, x_proj, y_proj]

##############
### Part 4 ###
##############
# various clustering algorithms
# see cluster.py
    

    
    
    
    
### Plotting

def plot_clstr_table(table, name="Data", min_count=0):
    fig = plt.figure(figsize=(10,5))
    plt.xlim([0, 360])
    plt.ylim([0, 50])
    ax = plt.gca()

    for index, row in table[table['count']>min_count].iterrows():
        clr  = det_color(row["class"])
        #alph = det_alph(row["count"])
        alph= 0.5
        plt.scatter(row['dir'], row['spd'], c=clr,  alpha = 1)
        spread = matplotlib.patches.Ellipse((row['dir'], row['spd']), 3*row['dir_std'], 3*row['spd_std'], angle=0, color=clr, alpha=alph)
        ax.add_patch(spread)

    plt.title(name +' Estimates, all')
    plt.ylabel('wind speed')
    plt.xlabel('wind dir')
    plt.legend()
    
    return fig

def det_color(c):
    if c == "NA":
        return 'black'
    elif c == "GL":
        return 'blue'
    elif c == "FL":
        return 'orange'