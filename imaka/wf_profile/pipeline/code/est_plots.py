## est_plots.py
## Eden McEwen
## September 25, 2020
# all plots used by the Estimator class

import math
import time
import os
import sys
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import statistics as stat
from scipy import signal
from astropy.io import fits
import pandas as pd
from astropy.stats import sigma_clipped_stats
from celluloid import Camera

# Using Clustering
import hdbscan
import seaborn as sns

# Self-witten code
from Estimator import *
from tpoint import TPoint

############################        
### Plot Naming Functs #####
############################ 

#check and make file paths
def path_check(path):
    if os.path.isdir(path):
        return True
    else:
        try:
            os.makedirs(path)
            print("made: " + path)
            return True
        except:
            print("FAILED at making: " +path)
            return False

# adds in the 4 params to the plot type to give you a file
def plot_file_gen(out_file, plot_type, file_type, xcor, avg_len, sig, thresh, c_min):
    """
    Generates a file name for a plot, renumbers if name already taken
        args: plot_type (varies by plot func), file_type ("png" or "gif"), Graphing specifications
        retuns: out_file (descriptive string)
    """
    # generating file structure more finely
    if xcor: out_file = out_file + "_xcor_"
    out_file = out_file + plot_type
    out_file = out_file + "_asub" + str(avg_len)
    out_file = out_file + "_sig" + str(sig)
    out_file = out_file + "_thresh" + str(thresh)
    out_file = out_file + "_cmin" + str(c_min)
    
    # Renumber if repeats 
    i = 0
    while os.path.exists(out_file + "_%s.%s" % (i, file_type)):
        i += 1
    out_f= out_file + "_%s.%s" % (i, file_type)
    return out_f

############################        
### Plotting functions #####
############################ 

# Plot spd dir
# Plot est centers
#
##################
def plot_speed_dir_est(out_dir, name, dirs, v, clsr, estimates, params):
    print(name + " plot_speed_dir_est")
    xcor, avg_len, sig, thresh, c_min = params    
    fig, ax = plt.subplots()
    fig.set_size_inches(12, 8)
    plt.xlim([0,360])
    plt.ylim([0,100])    
    ## plotting angle gridlines
    for d in range(0, 360, 45):
        plt.axvline(x=d, alpha=0.1, color='black')   
    ## plotting extimated centers
    n_est = len(estimates)
    for n in range(n_est):
            plt.errorbar(estimates[n][1][0], estimates[n][0][0], xerr= estimates[n][1][1], yerr=estimates[n][0][1], fmt='o', label="Peak {}".format(n))
    
    palette = sns.color_palette()
    cluster_colors = [sns.desaturate(palette[col], sat)
                  if col >= 0 else (0.5, 0.5, 0.5) for col, sat in
                  zip(clsr.labels_, clsr.probabilities_)]    
    ## plotting data    
    sc = plt.scatter(dirs, v, c=cluster_colors, alpha=0.75,  cmap=plt.cm.tab10)   
    plt.title('{} Scatter plot velocity and direction with estimates'.format(name) + "\n xcor %s, avg_len %s, sig %s, thresh %s, c_min %s"%(xcor, avg_len, sig, thresh, c_min))
    plt.xlabel('dir (degree)')
    plt.ylabel('vel (m/s)')
    plt.legend()   
    ## save file: 
    path = out_dir + "speed_vel_png/"
    ## check if path exists
    path_check(path)
    plot_f = path + name
    plot_type = "_spd_dir_clstr"
    file_type = "png"
    out_file = plot_file_gen(plot_f, plot_type, file_type, xcor, avg_len, sig, thresh, c_min)
    plt.savefig(out_file)
    plt.close('all') 

# Plot XY
# Clustered Peaks
################################
def plot_max_pts(out_dir, name, x, y, v, params):
    # translate to all of the points in that list, turn them into 
    #finding velocities and ploting?  
    xcor, avg_len, sig, thresh, c_min = params
    
    fig, ax = plt.subplots()
    fig.set_size_inches(10, 8) 
    ax.set_xlim([0,15])
    ax.set_ylim([0,15])
    
    # generating lines to graph
    x_0, y_0 = angle_line((7,7), 0, 10)
    x_90, y_90 = angle_line((7,7), 90, 10)
    plt.plot(x_0, y_0, color = "black", alpha=0.25, ls="-")
    plt.plot(x_90, y_90, color = "black", alpha=0.25, ls="-")
    
    # plotting Data
    sc = ax.scatter(x, y, c=v, s=80, edgecolor='', cmap=plt.cm.cool, vmin=0, vmax=60, alpha=0.3)
    fig.colorbar(sc)
    plt.title("{} Clustered peaks, in xy".format(name)+ "\n xcor %s, avg_len %s, sig %s, thresh %s, c_min %s"%(xcor, avg_len, sig, thresh, c_min))
    
    ## save file: 
    path = out_dir + "xy_max_png/"
    ## check if path exists
    path_check(path)
    plot_f = path + name
    plot_type = "_xy_max_clstr"
    file_type = "png"
    out_file = plot_file_gen(plot_f, plot_type, file_type, xcor, avg_len, sig, thresh, c_min)
    plt.savefig(out_file)
    plt.close('all') 
    
# Animate XY
# Clustered Peaks
################################
def gif_cor_peaks(mat, points, tmax, params, clstr = False, out_d = "", name = None, scale_v = True):
    print(name + " gif_cor_peaks")
    xcor, avg_len, sig, thresh, c_min = params
    
    fig, ax = plt.subplots(ncols=1,figsize=(6,6))
    fig.subplots_adjust(wspace=0.3)
    if clstr: 
        plt.suptitle('{} Used CLUSTERED Max peaks'.format(name) + "\n xcor %s, avg_len %s, sig %s, thresh %s, c_min %s"%(xcor, avg_len, sig, thresh, c_min))
    else:
        plt.suptitle('{} Used Max peaks'.format(name) + "\n xcor %s, avg_len %s, sig %s, thresh %s, c_min %s"%(xcor, avg_len, sig, thresh, c_min))
    ax.set_xticklabels([])
    
    max_val = np.amax(mat)
    min_val = np.amin(mat)

    camera = Camera(fig)
    for t in range(tmax):
        if scale_v:
            im  = ax.matshow(mat[t])
        else:
            im  = ax.matshow(mat[t], vmin=min_val, vmax=max_val)
        for pnt in points:
            if pnt.t == t: ax.scatter(pnt.x, pnt.y, c="r")
            if pnt.t < t: ax.scatter(pnt.x, pnt.y, c="w", alpha = 0.5)
        ax.text(0.05,0.95," t = " + str(t), bbox=dict(facecolor='white', alpha=0.5), transform=ax.transAxes)
        #fig.colorbar(im, cax=cax)
        camera.snap()
    animation = camera.animate()  
    
    # save file: 
    if clstr:
        path = out_d + "xy_clstr_pts_gif/"
        plot_type = "_xy_clstr_pts"
    else:
        path = out_d + "xy_max_pts_gif/"
        plot_type = "_xy_max_pts"
    #check if path exists
    path_check(path)
    plot_f = path + name
    file_type = "gif"
    out_file = plot_file_gen(plot_f, plot_type, file_type, xcor, avg_len, sig, thresh, c_min)
    
    animation.save(out_file, writer = 'imagemagick')


#####################################       
### Plotting functions BATCH FILES 
##################################### 
    
def plot_dir(df, name, upper=False):
    est_dirs = df["p_dir"].tolist()
    cfht_dirs = df["cfht_wdir"].tolist()
    fig, ax = plt.subplots()
    fig.set_size_inches(8, 6)
    plt.hist(est_dirs, bins=np.linspace(0, 360, 36), histtype=u'step', density=True, label='Estimates')
    plt.hist(cfht_dirs, bins=np.linspace(0, 360, 36), histtype=u'step', density=True, label='CFHT measurement')
    if upper:
        mb_dirs = df["250_wdir"].tolist()
        plt.hist(mb_dirs, bins=np.linspace(0, 360, 36), histtype=u'step', density=True, label='250mb measurement')
    plt.title(name + 'Wind direction population histogram')
    plt.xlabel('direction (degree)')
    plt.ylabel('Counts (%)')
    plt.legend()
    return plt
    
    
def plot_spd(df, name,  upper=False):
    est_spd = df["p_spd"].tolist()
    cfht_spd = df["cfht_wspd"].tolist()
    fig, ax = plt.subplots()
    fig.set_size_inches(8, 6)
    plt.hist(est_spd, bins=np.linspace(0, 60, 60), histtype=u'step', density=True, label='Estimates')
    plt.hist(cfht_spd, bins=np.linspace(0, 60, 60), histtype=u'step', density=True, label='Cfht', alpha =0.5)    
    if upper:
        mb_spd = df["250_wspd"].tolist()
        plt.hist(mb_spd, bins=np.linspace(0, 60, 60), histtype=u'step', density=True, label='250mb', alpha =0.5)
    plt.title(name + ' Wind speed population histogram')
    plt.xlabel('speed (m/s)')
    plt.ylabel('Counts (%)')
    plt.legend()
    return plt
    
    
def plot_spd_dir(df, name, upper=False):
    #TODO: subplots of lines
    est_spd = df["p_spd"].tolist()
    cfht_spd = df["cfht_wspd"].tolist()
    mb_spd = df["250_wspd"].tolist()
    est_dirs = df["p_dir"].tolist()
    cfht_dirs = df["cfht_wdir"].tolist()
    mb_dirs = df["250_wdir"].tolist()
    
    fig, (ax1, ax2) = plt.subplots(ncols=2, figsize=(10,4))
    fig.suptitle(name + ' Speed and direction population histogram')
    
    #AX1: SPD
    ax1.hist(est_spd, bins=np.linspace(0, 50, 50), histtype=u'step', density=True, label='Est')
    ax1.hist(cfht_spd, bins=np.linspace(0, 50, 50), histtype=u'step', density=True, label='CFHT', alpha =0.5)
    if upper: ax1.hist(mb_spd, bins=np.linspace(0, 50, 50), histtype=u'step', density=True, label='250mb', alpha =0.5)
    ax1.set_xlabel('speed (m/s)')
    ax1.legend()
    
    #AX2: DIR
    ax2.hist(est_dirs, bins=np.linspace(0, 360, 36), histtype=u'step', density=True, label='Est')
    ax2.hist(cfht_dirs, bins=np.linspace(0, 360, 36), histtype=u'step', density=True, label='CFHT')
    if upper: ax2.hist(mb_dirs, bins=np.linspace(0, 360, 36), histtype=u'step', density=True, label='250mb')
    ax2.set_xlabel('direction (degrees)')

    plt.legend()
    return plt

def hist2d_spd_dir(df, name, upper=False):
    est_spd = df["p_spd"].tolist()
    cfht_spd = df["cfht_wspd"].tolist()
    mb_spd = df["250_wspd"].tolist()
    est_dirs = df["p_dir"].tolist()
    cfht_dirs = df["cfht_wdir"].tolist()
    mb_dirs = df["250_wdir"].tolist()
    
    fig, (ax1, ax2) = plt.subplots(1, 2)
    fig.set_size_inches(13, 5)
    fig.suptitle("Directional counts")

    x = np.arange(360)
    y = x
    
    if upper:
        h1 = ax1.hist2d(mb_dirs, est_dirs, bins = 36, range=[[0, 360], [0, 360]])
        h2 = ax2.hist2d(mb_spd, est_spd, bins = 25, range=[[0, 50], [0, 50]])
        name = name + " 250mb"
    else:
        h1 = ax1.hist2d(cfht_dirs, est_dirs, bins = 36, range=[[0, 360], [0, 360]])
        h2 = ax2.hist2d(cfht_spd, est_spd, bins = 25, range=[[0, 50], [0, 50]])
        name = name + " CFHT"
    
    ax1.plot(x, y, color = "white")
    ax1.set_ylabel("Estimated direction (degree)")
    ax1.set_xlabel(name + " measured direction")
    fig.colorbar(h1[3], ax=ax1)

    ax2.plot(x, y, color = "white")
    ax2.set_ylabel("Estimated speed (m/s)")
    ax2.set_xlabel(name + " measured speed, shifted")
    fig.colorbar(h2[3], ax=ax2)

    return plt

#Plot 6 graps: step, histogram for both speed and direction
def plot_est_6(df, name=""):
    fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(12,6))
    fig.tight_layout()
    #fig.suptitle(str(name) + "Directional counts")
    
    #spd data
    est_spd = df["p_spd"].tolist()
    cfht_spd = df["cfht_wspd"].tolist()
    mb_spd = df["250_wspd"].tolist()
    #dir data
    est_dirs = df["p_dir"].tolist()
    cfht_dirs = df["cfht_wdir"].tolist()
    mb_dirs = df["250_wdir"].tolist()
    #prep
    x = np.arange(360)
    y = x
    
    ########### ROW 1 ############
    # Plot all speeds
    # Step hist
    ax1 = axes[0][0]
    ax1.hist(est_spd, bins=np.linspace(0, 50, 50), histtype=u'step', density=True, label='Est')
    ax1.hist(cfht_spd, bins=np.linspace(0, 50, 50), histtype=u'step', density=True, label='CFHT', alpha =0.5)
    ax1.hist(mb_spd, bins=np.linspace(0, 50, 50), histtype=u'step', density=True, label='250mb', alpha =0.5)
    ax1.set_xlabel('speed (m/s)')
    ax1.legend()
    
    # CFHT 2dhist
    ax = axes[0][1]
    h = ax.hist2d(cfht_spd, est_spd, bins = 36, range=[[0, 50], [0, 50]])
    ax.plot(x, y, color = "white")
    #ax.set_ylabel("Estimated direction (m/s)")
    #ax.set_xlabel("CFHT measured speed")
    ax.set_title("CFHT  measured v est speed")
    #fig.colorbar(h[3], ax=ax)
    
    # GFS 2dhist
    ax = axes[0][2]
    h = ax.hist2d(mb_spd, est_spd, bins = 36, range=[[0, 50], [0, 50]])
    ax.plot(x, y, color = "white")
    #ax.set_ylabel("Estimated direction (m/s)")
    #ax.set_xlabel("250mb measured speed")
    ax.set_title("250mb model vs est speed")
    #fig.colorbar(h[3], ax=ax)
    
    
    ########### ROW 2 ############
    # Plot all directions
    
    # Step hist
    ax2 = axes[1][0]
    ax2.hist(est_dirs, bins=np.linspace(0, 360, 36), histtype=u'step', density=True, label='Est')
    ax2.hist(cfht_dirs, bins=np.linspace(0, 360, 36), histtype=u'step', density=True, label='CFHT')
    ax2.hist(mb_dirs, bins=np.linspace(0, 360, 36), histtype=u'step', density=True, label='250mb')
    ax2.set_xlabel('direction (degrees)')
    
    # CFHT 2dhist
    ax = axes[1][1]
    h = ax.hist2d(cfht_dirs, est_dirs, bins = 36, range=[[0, 360], [0, 360]])
    ax.plot(x, y, color = "white")
    #ax.set_ylabel("Estimated direction (degree)")
    #ax.set_xlabel("CFHT measured direction")
    ax.set_title("CFHT measured v est dir")
    #fig.colorbar(h[3], ax=ax)
    
    # GFS 2dhist
    ax = axes[1][2]
    h = ax.hist2d(mb_dirs, est_dirs, bins = 36, range=[[0, 360], [0, 360]])
    ax.plot(x, y, color = "white")
    #ax.set_ylabel("Estimated direction (degree)")
    #ax.set_xlabel("250mb measured direction")
    ax.set_title("250mb model vs est dir")
    #fig.colorbar(h[3], ax=ax)
    
    return plt

###################
### PLOT HELPERS###
###################

def angle_line(point, angle, length):
    # unpack the first point
    x, y = point
    # find the start point
    starty = y - length * math.cos(math.radians(angle))
    startx = x - length * math.sin(math.radians(angle))
    # find the end point
    endy = y + length * math.cos(math.radians(angle))
    endx = x + length * math.sin(math.radians(angle))
    
    return ([startx, endx], [starty, endy])