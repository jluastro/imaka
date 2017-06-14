import csv
import numpy as np
import pylab as plt
import scipy.spatial
from astropy.table import Table
import subprocess
import os, pdb
from imaka.reduce import util
from datetime import datetime
from matplotlib import dates as mp_dates
from matplotlib import ticker
import glob
import matplotlib.pyplot as plt
import matplotlib
from imaka.analysis import add_data
from astropy.io import fits
from pandas import read_csv


def fetch_stats_from_onaga(dates, output_root):
    """
        SCP stats files (FITS tables) from onaga and place them in a mirrored
        directory structure on your local machine.
        Parameters
        ----------
        dates : list/array of strings
        List of date strings (e.g. '20170113') to process.
        You can also use 'all' or None and the entire available
        list of all dates will be used.
        output_root : str
        The root directory where the transferred files will be stored. Note
        that this is just the root and under this will be stored a structure
        that parallels that on onaga:
        <output_root>/<date>/fli/reduce/stats/*
        """
    all_dates = ['20170109', '20170110', '20170111', '20170112', '20170113']
    
    if dates == 'all' or dates == None:
        dates = all_dates
    
    for date in dates:
        remote_file = 'imaka@onaga.ifa.hawaii.edu:/Volumes/DATA4/imaka/{0:s}/fli/reduce/stats/stats*.fits'.format(date)
        destination = '{0:s}/{1:s}/fli/reduce/stats/'.format(output_root, date)
        
        util.mkdir(destination)
        
        print(remote_file)
        print(destination)
        p = subprocess.Popen(["scp", remote_file, destination])
        sts = os.waitpid(p.pid, 0)

    return

def load_stats_to_onaga(dates, input_root):
    """
        SCP stats files (FITS tables) from onaga and place them in a mirrored
        directory structure on your local machine.
        Parameters
        ----------
        dates : list/array of strings
        List of date strings (e.g. '20170113') to process.
        You can also use 'all' or None and the entire available
        list of all dates will be used.
        input_root : str
        The root directory where the files to be transferred are pulled from. Note
        that this is just the root and under this should be a structure
        that parallels that on onaga:
        <input_root>/<date>/fli/reduce/stats/*
        """
    all_dates = ['20170110', '20170111', '20170112', '20170113']
    
    if dates == 'all' or dates == None:
        dates = all_dates
    
    for date in dates:
        local_file = '{0:s}/{1:s}/fli/reduce/stats/stats*.fits'.format(input_root, date)
        destination = 'imaka@onaga.ifa.hawaii.edu:/Volumes/DATA4/imaka/{0:s}/fli/reduce/stats/'.format(date)
        
        cmd = "scp " + local_file + " " + destination
        print(cmd)
        p = subprocess.Popen(cmd, shell=True)
        sts = os.waitpid(p.pid, 0)

    return


def plot_stack_stats(date, suffixes=['open', 'ttf', 'closed'], root_dir='/Users/jlu/work/imaka/pleiades/'):
    """
        Make a suite of standard plots for the stats on a given night.
        Parameters
        ----------
        date : str
        The date string for which to plot up the stats (i.e. '20170113').
        Optional Parameters
        -------------------
        suffix : str
        stats files have the name stats_open<suffix>.fits, stats_ttf<suffix>.fits, and
        stats_closed<suffix>.fits.
        root_dir : str
        The root directory for the <date> observing run directories. The
        stats files will be searched for in:
        <root_dir>/<date>/fli/reduce/stats/
        """
    stats_dir = root_dir + date + '/fli/reduce/stats/'
    plots_dir = root_dir + date + '/fli/reduce/plots/'
    
    util.mkdir(plots_dir)
    
    stats = []
    for suffix in suffixes:
        stats.append(Table.read(stats_dir + 'stats_' + suffix + '.fits'))
    
    scale = 0.04
    colors = get_color_list()
    
    #
    # Plots for ratio of improvements. First we need to find a matching
    # closed loop image to go with each open (and TTF) image.
    #
    tree_indices = []
    tree_so = scipy.spatial.KDTree(np.array([stats[0]['Index']]).T)
    for ss in range(len(suffixes)):
        dist_ss, idx_ss = tree_so.query(np.array([stats[ss]['Index']]).T, 1)
        tree_indices.append(idx_ss)
    
    #####
    # FWHM vs. Frame
    #####
    plt.figure(1, figsize=(12, 6))
    plt.subplots_adjust(left=0.1, bottom=0.15)
    plt.clf()
    plt.subplot(121)
    for ss in range(len(suffixes)):
        plt.plot(stats[ss]['Index'], stats[ss]['FWHM']*stats[ss]['BINFAC']*scale, marker='o', linestyle='none', label=suffixes[ss])
    plt.xlabel('Frame Number')
    plt.ylabel('Gaussian-Fit FWHM (")')
    plt.legend(numpoints=1)
    plt.ylim(0, 1.5)
    plt.title(date)

    plt.subplot(122)
    for ss in range(1, len(suffixes)):
        idx = tree_indices[ss]
        label = suffixes[ss] + " / " + suffixes[0]
        plt.plot(stats[0]['Index'][idx], stats[ss]['FWHM'] / stats[0]['FWHM'][idx], marker='o', linestyle='none', label=label)
    plt.xlabel('Frame Number')
    plt.ylabel('Ratio of Gaussian-Fit FWHM')
    plt.legend(numpoints=1)
    plt.axhline(1.0, color='k', linestyle='--', linewidth=2)
    plt.ylim(0, 1.3)
    plt.savefig(plots_dir + 'fwhm_vs_frame' + suffix + '.png')
    
    #####
    # Empirical FWHM
    #####
    plt.figure(2, figsize=(12, 6))
    plt.subplots_adjust(left=0.1, bottom=0.15)
    plt.clf()
    plt.subplot(121)
    for ss in range(len(suffixes)):
        plt.plot(stats[ss]['Index'], stats[ss]['emp_fwhm']*stats[ss]['BINFAC']*scale, marker='o', linestyle='none', label=suffixes[ss])
    plt.xlabel('Frame Number')
    plt.ylabel('Empirical FWHM (")')
    plt.legend(numpoints=1)
    plt.ylim(0, 1.5)

    plt.subplot(122)
    for ss in range(1, len(suffixes)):
        idx = tree_indices[ss]
        label = suffixes[ss] + " / " + suffixes[0]
        plt.plot(stats[0]['Index'][idx], stats[ss]['emp_fwhm'] / stats[0]['emp_fwhm'][idx], marker='o', linestyle='none', label=label)
    plt.xlabel('Frame Number')
    plt.ylabel('Ratio of Empircal FWHM')
    plt.legend(numpoints=1)
    plt.axhline(1.0, color='k', linestyle='--', linewidth=2)
    plt.ylim(0, 1.3)
    plt.savefig(plots_dir + 'efwhm_vs_frame' + suffix + '.png')
    
    
    #####
    # EE 50
    #####
    plt.figure(3, figsize=(12, 6))
    plt.subplots_adjust(left=0.1, bottom=0.15)
    plt.clf()
    plt.subplot(121)
    for ss in range(len(suffixes)):
        plt.plot(stats[ss]['Index'], stats[ss]['EE50'], marker='o', linestyle='none', label=suffixes[ss])
    plt.xlabel('Frame Number')
    plt.ylabel('Radius of 50% Encircled Energy (")')
    plt.legend(numpoints=1)
    plt.ylim(0, 1.5)
    plt.title(date)
    
    plt.subplot(122)
    for ss in range(1, len(suffixes)):
        idx = tree_indices[ss]
        label = suffixes[ss] + " / " + suffixes[0]
        plt.plot(stats[0]['Index'][idx], stats[ss]['EE50'] / stats[0]['EE50'][idx], marker='o', linestyle='none', label=label)
    plt.xlabel('Frame Number')
    plt.ylabel('Ratio of 50% EE Radius')
    plt.legend(numpoints=1)
    plt.axhline(1.0, color='k', linestyle='--', linewidth=2)
    plt.ylim(0, 1.3)
    plt.savefig(plots_dir + 'ee50_vs_frame' + suffix + '.png')

#####
# EE 80
#####
    plt.figure(4, figsize=(12, 6))
    plt.subplots_adjust(left=0.1, bottom=0.15)
    plt.clf()
    plt.subplot(121)
    for ss in range(len(suffixes)):
        plt.plot(stats[ss]['Index'], stats[ss]['EE80'], marker='o', linestyle='none', label=suffixes[ss])
    plt.xlabel('Frame Number')
    plt.ylabel('Radius of 80% Encircled Energy (")')
    plt.legend(numpoints=1)
    plt.ylim(0, 1.5)
    plt.title(date)
    
    plt.subplot(122)
    for ss in range(1, len(suffixes)):
        idx = tree_indices[ss]
        label = suffixes[ss] + " / " + suffixes[0]
        plt.plot(stats[0]['Index'][idx], stats[ss]['EE80'] / stats[0]['EE80'][idx], marker='o', linestyle='none', label=label)
    plt.xlabel('Frame Number')
    plt.ylabel('Ratio of 80% EE Radius')
    plt.legend(numpoints=1)
    plt.axhline(1.0, color='k', linestyle='--', linewidth=2)
    plt.ylim(0, 1.3)
    plt.savefig(plots_dir + 'ee80_vs_frame' + suffix + '.png')
    
    #####
    # NEA
    #####
    plt.figure(5, figsize=(12, 6))
    plt.subplots_adjust(left=0.1, bottom=0.15)
    plt.clf()
    plt.subplot(121)
    for ss in range(len(suffixes)):
        plt.plot(stats[ss]['Index'], stats[ss]['NEA'], marker='o', linestyle='none', label=suffixes[ss])
    plt.xlabel('Frame Number')
    plt.ylabel('NEA (Sq. Arcsec)')
    plt.legend(numpoints=1)
    plt.ylim(0, 5)
    plt.title(date)
    
    plt.subplot(122)
    for ss in range(1, len(suffixes)):
        idx = tree_indices[ss]
        label = suffixes[ss] + " / " + suffixes[0]
        plt.plot(stats[0]['Index'][idx], stats[ss]['NEA'] / stats[0]['NEA'][idx], marker='o', linestyle='none', label=label)
    plt.xlabel('Frame Number')
    plt.ylabel('Ratio of NEA')
    plt.legend(numpoints=1)
    plt.axhline(1.0, color='k', linestyle='--', linewidth=2)
    plt.ylim(0, 1.5)
    plt.savefig(plots_dir + 'nea_vs_frame' + suffix + '.png')
    
    #####
    # NEA2
    #####
    plt.figure(6, figsize=(12, 6))
    plt.subplots_adjust(left=0.1, bottom=0.15)
    plt.clf()
    plt.subplot(121)
    for ss in range(len(suffixes)):
        plt.plot(stats[ss]['Index'], stats[ss]['NEA2'], marker='o', linestyle='none', label=suffixes[ss])
    plt.xlabel('Frame Number')
    plt.ylabel('NEA2 (Sq. Arcsec)')
    plt.legend(numpoints=1)
    plt.ylim(0, 5)
    plt.title(date)
    
    plt.subplot(122)
    for ss in range(1, len(suffixes)):
        idx = tree_indices[ss]
        label = suffixes[ss] + " / " + suffixes[0]
        plt.plot(stats[0]['Index'][idx], stats[ss]['NEA2'] / stats[0]['NEA2'][idx], marker='o', linestyle='none', label=label)
    plt.xlabel('Frame Number')
    plt.ylabel('Ratio of NEA2')
    plt.legend(numpoints=1)
    plt.axhline(1.0, color='k', linestyle='--', linewidth=2)
    plt.ylim(0, 1.5)
    plt.savefig(plots_dir + 'nea2_vs_frame' + suffix + '.png')

#####
# FWHM for each direction
#####
    plt.figure(7, figsize=(12, 6))
    plt.subplots_adjust(left=0.1, bottom=0.15)
    plt.clf()
    plt.subplot(121)
    for ss in range(len(suffixes)):
        c = np.take(colors, ss, mode='wrap')
        plt.plot(stats[ss]['Index'], stats[ss]['xFWHM']*stats[ss]['BINFAC']*scale, marker='o', color=c, linestyle='none', label='X ' + suffixes[ss])
        plt.plot(stats[ss]['Index'], stats[ss]['yFWHM']*stats[ss]['BINFAC']*scale, marker='^', color=c, linestyle='none', label='Y ' + suffixes[ss])
    plt.xlabel('Frame Number')
    plt.ylabel('Gaussian-Fit FWHM (")')
    plt.legend(numpoints=1, fontsize=10)
    plt.ylim(0, 1.5)
    plt.title(date)
    
    plt.subplot(122)
    for ss in range(1, len(suffixes)):
        idx = tree_indices[ss]
        label = suffixes[ss] + " / " + suffixes[0]
        plt.plot(stats[0]['Index'][idx], stats[ss]['xFWHM'] / stats[0]['xFWHM'][idx], marker='o', color=c, linestyle='none', label=label)
        plt.plot(stats[0]['Index'][idx], stats[ss]['yFWHM'] / stats[0]['yFWHM'][idx], marker='^', color=c, linestyle='none', label=label)
    plt.xlabel('Frame Number')
    plt.ylabel('Ratio of Gaussian-Fit FWHM')
    plt.legend(numpoints=1, fontsize=10)
    plt.axhline(1.0, color='k', linestyle='--', linewidth=2)
    plt.ylim(0, 1.3)
    plt.savefig(plots_dir + 'xyfwhm_vs_frame' + suffix + '.png')


##########
# All plots vs. time.
##########

    utcs = []
    for ss in range(len(suffixes)):
        utc_dt = [datetime.strptime(stats[ss]['TIME_UTC'][ii], '%I:%M:%S') for ii in range(len(stats[ss]))]
        utcs.append(utc_dt)

    time_fmt = mp_dates.DateFormatter('%H:%M')
    
    
    #####
    # FWHM
    #####
    plt.figure(8, figsize=(12, 6))
    plt.subplots_adjust(left=0.1, bottom=0.15)
    plt.clf()
    plt.subplot(121)
    for ss in range(len(suffixes)):
        plt.plot(utcs[ss], stats[ss]['FWHM']*stats[ss]['BINFAC']*scale, marker='o', linestyle='none', label=suffixes[ss])
    plt.gca().xaxis.set_major_formatter(time_fmt)
    plt.xticks(rotation=35)
    plt.xlabel('UTC Time (hr)')
    plt.ylabel('Gaussian-Fit FWHM (")')
    plt.legend(numpoints=1)
    plt.ylim(0, 1.5)
    plt.title(date)
    
    plt.subplot(122)
    for ss in range(1, len(suffixes)):
        idx = tree_indices[ss]
        label = suffixes[ss] + " / " + suffixes[0]
        utc = [utcs[0][ii] for ii in idx]
        plt.plot(utc, stats[ss]['FWHM'] / stats[0]['FWHM'][idx], marker='o', linestyle='none', label=label)
    plt.gca().xaxis.set_major_formatter(time_fmt)
    plt.xticks(rotation=35)
    plt.xlabel('UTC Time (hr)')
    plt.ylabel('Ratio of Gaussian-Fit FWHM')
    plt.legend(numpoints=1)
    plt.axhline(1.0, color='k', linestyle='--', linewidth=2)
    plt.ylim(0, 1.3)
    plt.savefig(plots_dir + 'fwhm_vs_time' + suffix + '.png')

#####
# Empirical FWHM
#####
    plt.figure(9, figsize=(12, 6))
    plt.subplots_adjust(left=0.1, bottom=0.15)
    plt.clf()
    plt.subplot(121)
    for ss in range(len(suffixes)):
        plt.plot(utcs[ss], stats[ss]['emp_fwhm']*stats[ss]['BINFAC']*scale, marker='o', linestyle='none', label=suffixes[ss])
    plt.gca().xaxis.set_major_formatter(time_fmt)
    plt.xticks(rotation=35)
    plt.xlabel('UTC Time (hr)')
    plt.ylabel('Empirical FWHM (")')
    plt.legend(numpoints=1)
    plt.ylim(0, 1.5)
    plt.title(date)
    
    plt.subplot(122)
    for ss in range(1, len(suffixes)):
        idx = tree_indices[ss]
        label = suffixes[ss] + " / " + suffixes[0]
        utc = [utcs[0][ii] for ii in idx]
        plt.plot(utc, stats[ss]['emp_fwhm'] / stats[0]['emp_fwhm'][idx], marker='o', linestyle='none', label=suffixes[ss])
    plt.gca().xaxis.set_major_formatter(time_fmt)
    plt.xticks(rotation=35)
    plt.xlabel('UTC Time (hr)')
    plt.ylabel('Ratio of Empircal FWHM')
    plt.legend(numpoints=1)
    plt.axhline(1.0, color='k', linestyle='--', linewidth=2)
    plt.ylim(0, 1.3)
    plt.savefig(plots_dir + 'efwhm_vs_time' + suffix + '.png')
    
    #####
    # EE 50
    #####
    plt.figure(10, figsize=(12, 6))
    plt.subplots_adjust(left=0.1, bottom=0.15)
    plt.clf()
    plt.subplot(121)
    for ss in range(len(suffixes)):
        plt.plot(utcs[ss], stats[ss]['EE50'], marker='o', linestyle='none', label=suffixes[ss])
    plt.gca().xaxis.set_major_formatter(time_fmt)
    plt.xticks(rotation=35)
    plt.xlabel('UTC Time (hr)')
    plt.ylabel('Radius of 50% Encircled Energy (")')
    plt.legend(numpoints=1)
    plt.ylim(0, 1.5)
    plt.title(date)
    
    plt.subplot(122)
    for ss in range(1, len(suffixes)):
        idx = tree_indices[ss]
        label = suffixes[ss] + " / " + suffixes[0]
        utc = [utcs[0][ii] for ii in idx]
        plt.plot(utc, stats[ss]['EE50'] / stats[0]['EE50'][idx], marker='o', linestyle='none', label=label)
    plt.gca().xaxis.set_major_formatter(time_fmt)
    plt.xticks(rotation=35)
    plt.xlabel('UTC Time (hr)')
    plt.ylabel('Ratio of 50% EE Radius')
    plt.legend(numpoints=1)
    plt.axhline(1.0, color='k', linestyle='--', linewidth=2)
    plt.ylim(0, 1.3)
    plt.savefig(plots_dir + 'ee50_vs_time' + suffix + '.png')

#####
# EE 80
#####
    plt.figure(11, figsize=(12, 6))
    plt.subplots_adjust(left=0.1, bottom=0.15)
    plt.clf()
    plt.subplot(121)
    for ss in range(len(suffixes)):
        plt.plot(utcs[ss], stats[ss]['EE80'], marker='o', linestyle='none', label=suffixes[ss])
    plt.gca().xaxis.set_major_formatter(time_fmt)
    plt.xticks(rotation=35)
    plt.xlabel('UTC Time (hr)')
    plt.ylabel('Radius of 80% Encircled Energy (")')
    plt.legend(numpoints=1)
    plt.ylim(0, 1.5)
    plt.title(date)
    
    plt.subplot(122)
    for ss in range(1, len(suffixes)):
        idx = tree_indices[ss]
        label = suffixes[ss] + " / " + suffixes[0]
        utc = [utcs[0][ii] for ii in idx]
        plt.plot(utc, stats[ss]['EE80'] / stats[0]['EE80'][idx], marker='o', linestyle='none', label=label)
    plt.gca().xaxis.set_major_formatter(time_fmt)
    plt.xticks(rotation=35)
    plt.xlabel('UTC Time (hr)')
    plt.ylabel('Ratio of 80% EE Radius')
    plt.legend(numpoints=1)
    plt.axhline(1.0, color='k', linestyle='--', linewidth=2)
    plt.ylim(0, 1.3)
    plt.savefig(plots_dir + 'ee80_vs_time' + suffix + '.png')

#####
# NEA
#####
    plt.figure(12, figsize=(12, 6))
    plt.subplots_adjust(left=0.1, bottom=0.15)
    plt.clf()
    plt.subplot(121)
    for ss in range(len(suffixes)):
        plt.plot(utcs[ss], stats[ss]['NEA'], marker='o', linestyle='none', label=suffixes[ss])
    plt.gca().xaxis.set_major_formatter(time_fmt)
    plt.xticks(rotation=35)
    plt.xlabel('UTC Time (hr)')
    plt.ylabel('NEA (Sq. Arcsec)')
    plt.legend(numpoints=1)
    plt.ylim(0, 5)
    plt.title(date)
    
    plt.subplot(122)
    for ss in range(1, len(suffixes)):
        idx = tree_indices[ss]
        label = suffixes[ss] + " / " + suffixes[0]
        utc = [utcs[0][ii] for ii in idx]
        plt.plot(utc, stats[ss]['NEA'] / stats[0]['NEA'][idx], marker='o', linestyle='none', label=label)
    plt.gca().xaxis.set_major_formatter(time_fmt)
    plt.xticks(rotation=35)
    plt.xlabel('UTC Time (hr)')
    plt.ylabel('Ratio of NEA')
    plt.legend(numpoints=1)
    plt.axhline(1.0, color='k', linestyle='--', linewidth=2)
    plt.ylim(0, 1.5)
    plt.savefig(plots_dir + 'nea_vs_time' + suffix + '.png')

#####
# NEA2
#####
    plt.figure(13, figsize=(12, 6))
    plt.subplots_adjust(left=0.1, bottom=0.15)
    plt.clf()
    plt.subplot(121)
    for ss in range(len(suffixes)):
        plt.plot(utcs[ss], stats[ss]['NEA2'], marker='o', linestyle='none', label=suffixes[ss])
    plt.gca().xaxis.set_major_formatter(time_fmt)
    plt.xticks(rotation=35)
    plt.xlabel('UTC Time (hr)')
    plt.ylabel('NEA2 (Sq. Arcsec)')
    plt.legend(numpoints=1)
    plt.ylim(0, 5)
    plt.title(date)
    
    plt.subplot(122)
    for ss in range(1, len(suffixes)):
        idx = tree_indices[ss]
        label = suffixes[ss] + " / " + suffixes[0]
        utc = [utcs[0][ii] for ii in idx]
        plt.plot(utc, stats[ss]['NEA2'] / stats[0]['NEA2'][idx], marker='o', linestyle='none', label=label)
    plt.gca().xaxis.set_major_formatter(time_fmt)
    plt.xticks(rotation=35)
    plt.xlabel('UTC Time (hr)')
    plt.ylabel('Ratio of NEA2')
    plt.legend(numpoints=1)
    plt.axhline(1.0, color='k', linestyle='--', linewidth=2)
    plt.ylim(0, 1.5)
    plt.savefig(plots_dir + 'nea2_vs_time' + suffix + '.png')
    
    #####
    # FWHM for each direction
    #####
    plt.figure(14, figsize=(12, 6))
    plt.subplots_adjust(left=0.1, bottom=0.15)
    plt.subplot(121)
    for ss in range(len(suffixes)):
        c = np.take(colors, ss, mode='wrap')
        plt.plot(utcs[ss], stats[ss]['xFWHM']*stats[ss]['BINFAC']*scale, marker='o', color=c, linestyle='none', label='X ' + suffixes[ss])
        plt.plot(utcs[ss], stats[ss]['yFWHM']*stats[ss]['BINFAC']*scale, marker='^', color=c, linestyle='none', label='Y ' + suffixes[ss])
    plt.gca().xaxis.set_major_formatter(time_fmt)
    plt.xticks(rotation=35)
    plt.xlabel('UTC Time (hr)')
    plt.ylabel('Gaussian-Fit FWHM (")')
    plt.legend(numpoints=1, fontsize=10)
    plt.ylim(0, 1.5)
    plt.title(date)
    
    plt.subplot(122)
    for ss in range(1, len(suffixes)):
        c = np.take(colors, ss, mode='wrap')
        idx = tree_indices[ss]
        label = suffixes[ss] + " / " + suffixes[0]
        utc = [utcs[0][ii] for ii in idx]
        plt.plot(utc, stats[ss]['xFWHM'] / stats[0]['xFWHM'][idx], marker='o', color=c, linestyle='none', label=label)
        plt.plot(utc, stats[ss]['yFWHM'] / stats[0]['yFWHM'][idx], marker='^', color=c, linestyle='none', label=label)
    plt.gca().xaxis.set_major_formatter(time_fmt)
    plt.xticks(rotation=35)
    plt.xlabel('UTC Time (hr)')
    plt.ylabel('Ratio of Gaussian-Fit FWHM')
    plt.legend(numpoints=1, fontsize=10)
    plt.axhline(1.0, color='k', linestyle='--', linewidth=2)
    plt.ylim(0, 1.3)
    plt.savefig(plots_dir + 'xyfwhm_vs_time' + suffix + '.png')
    
    return

def compare_fwhm(stats_files, out_dir):
    #Makes a figure of several plots comparing different measures of FWHM
    #(empirical, gaussian, NEA)
    
    for ss in range(len(stats_files)):
        stats = Table.read(stats_files[ss])
        
        FWHM = stats['FWHM']
        xFWHM = stats['xFWHM']
        yFWHM = stats['yFWHM']
        NEA_FWHM = nea_to_fwhm(stats['NEA'])
        emp_FWHM = stats['emp_fwhm']
        
        # Make plot
        
        # Get plot title from file name
        title_seg = stats_files[ss].split("_")[-1]
        title_seg = title_seg.split('.')[0]
        
        names = stats_files[ss].split("/")
        for directory in names:
            if "201" in directory:
                date = (directory[4:6]+"/"+directory[6:]+"/"+directory[0:4])
    
        plot_title = date + " FWHM Comparison: "+ title_seg
        
        # Make figure
        
        plt.figure(figsize=(15,10))
        plt.suptitle(plot_title, fontsize=24)
        
        max_val = np.amax([np.amax(FWHM), np.amax(xFWHM), np.amax(yFWHM), np.amax(NEA_FWHM), np.amax(emp_FWHM)])
        plot_edge = int(max_val+1)
        
        plt.subplot(2, 3, 1)
        plt.plot(FWHM, emp_FWHM, 'o', color='purple', alpha=0.5)
        plt.xlabel("Gaussian FWHM", fontsize=14)
        plt.ylabel("Empirical FWHM", fontsize=14)
        plt.title("Empirical v Gaussian", fontsize=18)
        plt.plot([0, plot_edge], [0, plot_edge], 'k-', alpha=0.5)
        plt.axis([0, plot_edge, 0, plot_edge])
        
        plt.subplot(2, 3, 2)
        plt.plot(xFWHM, emp_FWHM, 'bo', alpha=0.5)
        plt.plot(yFWHM, emp_FWHM, 'ro', alpha=0.5)
        plt.xlabel("Gaussian FWHM", fontsize=14)
        plt.ylabel("Empirical FWHM", fontsize=14)
        plt.title("Empirical v Gaussian", fontsize=18)
        plt.plot([0, plot_edge], [0, plot_edge], 'k-', alpha=0.5)
        plt.axis([0, plot_edge, 0, plot_edge])
        
        plt.subplot(2, 3, 3)
        plt.plot(NEA_FWHM, emp_FWHM, 'o', color='purple', alpha=0.5)
        plt.xlabel('NEA FWHM', fontsize=14)
        plt.ylabel('Empirical FWHM', fontsize=14)
        plt.title('Empirical v NEA', fontsize=18)
        plt.plot([0, plot_edge], [0, plot_edge], 'k-', alpha=0.5)
        plt.axis([0, plot_edge, 0, plot_edge])
        
        plt.subplot(2, 3, 4)
        plt.plot(FWHM, NEA_FWHM, 'o', color='purple', alpha=0.5)
        plt.xlabel("Gaussian FWHM", fontsize=14)
        plt.ylabel("NEA FWHM", fontsize=14)
        plt.title("NEA v Gaussian", fontsize=18)
        plt.plot([0, plot_edge], [0, plot_edge], 'k-', alpha=0.5)
        plt.axis([0, plot_edge, 0, plot_edge])
        
        plt.subplot(2, 3, 5)
        plt.plot(xFWHM, NEA_FWHM, 'bo', alpha=0.5)
        plt.plot(yFWHM, NEA_FWHM, 'ro', alpha=0.5)
        plt.xlabel("Gaussian FWHM", fontsize=14)
        plt.ylabel("NEA FWHM", fontsize=14)
        plt.title("NEA v Gaussian", fontsize=18)
        plt.plot([0, plot_edge], [0, plot_edge], 'k-', alpha=0.5)
        plt.axis([0, plot_edge, 0, plot_edge])
        
        plt.subplot(2, 3, 6)
        plt.plot(plot_edge+1, plot_edge+1, 'bo', markersize=10, label='x FWHM')
        plt.plot(plot_edge+1, plot_edge+1, 'ro', markersize=10, label='y FWHM')
        plt.plot(plot_edge+1, plot_edge+1, 'o', markersize=10, color='purple', label='Average')
        plt.legend(loc=10, frameon=False, fontsize=18)
        plt.axis([0, plot_edge, 0, plot_edge])
        plt.tick_params(axis='both', which='both', bottom='off', top='off', left='off', right='off',
                        labelbottom='off', labelleft='off') # labels along the bottom edge are off
            
        plt.subplots_adjust(hspace=0.25, wspace=0.25)
                        
        out_file = stats_files[ss].split("/")[-1].replace('.fits', '_fwhm_comp.png')
        plt.savefig(out_dir + out_file)

    return

def plot_stats_mdp(date, suffixes=['open', 'ttf', 'closed'], out_suffix='', root_dir='/Users/dorafohring/Desktop/imaka/data/'):
    """
        Make a suite of standard plots for the stats on a given night.
        Parameters
        ----------
        date : str
        The date string for which to plot up the stats (i.e. '20170113').
        Optional Parameters
        -------------------
        suffixes : numpy array of strings
        stats files have the name stats_<suffixes[0]>.fits, etc.
        root_dir : str
        The root directory for the <date> observing run directories. The
        stats files will be searched for in:
        <root_dir>/<date>/fli/reduce/stats/
        """
    stats_dir = root_dir + date + '/fli/reduce/stats/'
    plots_dir = root_dir + date + '/fli/reduce/plots/'
    
    util.mkdir(plots_dir)
    
    colors = get_color_list()
    
    stats = []
    utc = []
    all_utc = None
    all_dimm = None
    all_mass = None
    
    for ss in range(len(suffixes)):
        suffix = suffixes[ss]
        
        st = Table.read(stats_dir + 'stats_' + suffix + '.fits')
        stats.append(st)
        
        utc_dt = [datetime.strptime(st['TIME_UTC'][ii], '%I:%M:%S') for ii in range(len(st))]
        utc.append(utc_dt)
        
        # combine MASS/DIMM for all tables.
        if ss == 0:
            all_utc = utc_dt
            all_dimm = st['DIMM']
            all_mass = st['MASS']
        else:
            all_utc = all_utc + utc_dt
            all_dimm = np.concatenate([all_dimm, st['DIMM']])
            all_mass = np.concatenate([all_mass, st['MASS']])

    scale = 0.04
    
    time_fmt = mp_dates.DateFormatter('%H:%M')
    time_loc = ticker.MaxNLocator(nbins=6)
    
    plt.figure(1, figsize=(6, 6))
    plt.clf()
    
    for ii in range(len(suffixes)):
        c = np.take(colors, ii, mode='wrap')
        plt.plot(utc[ii], stats[ii]['emp_fwhm']*stats[ii]['BINFAC']*scale, marker='o', color=c, linestyle='none', label=suffixes[ii])
    
    plt.plot(all_utc, all_dimm, marker='x', color='fuchsia', linestyle='none', label='DIMM')
    plt.plot(all_utc, all_mass, marker='+', color='dodgerblue', linestyle='none', label='MASS')

    plt.gca().xaxis.set_major_formatter(time_fmt)
    plt.gca().xaxis.set_major_locator(time_loc)
    plt.xticks(rotation=35)
    plt.xlabel('UTC Time (hr)')
    plt.ylabel('Empirical FWHM (")')
    plt.legend(numpoints=1, fontsize=10)
    plt.ylim(0, 2.0)
    plt.title(date)
    plt.savefig(plots_dir + 'mdp_efwhm_vs_time' + out_suffix + '.png')
    
    plt.figure(2, figsize=(6, 6))
    plt.clf()
    for ii in range(len(suffixes)):
        plt.plot(utc[ii], stats[ii]['NEA'], marker='o', linestyle='none', label=suffixes[ii])
    plt.gca().xaxis.set_major_formatter(time_fmt)
    plt.gca().xaxis.set_major_locator(time_loc)
    plt.xticks(rotation=35)
    plt.xlabel('UTC Time (hr)')
    plt.ylabel('NEA (Sq. Arcsec)')
    plt.legend(numpoints=1, fontsize=10)
    plt.ylim(0, 5)
    plt.title(date)
    plt.savefig(plots_dir + 'mdp_nea_vs_time' + out_suffix + '.png')

    plt.figure(3, figsize=(6, 6))
    plt.clf()
    for ii in range(len(suffixes)):
        plt.plot(stats[ii]['MASS'], stats[ii]['emp_fwhm']*stats[ii]['BINFAC']*scale, marker='o', linestyle='none', label=suffixes[ii])
    plt.xlabel('MASS Seeing (")')
    plt.ylabel('Empirical FWHM (")')
    plt.legend()
    plt.ylim(0, 1.5)
    plt.title(date)
    plt.savefig(plots_dir + 'mdp_mass_vs_efwhm' + out_suffix + '.png')
    
    plt.figure(4, figsize=(6, 6))
    plt.clf()
    for ii in range(len(suffixes)):
        plt.plot(stats[ii]['MASS'], stats[ii]['NEA'], marker='o', linestyle='none', label=suffixes[ii])
    plt.xlabel('MASS Seeing (")')
    plt.ylabel('NEQ (Sq. Arcsec)')
    plt.legend()
    plt.ylim(0, 5)
    plt.title(date)
    plt.savefig(plots_dir + 'mdp_mass_vs_nea' + out_suffix + '.png')

    plt.figure(5, figsize=(6, 6))
    plt.clf()
    for ii in range(len(suffixes)):
        plt.plot(stats[ii]['DIMM'], stats[ii]['emp_fwhm']*stats[ii]['BINFAC']*scale, marker='o', linestyle='none', label=suffixes[ii])
    plt.xlabel('DIMM Seeing (")')
    plt.ylabel('Empirical FWHM (")')
    plt.legend()
    plt.ylim(0, 2.0)
    plt.title(date)
    plt.savefig(plots_dir + 'mdp_dimm_vs_efwhm' + out_suffix + '.png')
    
    plt.figure(6, figsize=(6, 6))
    plt.clf()
    for ii in range(len(suffixes)):
        plt.plot(stats[ii]['DIMM'], stats[ii]['NEA'], marker='o', linestyle='none', label=suffixes[ii])
    plt.xlabel('DIMM Seeing (")')
    plt.ylabel('NEQ (Sq. Arcsec)')
    plt.legend()
    plt.ylim(0, 5)
    plt.title(date)
    plt.savefig(plots_dir + 'mdp_dimm_vs_nea' + out_suffix + '.png')

    return

def plot_profile(date, suffixes=['open', 'ttf', 'closed'], out_suffix='', root_dir='/Users/dorafohring/Desktop/imaka/data/'):
    """
        Make a suite of standard plots for the stats on a given night.
        Parameters
        ----------
        date : str
        The date string for which to plot up the stats (i.e. '20170113').
        Optional Parameters
        -------------------
        suffixes : numpy array of strings
        stats files have the name stats_<suffixes[0]>.fits, etc.
        root_dir : str
        The root directory for the <date> observing run directories. The
        stats files will be searched for in:
        <root_dir>/<date>/fli/reduce/stats/
        """
    latexParams = {
        'figure.dpi'      :150,
            #'figure.figsize':[3.32,2.49],
            'figure.figsize'  : [3.3, 4.2],
                'axes.labelsize'  : 10,
                'xtick.labelsize' : 10,
                'ytick.labelsize' : 10,
                'font.family'     : 'serif',
                'font.serif'      : 'Computer Modern Roman',
                'text.usetex'     : True,
                'figure.subplot.top':0.90,
                'figure.subplot.right':0.95,
                'figure.subplot.bottom':0.15,
                'figure.subplot.left':0.15,
                'lines.linewidth':1.5,
                'lines.markersize':3
            }

    plt.rcParams.update(latexParams)
    
    stats_dir = root_dir + date + '/fli/reduce/stats/'
    plots_dir = root_dir + date + '/fli/reduce/plots/'
    
    util.mkdir(plots_dir)
    
    stats = []
    
    for suffix in suffixes:
        st = Table.read(stats_dir + 'stats_' + suffix + '_mdp.fits')
        stats.append(st)
    
    altitudes = ['Cn2dh_005', 'Cn2dh_010', 'Cn2dh_020', 'Cn2dh_040', 'Cn2dh_080', 'Cn2dh_160']
    
    allprofs = []
    for altitude in altitudes:
        combine = []
        for ii in range(len(suffixes)):
            combine.extend(stats[ii][altitude])
        allprofs.append(combine)
    gl = []
    for jj in range(len(suffixes)):
        gl.extend(stats[jj]['DIMM'] - stats[jj]['MASS'])

    allprofs = np.array(allprofs)
    mprofs = np.ma.masked_invalid(allprofs)
    profave = np.mean(mprofs, axis=1)

## sets zero and negative numbers to a very small value
    gl = np.array(gl)
    gindex = gl < 0.01
    gl[gindex] = 0.01
    
    cngl = (0.98*0.00000055*206265/gl)**(-5/3) / (16.7*(0.00000055)**(-2))
    profave = np.insert(profave, obj=0, values=np.mean(cngl))
    
    hlis = [0, 0.5, 1, 2, 4, 8, 16]
    
    plt.figure(1)
    plt.clf()
    plt.plot(profave, hlis)
    plt.xlabel(r'$C_n^2$ dh ($m^{1/3}$)')
    plt.ylabel(r'h (km)')
    plt.title(date, fontsize=12)
    plt.savefig(plots_dir + 'mass_profile' + out_suffix + '.png')
    plt.show()
    
    return

def plot_all_profiles(dates, root_dir='/Users/dorafohring/Desktop/imaka/data/'):
    """
        Make a suite of standard plots for the stats on a given night.
        Parameters
        ----------
        dates : str
        The date strings for which to plot up the stats (i.e. ['20170113', '20170114']).
        Optional Parameters
        -------------------
        root_dir : str
        The root directory for the <date> observing run directories. The
        stats files will be searched for in:
        <root_dir>/<date>/fli/reduce/stats/
    """
    savedatadir = '/Users/dorafohring/Desktop/imaka/data/'
    stats = []
    
    for date in dates:
        dir = root_dir + date + '/fli/reduce/stats/'
        mdpfiles = [item for item in os.listdir(dir) if item.endswith('mdp.fits')]
        for file in mdpfiles:
            st = Table.read(dir + file)
            stats.append(st)

    altitudes = ['Cn2dh_005', 'Cn2dh_010', 'Cn2dh_020', 'Cn2dh_040', 'Cn2dh_080', 'Cn2dh_160']
    
    allprofs = []
    for altitude in altitudes:
        combine = []
        for ii in range(len(st)):
            combine.extend(stats[ii][altitude])
        allprofs.append(combine)
    gl = []
    for jj in range(len(st)):
        gl.extend(stats[jj]['DIMM'] - stats[jj]['MASS'])
    allprofs = np.array(allprofs)
    
    ## sets zero and negative numbers to a very small value
    gl = np.array(gl)
    gindex = gl < 0.01
    gl[gindex] = 0.01
    
    cngl = (0.98*0.00000055*206265/gl)**(-5/3) / (16.7*(0.00000055)**(-2))

    glprofs = np.zeros([7,len(gl)])
    glprofs[0,:]  = cngl
    glprofs[1:,:] = allprofs

    mprofs = np.ma.masked_invalid(glprofs)

    np.savetxt(savedatadir+'allprofs.txt', glprofs)

    profave = np.mean(mprofs, axis=1)

    hlis = [0, 0.5, 1, 2, 4, 8, 16]
    plt.figure(1)
    plt.clf()
    plt.plot(profave, hlis)
    plt.xlabel(r'$C_n^2$ dh ($m^{1/3}$)')
    plt.ylabel(r'h (km)')
    plt.title('Average profile over all nights', fontsize=12)
    plt.savefig(root_dir + 'ave_profile.png')
    plt.show()
    
    return

def plot_best_stats(date, suffixes=['open', 'closed'], out_suffix='', root_dir=''):
    
    #plots different FWHMs as function of UTC time
    
    colors = ['b', 'g', 'r', 'c', 'm', 'k', 'yellow', 'purple', 'orange']
    
    stats_dir = root_dir + date + '/FLI/reduce/stats/'
    plots_dir = root_dir + date + '/FLI/reduce/plots/'
    
    scale = 0.08
    
    stats = []
    utcs = []
    max_fwhm = 0.0
    for suffix in suffixes:
        ss = Table.read(stats_dir + 'stats_' + suffix + '.fits')
        stats.append( ss )
        utc_dt = [datetime.strptime(ss['TIME_UTC'][ii], '%H:%M:%S') for ii in range(len(ss))]
        utcs.append(utc_dt)
        
        # Add NEA FWHM to table (temporarily)
        ss['NEA_FWHM'] = nea_to_fwhm(ss['NEA'])
        
        # Get the max value of all the FWHMs
        max_fwhm_tmp = np.max([ss['NEA_FWHM'], ss['emp_fwhm'], ss['xFWHM'], ss['yFWHM']])
        if max_fwhm_tmp > max_fwhm:
            max_fwhm = max_fwhm_tmp


    time_fmt = mp_dates.DateFormatter('%H:%M')
    time_loc = ticker.MaxNLocator(nbins=6)
    
    #####
    # NEA FWHM PLOT
    #####
    plt.figure()
    for ii in range(len(stats)):
        # Convert NEA to FWHM using NEA = pi * (FWHM*scale/2)**2
        nea_fwhm = nea_to_fwhm(stats[ii]['NEA'])
        plt.plot(utcs[ii], nea_fwhm, color=colors[ii], marker='o', label=suffixes[ii],
                 alpha=0.5, linestyle='none')

    plt.gca().xaxis.set_major_formatter(time_fmt)
    plt.gca().xaxis.set_major_locator(time_loc)
    plt.legend(loc=2, bbox_to_anchor=(1, 1), numpoints=1)
    plt.xticks(rotation=35)
    #plt.ylim(0, 2)
    plt.xlabel("UTC Time")
    plt.ylabel('NEA FWHM (")')
    plt.title(date)

#####
# EMPIRICAL FWHM PLOT
#####
    plt.figure()
    for ii in range(len(stats)):
        plt.plot(utcs[ii], stats[ii]['emp_fwhm']*scale, color=colors[ii], marker='o', label=suffixes[ii],
                 alpha=0.5, linestyle='none')

    plt.gca().xaxis.set_major_formatter(time_fmt)
    plt.gca().xaxis.set_major_locator(time_loc)
    plt.legend(loc=2, bbox_to_anchor=(1, 1), numpoints=1)
    plt.xticks(rotation=35)
    #plt.ylim(0, 2)
    plt.xlabel("UTC Time")
    plt.ylabel('Empirical FWHM (")')
    plt.title(date)


#####
# GAUSSIAN FWHM PLOT
#####
    plt.figure()
    for ii in range(len(stats)):
        plt.plot(utcs[ii], stats[ii]['xFWHM']*scale, color=colors[ii], label='X ' + suffixes[ii],
                 marker="o", alpha=0.5, linestyle='none')
        plt.plot(utcs[ii], stats[ii]['yFWHM']*scale, color=colors[ii], label='Y ' + suffixes[ii],
                          marker="^", alpha=0.5, linestyle='none')
    plt.gca().xaxis.set_major_formatter(time_fmt)
    plt.gca().xaxis.set_major_locator(time_loc)
    plt.legend(loc=2, bbox_to_anchor=(1, 1), numpoints=1)
    plt.xticks(rotation=35)
    #plt.ylim(0, 2)
    plt.xlabel("UTC Time")
    plt.ylabel('Gaussian FWHM (")')
    plt.title(date)
    plt.savefig(plots_dir + 'fwhm_vs_frame' + suffix + '.png')


    return


def get_color_list():
    rcParams = plt.matplotlib.rcParams
    prop_cycler = rcParams['axes.prop_cycle']
    if prop_cycler is None and 'axes.color_cycle' in rcParams:
        clist = rcParams['axes.color_cycle']
        prop_cycler = cycler('color', clist)
    
    colors = [item['color'] for item in prop_cycler]
    return colors


def nea_to_fwhm(nea):
    """
        Convert the Noise Equivalent Area using the following relation (King 1983, King 1971):
        \frac{b}{ \int \phi^2 dA } = 30.8 \sigma^2 b
        which is appropriate for a seeing-limited PSF (and probably GLAO as well. For a
        gaussian, the factor would be 4 * pi; but real PSFs have extended tails.
        The noise-equivalent area is:
        NEA = 1.0 / \int \phi^2 dA = 1.0 / \sum (f_i - b)^2
        """
    sigma = np.sqrt( nea / 30.8 )
    
    fwhm = 2.355 * sigma
    
    return fwhm

def plot_fwhmvt(open_file, closed_file, comp_col, title, plots_dir):
    
    #open_file and closed_file are the fits stats files
    #plots fwhm for open and closed and mass/dimm seeing vs time
    #'comp_col' is what column of data to compare (e.g. 'emp_fwhm')
    #obs_wav is the wavelength of that observation, to scale it to the 500 nm mass/dimm observations
    #title: title for generated plot
    #plots_dir: directory to put generated plot in
    
    #Read in data
    stats1 = Table.read(open_file)
    stats2 = Table.read(closed_file)
    
    #Match open and closed data in time
    time, date, data1, data2, err1, err2 = add_data.match_cols(open_file, closed_file, comp_col)
    calib = []
    
    #Get mass/dimm data
    if len(stats1) < len(stats2):
        mass = stats1['MASS']
        dimm = stats1['DIMM']
        for i in range(len(stats1)):
            wvln = filter2wv(stats1['FILTER'][i])
            scale = 0.04 * stats1['BINFAC'][i]
            factor = ((500/wvln)**0.2) * scale
            calib.append(factor)

    else:
        mass = stats2['MASS']
        dimm = stats2['DIMM']
        for i in range(len(stats2)):
            wvln = filter2wv(stats2['FILTER'][i])
            scale = 0.04 * stats2['BINFAC'][i]
            factor = ((500/wvln)**0.2) * scale
            calib.append(factor)

    open_err = err1 * calib
    closed_err = err2 * calib

    
    #Plot fwhm and seeing vs time
    times = []
    for i in range(len(time)):
        string = str(date[i])+str(time[i])
        dt_obj = datetime.strptime(string, "b'%Y-%m-%d'b'%H:%M:%S'")
        times.append(dt_obj)
    print("times", len(times)) #####

    plt.figure(1, figsize=(12, 6))
    plt.errorbar(times, data1*calib, yerr=open_err, fmt='o', label="Open")
    plt.errorbar(times, data2*calib, yerr=closed_err, fmt='ro', label="Closed")
    plt.plot(times, dimm, 'b-')
    plt.plot(times, mass, 'r-')
    plt.ylabel(comp_col)
    plt.xlabel("UTC Time")
    plt.xticks(rotation=35)
    plt.ylim(0, 2.5)
    plt.title(title)
    plt.gca().xaxis.set_major_formatter(mp_dates.DateFormatter('%H:%M'))
    plt.legend()
    
    plt.savefig(plots_dir + 'fwhm_v_time' + '.png')
    return



def plot_EE(labels, data_dir_root, stats_dir_end):
    
    plt.figure(1, figsize=(12, 4))
    plt.title('Main Title')
    plt.subplots_adjust(left=0.1, bottom=0.15)
    plt.subplot(131)
    
    all_data = []
    for item in labels:
        
        root_dir = data_dir_root + item[0] + stats_dir_end
        open_file = root_dir + "stats_open_mdp.fits"
        closed_file = root_dir + "stats_"+item[1]+"_mdp.fits"
        open_data = Table.read(open_file)
        closed_data = Table.read(closed_file)
        time, date, data1, data2, err1, err2 = add_data.match_cols(open_file, closed_file, 'EE80')
        
        if len(open_data) >= len(closed_data):
            co_rat = data2/data1 #closed to open ratio of 80% EE
            plt.plot(closed_data['DIMM']-closed_data['MASS'], co_rat, 'o', alpha=.5, label=item[0])
        else:
            co_rat = data2/data1 #closed to open ratio of 80% EE
            plt.plot(open_data['DIMM']-open_data['MASS'], co_rat, 'o',  alpha=.5, label=item[0])
        all_data.append(co_rat)
    con = np.concatenate(all_data)
    mean = np.mean(con)
    plt.plot([0,1.2],[1,1], 'k--')
    plt.xlabel('DIMM-MASS (as)')
    plt.ylabel('Closed EE / Open EE')
    plt.title('80EE')
    plt.axis([0, 1, 0, 1.25])
    plt.suptitle('Encircled Energy Profile')
    plt.plot([0,1.2],[mean,mean], '--', color='Gray')

#plt.legend(loc=4)

#####


    plt.subplot(132)
    all_data = []
    for item in labels:
        
        root_dir = data_dir_root + item[0] + stats_dir_end
        open_file = root_dir + "stats_open_mdp.fits"
        closed_file = root_dir + "stats_"+item[1]+"_mdp.fits"
        open_data = Table.read(open_file)
        closed_data = Table.read(closed_file)
        time, date, data1, data2, err1, err2= add_data.match_cols(open_file, closed_file, 'EE50')
        
        if len(open_data) >= len(closed_data):
            co_rat = data2/data1 #closed to open ratio of 80% EE
            plt.plot(closed_data['DIMM']-closed_data['MASS'], co_rat, 'o',  alpha=.5, label=item[0])
        else:
            co_rat = data2/data1 #closed to open ratio of 80% EE
            plt.plot(open_data['DIMM']-open_data['MASS'], co_rat, 'o',  alpha=.5, label=item[0])
        all_data.append(co_rat)
    con = np.concatenate(all_data)
    mean = np.mean(con)
    plt.plot([0,1.2],[1,1], 'k--')
    plt.xlabel('DIMM-MASS (as)')
    #plt.ylabel('Closed EE50 / Open EE50')
    plt.title('50EE')
    plt.axis([0, 1, 0, 1.25])
    plt.plot([0,1.2],[mean,mean], '--', color='Gray')
    
    ####
    
    plt.subplot(133)
    all_data = []
    for item in labels:
        
        root_dir = data_dir_root + item[0] + stats_dir_end
        open_file = root_dir + "stats_open_mdp.fits"
        closed_file = root_dir + "stats_"+item[1]+"_mdp.fits"
        open_data = Table.read(open_file)
        closed_data = Table.read(closed_file)
        time, date, data1, data2, err1, err2 = add_data.match_cols(open_file, closed_file, 'EE25')
        
        if len(open_data) >= len(closed_data):
            co_rat = data2/data1 #closed to open ratio of 80% EE
            plt.plot(closed_data['DIMM']-closed_data['MASS'], co_rat, 'o',  alpha=.5, label=item[0])
        else:
            co_rat = data2/data1 #closed to open ratio of 80% EE
            plt.plot(open_data['MASS']-open_data['MASS'], co_rat, 'o',  alpha=.5, label=item[0])
        all_data.append(co_rat)
    con = np.concatenate(all_data)
    mean = np.mean(con)
    plt.plot([0,1.2],[1,1], 'k--', Label = '1')
    plt.xlabel('DIMM-MASS (as)')
    #plt.ylabel('Closed EE25 / Open EE25')
    plt.title('25EE')
    plt.axis([0, 1, 0, 1.25])
    plt.plot([0,1.2],[mean,mean], '--', color='Gray', label='Mean')
    plt.legend(bbox_to_anchor=(1.5, 1))
    
    return


def plot_week_fwhm(labels, data_dir_root, stats_dir_end, title):
    
    plt.figure(1, figsize=(8, 8))
    scale = 0.04
    for day in labels:
        open_file = data_dir_root+day[0]+stats_dir_end + "stats_open_mdp.fits"
        open_data = Table.read(open_file)
        o_data = np.array(open_data['emp_fwhm'])
        o_binfac = np.array(open_data['BINFAC'])
        o_filt =  filter2wv(np.array(open_data['FILTER']))
        open_fin = o_data * scale* o_binfac * (500/o_filt)**(1/5)
        DIMM = np.array(open_data['DIMM'])
        
        closed_file = data_dir_root +day[0]+stats_dir_end + "stats_"+day[1]+"_mdp.fits"
        closed_data = Table.read(closed_file)
        c_data = np.array(closed_data['emp_fwhm'])
        c_binfac = np.array(closed_data['BINFAC'])
        c_filt =  filter2wv(np.array(closed_data['FILTER']))
        closed_fin = c_data * scale* c_binfac * (500/c_filt)**(1/5)
        MASS = np.array(closed_data['MASS'])
        
        if day == labels[0]:
            plt.plot(DIMM, open_fin, 'bo', label='Open vs DIMM')
            plt.plot(MASS, closed_fin, 'ro', label='Closed vs MASS')
        else:
            plt.plot(DIMM, open_fin, 'bo')
            plt.plot(MASS, closed_fin, 'ro')

    plt.plot([0, 2], [0, 2], 'k-')
    plt.xlabel('Seeing (as)')
    plt.ylabel('Empirical FWHM (as)')
    plt.title(title)
    plt.legend(loc=2)
    
    return


def plot_hist(labels, data_dir_root, stats_dir_end, title):
    
    scale=0.08
    
    open_files = []
    closed_files = []
    for night in labels:
        open_file = data_dir_root + night[0] + stats_dir_end + "stats_open_mdp.fits"
        closed_file = data_dir_root + night[0] + stats_dir_end + "stats_" + night[1] + "_mdp.fits"
        open_files.append(open_file)
        closed_files.append(closed_file)
    
    open_tables = []
    closed_tables = []
    for ii in range(len(labels)):
        open_table = np.array(Table.read(open_files[ii])['emp_fwhm']*((500/labels[ii][2])**(1/5))*scale)
        closed_table = np.array(Table.read(closed_files[ii])['emp_fwhm']*((500/labels[ii][2])**(1/5))*scale)
        open_tables.append(open_table)
        closed_tables.append(closed_table)



    all_fwhm_open = np.concatenate(open_tables, axis=0)
    all_fwhm_closed = np.concatenate(closed_tables, axis=0)
    max_val = np.amax(all_fwhm_open)
    bins_set = np.arange(0, max_val, max_val/30)


    plt.figure(2, figsize =(10, 6))
    plt.hist(all_fwhm_open, bins=bins_set, alpha=0.5, color='blue', label='Open Loop')
    plt.hist(all_fwhm_closed, alpha=0.5, color='red', label='Closed Loop')
    plt.legend(loc=1)
    plt.xlabel('Empirical FWHM')
    plt.title(title)
    
    return


def filter2wv(filter_label):
    
    #converts filter label from fits header into wavelength in nanometers
    
    if type(filter_label)==str:
        if filter_label == "I":
            return 806
        elif filter_label == "R":
            return 658
        elif filter_label == "1_micronlp":
            return 1000
        else:
            print("Filter not found: defaulting to 500 nm")
            return 500
    else:
        new_array = np.zeros(len(filter_label))
        for i in range(len(filter_label)):
            if filter_label[i] == "I":
                new_array[i] = 806
            elif filter_label[i] == "R":
                new_array[i] = 658
            elif filter_label == "1_micronlp":
                rew_array[i] = 1000
            else:
                print("Filter not found: defaulting to 500 nm")
                new_array[i] = 500
        return new_array




def plot_field_var(starlist):
    
    mode=starlist.split("_")[-2]
    x_cents, y_cents, fwhm, x_fwhm, y_fwhm, roundness, dist = add_data.read_starlist(starlist)
    plt.figure(1, figsize=(12, 8))
    plt.suptitle('Field Variability', fontsize=18)
    
    #calculate elongation of psfs
    elon = []
    for i in range(len(fwhm)):
        if x_fwhm[i]>y_fwhm[i]:
            elon.append(y_fwhm[i]/x_fwhm[i])
        else:
            elon.append(x_fwhm[i]/y_fwhm[i])

    #calculate mean and std of psf parameters
    mean_size = np.mean(fwhm)
    std_size = np.std(fwhm)
    mean_elon = np.mean(elon)
    std_elon = np.std(elon)

    plt.subplot(221)
    plt.plot(dist, y_fwhm,'o', alpha=0.5)
    plt.xlabel('Distance from center of image (pixels)')
    plt.ylabel('Gaussian fit FWHM (not scaled)')
    plt.title('PSF Size')
    
    plt.subplot(222)
    plt.plot(dist, elon,'o', alpha=0.5);
    plt.xlabel('Distance from center of image (pixels)')
    plt.ylabel('PSF Elongation (x_fwhm/y_fwhm)')
    plt.title('PSF Elongation')
    
    plt.subplot(223)
    plt.scatter(x_cents, y_cents, c=fwhm, vmin = mean_size-(std_size*2), vmax = mean_size+(std_size*2))
    plt.colorbar()
    plt.xlabel('x coordinate')
    plt.ylabel('y coordinate')
    
    plt.subplot(224)
    plt.scatter(x_cents, y_cents, c=elon, vmin = mean_elon-(std_elon*2), vmax = mean_elon+(std_elon*2))
    plt.colorbar()
    plt.xlabel('x coordinate')
    plt.ylabel('y coordinate')
    
    return



def plot_fwhmvt_nomatch(open_file, closed_file, comp_col, title, plots_dir):
    
    #open_file and closed_file are the fits stats files
    #plots fwhm for open and closed and mass/dimm seeing vs time
    #'comp_col' is what column of data to compare (e.g. 'emp_fwhm')
    #obs_wav is the wavelength of that observation, to scale it to the 500 nm mass/dimm observations
    #title: title for generated plot
    #plots_dir: directory to put generated plot in
    
    #Read in data
    comp_col = 'emp_fwhm'
    #Read in data
    stats1 = Table.read(open_file)
    stats2 = Table.read(closed_file)
    
    time1, date1, data1, data2, err1, err2 = add_data.match_cols(open_file, open_file, comp_col)
    time2, date2, data1, data2, err1, err2 = add_data.match_cols(closed_file, closed_file, comp_col)
    
    
    
    calib1 = []; calib2=[]
    
    #Get mass/dimm data
    mass = stats2['MASS']
    dimm = stats1['DIMM']
    for i in range(len(stats1)):
        wvln = filter2wv(stats1['FILTER'][i])
        scale = .04 * stats1['BINFAC'][i]
        factor = ((500/wvln)**0.2) * scale
        calib1.append(factor)
    
    for i in range(len(stats2)):
        wvln = filter2wv(stats2['FILTER'][i])
        scale = .04 * stats2['BINFAC'][i]
        factor = ((500/wvln)**0.2) * scale
        calib2.append(factor)


    open_err = np.array(stats1['emp_fwhm_std']) * np.array(calib1)[:,0]
    closed_err = np.array(stats2['emp_fwhm_std']) * np.array(calib2)[:,0]

#Plot fwhm and seeing vs time
    times1 = []
    for i in range(len(time1)):
        string = str(date1[i])+str(time1[i])
        dt_obj = datetime.strptime(string, "b'%Y-%m-%d'b'%H:%M:%S'")
        times1.append(dt_obj)

    times2 = []
    for i in range(len(time2)):
        string = str(date2[i])+str(time2[i])
        dt_obj = datetime.strptime(string, "b'%Y-%m-%d'b'%H:%M:%S'")
        times2.append(dt_obj)

    plt.figure(1, figsize=(12, 4))
    plt.errorbar(times1, np.array(stats1['emp_fwhm'])*np.array(calib1)[:,0], yerr=open_err, fmt='o', label="Open")
    plt.errorbar(times2, np.array(stats2['emp_fwhm'])*np.array(calib2)[:,0], yerr=closed_err, fmt='ro', label="Closed")
    plt.plot(times1, dimm, 'b-')
    plt.plot(times2, mass, 'r-')
    plt.ylabel(comp_col)
    plt.xlabel("UTC Time")
    plt.xticks(rotation=35)
    plt.ylim(0, 2)
    plt.title(title)
    plt.gca().xaxis.set_major_formatter(mp_dates.DateFormatter('%H:%M'))
    plt.legend()
    
    plt.savefig(plots_dir + 'fwhm_v_time' + '.png')
    return



def compare_seeing(date):
    
    #compares Olivier's seeing data for integrated seeing and free atmosphere to Mauna Kea's MASS/DIMM measurments.  only date string needed, but assumes organization of files like on onaga, with massdimm and stats files available.
    
    if date[5] == "5":
        run = "5"
    elif date[5] == "2":
        run = "4"
    elif date[5] == "1":
        run = "3"
    stat_dir = "/Users/fatimaabdurrahman/Desktop/Research/RUN"+run+"/"+date+"/FLI/reduce/stats/"
    md_dir = "/Users/fatimaabdurrahman/Desktop/Research/RUN"+run+"/"+date+"/FLI/reduce/massdimm/"
    alt_file = stat_dir + "profile-data_"+date+"-noTT.fits"
    dimm_file = md_dir +date+".dimm.dat"
    mass_file = md_dir +date+".mass.dat"

    table = read_csv(dimm_file, delim_whitespace=True, names= \
                 ['year', 'month', 'day', 'hour', 'minute', 'second', 'seeing'])
    
    year  = np.array(table['year'])
    month = np.array(table['month'])
    day   = np.array(table['day'])
    hour   = np.array(table['hour'])
    minute = np.array(table['minute'])
    second = np.array(table['second'])
    dimm_seeing = np.array(table['seeing'])
    hour += 10
    idx = np.where(hour >=24)[0]
    day[idx] += 1
    hour[idx] -= 24
    timeInHours_dimm = hour + (minute/60.0) + (second/3600.0)
    
    
    table = read_csv(mass_file, delim_whitespace=True, names= \
                     ['year', 'month', 'day', 'hour', 'minute', 'second', 'seeing'])
        
    year  = np.array(table['year'])
    month = np.array(table['month'])
    day   = np.array(table['day'])
    hour   = np.array(table['hour'])
    minute = np.array(table['minute'])
    second = np.array(table['second'])
    mass_seeing = np.array(table['seeing'])
    hour += 10
    idx = np.where(hour >=24)[0]
    day[idx] += 1
    hour[idx] -= 24
    timeInHours_mass = hour + (minute/60.0) + (second/3600.0)
                     
    data = fits.getdata(alt_file)
                     
                     #blues
    f, (ax1, ax2) = plt.subplots(2, figsize=(15,5), sharex=True, sharey=True)
    ax1.plot(timeInHours_dimm, dimm_seeing, '.-', color='#00bfff', markersize=10, label='DIMM')
    ax1.plot(data[:,0], data[:,1], '.-', color='Navy', markersize=10, label='Our Data')
    ax1.axis(plt.axis([np.amin(data[:,0])-((np.amax(data[:,0]-np.amin(data[:,0])))*0.1), np.amax(data[:,0])+((np.amax(data[:,0]-np.amin(data[:,0])))*0.1),0, np.amax(data[:,1])+0.2]))
    ax1.set_title(date, fontsize=15)
    ax1.set_ylabel('Integrated Seeing')
                     
                     #red
    ax2.plot(timeInHours_mass, mass_seeing, '.-', color='Tomato', markersize=10, label='MASS')
    ax2.plot(data[:,0], data[:,-1], '.-', color='Firebrick', markersize=10, label='Our Data')
    ax2.set_ylabel('Free Atmosphere Seeing')
    ax2.set_xlabel('UT Time')
                     
                     # Fine-tune figure; make subplots close to each other and hide x ticks for
                     # all but bottom plot.
    f.subplots_adjust(hspace=0)
    plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)
                     
    return
