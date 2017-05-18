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
            scale = 4 * stats1['BINFAC'][i]
            factor = ((500/wvln)**0.2) / scale
            calib.append(factor)

    else:
        mass = stats2['MASS']
        dimm = stats2['DIMM']
        for i in range(len(stats2)):
            wvln = filter2wv(stats2['FILTER'][i])
            scale = 4 * stats2['BINFAC'][i]
            factor = ((500/wvln)**0.2) / scale
            calib.append(factor)
                
    open_err = err1 * calib
    closed_err = err2 * calib

    #Plot fwhm and seeing vs time
    times = []
    for i in range(len(time)):
        string = str(date[i])+str(time[i])
        dt_obj = datetime.strptime(string, "b'%Y-%m-%d'b'%H:%M:%S'")
        times.append(dt_obj)

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

    for item in labels:

        root_dir = data_dir_root + item[0] + stats_dir_end
        open_file = root_dir + "stats_open_mdp.fits"
        closed_file = root_dir + "stats_"+item[1]+"_mdp.fits"
        open_data = Table.read(open_file)
        closed_data = Table.read(closed_file)
        time, date, data1, data2, err1, err2 = add_data.match_cols(open_file, closed_file, 'EE80')

        if len(open_data) >= len(closed_data):
            co_rat = data2/data1 #closed to open ratio of 80% EE
            plt.plot(closed_data['DIMM']-closed_data['MASS'], co_rat, 'o', label=item[0])
        else:
            co_rat = data2/data1 #closed to open ratio of 80% EE
            plt.plot(open_data['DIMM']-open_data['MASS'], co_rat, 'o', label=item[0])

    plt.plot([0,1.2],[1,1], 'k--')
    plt.xlabel('DIMM-MASS (as)')
    plt.ylabel('Closed EE / Open EE')
    plt.title('80EE')
    plt.axis([0, 2, 0, 2])
    #plt.legend(loc=4)

    #####


    plt.subplot(132)

    for item in labels:

        root_dir = data_dir_root + item[0] + stats_dir_end
        open_file = root_dir + "stats_open_mdp.fits"
        closed_file = root_dir + "stats_"+item[1]+"_mdp.fits"
        open_data = Table.read(open_file)
        closed_data = Table.read(closed_file)
        time, date, data1, data2, err1, err2= add_data.match_cols(open_file, closed_file, 'EE50')

        if len(open_data) >= len(closed_data):
            co_rat = data2/data1 #closed to open ratio of 80% EE
            plt.plot(closed_data['DIMM']-closed_data['MASS'], co_rat, 'o', label=item[0])
        else:
            co_rat = data2/data1 #closed to open ratio of 80% EE
            plt.plot(open_data['DIMM']-open_data['MASS'], co_rat, 'o', label=item[0])

    plt.plot([0,1.2],[1,1], 'k--')
    plt.xlabel('DIMM-MASS (as)')
    #plt.ylabel('Closed EE50 / Open EE50')
    plt.title('50EE')
    plt.axis([0, 2, 0, 2])
    #plt.legend(loc=1)

    ####

    plt.subplot(133)

    for item in labels:

        root_dir = data_dir_root + item[0] + stats_dir_end
        open_file = root_dir + "stats_open_mdp.fits"
        closed_file = root_dir + "stats_"+item[1]+"_mdp.fits"
        open_data = Table.read(open_file)
        closed_data = Table.read(closed_file)
        time, date, data1, data2, err1, err2 = add_data.match_cols(open_file, closed_file, 'EE25')

        if len(open_data) >= len(closed_data):
            co_rat = data2/data1 #closed to open ratio of 80% EE
            plt.plot(closed_data['DIMM']-closed_data['MASS'], co_rat, 'o', label=item[0])
        else:
            co_rat = data2/data1 #closed to open ratio of 80% EE
            plt.plot(open_data['MASS']-open_data['MASS'], co_rat, 'o', label=item[0])

    plt.plot([0,1.2],[1,1], 'k--')
    plt.xlabel('DIMM-MASS (as)')
    #plt.ylabel('Closed EE25 / Open EE25')
    plt.title('25EE')
    plt.axis([0, 2, 0, 2])
    plt.legend(loc=1)

    return



def plot_week_fwhm(labels, data_dir_root, stats_dir_end, title):
    
    #plots fwhm vs seeing (open v dimm and closed v mass) for a week of observing
    #inputs:
        #labels: a list of entries of the form [date, closed stat ending, wavelength]
                #e.g ['20170217', 'closeda', 658]
        #data_dir_root: directory before data, e.g. "/Users/astrouser/Desktop/Reserach/imaka/RUN4/"
        #stats_dir_end: location of stats file after date dir, e.g.: "/FLI/reduce/stats/"
        #title: a string for the title of the plot
    
    
    open_file1 = data_dir_root+labels[0][0]+stats_dir_end + "stats_open_mdp.fits"
    open_file2 = data_dir_root+labels[1][0]+stats_dir_end + "stats_open_mdp.fits"
    open_file3 = data_dir_root+labels[2][0]+stats_dir_end + "stats_open_mdp.fits"
    if len(labels) > 3:
        open_file4 = data_dir_root+labels[3][0]+stats_dir_end + "stats_open_mdp.fits"
    if len(labels) > 4:
        open_file5 = data_dir_root+labels[4][0]+stats_dir_end + "stats_open_mdp.fits"

    o1 = Table.read(open_file1)
    o2 = Table.read(open_file2)
    o3 = Table.read(open_file3)
    if len(labels) > 3:
        o4 = Table.read(open_file4)
    if len(labels) > 4:
        o5 = Table.read(open_file5)

    closed_file1 = data_dir_root +labels[0][0]+stats_dir_end + "stats_"+labels[0][1]+"_mdp.fits"
    closed_file2 = data_dir_root +labels[1][0]+stats_dir_end + "stats_"+labels[1][1]+"_mdp.fits"
    closed_file3 = data_dir_root +labels[2][0]+stats_dir_end + "stats_"+labels[2][1]+"_mdp.fits"
    if len(labels) > 3:
        closed_file4 = data_dir_root +labels[3][0]+stats_dir_end + "stats_"+labels[3][1]+"_mdp.fits"
    if len(labels) > 4:
        closed_file5 = data_dir_root +labels[4][0]+stats_dir_end + "stats_"+labels[4][1]+"_mdp.fits"

    c1 = Table.read(closed_file1)
    c2 = Table.read(closed_file2)
    c3 = Table.read(closed_file3)
    if len(labels) > 3:
        c4 = Table.read(closed_file4)
    if len(labels) > 4:
        c5 = Table.read(closed_file5)

    scale=12

    col = 'emp_fwhm'

    plt.figure(1, figsize=(8, 8))
    plt.plot(o1['DIMM'], o1[col]*((500/labels[0][2])**(1/5))/scale, 'bo', label='Open vs DIMM')
    plt.plot(o2['DIMM'], o2[col]*((500/labels[1][2])**(1/5))/scale, 'bo', label='_nolegend_')
    plt.plot(o3['DIMM'], o3[col]*((500/labels[2][2])**(1/5))/scale, 'bo', label='_nolegend_')
    if len(labels) > 3:
        plt.plot(o4['DIMM'], o4[col]*((500/labels[3][2])**(1/5))/scale, 'bo', label='_nolegend_')
    if len(labels) > 4:
        plt.plot(o5['DIMM'], o5[col]*((500/labels[4][2])**(1/5))/scale, 'bo', label='_nolegend_')

    plt.plot(c1['MASS'], c1[col]*((500/labels[0][2])**(1/5))/scale, 'ro', label='Closed vs MASS')
    plt.plot(c2['MASS'], c2[col]*((500/labels[1][2])**(1/5))/scale, 'ro', label='_nolegend_')
    plt.plot(c3['MASS'], c3[col]*((500/labels[2][2])**(1/5))/scale, 'ro', label='_nolegend_')
    if len(labels) > 3:
        plt.plot(c4['MASS'], c4[col]*((500/labels[3][2])**(1/5))/scale, 'ro', label='_nolegend_')
    if len(labels) > 4:
        plt.plot(c5['MASS'], c5[col]*((500/labels[4][2])**(1/5))/scale, 'ro', label='_nolegend_')

    plt.plot([0, 2], [0, 2], 'k-')
    plt.xlabel('Seeing (as)')
    plt.ylabel('Empirical FWHM (as)')
    plt.title(title)
    plt.legend(loc=2)

    return


def plot_hist(labels, data_dir_root, stats_dir_end, title):

    scale=12

    open_file1 = data_dir_root + labels[0][0] + stats_dir_end + "stats_open_mdp.fits"
    open_file2 = data_dir_root + labels[1][0] + stats_dir_end + "stats_open_mdp.fits"
    open_file3 = data_dir_root + labels[2][0] + stats_dir_end + "stats_open_mdp.fits"

    o1 = np.array(Table.read(open_file1)['emp_fwhm']*((500/labels[0][2])**(1/5))/scale)
    o2 = np.array(Table.read(open_file2)['emp_fwhm']*((500/labels[1][2])**(1/5))/scale)
    o3 = np.array(Table.read(open_file3)['emp_fwhm']*((500/labels[2][2])**(1/5))/scale)

    closed_file1 = data_dir_root +labels[0][0]+stats_dir_end + "stats_"+labels[0][1]+"_mdp.fits"
    closed_file2 = data_dir_root +labels[1][0]+stats_dir_end + "stats_"+labels[1][1]+"_mdp.fits"
    closed_file3 = data_dir_root +labels[2][0]+stats_dir_end + "stats_"+labels[2][1]+"_mdp.fits"

    c1 = np.array(Table.read(closed_file1)['emp_fwhm']*((500/labels[0][2])**(1/5))/scale)
    c2 = np.array(Table.read(closed_file2)['emp_fwhm']*((500/labels[1][2])**(1/5))/scale)
    c3 = np.array(Table.read(closed_file3)['emp_fwhm']*((500/labels[2][2])**(1/5))/scale)

    all_fwhm_open = np.concatenate((o1, o2, o3), axis=0)
    all_fwhm_closed = np.concatenate((c1, c2, c3), axis=0)
        
    if len(labels) > 3:
        open_file4 = data_dir_root + labels[3][0] + stats_dir_end + "stats_open_mdp.fits"
        o4 = np.array(Table.read(open_file4)['emp_fwhm']*((500/labels[3][2])**(1/5))/scale)
        closed_file4 = data_dir_root +labels[3][0]+stats_dir_end + "stats_"+labels[3][1]+"_mdp.fits"
        c4 = np.array(Table.read(closed_file4)['emp_fwhm']*((500/labels[3][2])**(1/5))/scale)
        all_fwhm_open = np.concatenate((o1, o2, o3, o4), axis=0)
        all_fwhm_closed = np.concatenate((c1, c2, c3, c4), axis=0)
    
    if len(labels) > 4:
        open_file5 = data_dir_root + labels[4][0] + stats_dir_end + "stats_open_mdp.fits"
        o5 = np.array(Table.read(open_file5)['emp_fwhm']*((500/labels[4][2])**(1/5))/scale)
        closed_file5 = data_dir_root +labels[4][0]+stats_dir_end + "stats_"+labels[4][1]+"_mdp.fits"
        c5 = np.array(Table.read(closed_file5)['emp_fwhm']*((500/labels[4][2])**(1/5))/scale)
        all_fwhm_open = np.concatenate((o1, o2, o3, o4, o5), axis=0)
        all_fwhm_closed = np.concatenate((c1, c2, c3, c4, c5), axis=0) 
        
    plt.figure(2, figsize =(10, 6))
    plt.hist(all_fwhm_open, bins=15, alpha=0.5, color='blue', label='Open Loop')
    plt.hist(all_fwhm_closed, bins=15, alpha=0.5, color='red', label='Closed Loop')
    plt.legend(loc=1)
    plt.xlabel('Empirical FWHM')
    plt.title(title)
    
    return


def filter2wv(filter_label):

    #converts filter label from fits header into wavelength in nanometers 
    
    if filter_label == "I":
        return 806
    elif filter_label == "R":
        return 658
    elif filter_label == "1_micronlp":
        return 1000
    else: 
        print("Filter not found: defaulting to 500 nm")
        return 500


def plot_field_var(starlist):

    mode=starlist.split("_")[-2]
    x_cents, y_cents, fwhm, x_fwhm, y_fwhm, roundness, dist = add_data.read_starlist(starlist)
    plt.figure(1, figsize=(12, 8))
    plt.suptitle('Field Variability', fontsize=18)

    plt.subplot(221)
    plt.plot(dist, y_fwhm,'o', alpha=0.5)
    plt.xlabel('Distance from center of image (pixels)')
    plt.ylabel('Gaussian fit FWHM (not scaled)')
    plt.title('PSF Size')

    plt.subplot(222)
    plt.plot(dist, x_fwhm/y_fwhm,'o', alpha=0.5); 
    plt.xlabel('Distance from center of image (pixels)')
    plt.ylabel('PSF Elongation (x_fwhm/y_fwhm)')
    plt.title('PSF Elongation')

    plt.subplot(223)
    plt.scatter(x_cents, y_cents, c=fwhm)
    plt.colorbar()
    plt.xlabel('x coordinate')
    plt.ylabel('y coordinate')

    plt.subplot(224)
    plt.scatter(x_cents, y_cents, c=x_fwhm/y_fwhm)
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

                
    open_err = stats1['emp_fwhm_std'] * calib1
    closed_err = stats2['emp_fwhm_std'] * calib2

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

    plt.figure(1, figsize=(12, 6))
    plt.errorbar(times1, stats1['emp_fwhm']*calib1, yerr=open_err, fmt='o', label="Open")
    plt.errorbar(times2, stats2['emp_fwhm']*calib2, yerr=closed_err, fmt='ro', label="Closed")
    plt.plot(times1, dimm, 'b-')
    plt.plot(times2, mass, 'r-')
    plt.ylabel(comp_col)
    plt.xlabel("UTC Time")
    plt.xticks(rotation=35)
    plt.ylim(0, 2.5)
    plt.title(title)
    plt.gca().xaxis.set_major_formatter(mp_dates.DateFormatter('%H:%M'))
    plt.legend()

    plt.savefig(plots_dir + 'fwhm_v_time' + '.png')
    return
