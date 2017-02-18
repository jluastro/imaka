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

    scale = 0.12
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
        plt.plot(stats[ss]['Index'], stats[ss]['FWHM']*scale, marker='o', linestyle='none', label=suffixes[ss])
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
        plt.plot(stats[ss]['Index'], stats[ss]['emp_fwhm']*scale, marker='o', linestyle='none', label=suffixes[ss])
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
        plt.plot(stats[ss]['Index'], stats[ss]['xFWHM']*scale, marker='o', color=c, linestyle='none', label='X ' + suffixes[ss])
        plt.plot(stats[ss]['Index'], stats[ss]['yFWHM']*scale, marker='^', color=c, linestyle='none', label='Y ' + suffixes[ss])
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
        plt.plot(utcs[ss], stats[ss]['FWHM']*scale, marker='o', linestyle='none', label=suffixes[ss])
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
        plt.plot(utcs[ss], stats[ss]['emp_fwhm']*scale, marker='o', linestyle='none', label=suffixes[ss])
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
        plt.plot(utcs[ss], stats[ss]['xFWHM']*scale, marker='o', color=c, linestyle='none', label='X ' + suffixes[ss])
        plt.plot(utcs[ss], stats[ss]['yFWHM']*scale, marker='^', color=c, linestyle='none', label='Y ' + suffixes[ss])
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

    scale = 0.12
    
    time_fmt = mp_dates.DateFormatter('%H:%M')    
    time_loc = ticker.MaxNLocator(nbins=6)

    plt.figure(1, figsize=(6, 6))
    plt.clf()

    for ii in range(len(suffixes)):
        c = np.take(colors, ii, mode='wrap')
        plt.plot(utc[ii], stats[ii]['emp_fwhm']*scale, marker='o', color=c, linestyle='none', label=suffixes[ii])
        
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
        plt.plot(stats[ii]['MASS'], stats[ii]['emp_fwhm']*scale, marker='o', linestyle='none', label=suffixes[ii])
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
        plt.plot(stats[ii]['DIMM'], stats[ii]['emp_fwhm']*scale, marker='o', linestyle='none', label=suffixes[ii])
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

    stats_dir = root_dir + date + '/fli/reduce/stats/'
    plots_dir = root_dir + date + '/fli/reduce/plots/'

    scale = 0.12
    
    stats = []
    utcs = []
    max_fwhm = 0.0
    for suffix in suffixes:
        ss = Table.read(stats_dir + 'stats_' + suffix + '.fits')
        stats.append( ss )
        utc_dt = [datetime.strptime(ss['TIME_UTC'][ii], '%I:%M:%S') for ii in range(len(ss))]
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
    plt.ylim(0, 2)
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
    plt.ylim(0, 2)
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
    plt.ylim(0, 2)
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
