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
import glob
import matplotlib.pyplot as plt
import matplotlib
import datetime

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

        print(local_file)
        print(destination)
        p = subprocess.Popen("scp " + local_file + " " + destination, shell=True)
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
    
    for file in stats_files:
        
        #extract/calculate FWHMs
    
        FWHM = []
        xFWHM = []
        yFWHM = []
        NEA_FWHM = []
        emp_FWHM = []
        with open(file, "r") as f:
            reader = csv.reader(f)
            for row in reader:
                if row[0] == 'Image':
                    labels = [row[6], row[10], row[12], row[13], row[15]]
                    full_labels = row
                else:
                    FWHM.append(float(row[6]))
                    NEA_FWHM.append((((float(row[10])/0.12)/np.pi)**0.5)*2)
                    xFWHM.append(float(row[12]))
                    yFWHM.append(float(row[13]))
                    emp_FWHM.append(float(row[15]))

        #Make plot

        #Get plot title from file name
        title_seg = file.split("_")[-1]
        if title_seg == 'closed.csv':
            loop_type = "Closed Loop"
        elif title_seg == 'open.csv':
            loop_type = "Open Loop"
        elif title_seg == "tt.csv":
            loop_type = "TT"
        elif title_seg == "ttf.csv":
            loop_type = "TTF"

        names = file.split("/")
        for directory in names:
            if "201" in directory:
                date = (directory[4:6]+"/"+directory[6:]+"/"+directory[0:4])

        plot_title = date + " FWHM Comparison: "+loop_type

        #Make figure
        
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
        plt.tick_params(axis='both', which='both', bottom='off', top='off', left='off', right='off', labelbottom='off', labelleft='off') # labels along the bottom edge are off

        plt.subplots_adjust(hspace=0.25, wspace=0.25)

        plt.savefig(out_dir + file.split("/")[-1].replace('.csv', '_plot.png'))

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

    stats = []
    utc = []

    for suffix in suffixes:
        st = Table.read(stats_dir + 'stats_' + suffix + '.fits')
        stats.append(st)

        utc_dt = [datetime.strptime(st['TIME_UTC'][ii], '%I:%M:%S') for ii in range(len(st))]
        utc.append(utc_dt)

    scale = 0.12
    
    time_fmt = mp_dates.DateFormatter('%H:%M')    

    plt.figure(1, figsize=(6, 6))
    plt.clf()
    for ii in range(len(suffixes)):
        plt.plot(utc[ii], stats[ii]['emp_fwhm']*scale, marker='o', linestyle='none', label=suffixes[ii])
    plt.gca().xaxis.set_major_formatter(time_fmt)
    plt.xticks(rotation=35)
    plt.xlabel('UTC Time (hr)')
    plt.ylabel('Empirical FWHM (")')
    plt.legend(numpoints=1)
    plt.ylim(0, 1.5)
    plt.title(date)
    plt.savefig(plots_dir + 'mdp_efwhm_vs_time' + out_suffix + '.png')

    plt.figure(2, figsize=(6, 6))
    plt.clf()
    for ii in range(len(suffixes)):
        plt.plot(utc[ii], stats[ii]['NEA'], marker='o', linestyle='none', label=suffixes[ii])
    plt.gca().xaxis.set_major_formatter(time_fmt)
    plt.xticks(rotation=35)
    plt.xlabel('UTC Time (hr)')
    plt.ylabel('NEA (Sq. Arcsec)')
    plt.legend(numpoints=1)
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
    plt.ylim(0, 1.5)
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
    utc = []
    
    for suffix in suffixes:
        st = Table.read(stats_dir + 'stats_' + suffix + '_mdp.fits')
        stats.append(st)

        utc_dt = [datetime.strptime(st['TIME_UTC'][ii], '%I:%M:%S') for ii in range(len(st))]
        utc.append(utc_dt)

    scale = 0.12
    
    time_fmt = mp_dates.DateFormatter('%H:%M')    

    hlis = [0, 0.5, 1, 2, 4, 8, 16]

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
    
    cngl = 0.98*0.00000055*206265/gl**(-5/3) / (16.7*(0.00000055)**(-2))
    profave = np.insert(profave, obj=0, values=np.mean(cngl))
    plt.figure(1)
    plt.clf()
    plt.plot(profave, hlis)
    plt.xlabel(r'C_n^2 dh')
    plt.ylabel(r'h (km)')
    plt.title(date, fontsize=12)
    plt.savefig(plots_dir + 'mass_profile' + out_suffix + '.png')
    plt.show()

    return

def plot_abc(date, suffixes=['open', 'closed'], out_suffix='', root_dir=''):

    #plots different FWHMs as function of UTC time

    colors = ['b', 'g', 'r', 'c', 'm', 'k', 'yellow', 'purple', 'orange']

    stats_file_root = root_dir + date + '/FLI/reduce/stats/'
    stats_files = []
    for suffix in suffixes:
        name = stats_file_root + "stats_" + suffix + ".csv"
        stats_files.append(name)
        
    ###NEA FWHM PLOT###

    plt.figure()
    for ii in range(len(stats_files)):
        data_label = stats_files[ii].split("/")[-1].split("_")[1].split(".")[0]
        TIME_UTC = []
        NEA_FWHM = []
        with open(stats_files[ii], "r") as f:
            reader = csv.reader(f)
            for row in reader:
                if row[0] == 'Image':
                    labels = [row[5], row[10], row[15]]
                    full_labels = row
                else:
                    TIME_UTC.append(row[5])
                    NEA_FWHM.append((((float(row[10])/0.12)/np.pi)**0.5)*2)

        date_times = []
        for time in TIME_UTC:
            dt_obj = datetime.strptime(time, '%H:%M:%S')
            dt_obj.strftime("%I:%M:%S %p")
            date_times.append(dt_obj)
        dates = matplotlib.dates.date2num(date_times)
        plt.plot_date(dates, NEA_FWHM, colors[ii]+'o', alpha=0.5, label=data_label)

    plt.legend(loc=2, bbox_to_anchor=(1, 1))
    plt.xticks(rotation=45)
    plt.xlabel("UTC Time")
    plt.ylabel("NEA FWHM")
    plt.title("NEA FWHM over Time", fontsize=18)

    ###EMPIRICAL FWHM PLOT###
    
    plt.figure()
    for ii in range(len(stats_files)):
        data_label = stats_files[ii].split("/")[-1].split("_")[1].split(".")[0]
        TIME_UTC = []
        emp_FWHM = []
        with open(stats_files[ii], "r") as f:
            reader = csv.reader(f)
            for row in reader:
                if row[0] == 'Image':
                    labels = [row[5], row[10], row[15]]
                    full_labels = row
                else:
                    TIME_UTC.append(row[5])
                    emp_FWHM.append(float(row[15]))

        date_times = []
        for time in TIME_UTC:
            dt_obj = datetime.strptime(time, '%H:%M:%S')
            dt_obj.strftime("%I:%M:%S %p")
            date_times.append(dt_obj)
        dates = matplotlib.dates.date2num(date_times)
        plt.plot_date(dates, emp_FWHM, colors[ii]+'o', alpha=0.5, label=data_label)

    plt.legend(loc=2, bbox_to_anchor=(1, 1))
    plt.xticks(rotation=45)
    plt.xlabel("UTC Time")
    plt.ylabel("Empirical FWHM")
    plt.title("Empirical FWHM over Time", fontsize=18)    


    ###GAUSSIAN FWHM PLOT###

    plt.figure()
    for ii in range(len(stats_files)):
        data_label = stats_files[ii].split("/")[-1].split("_")[1].split(".")[0]
        TIME_UTC = []
        x_FWHM = []
        y_FWHM = []
        with open(stats_files[ii], "r") as f:
            reader = csv.reader(f)
            for row in reader:
                if row[0] == 'Image':
                    labels = [row[5], row[10], row[15]]
                    full_labels = row
                else:
                    TIME_UTC.append(row[5])
                    x_FWHM.append(float(row[12]))
                    y_FWHM.append(float(row[13]))

        date_times = []
        for time in TIME_UTC:
            dt_obj = datetime.strptime(time, '%H:%M:%S')
            dt_obj.strftime("%I:%M:%S %p")
            date_times.append(dt_obj)
        dates = matplotlib.dates.date2num(date_times)
        plt.plot_date(dates, x_FWHM, colors[ii]+"o", alpha=0.5, label=data_label)
        plt.plot_date(dates, y_FWHM, colors[ii]+"s", alpha=0.5, label=data_label)

    plt.legend(loc=2, bbox_to_anchor=(1, 1))
    plt.xticks(rotation=45)
    plt.xlabel("UTC Time")
    plt.ylabel("Gaussian FWHM")
    plt.title("Gaussian FWHM over Time (Circle = x, Square = y)", fontsize=18)

    return


def get_color_list():
    rcParams = plt.matplotlib.rcParams
    prop_cycler = rcParams['axes.prop_cycle']
    if prop_cycler is None and 'axes.color_cycle' in rcParams:
        clist = rcParams['axes.color_cycle']
        prop_cycler = cycler('color', clist)

    colors = [item['color'] for item in prop_cycler]
    return colors

