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
    all_dates = ['20170110', '20170111', '20170112', '20170113']

    if dates == 'all' or dates == None:
        dates = all_dates

    for date in dates:
        remote_file = 'imaka@onaga.ifa.hawaii.edu:/Volumes/DATA/imaka/{0:s}/fli/reduce/stats/stats*.fits'.format(date)
        destination = '{0:s}/{1:s}/fli/reduce/stats/'.format(output_root, date)

        util.mkdir(destination)

        print(remote_file)
        print(destination)
        p = subprocess.Popen(["scp", remote_file, destination])
        sts = os.waitpid(p.pid, 0)

    return
    

def plot_stack_stats(date, suffix='', root_dir='/Users/jlu/work/imaka/pleiades/'):
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
    
    so = Table.read(stats_dir + 'stats_open' + suffix + '.fits')   # Open Loop
    st = Table.read(stats_dir + 'stats_ttf' + suffix + '.fits')    # TTF-only loop
    sc = Table.read(stats_dir + 'stats_closed' + suffix + '.fits') # Closed loop

    scale = 0.12

    # 
    # Plots for ratio of improvements. First we need to find a matching
    # closed loop image to go with each open (and TTF) image.
    #
    tree_sc = scipy.spatial.KDTree(np.array([sc['Index']]).T)
    dist_oc, idx_in_c_of_o = tree_sc.query(np.array([so['Index']]).T, 1)
    dist_tc, idx_in_c_of_t = tree_sc.query(np.array([st['Index']]).T, 1)

    #####
    # FWHM vs. Frame
    #####
    plt.figure(1, figsize=(12, 6))
    plt.subplots_adjust(left=0.1, bottom=0.15)
    plt.clf()
    plt.subplot(121)
    plt.plot(so['Index'], so['FWHM']*scale, 'bo', label='Open')
    plt.plot(st['Index'], st['FWHM']*scale, 'go', label='TTF')
    plt.plot(sc['Index'], sc['FWHM']*scale, 'ro', label='Closed')
    plt.xlabel('Frame Number')
    plt.ylabel('Gaussian-Fit FWHM (")')
    plt.legend(numpoints=1)
    plt.ylim(0, 1.5)
    plt.title(date)

    plt.subplot(122)
    plt.plot(so['Index'], sc['FWHM'][idx_in_c_of_o] / so['FWHM'], 'bo', label='Closed / Open')
    plt.plot(st['Index'], sc['FWHM'][idx_in_c_of_t] / st['FWHM'], 'ro', label='Closed / TTF')
    plt.xlabel('Frame Number')
    plt.ylabel('Ratio of Gaussian-Fit FWHM')
    plt.legend(numpoints=1)
    plt.axhline(1.0, color='k', linestyle='--', linewidth=2)
    plt.ylim(0, 1.3)
    plt.savefig(plots_dir + 'fwhm_vs_frame.png')

    #####
    # Empirical FWHM
    #####
    plt.figure(2, figsize=(12, 6))
    plt.subplots_adjust(left=0.1, bottom=0.15)
    plt.clf()
    plt.subplot(121)
    plt.plot(so['Index'], so['emp_fwhm']*scale, 'bo', label='Open')
    plt.plot(st['Index'], st['emp_fwhm']*scale, 'go', label='TTF')
    plt.plot(sc['Index'], sc['emp_fwhm']*scale, 'ro', label='Closed')
    # plt.errorbar(so['Index'], so['emp_fwhm']*scale, fmt='bo', yerr=so['emp_fwhm_std'], label='Open')
    # plt.errorbar(st['Index'], st['emp_fwhm']*scale, fmt='go', yerr=st['emp_fwhm_std'], label='TTF')
    # plt.errorbar(sc['Index'], sc['emp_fwhm']*scale, fmt='ro', yerr=sc['emp_fwhm_std'], label='Closed')
    plt.xlabel('Frame Number')
    plt.ylabel('Empirical FWHM (")')
    plt.legend(numpoints=1)
    plt.ylim(0, 1.5)

    plt.subplot(122)
    plt.plot(so['Index'], sc['emp_fwhm'][idx_in_c_of_o] / so['emp_fwhm'], 'bo', label='Closed / Open')
    plt.plot(st['Index'], sc['emp_fwhm'][idx_in_c_of_t] / st['emp_fwhm'], 'ro', label='Closed / TTF')
    plt.xlabel('Frame Number')
    plt.ylabel('Ratio of Empircal FWHM')
    plt.legend(numpoints=1)
    plt.axhline(1.0, color='k', linestyle='--', linewidth=2)
    plt.ylim(0, 1.3)
    plt.savefig(plots_dir + 'efwhm_vs_frame.png')
    
    
    #####
    # EE 50
    #####
    plt.figure(3, figsize=(12, 6))
    plt.subplots_adjust(left=0.1, bottom=0.15)
    plt.clf()
    plt.subplot(121)
    plt.plot(so['Index'], so['EE50'], 'bo', label='Open')
    plt.plot(st['Index'], st['EE50'], 'go', label='TTF')
    plt.plot(sc['Index'], sc['EE50'], 'ro', label='Closed')
    plt.xlabel('Frame Number')
    plt.ylabel('Radius of 50% Encircled Energy (")')
    plt.legend(numpoints=1)
    plt.ylim(0, 1.5)
    plt.title(date)

    plt.subplot(122)
    plt.plot(so['Index'], sc['EE50'][idx_in_c_of_o] / so['EE50'], 'bo', label='Closed / Open')
    plt.plot(st['Index'], sc['EE50'][idx_in_c_of_t] / st['EE50'], 'ro', label='Closed / TTF')
    plt.xlabel('Frame Number')
    plt.ylabel('Ratio of 50% EE Radius')
    plt.legend(numpoints=1)
    plt.axhline(1.0, color='k', linestyle='--', linewidth=2)
    plt.ylim(0, 1.3)
    plt.savefig(plots_dir + 'ee50_vs_frame.png')

    #####
    # EE 80
    #####
    plt.figure(4, figsize=(12, 6))
    plt.subplots_adjust(left=0.1, bottom=0.15)
    plt.clf()
    plt.subplot(121)
    plt.plot(so['Index'], so['EE80'], 'bo', label='Open')
    plt.plot(st['Index'], st['EE80'], 'go', label='TTF')
    plt.plot(sc['Index'], sc['EE80'], 'ro', label='Closed')
    plt.xlabel('Frame Number')
    plt.ylabel('Radius of 80% Encircled Energy (")')
    plt.legend(numpoints=1)
    plt.ylim(0, 1.5)
    plt.title(date)
    
    plt.subplot(122)
    plt.plot(so['Index'], sc['EE80'][idx_in_c_of_o] / so['EE80'], 'bo', label='Closed / Open')
    plt.plot(st['Index'], sc['EE80'][idx_in_c_of_t] / st['EE80'], 'ro', label='Closed / TTF')
    plt.xlabel('Frame Number')
    plt.ylabel('Ratio of 80% EE Radius')
    plt.legend(numpoints=1)
    plt.axhline(1.0, color='k', linestyle='--', linewidth=2)
    plt.ylim(0, 1.3)
    plt.savefig(plots_dir + 'ee80_vs_frame.png')

    #####
    # NEA
    #####
    plt.figure(5, figsize=(12, 6))
    plt.subplots_adjust(left=0.1, bottom=0.15)
    plt.clf()
    plt.subplot(121)
    plt.plot(so['Index'], so['NEA'], 'bo', label='Open')
    plt.plot(st['Index'], st['NEA'], 'go', label='TTF')
    plt.plot(sc['Index'], sc['NEA'], 'ro', label='Closed')
    plt.xlabel('Frame Number')
    plt.ylabel('Noise Equivalent Area (Sq. Arcsec)')
    plt.legend(numpoints=1)
    plt.ylim(0, 3)
    plt.title(date)
    
    plt.subplot(122)
    plt.plot(so['Index'], sc['NEA'][idx_in_c_of_o] / so['NEA'], 'bo', label='Closed / Open')
    plt.plot(st['Index'], sc['NEA'][idx_in_c_of_t] / st['NEA'], 'ro', label='Closed / TTF')
    plt.xlabel('Frame Number')
    plt.ylabel('Ratio of NEA')
    plt.legend(numpoints=1)
    plt.axhline(1.0, color='k', linestyle='--', linewidth=2)
    plt.ylim(0, 1.3)
    plt.savefig(plots_dir + 'nea_vs_frame.png')
    
    #####
    # NEA2
    #####
    plt.figure(6, figsize=(12, 6))
    plt.subplots_adjust(left=0.1, bottom=0.15)
    plt.clf()
    plt.subplot(121)
    plt.plot(so['Index'], so['NEA2'], 'bo', label='Open')
    plt.plot(st['Index'], st['NEA2'], 'go', label='TTF')
    plt.plot(sc['Index'], sc['NEA2'], 'ro', label='Closed')
    plt.xlabel('Frame Number')
    plt.ylabel('Noise Equivalent Area (Sq. Arcsec)')
    plt.legend(numpoints=1)
    plt.ylim(0, 3)
    plt.title(date)

    plt.subplot(122)
    plt.plot(so['Index'], sc['NEA2'][idx_in_c_of_o] / so['NEA2'], 'bo', label='Closed / Open')
    plt.plot(st['Index'], sc['NEA2'][idx_in_c_of_t] / st['NEA2'], 'ro', label='Closed / TTF')
    plt.xlabel('Frame Number')
    plt.ylabel('Ratio of NEA2')
    plt.legend(numpoints=1)
    plt.axhline(1.0, color='k', linestyle='--', linewidth=2)
    plt.ylim(0, 1.3)
    plt.savefig(plots_dir + 'nea2_vs_frame.png')

    #####
    # FWHM for each direction
    #####
    plt.figure(7, figsize=(12, 6))
    plt.subplots_adjust(left=0.1, bottom=0.15)
    plt.clf()
    plt.subplot(121)
    plt.plot(so['Index'], so['xFWHM']*scale, 'bo', label='Open X')
    plt.plot(st['Index'], st['xFWHM']*scale, 'go', label='TTF X')
    plt.plot(sc['Index'], sc['xFWHM']*scale, 'ro', label='Closed X')
    plt.plot(so['Index'], so['yFWHM']*scale, 'b^', label='Open Y')
    plt.plot(st['Index'], st['yFWHM']*scale, 'g^', label='TTF Y')
    plt.plot(sc['Index'], sc['yFWHM']*scale, 'r^', label='Closed Y')
    plt.xlabel('Frame Number')
    plt.ylabel('Gaussian-Fit FWHM (")')
    plt.legend(numpoints=1, fontsize=10)
    plt.ylim(0, 1.5)
    plt.title(date)

    plt.subplot(122)
    plt.plot(so['Index'], sc['xFWHM'][idx_in_c_of_o] / so['xFWHM'], 'bo', label='X Closed / Open')
    plt.plot(st['Index'], sc['xFWHM'][idx_in_c_of_t] / st['xFWHM'], 'ro', label='X Closed / TTF')
    plt.plot(so['Index'], sc['yFWHM'][idx_in_c_of_o] / so['yFWHM'], 'b^', label='Y Closed / Open')
    plt.plot(st['Index'], sc['yFWHM'][idx_in_c_of_t] / st['yFWHM'], 'r^', label='Y Closed / TTF')
    plt.xlabel('Frame Number')
    plt.ylabel('Ratio of Gaussian-Fit FWHM')
    plt.legend(numpoints=1, fontsize=10)
    plt.axhline(1.0, color='k', linestyle='--', linewidth=2)
    plt.ylim(0, 1.3)
    plt.savefig(plots_dir + 'xyfwhm_vs_frame.png')

    
    ##########
    # All plots vs. time.
    ##########

    utc_o_dt = [datetime.strptime(so['TIME_UTC'][ii], '%I:%M:%S') for ii in range(len(so))]
    utc_t_dt = [datetime.strptime(st['TIME_UTC'][ii], '%I:%M:%S') for ii in range(len(st))]
    utc_c_dt = [datetime.strptime(sc['TIME_UTC'][ii], '%I:%M:%S') for ii in range(len(sc))]

    time_fmt = mp_dates.DateFormatter('%H:%M')    
    
    
    #####
    # FWHM
    #####
    plt.figure(8, figsize=(12, 6))
    plt.subplots_adjust(left=0.1, bottom=0.15)
    plt.clf()
    plt.subplot(121)
    plt.plot(utc_o_dt, so['FWHM']*scale, 'bo', label='Open')
    plt.plot(utc_t_dt, st['FWHM']*scale, 'go', label='TTF')
    plt.plot(utc_c_dt, sc['FWHM']*scale, 'ro', label='Closed')
    plt.gca().xaxis.set_major_formatter(time_fmt)
    plt.xticks(rotation=35)
    plt.xlabel('UTC Time (hr)')
    plt.ylabel('Gaussian-Fit FWHM (")')
    plt.legend(numpoints=1)
    plt.ylim(0, 1.5)
    plt.title(date)

    plt.subplot(122)
    plt.plot(utc_o_dt, sc['FWHM'][idx_in_c_of_o] / so['FWHM'], 'bo', label='Closed / Open')
    plt.plot(utc_t_dt, sc['FWHM'][idx_in_c_of_t] / st['FWHM'], 'ro', label='Closed / TTF')
    plt.gca().xaxis.set_major_formatter(time_fmt)
    plt.xticks(rotation=35)
    plt.xlabel('UTC Time (hr)')
    plt.ylabel('Ratio of Gaussian-Fit FWHM')
    plt.legend(numpoints=1)
    plt.axhline(1.0, color='k', linestyle='--', linewidth=2)
    plt.ylim(0, 1.3)
    plt.savefig(plots_dir + 'fwhm_vs_time.png')

    #####
    # Empirical FWHM
    #####
    plt.figure(9, figsize=(12, 6))
    plt.subplots_adjust(left=0.1, bottom=0.15)
    plt.clf()
    plt.subplot(121)
    plt.plot(utc_o_dt, so['emp_fwhm']*scale, 'bo', label='Open')
    plt.plot(utc_t_dt, st['emp_fwhm']*scale, 'go', label='TTF')
    plt.plot(utc_c_dt, sc['emp_fwhm']*scale, 'ro', label='Closed')
    plt.gca().xaxis.set_major_formatter(time_fmt)
    plt.xticks(rotation=35)
    plt.xlabel('UTC Time (hr)')
    plt.ylabel('Empirical FWHM (")')
    plt.legend(numpoints=1)
    plt.ylim(0, 1.5)
    plt.title(date)

    plt.subplot(122)
    plt.plot(utc_o_dt, sc['emp_fwhm'][idx_in_c_of_o] / so['emp_fwhm'], 'bo', label='Closed / Open')
    plt.plot(utc_t_dt, sc['emp_fwhm'][idx_in_c_of_t] / st['emp_fwhm'], 'ro', label='Closed / TTF')
    plt.gca().xaxis.set_major_formatter(time_fmt)
    plt.xticks(rotation=35)
    plt.xlabel('UTC Time (hr)')
    plt.ylabel('Ratio of Empircal FWHM')
    plt.legend(numpoints=1)
    plt.axhline(1.0, color='k', linestyle='--', linewidth=2)
    plt.ylim(0, 1.3)
    plt.savefig(plots_dir + 'efwhm_vs_time.png')
    
    #####
    # EE 50
    #####
    plt.figure(10, figsize=(12, 6))
    plt.subplots_adjust(left=0.1, bottom=0.15)
    plt.clf()
    plt.subplot(121)
    plt.plot(utc_o_dt, so['EE50'], 'bo', label='Open')
    plt.plot(utc_t_dt, st['EE50'], 'go', label='TTF')
    plt.plot(utc_c_dt, sc['EE50'], 'ro', label='Closed')
    plt.gca().xaxis.set_major_formatter(time_fmt)
    plt.xticks(rotation=35)
    plt.xlabel('UTC Time (hr)')
    plt.ylabel('Radius of 50% Encircled Energy (")')
    plt.legend(numpoints=1)
    plt.ylim(0, 1.5)
    plt.title(date)

    plt.subplot(122)
    plt.plot(utc_o_dt, sc['EE50'][idx_in_c_of_o] / so['EE50'], 'bo', label='Closed / Open')
    plt.plot(utc_t_dt, sc['EE50'][idx_in_c_of_t] / st['EE50'], 'ro', label='Closed / TTF')
    plt.gca().xaxis.set_major_formatter(time_fmt)
    plt.xticks(rotation=35)
    plt.xlabel('UTC Time (hr)')
    plt.ylabel('Ratio of 50% EE Radius')
    plt.legend(numpoints=1)
    plt.axhline(1.0, color='k', linestyle='--', linewidth=2)
    plt.ylim(0, 1.3)
    plt.savefig(plots_dir + 'ee50_vs_time.png')

    #####
    # EE 80
    #####
    plt.figure(11, figsize=(12, 6))
    plt.subplots_adjust(left=0.1, bottom=0.15)
    plt.clf()
    plt.subplot(121)
    plt.plot(utc_o_dt, so['EE80'], 'bo', label='Open')
    plt.plot(utc_t_dt, st['EE80'], 'go', label='TTF')
    plt.plot(utc_c_dt, sc['EE80'], 'ro', label='Closed')
    plt.gca().xaxis.set_major_formatter(time_fmt)
    plt.xticks(rotation=35)
    plt.xlabel('UTC Time (hr)')
    plt.ylabel('Radius of 80% Encircled Energy (")')
    plt.legend(numpoints=1)
    plt.ylim(0, 1.5)
    plt.title(date)

    plt.subplot(122)
    plt.plot(utc_o_dt, sc['EE80'][idx_in_c_of_o] / so['EE80'], 'bo', label='Closed / Open')
    plt.plot(utc_t_dt, sc['EE80'][idx_in_c_of_t] / st['EE80'], 'ro', label='Closed / TTF')
    plt.gca().xaxis.set_major_formatter(time_fmt)
    plt.xticks(rotation=35)
    plt.xlabel('UTC Time (hr)')
    plt.ylabel('Ratio of 80% EE Radius')
    plt.legend(numpoints=1)
    plt.axhline(1.0, color='k', linestyle='--', linewidth=2)
    plt.ylim(0, 1.3)
    plt.savefig(plots_dir + 'ee80_vs_time.png')

    #####
    # NEA
    #####
    plt.figure(12, figsize=(12, 6))
    plt.subplots_adjust(left=0.1, bottom=0.15)
    plt.clf()
    plt.subplot(121)
    plt.plot(utc_o_dt, so['NEA'], 'bo', label='Open')
    plt.plot(utc_t_dt, st['NEA'], 'go', label='TTF')
    plt.plot(utc_c_dt, sc['NEA'], 'ro', label='Closed')
    plt.gca().xaxis.set_major_formatter(time_fmt)
    plt.xticks(rotation=35)
    plt.xlabel('UTC Time (hr)')
    plt.ylabel('Noise Equivalent Area (Sq. Arcsec)')
    plt.legend(numpoints=1)
    plt.ylim(0, 3)
    plt.title(date)

    plt.subplot(122)
    plt.plot(utc_o_dt, sc['NEA'][idx_in_c_of_o] / so['NEA'], 'bo', label='Closed / Open')
    plt.plot(utc_t_dt, sc['NEA'][idx_in_c_of_t] / st['NEA'], 'ro', label='Closed / TTF')
    plt.gca().xaxis.set_major_formatter(time_fmt)
    plt.xticks(rotation=35)
    plt.xlabel('UTC Time (hr)')
    plt.ylabel('Ratio of NEA')
    plt.legend(numpoints=1)
    plt.axhline(1.0, color='k', linestyle='--', linewidth=2)
    plt.ylim(0, 1.3)
    plt.savefig(plots_dir + 'nea_vs_time.png')

    #####
    # NEA2
    #####
    plt.figure(13, figsize=(12, 6))
    plt.subplots_adjust(left=0.1, bottom=0.15)
    plt.clf()
    plt.subplot(121)
    plt.plot(utc_o_dt, so['NEA2'], 'bo', label='Open')
    plt.plot(utc_t_dt, st['NEA2'], 'go', label='TTF')
    plt.plot(utc_c_dt, sc['NEA2'], 'ro', label='Closed')
    plt.gca().xaxis.set_major_formatter(time_fmt)
    plt.xticks(rotation=35)
    plt.xlabel('UTC Time (hr)')
    plt.ylabel('Noise Equivalent Area (Sq. Arcsec)')
    plt.legend(numpoints=1)
    plt.ylim(0, 3)
    plt.title(date)

    plt.subplot(122)
    plt.plot(utc_o_dt, sc['NEA2'][idx_in_c_of_o] / so['NEA2'], 'bo', label='Closed / Open')
    plt.plot(utc_t_dt, sc['NEA2'][idx_in_c_of_t] / st['NEA2'], 'ro', label='Closed / TTF')
    plt.gca().xaxis.set_major_formatter(time_fmt)
    plt.xticks(rotation=35)
    plt.xlabel('UTC Time (hr)')
    plt.ylabel('Ratio of NEA2')
    plt.legend(numpoints=1)
    plt.axhline(1.0, color='k', linestyle='--', linewidth=2)
    plt.ylim(0, 1.3)
    plt.savefig(plots_dir + 'nea2_vs_time.png')
    
    #####
    # FWHM for each direction
    #####
    plt.figure(14, figsize=(12, 6))
    plt.subplots_adjust(left=0.1, bottom=0.15)
    plt.subplot(121)
    plt.plot(utc_o_dt, so['xFWHM']*scale, 'bo', label='Open X')
    plt.plot(utc_t_dt, st['xFWHM']*scale, 'go', label='TTF X')
    plt.plot(utc_c_dt, sc['xFWHM']*scale, 'ro', label='Closed X')
    plt.plot(utc_o_dt, so['yFWHM']*scale, 'b^', label='Open Y')
    plt.plot(utc_t_dt, st['yFWHM']*scale, 'g^', label='TTF Y')
    plt.plot(utc_c_dt, sc['yFWHM']*scale, 'r^', label='Closed Y')
    plt.gca().xaxis.set_major_formatter(time_fmt)
    plt.xticks(rotation=35)
    plt.xlabel('UTC Time (hr)')
    plt.ylabel('Gaussian-Fit FWHM (")')
    plt.legend(numpoints=1, fontsize=10)
    plt.ylim(0, 1.5)
    plt.title(date)
    
    plt.subplot(122)
    plt.plot(utc_o_dt, sc['xFWHM'][idx_in_c_of_o] / so['xFWHM'], 'bo', label='X Closed / Open')
    plt.plot(utc_t_dt, sc['xFWHM'][idx_in_c_of_t] / st['xFWHM'], 'ro', label='X Closed / TTF')
    plt.plot(utc_o_dt, sc['yFWHM'][idx_in_c_of_o] / so['yFWHM'], 'b^', label='Y Closed / Open')
    plt.plot(utc_t_dt, sc['yFWHM'][idx_in_c_of_t] / st['yFWHM'], 'r^', label='Y Closed / TTF')
    plt.gca().xaxis.set_major_formatter(time_fmt)
    plt.xticks(rotation=35)
    plt.xlabel('UTC Time (hr)')
    plt.ylabel('Ratio of Gaussian-Fit FWHM')
    plt.legend(numpoints=1, fontsize=10)
    plt.axhline(1.0, color='k', linestyle='--', linewidth=2)
    plt.ylim(0, 1.3)
    plt.savefig(plots_dir + 'xyfwhm_vs_time.png')
    
    return

def compare_fwhm(stats_csv_file, plot_title):
    
    #Makes a figure of several plots comparing different measures of FWHM
    #(empirical, gaussian, NEA)
    
    FWHM = []
    xFWHM = []
    yFWHM = []
    NEA_FWHM = []
    emp_FWHM = []
    stat_list = stats_csv_file
    with open(stat_list, "r") as f:
        reader = csv.reader(f)
        for row in reader:
            if row[0] == 'Image':
                labels = [row[6], row[10], row[12], row[13], row[15]]
                full_labels = row
            else:
                FWHM.append(float(row[6]))
                NEA_FWHM.append((2*(float(row[10])/np.pi)**0.5))
                xFWHM.append(float(row[12]))
                yFWHM.append(float(row[13]))
                emp_FWHM.append(float(row[15]))

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

    plt.savefig(stats_csv_file.replace('.csv', '_plot.png'))
    
    return
