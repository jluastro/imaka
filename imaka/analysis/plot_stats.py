from astropy.table import Table
import subprocess
import os
from imaka.reduce import util


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
    

def plot_stack_stats():
    
    return

import csv
import numpy as np
import matplotlib.pyplot as plt

def make_stats_plot(stats_csv_file, plot_title):
    
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
