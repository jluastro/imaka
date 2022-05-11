## imaka/analysis/plots.py
## Eden McEwen
# A series of plotting codes for analysis of imaka

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
from astropy.io import fits
from imaka.analysis import moffat
from astropy.modeling import fitting
from astropy.stats import sigma_clip
from astropy.table import Column, Table
import scipy.linalg
from matplotlib import cm
from astropy.stats import sigma_clipped_stats
import matplotlib.lines as mlines
import matplotlib.patches as mpatches


# for plotting and checking one files starlist
def plot_starlist_stats(test_img_base, root_dir, fld, night, rebinned=False):

    # creating grid for subplots
    fig, ax = plt.subplots(ncols=2, nrows=1, figsize=(20,10))

    ax1 = plt.subplot2grid((3, 3),(0,0))
    ax2 = plt.subplot2grid((3, 3), (1, 0))
    ax3 = plt.subplot2grid((3, 3), (0, 1), rowspan=2, colspan=1)
    ax4 = plt.subplot2grid((3, 3), (0, 2))
    ax5 = plt.subplot2grid((3, 3), (1, 2))

    # Plotting star FWHM vs. Mag
    b2 = "bin2" if rebinned else ""
    b2_suff = "_bin2" if rebinned else ""
    img, hdr = fits.getdata(root_dir + f'reduce/{fld}/{b2}/{test_img_base}_clean{b2_suff}.fits', header=True)
    scale = util.get_plate_scale(img, hdr)
    del img
    del hdr
    stars = Table.read(root_dir + f'reduce/{fld}/{b2}/{test_img_base}_clean{b2_suff}_stars.txt', format='ascii')
    stars2 = Table.read(root_dir + f'reduce/{fld}/{b2}/{test_img_base}_clean{b2_suff}_stars_stats.fits')
    img, hdr = fits.getdata(root_dir + f'reduce/{fld}/{b2}/{test_img_base}_clean{b2_suff}.fits', header=True)
    print(root_dir + f'reduce/{fld}/{b2}/{test_img_base}_clean{b2_suff}.fits')

    ax1.plot(stars['mag'], stars['x_fwhm'] * scale, 'r.', label='X', alpha=0.5)
    ax1.plot(stars['mag'], stars['y_fwhm'] * scale, 'b.', label='Y', alpha=0.5)
    ax1.set_xlabel('Mag')
    ax1.set_ylabel('FWHM (")')
    ax1.legend()
    ax1.set_ylim(0, 2)

    # Plotting star FWHM vs. Rad
    r = np.hypot(stars['xcentroid'] - (stars['xcentroid'].max() / 2.0), 
                 stars['ycentroid'] - (stars['ycentroid'].max() / 2.0)) * scale / 60.0
    ax2.plot(r, stars['x_fwhm'] * scale, 'r.', label='X', alpha=0.5)
    ax2.plot(r, stars['y_fwhm'] * scale, 'b.', label='Y', alpha=0.5)
    ax2.set_xlabel('Radius (amin)')
    ax2.set_ylabel('FWHM (")')
    ax2.legend()
    ax2.set_ylim(0, 2)

    #plotting starlist
    xscale = scale
    yscale = 1.073 * scale
    xref = 2525
    yref = 2810
    #img, hdr = fits.getdata(root_dir + f'reduce/{fld}/' + test_img_base + '_clean.fits', header=True)
    x_arc = (stars['xcentroid'] - xref) * xscale
    y_arc = (stars['ycentroid'] - yref) * yscale
    img_x_arc = (np.arange(img.shape[1]) - xref) * xscale
    img_y_arc = (np.arange(img.shape[0]) - yref) * yscale
    extent = [img_x_arc.min(), img_x_arc.max(), img_y_arc.min(), img_y_arc.max()]
    ax3.imshow(img, vmin=0, vmax=500, origin='lower', cmap=plt.cm.gray, extent=extent) 
    idx = np.where((stars['mag'] > -8) & (stars['mag'] < -2))[0]
    ax3.plot(x_arc[idx], y_arc[idx], 'ro', mec='red', mfc='none', label='Starlist')
    #ax3.gca().invert_xaxis()
    #plt.xlim(5, -5)
    #plt.ylim(-5, 5)
    ax3.legend()

    # plotting emperical FWHM
    #stars2 = Table.read(root_dir + f'reduce/{fld}/{test_img_base}_clean_stars_stats.fits')
    ax4.plot(stars2['mag'], stars2['fwhm_emp'] * scale, 'k.', alpha=0.4)
    ax4.set_ylim(0, 1)
    ax4.set_xlabel('Mag')
    ax4.set_ylabel('Emp. FWHM (")')

    r = np.hypot(stars2['xcentroid'] - (stars2['xcentroid'].max() / 2.0), 
                 stars2['ycentroid'] - (stars2['ycentroid'].max() / 2.0)) * scale / 60.0
    ax5.plot(r, stars2['fwhm_emp'] * scale, 'k.', alpha=0.4)
    ax5.set_ylim(0, 1)
    ax5.set_xlabel('Radius (amin)')
    ax5.set_ylabel('Emp. FWHM (")')

    # automatically adjust padding horizontally
    # as well as vertically.
    ax3.set_title(night + " " + test_img_base + ' Starlist')
    plt.tight_layout()

    # display plot
    plt.show()

    fwhm_x_avg, fwhm_x_med, fwhm_x_std = sigma_clipped_stats(stars['x_fwhm'] * scale)
    fwhm_y_avg, fwhm_y_med, fwhm_y_std = sigma_clipped_stats(stars['y_fwhm'] * scale)
    print(f'x: fwhm_x_avg = {fwhm_x_avg:.2f}" fwhm_x_med = {fwhm_x_med:.2f}" fwhm_x_std = {fwhm_x_std:.2f}"')
    print(f'x: fwhm_x_avg = {fwhm_x_avg:.2f}" fwhm_x_med = {fwhm_x_med:.2f}" fwhm_x_std = {fwhm_x_std:.2f}"')
    
    return


def plot_ee_50(test_img_base, root_dir, fld, night):
    img, hdr = fits.getdata(root_dir + f'reduce/{fld}/' + test_img_base + '_clean.fits', header=True)
    scale = util.get_plate_scale(img, hdr)
    
    ee = Table.read(root_dir + f'reduce/{fld}/ee/' + test_img_base + '_clean_ee.txt', format='ascii')

    plt.plot(ee['Radius'], ee['EE'])
    plt.xlabel('Radius (arcsec)')
    plt.ylabel('Encircled Energy')

    fdx = np.where(ee['EE'] < 0.5)[0][-1]
    fwhm = ee['Radius'][fdx]

    plt.axhline(0.5, alpha=0.5)
    plt.axvline(fwhm, ls ='--',alpha=0.5)
    plt.title(f"EE for {night}, {test_img_base} \n 50% EE radius = {fwhm:.3f}")

    print(f'50% EE radius = {fwhm:.3f}"')
    return

def plot_moffat_fit(root_dir, c_key, fld, night, o_key = "_o"):
    stats_c = Table.read(f'{root_dir}reduce/stats/stats_{c_key}_mdp.fits')
    stats_o = Table.read(f'{root_dir}reduce/stats/stats_{o_key}_mdp.fits')
    
    plt.figure(figsize=(10,5))
    plt.clf()
    plt.subplot(121)
    plt.errorbar(stats_c['Minor Alpha'], stats_c['Beta'], xerr=stats_c['Minor Alpha std'], yerr=stats_c['Beta std'], fmt='r.', label='Closed (LS)')
    plt.errorbar(stats_o['Minor Alpha'], stats_o['Beta'], xerr=stats_o['Minor Alpha std'], yerr=stats_o['Beta std'], fmt='k.', label='Open')
    plt.xlabel(r'Minor $\alpha$')
    plt.ylabel(r'$\beta$')
    plt.legend()

    plt.subplot(122)
    plt.errorbar(stats_c['Phi'], stats_c['Beta'], xerr=stats_c['Phi std'], yerr=stats_c['Beta std'], fmt='r.', label='Closed (LS)')
    plt.errorbar(stats_o['Phi'], stats_o['Beta'], xerr=stats_o['Phi std'], yerr=stats_o['Beta std'], fmt='k.', label='Open')
    plt.xlabel(r'$\phi$')
    plt.ylabel(r'$\beta$')
    plt.legend()

    plt.suptitle(f"Moffat Fit Open vs. {c_key} Closed, \n summary {night}")
    return


def plot_moffat_fit_4F(fil_band, root_dir, c_key, fld, night, o_key = "_o"):
    stats_c = Table.read(f'{root_dir}reduce/stats/stats_{c_key}_{fil_band}_mdp.fits')
    stats_o = Table.read(f'{root_dir}reduce/stats/stats_{o_key}_{fil_band}_mdp.fits')
    
    plt.figure(figsize=(10,5))
    plt.clf()
    plt.subplot(121)
    plt.errorbar(stats_c['Minor Alpha'], stats_c['Beta'], xerr=stats_c['Minor Alpha std'], yerr=stats_c['Beta std'], fmt='r.', label='Closed (LS)')
    plt.errorbar(stats_o['Minor Alpha'], stats_o['Beta'], xerr=stats_o['Minor Alpha std'], yerr=stats_o['Beta std'], fmt='k.', label='Open')
    plt.xlabel(r'Minor $\alpha$')
    plt.ylabel(r'$\beta$')
    plt.ylim(1,10)
    plt.xlim(1,25)
    plt.legend()

    plt.subplot(122)
    plt.errorbar(stats_c['emp_fwhm'], stats_c['Beta'], xerr=stats_c['emp_fwhm_std'], yerr=stats_c['Beta std'], fmt='r.', label='Closed (LS)')
    plt.errorbar(stats_o['emp_fwhm'], stats_o['Beta'], xerr=stats_o['emp_fwhm_std'], yerr=stats_o['Beta std'], fmt='k.', label='Open')
    plt.xlabel(r'emp_fwhm')
    plt.ylabel(r'$\beta$')
    plt.ylim(1,10)
    plt.xlim(1,20)
    plt.legend()

    plt.suptitle(f"Moffat Fit Open vs. {c_key} Closed, \n filter {fil_band} summary {night}")
    return


def plot_stack_stats_4F_frame(date, suffixes=['open', 'ttf', 'closed'], root_dir='/Users/jlu/work/imaka/pleiades/', reduce_dir='fli/reduce/'):
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
    stats_dir = root_dir + date + '/' + reduce_dir + 'stats/'
    plots_dir = root_dir + date + '/' + reduce_dir + 'plots/'
    
    labels = ['B', 'V', 'R', 'I']
    dict_filt = {"B":"b", "V":"c", "R":"r", "I":"m"}
    f_fmt = {'LS_c':'x', 'docz2_c':'+', 'tt_c':'^', "_o":'o'}
    f_name = {"_o":'Open', 'LS_c':'Closed (LS)', 'docz2_c':'Closed (DOC)', 'tt_c':'Closed (tt)'}
        
    util.mkdir(plots_dir)
    
    suff_filt = []
    stats = []
    for suffix in suffixes:
        for f in labels:
            stats_file = f'{stats_dir}stats_{suffix}_{f}.fits'
            table_tmp = Table.read(stats_file)
            scale_tmp = table_tmp[0].meta['SCALE']
            
            FWHM_min, sig_FWHM_min, FWHM_maj, sig_FWHM_maj = moffat.calc_mof_fwhm(stats_file, filt=False, plate_scale=scale_tmp)
            
            N_frames = len(table_tmp)
            table_tmp.add_column(Column(FWHM_maj, name='FWHM_maj'))
            table_tmp.add_column(Column(FWHM_min, name='FWHM_min'))
            table_tmp.add_column(Column(np.repeat(f, N_frames), name='filter'))

            stats.append(table_tmp)
            suff_filt.append(f"{suffix}_{f}")

    scale = stats[0].meta['SCALE']
    colors = get_color_list()
    
    #
    # Plots for ratio of improvements. First we need to find a matching
    # closed loop image to go with each open (and TTF) image.
    #
    tree_indices = []
    tree_so = scipy.spatial.KDTree(np.array([stats[0]['Index']]).T)
    for ss in range(len(stats)):
        dist_ss, idx_ss = tree_so.query(np.array([stats[ss]['Index']]).T, 1)
        tree_indices.append(idx_ss)
    
    #####
    # MOFFAT FWHM vs. Frame
    #####
    
    plt.figure(1, figsize=(12, 4))
    plt.subplots_adjust(left=0.1, bottom=0.15)
    plt.clf()
    plt.subplot(121)
    for ss in range(len(suff_filt)):   
        plt.plot(stats[ss]['Index'], stats[ss]['FWHM_min'], marker=f_fmt[suffixes[ss//4]], color=dict_filt[labels[ss%4]], linestyle='none', label=suff_filt[ss], alpha=0.7)
    plt.xlabel('Frame Number')
    plt.ylabel('Moffat-Fit FWHM_min (")')
    #plt.legend(numpoints=1)
    plt.legend(handles=[mlines.Line2D([0], [0], marker=value, label=f_name[key], color='grey',  lw=0) for key, value in f_fmt.items()])
    plt.ylim(0, 1.5)
    plt.title(date)

    plt.subplot(122)
    for ss in range(4, len(suff_filt)):
        idx = tree_indices[ss]
        label = suff_filt[ss] + " / " + suff_filt[ss%4]
        plt.plot(stats[0]['Index'][idx], stats[ss]['FWHM_min'] / stats[ss%4]['FWHM_min'][idx],  linestyle='none', marker=f_fmt[suffixes[ss//4]], color=dict_filt[labels[ss%4]], label=label, alpha=0.7)
    plt.legend(handles=[mpatches.Patch(color=value, label=key) for key, value in dict_filt.items()], loc=2)
    plt.xlabel('Frame Number')
    plt.ylabel('Ratio of Moffat-Fit FWHM_min')
    #plt.legend(numpoints=1)
    plt.axhline(1.0, color='k', linestyle='--', linewidth=2)
    plt.ylim(0, 1.3)
    plt.savefig(plots_dir + 'fwhm_vs_frame' + suffix + '.png')
    
    #####
    # Empirical FWHM
    #####
    plt.figure(2, figsize=(12, 4))
    plt.subplots_adjust(left=0.1, bottom=0.15)
    plt.clf()
    plt.subplot(121)
    for ss in range(len(suff_filt)):
        plt.plot(stats[ss]['Index'], stats[ss]['emp_fwhm']*scale, marker=f_fmt[suffixes[ss//4]], color=dict_filt[labels[ss%4]], linestyle='none', label=suff_filt[ss], alpha=0.7)
    plt.xlabel('Frame Number')
    plt.ylabel('Empirical FWHM (")')
    #plt.legend(numpoints=1)
    plt.legend(handles=[mlines.Line2D([0], [0], marker=value, label=f_name[key], color='grey',  lw=0) for key, value in f_fmt.items()])
    plt.ylim(0, 1.5)

    plt.subplot(122)
    for ss in range(4, len(suff_filt)):
        idx = tree_indices[ss]
        label = suff_filt[ss] + " / " + suff_filt[ss%4]
        plt.plot(stats[0]['Index'][idx], stats[ss]['emp_fwhm'] / stats[ss%4]['emp_fwhm'][idx], linestyle='none', 
                 marker=f_fmt[suffixes[ss//4]], color=dict_filt[labels[ss%4]], label=label, alpha=0.7)
    plt.xlabel('Frame Number')
    plt.ylabel('Ratio of Empircal FWHM')
    #plt.legend(numpoints=1)
    plt.legend(handles=[mlines.Line2D([0], [0], marker=value, label=f_name[key]+" / Open", color='grey',  lw=0) for key, value in f_fmt.items()][:ss//4])
    plt.axhline(1.0, color='k', linestyle='--', linewidth=2)
    plt.ylim(0, 1.3)
    plt.savefig(plots_dir + 'efwhm_vs_frame' + suffix + '.png')
    
    
    #####
    # EE 50
    #####
    plt.figure(3, figsize=(12, 4))
    plt.subplots_adjust(left=0.1, bottom=0.15)
    plt.clf()
    plt.subplot(121)
    for ss in range(len(suff_filt)):
        plt.plot(stats[ss]['Index'], stats[ss]['EE50'], marker=f_fmt[suffixes[ss//4]], 
                 color=dict_filt[labels[ss%4]], linestyle='none', label=suff_filt[ss], alpha=0.7)
    plt.xlabel('Frame Number')
    plt.ylabel('Radius of 50% Encircled Energy (")')
    #plt.legend(numpoints=1)
    plt.legend(handles=[mlines.Line2D([0], [0], marker=value, label=f_name[key], color='grey',  lw=0) 
                        for key, value in f_fmt.items()])
    plt.ylim(0, 1.5)
    plt.title(date)
    
    plt.subplot(122)
    for ss in range(4, len(suff_filt)):
        idx = tree_indices[ss]
        label = suff_filt[ss] + " / " + suff_filt[ss%4]
        plt.plot(stats[0]['Index'][idx], stats[ss]['EE50'] / stats[ss%4]['EE50'][idx], linestyle='none', 
                 marker=f_fmt[suffixes[ss//4]], color=dict_filt[labels[ss%4]], label=label, alpha=0.7)
    plt.xlabel('Frame Number')
    plt.ylabel('Ratio of 50% EE Radius')
    #plt.legend(numpoints=1)
    plt.legend(handles=[mlines.Line2D([0], [0], marker=value, label=f_name[key]+" / Open", color='grey',  lw=0) for key, value in f_fmt.items()][:ss//4])
    plt.axhline(1.0, color='k', linestyle='--', linewidth=2)
    plt.ylim(0, 1.3)
    plt.savefig(plots_dir + 'ee50_vs_frame' + suffix + '.png')

#####
# EE 80
#####
    plt.figure(4, figsize=(12, 4))
    plt.subplots_adjust(left=0.1, bottom=0.15)
    plt.clf()
    plt.subplot(121)
    for ss in range(len(suff_filt)):
        plt.plot(stats[ss]['Index'], stats[ss]['EE80'], marker=f_fmt[suffixes[ss//4]], 
                 color=dict_filt[labels[ss%4]], linestyle='none', label=suff_filt[ss], alpha=0.7)
    plt.xlabel('Frame Number')
    plt.ylabel('Radius of 80% Encircled Energy (")')
    #plt.legend(numpoints=1)
    plt.legend(handles=[mlines.Line2D([0], [0], marker=value, label=f_name[key], color='grey',  lw=0) 
                        for key, value in f_fmt.items()])
    plt.ylim(0, 1.5)
    plt.title(date)
    
    plt.subplot(122)
    for ss in range(4, len(suff_filt)):
        idx = tree_indices[ss]
        label = suff_filt[ss] + " / " + suff_filt[ss%4]
        plt.plot(stats[0]['Index'][idx], stats[ss]['EE80'] / stats[ss%4]['EE80'][idx], linestyle='none', 
                 marker=f_fmt[suffixes[ss//4]], color=dict_filt[labels[ss%4]], label=label, alpha=0.7)
    plt.xlabel('Frame Number')
    plt.ylabel('Ratio of 80% EE Radius')
    #plt.legend(numpoints=1)
    plt.legend(handles=[mlines.Line2D([0], [0], marker=value, label=f_name[key]+" / Open", color='grey',  lw=0) for key, value in f_fmt.items()][:ss//4])
    plt.axhline(1.0, color='k', linestyle='--', linewidth=2)
    plt.ylim(0, 1.3)
    plt.savefig(plots_dir + 'ee80_vs_frame' + suffix + '.png')
    
    #####
    # NEA
    #####
    plt.figure(5, figsize=(12, 4))
    plt.subplots_adjust(left=0.1, bottom=0.15)
    plt.clf()
    plt.subplot(121)
    for ss in range(len(suff_filt)):
        plt.plot(stats[ss]['Index'], stats[ss]['NEA'] * scale**2, marker=f_fmt[suffixes[ss//4]], 
                 color=dict_filt[labels[ss%4]], linestyle='none', label=suff_filt[ss], alpha=0.7)
    plt.xlabel('Frame Number')
    plt.ylabel('NEA (Sq. Arcsec)')
    #plt.legend(numpoints=1)
    plt.legend(handles=[mlines.Line2D([0], [0], marker=value, label=f_name[key], color='grey',  lw=0) 
                        for key, value in f_fmt.items()])
    plt.ylim(0, 5)
    plt.title(date)
    
    plt.subplot(122)
    for ss in range(4, len(suff_filt)):
        idx = tree_indices[ss]
        label = suff_filt[ss] + " / " + suff_filt[ss%4]
        plt.plot(stats[0]['Index'][idx], stats[ss]['NEA'] / stats[ss%4]['NEA'][idx], linestyle='none', 
                 marker=f_fmt[suffixes[ss//4]], color=dict_filt[labels[ss%4]], label=label, alpha=0.7)
    plt.xlabel('Frame Number')
    plt.ylabel('Ratio of NEA')
    #plt.legend(numpoints=1)
    plt.legend(handles=[mlines.Line2D([0], [0], marker=value, label=f_name[key]+" / Open", color='grey',  lw=0) for key, value in f_fmt.items()][:ss//4])
    plt.axhline(1.0, color='k', linestyle='--', linewidth=2)
    plt.ylim(0, 1.5)
    plt.savefig(plots_dir + 'nea_vs_frame' + suffix + '.png')
    
    #####
    # NEA2
    #####
    plt.figure(6, figsize=(12, 4))
    plt.subplots_adjust(left=0.1, bottom=0.15)
    plt.clf()
    plt.subplot(121)
    for ss in range(len(suff_filt)):
        plt.plot(stats[ss]['Index'], stats[ss]['NEA2']*scale**2, marker=f_fmt[suffixes[ss//4]], 
                 color=dict_filt[labels[ss%4]], linestyle='none', label=suff_filt[ss], alpha=0.7)
    plt.xlabel('Frame Number')
    plt.ylabel('NEA2 (Sq. Arcsec)')
    #plt.legend(numpoints=1)
    plt.legend(handles=[mlines.Line2D([0], [0], marker=value, label=f_name[key], color='grey',  lw=0) 
                        for key, value in f_fmt.items()])
    plt.ylim(0, 5)
    plt.title(date)
    
    plt.subplot(122)
    for ss in range(4, len(suff_filt)):
        idx = tree_indices[ss]
        label = suff_filt[ss] + " / " + suff_filt[ss%4]
        plt.plot(stats[0]['Index'][idx], stats[ss]['NEA2'] / stats[ss%4]['NEA2'][idx], linestyle='none', 
                 marker=f_fmt[suffixes[ss//4]], color=dict_filt[labels[ss%4]], label=label, alpha=0.7)
    plt.xlabel('Frame Number')
    plt.ylabel('Ratio of NEA2')
    #plt.legend(numpoints=1)
    plt.legend(handles=[mlines.Line2D([0], [0], marker=value, label=f_name[key]+" / Open", color='grey',  lw=0) for key, value in f_fmt.items()][:ss//4])
    plt.axhline(1.0, color='k', linestyle='--', linewidth=2)
    plt.ylim(0, 1.5)
    plt.savefig(plots_dir + 'nea2_vs_frame' + suffix + '.png')

#####
# FWHM for each direction
#####
    plt.figure(7, figsize=(12, 8))
    plt.subplots_adjust(left=0.1, bottom=0.15)
    plt.clf()
    plt.subplot(121)
    for ss in range(len(suff_filt)):
        c = np.take(colors, ss, mode='wrap')
        plt.plot(stats[ss]['Index'], stats[ss]['xFWHM']*scale, marker=f_fmt[suffixes[ss//4]],   
                 color=dict_filt[labels[ss%4]], linestyle='none', label='X ' + suff_filt[ss])
        plt.plot(stats[ss]['Index'], stats[ss]['yFWHM']*scale, marker=f_fmt[suffixes[ss//4]], 
                 color=dict_filt[labels[ss%4]], linestyle='none', label='Y ' + suff_filt[ss], alpha = 0.5)
    plt.xlabel('Frame Number')
    plt.ylabel('Gaussian-Fit FWHM (")')
    #plt.legend(numpoints=1, fontsize=10)
    handles_xy = [mlines.Line2D([0], [0], marker=value, label="x "+f_name[key], color='grey',  lw=0) for key, value in f_fmt.items()] + [mlines.Line2D([0], [0], marker=value, label="y "+f_name[key], color='grey',  lw=0, alpha=0.3) for key, value in f_fmt.items()]
    plt.legend(handles=handles_xy)
    plt.ylim(0, 1.5)
    plt.title(date)
    
    plt.subplot(122)
    for ss in range(4, len(suff_filt)):
        idx = tree_indices[ss]
        label = suff_filt[ss] + " / " + suff_filt[ss%4]
        plt.plot(stats[0]['Index'][idx], stats[ss]['xFWHM'] / stats[ss%4]['xFWHM'][idx], 
                 marker=f_fmt[suffixes[ss//4]], color=dict_filt[labels[ss%4]], linestyle='none', label=label)
        plt.plot(stats[0]['Index'][idx], stats[ss]['yFWHM'] / stats[ss%4]['yFWHM'][idx], 
                 marker=f_fmt[suffixes[ss//4]], color=dict_filt[labels[ss%4]],linestyle='none', label=label, alpha = 0.5)
    plt.xlabel('Frame Number')
    plt.ylabel('Ratio of Gaussian-Fit FWHM')
    #plt.legend(numpoints=1, fontsize=10)
    handles_xy = [mlines.Line2D([0], [0], marker=value, label="x "+f_name[key]+" / Open", color='grey',  lw=0) for key, value in f_fmt.items()] + [mlines.Line2D([0], [0], marker=value, label="y "+f_name[key]+" / Open", color='grey',  lw=0, alpha=0.5) for key, value in f_fmt.items()]
    plt.legend(handles=handles_xy)
    plt.axhline(1.0, color='k', linestyle='--', linewidth=2)
    plt.ylim(0, 1.3)
    plt.savefig(plots_dir + 'xyfwhm_vs_frame' + suffix + '.png')
    
    return



def plot_stack_stats_4F_time(date, suffixes=['open', 'ttf', 'closed'], root_dir='/Users/jlu/work/imaka/pleiades/', reduce_dir='fli/reduce/'):
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
    stats_dir = root_dir + date + '/' + reduce_dir + 'stats/'
    plots_dir = root_dir + date + '/' + reduce_dir + 'plots/'
    
    labels = ['B', 'V', 'R', 'I']
    dict_filt = {"B":"b", "V":"c", "R":"r", "I":"m"}
    f_fmt = {'LS_c':'x', 'docz2_c':'+', 'tt_c':'^', "_o":'o'}
    f_name = {"_o":'Open', 'LS_c':'Closed (LS)', 'docz2_c':'Closed (DOC)', 'tt_c':'Closed (tt)'}
        
    util.mkdir(plots_dir)
    
    suff_filt = []
    stats = []
    for suffix in suffixes:
        for f in labels:
            stats_file = f'{stats_dir}stats_{suffix}_{f}.fits'
            table_tmp = Table.read(stats_file)
            scale_tmp = table_tmp[0].meta['SCALE']
            
            FWHM_min, sig_FWHM_min, FWHM_maj, sig_FWHM_maj = moffat.calc_mof_fwhm(stats_file, filt=False, plate_scale=scale_tmp)
            
            N_frames = len(table_tmp)
            table_tmp.add_column(Column(FWHM_maj, name='FWHM_maj'))
            table_tmp.add_column(Column(FWHM_min, name='FWHM_min'))
            table_tmp.add_column(Column(np.repeat(f, N_frames), name='filter'))

            stats.append(table_tmp)
            suff_filt.append(f"{suffix}_{f}")

    scale = stats[0].meta['SCALE']
    colors = get_color_list()
    
    #
    # Plots for ratio of improvements. First we need to find a matching
    # closed loop image to go with each open (and TTF) image.
    #
    tree_indices = []
    tree_so = scipy.spatial.KDTree(np.array([stats[0]['Index']]).T)
    for ss in range(len(stats)):
        dist_ss, idx_ss = tree_so.query(np.array([stats[ss]['Index']]).T, 1)
        tree_indices.append(idx_ss)

##########
# All plots vs. time.
##########

    utcs = []
    for ss in range(len(suff_filt)):
        utc_dt = [datetime.strptime(stats[ss]['TIME_UTC'][ii], '%H:%M:%S') for ii in range(len(stats[ss]))]
        utcs.append(utc_dt)

    time_fmt = mp_dates.DateFormatter('%H:%M')
    
    
    #####
    # FWHM
    #####
    plt.figure(8, figsize=(12, 6))
    plt.subplots_adjust(left=0.1, bottom=0.15)
    plt.clf()
    plt.subplot(121)
    for ss in range(len(suff_filt)):
        plt.plot(utcs[ss], stats[ss]['FWHM']*scale, marker='o', linestyle='none', label=suff_filt[ss])
    plt.gca().xaxis.set_major_formatter(time_fmt)
    plt.xticks(rotation=35)
    plt.xlabel('UTC Time (hr)')
    plt.ylabel('Gaussian-Fit FWHM (")')
    plt.legend(numpoints=1)
    plt.ylim(0, 1.5)
    plt.title(date)
    
    plt.subplot(122)
    for ss in range(4, len(suff_filt)):
        idx = tree_indices[ss]
        label = suff_filt[ss] + " / " + suff_filt[ss%4]
        utc = [utcs[0][ii] for ii in idx]
        plt.plot(utc, stats[ss]['FWHM'] / stats[ss%4]['FWHM'][idx], marker='o', linestyle='none', label=label)
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
    for ss in range(len(suff_filt)):
        plt.plot(utcs[ss], stats[ss]['emp_fwhm']*scale, marker='o', linestyle='none', label=suff_filt[ss])
    plt.gca().xaxis.set_major_formatter(time_fmt)
    plt.xticks(rotation=35)
    plt.xlabel('UTC Time (hr)')
    plt.ylabel('Empirical FWHM (")')
    plt.legend(numpoints=1)
    plt.ylim(0, 1.5)
    plt.title(date)
    
    plt.subplot(122)
    for ss in range(4, len(suff_filt)):
        idx = tree_indices[ss]
        label = suff_filt[ss] + " / " + suff_filt[ss%4]
        utc = [utcs[0][ii] for ii in idx]
        plt.plot(utc, stats[ss]['emp_fwhm'] / stats[ss%4]['emp_fwhm'][idx], marker='o', linestyle='none', label=label)
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
    for ss in range(len(suff_filt)):
        plt.plot(utcs[ss], stats[ss]['EE50'], marker='o', linestyle='none', label=suff_filt[ss])
    plt.gca().xaxis.set_major_formatter(time_fmt)
    plt.xticks(rotation=35)
    plt.xlabel('UTC Time (hr)')
    plt.ylabel('Radius of 50% Encircled Energy (")')
    plt.legend(numpoints=1)
    plt.ylim(0, 1.5)
    plt.title(date)
    
    plt.subplot(122)
    for ss in range(4, len(suff_filt)):
        idx = tree_indices[ss]
        label = suff_filt[ss] + " / " + suff_filt[ss%4]
        utc = [utcs[0][ii] for ii in idx]
        plt.plot(utc, stats[ss]['EE50'] / stats[ss%4]['EE50'][idx], marker='o', linestyle='none', label=label)
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
    for ss in range(len(suff_filt)):
        plt.plot(utcs[ss], stats[ss]['EE80'], marker='o', linestyle='none', label=suff_filt[ss])
    plt.gca().xaxis.set_major_formatter(time_fmt)
    plt.xticks(rotation=35)
    plt.xlabel('UTC Time (hr)')
    plt.ylabel('Radius of 80% Encircled Energy (")')
    plt.legend(numpoints=1)
    plt.ylim(0, 1.5)
    plt.title(date)
    
    plt.subplot(122)
    for ss in range(4, len(suff_filt)):
        idx = tree_indices[ss]
        label = suff_filt[ss] + " / " + suff_filt[ss%4]
        utc = [utcs[0][ii] for ii in idx]
        plt.plot(utc, stats[ss]['EE80'] / stats[ss%4]['EE80'][idx], marker='o', linestyle='none', label=label)
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
    for ss in range(len(suff_filt)):
        plt.plot(utcs[ss], stats[ss]['NEA']*scale**2, marker='o', linestyle='none', label=suff_filt[ss])
    plt.gca().xaxis.set_major_formatter(time_fmt)
    plt.xticks(rotation=35)
    plt.xlabel('UTC Time (hr)')
    plt.ylabel('NEA (Sq. Arcsec)')
    plt.legend(numpoints=1)
    plt.ylim(0, 5)
    plt.title(date)
    
    plt.subplot(122)
    for ss in range(4, len(suff_filt)):
        idx = tree_indices[ss]
        label = suff_filt[ss] + " / " + suff_filt[ss%4]
        utc = [utcs[0][ii] for ii in idx]
        plt.plot(utc, stats[ss]['NEA'] / stats[ss%4]['NEA'][idx], marker='o', linestyle='none', label=label)
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
    for ss in range(len(suff_filt)):
        plt.plot(utcs[ss], stats[ss]['NEA2']*scale**2, marker='o', linestyle='none', label=suff_filt[ss])
    plt.gca().xaxis.set_major_formatter(time_fmt)
    plt.xticks(rotation=35)
    plt.xlabel('UTC Time (hr)')
    plt.ylabel('NEA2 (Sq. Arcsec)')
    plt.legend(numpoints=1)
    plt.ylim(0, 5)
    plt.title(date)
    
    plt.subplot(122)
    for ss in range(4, len(suff_filt)):
        idx = tree_indices[ss]
        label = suff_filt[ss] + " / " + suff_filt[ss%4]
        utc = [utcs[0][ii] for ii in idx]
        plt.plot(utc, stats[ss]['NEA2'] / stats[ss%4]['NEA2'][idx], marker='o', linestyle='none', label=label)
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
    for ss in range(len(suff_filt)):
        c = np.take(colors, ss, mode='wrap')
        plt.plot(utcs[ss], stats[ss]['xFWHM']*scale, marker='o', color=c, linestyle='none', label='X ' + suff_filt[ss])
        plt.plot(utcs[ss], stats[ss]['yFWHM']*scale, marker='^', color=c, linestyle='none', label='Y ' + suff_filt[ss])
    plt.gca().xaxis.set_major_formatter(time_fmt)
    plt.xticks(rotation=35)
    plt.xlabel('UTC Time (hr)')
    plt.ylabel('Gaussian-Fit FWHM (")')
    plt.legend(numpoints=1, fontsize=10)
    plt.ylim(0, 1.5)
    plt.title(date)
    
    plt.subplot(122)
    for ss in range(4, len(suff_filt)):
        c = np.take(colors, ss, mode='wrap')
        idx = tree_indices[ss]
        label = suff_filt[ss] + " / " + suff_filt[ss%4]
        utc = [utcs[0][ii] for ii in idx]
        plt.plot(utc, stats[ss]['xFWHM'] / stats[ss%4]['xFWHM'][idx], marker='o', color=c, linestyle='none', label=label)
        plt.plot(utc, stats[ss]['yFWHM'] / stats[ss%4]['yFWHM'][idx], marker='^', color=c, linestyle='none', label=label)
    plt.gca().xaxis.set_major_formatter(time_fmt)
    plt.xticks(rotation=35)
    plt.xlabel('UTC Time (hr)')
    plt.ylabel('Ratio of Gaussian-Fit FWHM')
    plt.legend(numpoints=1, fontsize=10)
    plt.axhline(1.0, color='k', linestyle='--', linewidth=2)
    plt.ylim(0, 1.3)
    plt.savefig(plots_dir + 'xyfwhm_vs_time' + suffix + '.png')
    
    return




def plot_fwhmvt_nomatch_filt(open_file, closed_file, comp_col, title, plots_dir, wav_given):
    """
    Plot FWHM vs. time for both open, closed, MASS, and DIMM
    
    Inputs
    ----------
    open_file and closed_file: str
        The fits table stats files
        
    comp_col: str
        what column of data to compare (e.g. 'emp_fwhm')

    title: str
        title for generated plot
        
    plots_dir: str
        directory to put generated plot in
    """
    
    #Read in data
    stats1 = Table.read(open_file)
    stats2 = Table.read(closed_file)
    
    time1, date1, data1, data2, err1, err2 = add_data.match_cols(open_file, open_file, comp_col)
    time2, date2, data1, data2, err1, err2 = add_data.match_cols(closed_file, closed_file, comp_col)
    
    calib1 = []
    calib2 = []
    
    # Get mass/dimm data
    mass = stats2['MASS']
    dimm = stats1['DIMM']
    wvln = wav_given

    # Shift to match wavelengths
    for i in range(len(stats1)):
       # wvln = filter2wv(stats1['FILTER'][i])
        scale = stats1.meta['SCALE']
        factor = ((500/wvln)**0.2) * scale
        calib1.append(factor)
    
    for i in range(len(stats2)):
        #wvln = filter2wv(stats2['FILTER'][i])
        scale = stats2.meta['SCALE']
        factor = ((500/wvln)**0.2) * scale
        calib2.append(factor)

    calib1 = np.array(calib1)
    calib2 = np.array(calib2)

    comp_col_err = comp_col + '_std'
    if comp_col_err in stats1.colnames:
        open_err = stats1[comp_col_err] * calib1
        closed_err = stats2[comp_col_err] * calib2
    else:
        open_err = np.zeros(len(stats1))
        closed_err = np.zeros(len(stats2))

    # Plot fwhm and seeing vs time
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
    plt.errorbar(times1, stats1[comp_col]*calib1, yerr=open_err, fmt='o', label="Open")
    plt.errorbar(times2, stats2[comp_col]*calib2, yerr=closed_err, fmt='ro', label="Closed")
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

def get_color_list():
    rcParams = plt.matplotlib.rcParams
    prop_cycler = rcParams['axes.prop_cycle']
    if prop_cycler is None and 'axes.color_cycle' in rcParams:
        clist = rcParams['axes.color_cycle']
        prop_cycler = cycler('color', clist)
    
    colors = [item['color'] for item in prop_cycler]
    return colors