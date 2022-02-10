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
import scipy.linalg
from matplotlib import cm
from astropy.stats import sigma_clipped_stats


# for plotting and checking one files starlist
def plot_starlist_stats(test_img_base, root_dir, fld, night):

    # creating grid for subplots
    fig, ax = plt.subplots(ncols=2, nrows=1, figsize=(20,10))

    ax1 = plt.subplot2grid((3, 3),(0,0))
    ax2 = plt.subplot2grid((3, 3), (1, 0))
    ax3 = plt.subplot2grid((3, 3), (0, 1), rowspan=2, colspan=1)
    ax4 = plt.subplot2grid((3, 3), (0, 2))
    ax5 = plt.subplot2grid((3, 3), (1, 2))

    # Plotting star FWHM vs. Mag
    img, hdr = fits.getdata(root_dir + f'{fld}/' + test_img_base + '.fits', header=True)
    scale = util.get_plate_scale(img, hdr)
    del img
    del hdr
    stars = Table.read(root_dir + f'reduce/{fld}/' + test_img_base + '_clean_stars.txt', format='ascii')

    ax1.plot(stars['mag'], stars['x_fwhm'] * scale, 'r.', label='X', alpha=0.5)
    ax1.plot(stars['mag'], stars['y_fwhm'] * scale, 'b.', label='Y', alpha=0.5)
    ax1.set_xlabel('Mag')
    ax1.set_ylabel('FWHM (")')
    ax1.legend()
    ax1.set_ylim(0, 1.5)

    # Plotting star FWHM vs. Rad
    r = np.hypot(stars['xcentroid'] - (stars['xcentroid'].max() / 2.0), 
                 stars['ycentroid'] - (stars['ycentroid'].max() / 2.0)) * scale / 60.0
    ax2.plot(r, stars['x_fwhm'] * scale, 'r.', label='X', alpha=0.5)
    ax2.plot(r, stars['y_fwhm'] * scale, 'b.', label='Y', alpha=0.5)
    ax2.set_xlabel('Radius (amin)')
    ax2.set_ylabel('FWHM (")')
    ax2.legend()
    ax2.set_ylim(0, 1.5)

    #plotting starlist
    xscale = scale
    yscale = 1.073 * scale
    xref = 2525
    yref = 2810
    img, hdr = fits.getdata(root_dir + f'reduce/{fld}/' + test_img_base + '_clean.fits', header=True)
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
    stars2 = Table.read(root_dir + f'reduce/{fld}/{test_img_base}_clean_stars_stats.fits')
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