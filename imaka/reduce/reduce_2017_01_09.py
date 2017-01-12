import pylab as plt
import matplotlib.patches as mpatches
import numpy as np
from astropy.io import fits
from astropy import table
from astropy import units
from astropy.modeling.models import Ellipse2D
import glob
from imaka.reduce import reduce_fli
from imaka.reduce import calib
import pdb
import os
from flystar import match
import shutil

def make_sky():
    # Didn't take skies, so just copy over the one from Tuesday night.
    old_sky = '/Users/jlu/data/imaka/2017_01_10/fli/Pleiades/pleiades_sky.fits'
    new_sky = '/Users/jlu/data/imaka/2017_01_09/fli/Pleiades_E/pleiades_sky.fits'

    shutil.copyfile(old_sky, new_sky)
    
        
def reduce_pleiades_binned_open():
    sky_dir = '/Users/jlu/data/imaka/2017_01_09/fli/Pleiades_E/'
    data_dir = '/Users/jlu/data/imaka/2017_01_09/fli/Pleiades_E/'
    os.chdir(data_dir)

    fnum = np.arange(57, 67)
    img_files = ['obj{0:03d}.fits'.format(ii) for ii in fnum]

    reduce_fli.clean_images(img_files, rebin=1, sky_frame=sky_dir + 'pleiades_sky.fits')

    return
    

def reduce_pleiades_binned_closed():
    sky_dir = '/Users/jlu/data/imaka/2017_01_09/fli/Pleiades_E/'
    data_dir = '/Users/jlu/data/imaka/2017_01_09/fli/Pleiades_E/'
    os.chdir(data_dir)

    fnum = np.arange(47, 57)
    img_files = ['obj{0:03d}.fits'.format(ii) for ii in fnum]

    reduce_fli.clean_images(img_files, rebin=1, sky_frame=sky_dir + 'pleiades_sky.fits')

    return
    
def find_stars_pleiades_binned_open():
    data_dir = '/Users/jlu/data/imaka/2017_01_09/fli/Pleiades_E/'
    os.chdir(data_dir)
    
    fnum = np.arange(57, 67)
    img_files = ['obj{0:03d}_bin_nobkg.fits'.format(ii) for ii in fnum]

    reduce_fli.find_stars_bin(img_files, fwhm=2, threshold=6)

    return
    
def find_stars_pleiades_binned_closed():
    data_dir = '/Users/jlu/data/imaka/2017_01_09/fli/Pleiades_E/'
    os.chdir(data_dir)
    
    fnum = np.arange(47, 57)
    img_files = ['obj{0:03d}_bin_nobkg.fits'.format(ii) for ii in fnum]

    reduce_fli.find_stars_bin(img_files, fwhm=2, threshold=6)

    return


def compare_fwhm_list():
    data_dir = '/Users/jlu/data/imaka/2017_01_09/fli/Pleiades_E/'
    os.chdir(data_dir)
    
    o_list = np.arange(57, 67) # Open
    c_list = np.arange(47, 57) # Closed
    
    plt.ion()

    for ii in range(len(o_list)):
        open_list = 'obj{0:03d}_bin_nobkg_stars.txt'.format(o_list[ii])
        closed_list = 'obj{0:03d}_bin_nobkg_stars.txt'.format(c_list[ii])

        compare_fwhm(open_list, closed_list, scale=3*0.04, flux_min=4)
        pdb.set_trace()

def compare_fwhm(open_list, closed_list, scale=0.04, flux_min=2.0):
    topen = table.Table.read(open_list, format='ascii')
    tclose = table.Table.read(closed_list, format='ascii')

    print('')
    print('')
    print('Comparing:')
    print('    Open List:   ' + open_list)
    print('    Closed List: ' + closed_list)
    print("N_open = {0:3d}  N_close = {1:3d} in original lists".format(len(topen), len(tclose)))

    # Trim out any stars with FWHM > 20
    idx_o = np.where(topen['x_fwhm'] < 20)[0]
    idx_c = np.where(tclose['x_fwhm'] < 20)[0]

    topen = topen[idx_o]
    tclose = tclose[idx_c]
    
    print("N_open = {:3d}  N_close = {:3d} after trimming fwhm<20".format(len(topen), len(tclose)))

    # Trim out any stars with flux < flux_min
    idx_o = np.where(topen['flux'] > flux_min)[0]
    idx_c = np.where(tclose['flux'] > flux_min)[0]

    topen = topen[idx_o]
    tclose = tclose[idx_c]
    
    print("N_open = {:3d}  N_close = {:3d} after trimming low flux sources".format(len(topen), len(tclose)))
    
    m_c = np.ones(len(tclose))
    m_o = np.ones(len(topen))

    idx_c, idx_o, dr, dm = match.match(tclose['xcentroid'], tclose['ycentroid'], m_c,
                                       topen['xcentroid'], topen['ycentroid'], m_o,
                                       dr_tol=10)
    # Matched catalogs
    to_match = topen[idx_o]
    tc_match = tclose[idx_c]
    
    print("N_open = {:3d}  N_close = {:3d} after matching".format(len(to_match), len(tc_match)))

    # Plot
    plt.figure(1)
    plt.clf()
    plt.plot(to_match['x_fwhm']*scale, tc_match['x_fwhm']*scale, 'r.', label='Minor')
    plt.plot(to_match['y_fwhm']*scale, tc_match['y_fwhm']*scale, 'b.', label='Major')
    plt.plot([0, 10], [0, 10], 'k--')
    plt.xlabel('FWHM in Open Loop (")')
    plt.ylabel('FWHM in Closed Loop (")')
    plt.axis('equal')

    max_fwhm = 1.5*np.mean([to_match['x_fwhm'], to_match['y_fwhm'], tc_match['x_fwhm'], tc_match['y_fwhm']])
    plt.ylim(0, max_fwhm*scale)
    plt.xlim(0, max_fwhm*scale)
    plt.legend(numpoints=1, loc='upper left')
    plt.pause(0.05)

    plt.figure(2)
    plt.clf()
    x_ratio = tc_match['x_fwhm'] / to_match['x_fwhm']
    y_ratio = tc_match['y_fwhm'] / to_match['y_fwhm']
    
    plt.plot(to_match['x_fwhm']*scale, x_ratio, 'r.', label='Minor')
    plt.plot(to_match['y_fwhm']*scale, y_ratio, 'b.', label='Major')
    plt.axhline(1, linestyle='--', color='black')
    plt.xlabel('FWHM in Open Loop (")')
    plt.ylabel('Closed / Open FWHM')

    max_fwhm = 1.5*np.mean([to_match['x_fwhm'], to_match['y_fwhm'], tc_match['x_fwhm'], tc_match['y_fwhm']])
    plt.xlim(0, max_fwhm*scale)
    plt.ylim(0, 1.5)
    plt.legend(numpoints=1, loc='upper left')
    plt.pause(0.05)

    x_fwhm_o_med = np.median(to_match['x_fwhm'])    
    y_fwhm_o_med = np.median(to_match['y_fwhm'])    
    x_fwhm_c_med = np.median(tc_match['x_fwhm'])    
    y_fwhm_c_med = np.median(tc_match['y_fwhm'])
    theta_o_med = np.median(to_match['theta'])
    theta_c_med = np.median(tc_match['theta'])

    ### Plot Gaussian ellipses to show the closed and open loop FWHM.
    plt.figure(3)
    plt.clf()
    ec = mpatches.Ellipse((0, 0), x_fwhm_c_med * scale, y_fwhm_c_med * scale, theta_c_med,
                              edgecolor='red', facecolor='none', lw=4, label='Closed')
    eo = mpatches.Ellipse((0, 0), x_fwhm_o_med * scale, y_fwhm_o_med * scale, theta_o_med,
                              edgecolor='blue', facecolor='none', lw=4, label='Open')
    ax = plt.gca()
    ax.add_patch(ec)
    ax.add_patch(eo)
    plt.xlim(-1, 1)
    plt.ylim(-1, 1)
    plt.legend()
    plt.xlabel('X (")')
    plt.ylabel('Y (")')
    plt.pause(0.05)

    print('Open Loop Stats:')
    print('\t Median x_fwhm = {0:.2f} +/- {1:.2f}'.format(x_fwhm_o_med * scale,
                                                          to_match['x_fwhm'].std() * scale))
    print('\t Median y_fwhm = {0:.2f} +/- {1:.2f}'.format(y_fwhm_o_med * scale,
                                                          to_match['y_fwhm'].std() * scale))

    print('Closed Loop Stats:')
    print('\t Median x_fwhm = {0:.2f} +/- {1:.2f}'.format(x_fwhm_c_med * scale,
                                                          tc_match['x_fwhm'].std() * scale))
    print('\t Median y_fwhm = {0:.2f} +/- {1:.2f}'.format(y_fwhm_c_med * scale,
                                                          tc_match['y_fwhm'].std() * scale))

    print('Fractional Improvement')
    print('\t x_fwhm closed / open = {0:.2f}'.format(x_fwhm_c_med / x_fwhm_o_med))
    print('\t y_fwhm closed / open = {0:.2f}'.format(y_fwhm_c_med / y_fwhm_o_med))
    
