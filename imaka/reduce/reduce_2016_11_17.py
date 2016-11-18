import pylab as plt
import numpy as np
from astropy.io import fits
from astropy import table
from astropy import units
import glob
from imaka.reduce import reduce_fli
import pdb
import os
from flystar import match

def reduce_pleiades_unbinned_closed():
    data_dir = '/Users/jlu/data/imaka/2016_11_18/Pleiades_E/closed_loop/'
    os.chdir(data_dir)
    
    # img_files = ['obj_{0:03d}.fits'.format(ii) for ii in range(45, 55)]
    img_files = ['obj_050.fits']

    reduce_fli.clean_images(img_files, rebin=10)

    return
    
def find_stars_pleiades_unbinned_closed():
    data_dir = '/Users/jlu/data/imaka/2016_11_18/Pleiades_E/closed_loop/'
    os.chdir(data_dir)
    
    # img_files = ['obj_{0:03d}.fits'.format(ii) for ii in range(45, 55)]
    img_files = ['obj_050_bin_nobkg.fits']

    reduce_fli.find_stars_bin10(img_files, fwhm=5)

    return
    
def reduce_pleiades_unbinned_open():
    data_dir = '/Users/jlu/data/imaka/2016_11_18/Pleiades_E/open_loop/'
    os.chdir(data_dir)
    
    # img_files = ['obj_{0:03d}.fits'.format(ii) for ii in range(55, 57)]
    img_files = ['obj_056.fits']

    reduce_fli.clean_images(img_files, rebin=10)

    return
    
def find_stars_pleiades_unbinned_open():
    data_dir = '/Users/jlu/data/imaka/2016_11_18/Pleiades_E/open_loop/'
    os.chdir(data_dir)
    
    # img_files = ['obj_{0:03d}.fits'.format(ii) for ii in range(55, 57)]
    img_files = ['obj_056_bin_nobkg.fits']

    reduce_fli.find_stars_bin10(img_files, fwhm=5)

    return

def reduce_pleiades_binned_open():
    data_dir = '/Users/jlu/data/imaka/2016_11_18/Pleiades_E/open_loop/'
    os.chdir(data_dir)
    
    # img_files = ['obj_{0:03d}.fits'.format(ii) for ii in range(82, 91)]
    img_files = ['obj_{0:03d}.fits'.format(ii) for ii in range(163, 173)]

    reduce_fli.clean_images(img_files, rebin=2)

    return
    

def reduce_pleiades_binned_closed():
    data_dir = '/Users/jlu/data/imaka/2016_11_18/Pleiades_E/closed_loop/'
    os.chdir(data_dir)
    
    # img_files = ['obj_{0:03d}.fits'.format(ii) for ii in range(92, 101)]
    img_files = ['obj_{0:03d}.fits'.format(ii) for ii in range(153, 163)]

    reduce_fli.clean_images(img_files, rebin=2)

    return
    
def find_stars_pleiades_binned_open():
    data_dir = '/Users/jlu/data/imaka/2016_11_18/Pleiades_E/open_loop/'
    os.chdir(data_dir)
    
    # img_files = ['obj_{0:03d}_bin_nobkg.fits'.format(ii) for ii in range(82, 92)]
    # img_files = ['obj_{0:03d}_bin_nobkg.fits'.format(ii) for ii in range(163, 173)]
    img_files = ['obj_{0:03d}_bin_nobkg.fits'.format(ii) for ii in range(167, 173)]

    reduce_fli.find_stars_bin(img_files, fwhm=4, threshold=6)

    return
    
def find_stars_pleiades_binned_closed():
    data_dir = '/Users/jlu/data/imaka/2016_11_18/Pleiades_E/closed_loop/'
    os.chdir(data_dir)
    
    # img_files = ['obj_{0:03d}_bin_nobkg.fits'.format(ii) for ii in range(92, 102)]
    img_files = ['obj_{0:03d}_bin_nobkg.fits'.format(ii) for ii in range(153, 163)]

    reduce_fli.find_stars_bin(img_files, fwhm=4, threshold=6)

    return

def compare_fwhm_list():
    data_dir = '/Users/jlu/data/imaka/2016_11_18/Pleiades_E/'
    os.chdir(data_dir)
    
    # o_list = np.arange(163, 173) # Open loop star lists
    # c_list = np.arange(153, 163) # Closed loop star lists
    o_list = np.arange(163, 173) # Open loop star lists
    c_list = np.arange(153, 163) # Closed loop star lists

    plt.ion()

    for ii in range(len(o_list)):
        open_list = 'open_loop/obj_{0:03d}_bin_nobkg_stars.txt'.format(o_list[ii])
        closed_list = 'closed_loop/obj_{0:03d}_bin_nobkg_stars.txt'.format(c_list[ii])

        compare_fwhm(open_list, closed_list, scale=6*0.04)
        pdb.set_trace()

def compare_fwhm(open_list, closed_list, scale=0.04):
    topen = table.Table.read(open_list, format='ascii')
    tclose = table.Table.read(closed_list, format='ascii')

    # Trim out any stars with FWHM > 10
    idx_o = np.where(topen['x_fwhm'] < 10)[0]
    idx_c = np.where(tclose['x_fwhm'] < 10)[0]

    topen = topen[idx_o]
    tclose = tclose[idx_c]

    m_c = np.ones(len(tclose))
    m_o = np.ones(len(topen))

    idx_c, idx_o, dr, dm = match.match(tclose['xcentroid'], tclose['ycentroid'], m_c,
                                       topen['xcentroid'], topen['ycentroid'], m_o,
                                       dr_tol=10)
    # Matched catalogs
    to_match = topen[idx_o]
    tc_match = tclose[idx_c]

    # Plot
    plt.figure(1)
    plt.clf()
    plt.plot(to_match['x_fwhm']*scale, tc_match['x_fwhm']*scale, 'r.', label='X')
    plt.plot(to_match['y_fwhm']*scale, tc_match['y_fwhm']*scale, 'b.', label='Y')
    plt.plot([0, 10], [0, 10], 'k--')
    plt.xlabel('FWHM in Open Loop (")')
    plt.ylabel('FWHM in Closed Loop (")')
    plt.axis('equal')

    max_fwhm = 1.5*np.median([to_match['x_fwhm'], to_match['y_fwhm'], tc_match['x_fwhm'], tc_match['y_fwhm']])
    plt.ylim(0, max_fwhm*scale)
    plt.xlim(0, max_fwhm*scale)
    plt.legend(numpoints=1, loc='upper left')
    plt.pause(0.05)

    plt.figure(2)
    plt.clf()
    x_ratio = tc_match['x_fwhm'] / to_match['x_fwhm']
    y_ratio = tc_match['y_fwhm'] / to_match['y_fwhm']
    
    plt.plot(to_match['x_fwhm']*scale, x_ratio, 'r.', label='X')
    plt.plot(to_match['y_fwhm']*scale, y_ratio, 'b.', label='Y')
    plt.axhline(1, linestyle='--', color='black')
    plt.xlabel('FWHM in Open Loop (")')
    plt.ylabel('Closed / Open FWHM')

    max_fwhm = 1.5*np.median([to_match['x_fwhm'], to_match['y_fwhm'], tc_match['x_fwhm'], tc_match['y_fwhm']])
    plt.xlim(0, max_fwhm*scale)
    plt.ylim(0, 1.5)
    plt.legend(numpoints=1, loc='upper left')
    plt.pause(0.05)

    x_fwhm_o_med = np.median(to_match['x_fwhm'])    
    y_fwhm_o_med = np.median(to_match['y_fwhm'])    
    x_fwhm_c_med = np.median(tc_match['x_fwhm'])    
    y_fwhm_c_med = np.median(tc_match['y_fwhm'])
    
    print('Open Loop Stats:')
    print('\t Median x_fwhm = {0:.1f} +/- {1:.1f}'.format(x_fwhm_o_med,
                                                          to_match['x_fwhm'].std()))
    print('\t Median y_fwhm = {0:.1f} +/- {1:.1f}'.format(y_fwhm_o_med,
                                                          to_match['y_fwhm'].std()))

    print('Closed Loop Stats:')
    print('\t Median x_fwhm = {0:.1f} +/- {1:.1f}'.format(x_fwhm_c_med,
                                                          tc_match['x_fwhm'].std()))
    print('\t Median y_fwhm = {0:.1f} +/- {1:.1f}'.format(y_fwhm_c_med,
                                                          tc_match['y_fwhm'].std()))

    print('Fractional Improvement')
    print('\t x_fwhm closed / open = {0:.2f}'.format(x_fwhm_c_med / x_fwhm_o_med))
    print('\t y_fwhm closed / open = {0:.2f}'.format(y_fwhm_c_med / y_fwhm_o_med))
    
    
    
