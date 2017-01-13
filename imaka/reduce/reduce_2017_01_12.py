import pylab as plt
import numpy as np
from astropy.io import fits
from astropy import table
from astropy import units
import glob
from imaka.reduce import reduce_fli
from imaka.reduce import calib
import pdb
import os
from flystar import match

def make_sky():
    sky_dir = '/Users/jlu/data/imaka/2017_01_12/fli/Pleiades/'

    sky_num = [41, 42, 43, 71, 72, 73]
    sky_frames = ['{0:s}sky{1:03d}.fits'.format(sky_dir, ss) for ss in sky_num]
    calib.makedark(sky_frames, 'pleiades_sky.fits')

    return
    
        
def reduce_pleiades_binned():
    sky_dir = '/Users/jlu/data/imaka/2017_01_12/fli/Pleiades/'
    data_dir = '/Users/jlu/data/imaka/2017_01_12/fli/Pleiades/'
    os.chdir(data_dir)

    # # Open Loop
    # fnum = [56, 57, 58, 65, 66, 67]
    # img_files = ['obj_o{0:03d}.fits'.format(ii) for ii in fnum]
    # reduce_fli.clean_images(img_files, rebin=1, sky_frame=sky_dir + 'pleiades_sky.fits')

    # # TTF Closed Loop
    # fnum = [32, 33, 34, 38, 39, 40, 50, 51, 52, 59, 60, 61]
    # img_files = ['obj_ttf{0:03d}.fits'.format(ii) for ii in fnum]
    # reduce_fli.clean_images(img_files, rebin=1, sky_frame=sky_dir + 'pleiades_sky.fits')

    # Closed Loop
    fnum = [27, 28, 29, 30, 31, 35, 36, 37, 44, 45, 46, 47, 48, 49, 53, 54, 55, 62, 63, 64]
    img_files = ['obj_c{0:03d}.fits'.format(ii) for ii in fnum]
    reduce_fli.clean_images(img_files, rebin=1, sky_frame=sky_dir + 'pleiades_sky.fits')
    
    return
    
    
def find_stars_pleiades_binned_open():
    data_dir = '/Users/jlu/data/imaka/2017_01_12/fli/Pleiades/'
    os.chdir(data_dir)
    
    fnum = [56, 57, 58, 65, 66, 67]
    img_files = ['obj_o{0:03d}_bin_nobkg.fits'.format(ii) for ii in fnum]
    reduce_fli.find_stars_bin(img_files, fwhm=5, threshold=6)

    return

def find_stars_pleiades_binned_ttf():
    data_dir = '/Users/jlu/data/imaka/2017_01_12/fli/Pleiades/'
    os.chdir(data_dir)
    
    fnum = [32, 33, 34, 38, 39, 40, 50, 51, 52, 59, 60, 61]
    img_files = ['obj_ttf{0:03d}_bin_nobkg.fits'.format(ii) for ii in fnum]
    reduce_fli.find_stars_bin(img_files, fwhm=5, threshold=6)

    return

def find_stars_pleiades_binned_closed():
    data_dir = '/Users/jlu/data/imaka/2017_01_12/fli/Pleiades/'
    os.chdir(data_dir)
    
    fnum = [27, 28, 29, 30, 31, 35, 36, 37, 44, 45, 46, 47, 48, 49, 53, 54, 55, 62, 63, 64]
    img_files = ['obj_c{0:03d}_bin_nobkg.fits'.format(ii) for ii in fnum]
    reduce_fli.find_stars_bin(img_files, fwhm=3, threshold=6)

    return

def calc_star_stats_open():
    data_dir = '/Users/jlu/data/imaka/2017_01_10/fli/Pleiades/'
    os.chdir(data_dir)
    
    fnum = [10, 11]
    img_files = ['obj{0:03d}_bin_nobkg.fits'.format(ii) for ii in fnum]
    reduce_fli.calc_star_stats(img_files)

    return


def compare_fwhm_list():
    data_dir = '/Users/jlu/data/imaka/2017_01_10/fli/Pleiades/'
    os.chdir(data_dir)
    
    # o_list = np.arange(163, 173) # Open loop star lists
    # c_list = np.arange(153, 163) # Closed loop star lists
    # o_list = np.arange(163, 173) # Open loop star lists
    # c_list = np.arange(153, 163) # Closed loop star lists
    o_list = [10, 11, 14, 15] # open
    c_list = [8, 9, 12, 13]   # closed
    
    plt.ion()

    for ii in range(len(o_list)):
        open_list = 'obj{0:03d}_bin_nobkg_stars.txt'.format(o_list[ii])
        closed_list = 'obj{0:03d}_bin_nobkg_stars.txt'.format(c_list[ii])

        compare_fwhm(open_list, closed_list, scale=3*0.04, flux_min=1)
        pdb.set_trace()

def compare_fwhm(open_list, closed_list, scale=0.04, flux_min=2.0):
    topen = table.Table.read(open_list, format='ascii')
    tclose = table.Table.read(closed_list, format='ascii')

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
    plt.plot(to_match['x_fwhm']*scale, tc_match['x_fwhm']*scale, 'r.', label='X')
    plt.plot(to_match['y_fwhm']*scale, tc_match['y_fwhm']*scale, 'b.', label='Y')
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
    
    plt.plot(to_match['x_fwhm']*scale, x_ratio, 'r.', label='X')
    plt.plot(to_match['y_fwhm']*scale, y_ratio, 'b.', label='Y')
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
    
    print('Open Loop Stats:')
    print('\t Median x_fwhm = {0:.1f} +/- {1:.1f}'.format(x_fwhm_o_med * scale,
                                                          to_match['x_fwhm'].std() * scale))
    print('\t Median y_fwhm = {0:.1f} +/- {1:.1f}'.format(y_fwhm_o_med * scale,
                                                          to_match['y_fwhm'].std() * scale))

    print('Closed Loop Stats:')
    print('\t Median x_fwhm = {0:.1f} +/- {1:.1f}'.format(x_fwhm_c_med * scale,
                                                          tc_match['x_fwhm'].std() * scale))
    print('\t Median y_fwhm = {0:.1f} +/- {1:.1f}'.format(y_fwhm_c_med * scale,
                                                          tc_match['y_fwhm'].std() * scale))

    print('Fractional Improvement')
    print('\t x_fwhm closed / open = {0:.2f}'.format(x_fwhm_c_med / x_fwhm_o_med))
    print('\t y_fwhm closed / open = {0:.2f}'.format(y_fwhm_c_med / y_fwhm_o_med))
    
