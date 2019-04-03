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
    sky_dir = '/Users/jlu/data/imaka/2016_11_19/20161118/skies/binned/'
    sky_frames = [133, 145, 146, 163, 172, 179, 188]

    for ss in range(len(sky_frames)):
        sky_frames[ss] = '{0:s}sky_{1:03d}.fits'.format(sky_dir, sky_frames[ss])

    calib.makedark(sky_frames, 'pleiades_sky.fits')
    
def make_twilight_flats():
    twilights = ['twi_1001.fits', 'twi_2002.fits', 'twi_3003.fits',
                     'twi_4004.fits', 'twi_5005.fits', 'twi_6006.fits', 'twi_7007.fits',
                     'twi_8008.fits', 'twi_9009.fits', 'twi_10010.fits', 'twi_11011.fits',
                     'twi_12012.fits', 'twi_13013.fits', 'twi_14014.fits', 'twi_15015.fits',
                     'twi_16016.fits', 'twi_17017.fits', 'twi_18018.fits', 'twi_19019.fits',
                     'twi_20020.fits', 'twi_21021.fits', 'twi_22022.fits', 'twi_23023.fits',
                     'twi_24024.fits', 'twi_25025.fits', 'twi_26026.fits', 'twi_27027.fits',
                     'twi_28028.fits', 'twi_29029.fits', 'twi_30030.fits', 'twi_31031.fits']

    twi_darks = ['dark_1001002.fits', 'dark_2002013.fits', 'dark_3003024.fits',
                     'dark_4004035.fits', 'dark_5005044.fits', 'dark_6006045.fits', 'dark_7007046.fits',
                     'dark_8008047.fits',  'dark_9009048.fits', 'dark_10010001.fits', 'dark_11011003.fits',
                     'dark_12012004.fits', 'dark_13013005.fits', 'dark_14014006.fits', 'dark_15015007.fits',
                     'dark_16016008.fits', 'dark_17017009.fits', 'dark_18018010.fits', 'dark_19019011.fits',
                     'dark_20020012.fits', 'dark_21021014.fits', 'dark_22022015.fits', 'dark_23023016.fits',
                     'dark_24024017.fits', 'dark_25025018.fits', 'dark_26026019.fits', 'dark_27027020.fits',
                     'dark_28028021.fits', 'dark_29029022.fits', 'dark_30030023.fits', 'dark_31031025.fits']

    twilight_root = '/Users/jlu/data/imaka/2016_11_19/20161118/twilight/'
    twi_dark_root = '/Users/jlu/data/imaka/2016_11_19/20161118/darks/darks_for_twis/'

    for tt in range(len(twilights)):
        twilights[tt] = twilight_root + twilights[tt]
        twi_darks[tt] = twi_dark_root + twi_darks[tt]
        
    calib.makeflat(twilights, twi_darks, twilight_root + 'flat_r.fits')
    
    return

        
def reduce_pleiades_binned_open():
    sky_dir = '/Users/jlu/data/imaka/2016_11_19/20161118/skies/binned/'
    data_dir = '/Users/jlu/data/imaka/2016_11_19/20161118/Pleiades_E/open_loop/'
    os.chdir(data_dir)

    fnum = [137, 139, 142, 144, 148, 150, 152]
    img_files = ['obj_{0:03d}.fits'.format(ii) for ii in fnum]

    reduce_fli.clean_images(img_files, rebin=1, sky_frame=sky_dir + 'pleiades_sky.fits')

    return
    

def reduce_pleiades_binned_closed():
    sky_dir = '/Users/jlu/data/imaka/2016_11_19/20161118/skies/binned/'
    data_dir = '/Users/jlu/data/imaka/2016_11_19/20161118/Pleiades_E/closed_loop/'
    os.chdir(data_dir)

    fnum = [138, 141, 143, 147, 149, 151]
    img_files = ['obj_{0:03d}.fits'.format(ii) for ii in fnum]

    reduce_fli.clean_images(img_files, rebin=1, sky_frame=sky_dir + 'pleiades_sky.fits')

    return
    
def find_stars_pleiades_binned_open():
    data_dir = '/Users/jlu/data/imaka/2016_11_19/20161118/Pleiades_E/open_loop/'
    os.chdir(data_dir)
    
    fnum = [137, 139, 142, 144, 148, 150, 152]
    img_files = ['obj_{0:03d}_bin_nobkg.fits'.format(ii) for ii in fnum]

    reduce_fli.find_stars_bin(img_files, fwhm=2, threshold=6)

    return
    
def find_stars_pleiades_binned_closed():
    data_dir = '/Users/jlu/data/imaka/2016_11_19/20161118/Pleiades_E/closed_loop/'
    os.chdir(data_dir)
    
    fnum = [138, 141, 143, 147, 149, 151]
    img_files = ['obj_{0:03d}_bin_nobkg.fits'.format(ii) for ii in fnum]

    reduce_fli.find_stars_bin(img_files, fwhm=2, threshold=6)

    return


def compare_fwhm_list():
    data_dir = '/Users/jlu/data/imaka/2016_11_19/20161118/Pleiades_E/'
    os.chdir(data_dir)
    
    # o_list = np.arange(163, 173) # Open loop star lists
    # c_list = np.arange(153, 163) # Closed loop star lists
    # o_list = np.arange(163, 173) # Open loop star lists
    # c_list = np.arange(153, 163) # Closed loop star lists
    o_list = [137, 139, 142, 144, 148, 150]
    c_list = [138, 141, 143, 147, 149, 151]
    
    plt.ion()

    for ii in range(len(o_list)):
        open_list = 'open_loop/obj_{0:03d}_bin_nobkg_stars.txt'.format(o_list[ii])
        closed_list = 'closed_loop/obj_{0:03d}_bin_nobkg_stars.txt'.format(c_list[ii])

        compare_fwhm(open_list, closed_list, scale=3*0.04, flux_min=5)
        pdb.set_trace()

def compare_fwhm(open_list, closed_list, scale=0.04, flux_min=2.0):
    topen = table.Table.read(open_list, format='ascii')
    tclose = table.Table.read(closed_list, format='ascii')

    # Trim out any stars with FWHM > 10
    idx_o = np.where(topen['x_fwhm'] < 10)[0]
    idx_c = np.where(tclose['x_fwhm'] < 10)[0]

    topen = topen[idx_o]
    tclose = tclose[idx_c]

    # Trim out any stars with flux < flux_min
    idx_o = np.where(topen['flux'] > flux_min)[0]
    idx_c = np.where(tclose['flux'] > flux_min)[0]

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
    
