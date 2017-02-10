import pylab as plt
import numpy as np
from astropy.io import fits
from astropy import table
from astropy import units
import scipy
import glob
from imaka.reduce import reduce_fli
from imaka.reduce import calib
import pdb
import os
from flystar import match

def make_sky():
    sky_dir = '/Users/jlu/data/imaka/2017_01_13/fli/Pleiades/'

    sky_num = [46, 47, 48]
    sky_frames = ['{0:s}sky{1:03d}.fits'.format(sky_dir, ss) for ss in sky_num]
    calib.makedark(sky_frames, 'pleiades_sky.fits')

    return
    
        
def reduce_pleiades_binned():
    sky_dir = '/Users/jlu/data/imaka/2017_01_13/fli/Pleiades/'
    data_dir = '/Users/jlu/data/imaka/2017_01_13/fli/Pleiades/'
    os.chdir(data_dir)

    # Open Loop
    # fnum = [34, 35, 36, 43, 44, 45]
    fnum = [56, 57, 58, 59]
    img_files = ['obj_o{0:03d}.fits'.format(ii) for ii in fnum]
    reduce_fli.clean_images(img_files, rebin=1, sky_frame=sky_dir + 'pleiades_sky.fits')

    # TTF Closed Loop
    # fnum = [37, 38, 39]
    fnum = [49, 50, 51, 60, 61, 62, 63]
    img_files = ['obj_ttf{0:03d}.fits'.format(ii) for ii in fnum]
    reduce_fli.clean_images(img_files, rebin=1, sky_frame=sky_dir + 'pleiades_sky.fits')

    # Closed Loop
    # fnum = [31, 32, 33, 40, 41, 42]
    fnum = [52, 53, 54, 55]
    img_files = ['obj_c{0:03d}.fits'.format(ii) for ii in fnum]
    reduce_fli.clean_images(img_files, rebin=1, sky_frame=sky_dir + 'pleiades_sky.fits')
    
    return
    
    
def find_stars_pleiades_binned_open():
    data_dir = '/Users/jlu/data/imaka/2017_01_13/fli/Pleiades/'
    os.chdir(data_dir)
    
    # fnum = [34, 35, 36, 43, 44, 45]
    fnum = [56, 57, 58, 59]
    img_files = ['obj_o{0:03d}_bin_nobkg.fits'.format(ii) for ii in fnum]
    
    reduce_fli.find_stars_bin(img_files, fwhm=5, threshold=6)

    return

def find_stars_pleiades_binned_ttf():
    data_dir = '/Users/jlu/data/imaka/2017_01_13/fli/Pleiades/'
    os.chdir(data_dir)
    
    # fnum = [37, 38, 39]
    fnum = [49, 50, 51, 60, 61, 62, 63]
    img_files = ['obj_ttf{0:03d}_bin_nobkg.fits'.format(ii) for ii in fnum]
    
    reduce_fli.find_stars_bin(img_files, fwhm=5, threshold=6)

    return

def find_stars_pleiades_binned_closed():
    data_dir = '/Users/jlu/data/imaka/2017_01_13/fli/Pleiades/'
    os.chdir(data_dir)
    
    # fnum = [31, 32, 33, 40, 41, 42]
    fnum = [52, 53, 54, 55]
    img_files = ['obj_c{0:03d}_bin_nobkg.fits'.format(ii) for ii in fnum]
    
    reduce_fli.find_stars_bin(img_files, fwhm=3, threshold=6)

    return

def calc_star_stats():
    fnum = [34, 35, 36, 43, 44, 45, 56, 57, 58, 59]
    img_files = ['obj_o{0:03d}_bin_nobkg.fits'.format(ii) for ii in fnum]
    reduce_fli.calc_star_stats(img_files, output_stats='stats_open.fits')

    fnum = [37, 38, 39, 49, 50, 51, 60, 61, 62, 63]
    img_files = ['obj_ttf{0:03d}_bin_nobkg.fits'.format(ii) for ii in fnum]
    reduce_fli.calc_star_stats(img_files, output_stats='stats_ttf.fits')

    fnum = [31, 32, 33, 40, 41, 42, 52, 53, 54, 55]
    img_files = ['obj_c{0:03d}_bin_nobkg.fits'.format(ii) for ii in fnum]
    reduce_fli.calc_star_stats(img_files, output_stats='stats_closed.fits')

    plot_stats()
    
    return
    

def compare_fwhm_list():
    data_dir = '/Users/jlu/data/imaka/2017_01_13/fli/Pleiades/'
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
    

def plot_stats():
    so = table.Table.read('stats_open.fits')
    st = table.Table.read('stats_ttf.fits')
    sc = table.Table.read('stats_closed.fits')

    scale = 0.12

    # FWHM
    plt.clf()
    plt.plot(so['Index'], so['FWHM']*scale, 'bo', label='Open')
    plt.plot(st['Index'], st['FWHM']*scale, 'go', label='TTF')
    plt.plot(sc['Index'], sc['FWHM']*scale, 'ro', label='Closed')
    plt.xlabel('Frame Number')
    plt.ylabel('FWHM (")')
    plt.legend(numpoints=1)
    plt.ylim(0, 1.5)
    plt.savefig('plots/fwhm_vs_frame.png')

    # EE 50
    plt.clf()
    plt.plot(so['Index'], so['EE50'], 'bo', label='Open')
    plt.plot(st['Index'], st['EE50'], 'go', label='TTF')
    plt.plot(sc['Index'], sc['EE50'], 'ro', label='Closed')
    plt.xlabel('Frame Number')
    plt.ylabel('Radius of 50% Encircled Energy (")')
    plt.legend(numpoints=1)
    plt.ylim(0, 1.5)
    plt.savefig('plots/ee50_vs_frame.png')

    # EE 80
    plt.clf()
    plt.plot(so['Index'], so['EE80'], 'bo', label='Open')
    plt.plot(st['Index'], st['EE80'], 'go', label='TTF')
    plt.plot(sc['Index'], sc['EE80'], 'ro', label='Closed')
    plt.xlabel('Frame Number')
    plt.ylabel('Radius of 80% Encircled Energy (")')
    plt.legend(numpoints=1)
    plt.ylim(0, 1.5)
    plt.savefig('plots/ee80_vs_frame.png')

    # NEA
    plt.clf()
    plt.plot(so['Index'], so['NEA'], 'bo', label='Open')
    plt.plot(st['Index'], st['NEA'], 'go', label='TTF')
    plt.plot(sc['Index'], sc['NEA'], 'ro', label='Closed')
    plt.xlabel('Frame Number')
    plt.ylabel('Noise Equivalent Area (Sq. Arcsec)')
    plt.legend(numpoints=1)
    plt.ylim(0, 3)
    plt.savefig('plots/nea_vs_frame.png')

    # NEA2
    plt.clf()
    plt.plot(so['Index'], so['NEA2'], 'bo', label='Open')
    plt.plot(st['Index'], st['NEA2'], 'go', label='TTF')
    plt.plot(sc['Index'], sc['NEA2'], 'ro', label='Closed')
    plt.xlabel('Frame Number')
    plt.ylabel('Noise Equivalent Area (Sq. Arcsec)')
    plt.legend(numpoints=1)
    plt.ylim(0, 3)
    plt.savefig('plots/nea2_vs_frame.png')
    
    # FWHM for each direction
    plt.clf()
    plt.plot(so['Index'], so['xFWHM']*scale, 'bo', label='Open X')
    plt.plot(st['Index'], st['xFWHM']*scale, 'go', label='TTF X')
    plt.plot(sc['Index'], sc['xFWHM']*scale, 'ro', label='Closed X')
    plt.plot(so['Index'], so['yFWHM']*scale, 'b^', label='Open Y')
    plt.plot(st['Index'], st['yFWHM']*scale, 'g^', label='TTF Y')
    plt.plot(sc['Index'], sc['yFWHM']*scale, 'r^', label='Closed Y')
    plt.xlabel('Frame Number')
    plt.ylabel('FWHM (")')
    plt.legend(numpoints=1, fontsize=10)
    plt.ylim(0, 1.5)
    plt.savefig('plots/xyfwhm_vs_frame.png')

    
    # 
    # Plot ratio of improvements. First we need to find a matching
    # closed loop image to go with each open (and TTF) image.
    #
    tree_sc = scipy.spatial.KDTree(np.array([sc['Index']]).T)
    dist_oc, idx_in_c_of_o = tree_sc.query(np.array([so['Index']]).T, 1)
    dist_tc, idx_in_c_of_t = tree_sc.query(np.array([st['Index']]).T, 1)

    # FWHM
    plt.clf()
    plt.plot(so['Index'], sc['FWHM'][idx_in_c_of_o] / so['FWHM'], 'bo', label='Closed / Open')
    plt.plot(st['Index'], sc['FWHM'][idx_in_c_of_t] / st['FWHM'], 'ro', label='Closed / TTF')
    plt.xlabel('Frame Number')
    plt.ylabel('Ratio of FWHM')
    plt.legend(numpoints=1)
    plt.axhline(1.0, color='k', linestyle='--', linewidth=2)
    plt.ylim(0, 1.3)
    plt.savefig('plots/fwhm_ratio_vs_frame.png')

    # EE 50
    plt.clf()
    plt.plot(so['Index'], sc['EE50'][idx_in_c_of_o] / so['EE50'], 'bo', label='Closed / Open')
    plt.plot(st['Index'], sc['EE50'][idx_in_c_of_t] / st['EE50'], 'ro', label='Closed / TTF')
    plt.xlabel('Frame Number')
    plt.ylabel('Ratio of 50% EE Radius')
    plt.legend(numpoints=1)
    plt.axhline(1.0, color='k', linestyle='--', linewidth=2)
    plt.ylim(0, 1.3)
    plt.savefig('plots/ee50_ratio_vs_frame.png')

    # EE 80
    plt.clf()
    plt.plot(so['Index'], sc['EE80'][idx_in_c_of_o] / so['EE80'], 'bo', label='Closed / Open')
    plt.plot(st['Index'], sc['EE80'][idx_in_c_of_t] / st['EE80'], 'ro', label='Closed / TTF')
    plt.xlabel('Frame Number')
    plt.ylabel('Ratio of 80% EE Radius')
    plt.legend(numpoints=1)
    plt.axhline(1.0, color='k', linestyle='--', linewidth=2)
    plt.ylim(0, 1.3)
    plt.savefig('plots/ee80_ratio_vs_frame.png')

    # NEA
    plt.clf()
    plt.plot(so['Index'], sc['NEA'][idx_in_c_of_o] / so['NEA'], 'bo', label='Closed / Open')
    plt.plot(st['Index'], sc['NEA'][idx_in_c_of_t] / st['NEA'], 'ro', label='Closed / TTF')
    plt.xlabel('Frame Number')
    plt.ylabel('Ratio of NEA')
    plt.legend(numpoints=1)
    plt.axhline(1.0, color='k', linestyle='--', linewidth=2)
    plt.ylim(0, 1.3)
    plt.savefig('plots/nea_ratio_vs_frame.png')

    # NEA2
    plt.clf()
    plt.plot(so['Index'], sc['NEA2'][idx_in_c_of_o] / so['NEA2'], 'bo', label='Closed / Open')
    plt.plot(st['Index'], sc['NEA2'][idx_in_c_of_t] / st['NEA2'], 'ro', label='Closed / TTF')
    plt.xlabel('Frame Number')
    plt.ylabel('Ratio of NEA2')
    plt.legend(numpoints=1)
    plt.axhline(1.0, color='k', linestyle='--', linewidth=2)
    plt.ylim(0, 1.3)
    plt.savefig('plots/nea2_ratio_vs_frame.png')
    
    # FWHM for each direction
    plt.clf()
    plt.plot(so['Index'], sc['xFWHM'][idx_in_c_of_o] / so['xFWHM'], 'bo', label='X Closed / Open')
    plt.plot(st['Index'], sc['xFWHM'][idx_in_c_of_t] / st['xFWHM'], 'ro', label='X Closed / TTF')
    plt.plot(so['Index'], sc['yFWHM'][idx_in_c_of_o] / so['yFWHM'], 'b^', label='Y Closed / Open')
    plt.plot(st['Index'], sc['yFWHM'][idx_in_c_of_t] / st['yFWHM'], 'r^', label='Y Closed / TTF')
    plt.xlabel('Frame Number')
    plt.ylabel('Ratio of FWHM')
    plt.legend(numpoints=1, fontsize=10)
    plt.axhline(1.0, color='k', linestyle='--', linewidth=2)
    plt.ylim(0, 1.3)
    plt.savefig('plots/xyfwhm_ratio_vs_frame.png')
    
    return
