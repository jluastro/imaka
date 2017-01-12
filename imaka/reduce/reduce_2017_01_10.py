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
    sky_dir = '/Users/jlu/data/imaka/2017_01_10/fli/Pleiades/'

    sky_num = [20, 21, 30, 31, 46, 47, 56, 57, 70, 71, 84, 85, 98, 99, 112, 113,
                   126, 127, 140, 141, 140, 141, 154, 155, 168, 169, 186, 187]
    sky_frames = ['{0:s}sky{1:03d}.fits'.format(sky_dir, ss) for ss in sky_num]
    calib.makedark(sky_frames, 'pleiades_sky.fits')
    
    sky_num = [20, 21, 30, 31]
    sky_frames = ['{0:s}sky{1:03d}.fits'.format(sky_dir, ss) for ss in sky_num]
    calib.makedark(sky_frames, 'pleiades_sky_020.fits')

    sky_num = [46, 47, 56, 57]
    sky_frames = ['{0:s}sky{1:03d}.fits'.format(sky_dir, ss) for ss in sky_num]
    calib.makedark(sky_frames, 'pleiades_sky_046.fits')

    sky_num = [70, 71, 84, 85]
    sky_frames = ['{0:s}sky{1:03d}.fits'.format(sky_dir, ss) for ss in sky_num]
    calib.makedark(sky_frames, 'pleiades_sky_070.fits')

    sky_num = [98, 99, 112, 113]
    sky_frames = ['{0:s}sky{1:03d}.fits'.format(sky_dir, ss) for ss in sky_num]
    calib.makedark(sky_frames, 'pleiades_sky_098.fits')

    sky_num = [126, 127, 140, 141]
    sky_frames = ['{0:s}sky{1:03d}.fits'.format(sky_dir, ss) for ss in sky_num]
    calib.makedark(sky_frames, 'pleiades_sky_126.fits')

    sky_num = [140, 141, 154, 155]
    sky_frames = ['{0:s}sky{1:03d}.fits'.format(sky_dir, ss) for ss in sky_num]
    calib.makedark(sky_frames, 'pleiades_sky_140.fits')

    sky_num = [168, 169, 186, 187]
    sky_frames = ['{0:s}sky{1:03d}.fits'.format(sky_dir, ss) for ss in sky_num]
    calib.makedark(sky_frames, 'pleiades_sky_168.fits')

    # Note this is the 45 sec integration (the rest are the 30 sec integration)
    sky_num = [204, 205, 222, 240, 241]
    sky_frames = ['{0:s}sky{1:03d}.fits'.format(sky_dir, ss) for ss in sky_num]
    calib.makedark(sky_frames, 'pleiades_sky_204.fits')

    return
    
        
def reduce_pleiades_binned_open():
    sky_dir = '/Users/jlu/data/imaka/2017_01_10/fli/Pleiades/'
    data_dir = '/Users/jlu/data/imaka/2017_01_10/fli/Pleiades/'
    os.chdir(data_dir)

    fnum = [10, 11, 14, 15, 18, 19, 24, 25, 28, 29, 34, 35]
    img_files = ['obj{0:03d}.fits'.format(ii) for ii in fnum]
    reduce_fli.clean_images(img_files, rebin=1, sky_frame=sky_dir + 'pleiades_sky_020.fits')

    fnum = [38, 39, 40, 41, 44, 45, 50, 51, 60, 61]
    img_files = ['obj{0:03d}.fits'.format(ii) for ii in fnum]
    reduce_fli.clean_images(img_files, rebin=1, sky_frame=sky_dir + 'pleiades_sky_046.fits')

    fnum = [64, 65, 68, 69, 74, 75, 78, 79, 82, 83, 88, 89]
    img_files = ['obj{0:03d}.fits'.format(ii) for ii in fnum]
    reduce_fli.clean_images(img_files, rebin=1, sky_frame=sky_dir + 'pleiades_sky_070.fits')

    fnum_1 = [92, 93, 96, 97, 102, 103]
    fnum_2 = [106, 107, 110, 111, 116, 117]
    img_files_1 = ['obj{0:03d}.fits'.format(ii) for ii in fnum_1]
    img_files_2 = ['obj_o{0:03d}.fits'.format(ii) for ii in fnum_2]
    img_files = img_files_1 + img_files_2
    reduce_fli.clean_images(img_files, rebin=1, sky_frame=sky_dir + 'pleiades_sky_098.fits')

    fnum = [120, 121, 124, 125, 130, 131]
    img_files = ['obj_o{0:03d}.fits'.format(ii) for ii in fnum]
    reduce_fli.clean_images(img_files, rebin=1, sky_frame=sky_dir + 'pleiades_sky_126.fits')

    fnum = [134, 135, 138, 139, 144, 145, 148, 149, 152, 153, 158, 159]
    img_files = ['obj_o{0:03d}.fits'.format(ii) for ii in fnum]
    reduce_fli.clean_images(img_files, rebin=1, sky_frame=sky_dir + 'pleiades_sky_140.fits')

    fnum = [162, 163, 166, 167, 172, 173, 176, 177, 180, 181, 184, 185, 190, 191, 194, 195]
    img_files = ['obj_o{0:03d}.fits'.format(ii) for ii in fnum]
    reduce_fli.clean_images(img_files, rebin=1, sky_frame=sky_dir + 'pleiades_sky_168.fits')

    fnum = [198, 199, 202, 203, 208, 209, 212, 213, 216, 217, 220,
                221, 226, 227, 230, 231, 234, 235, 238, 239, 244, 245, 248, 249, 252, 253]
    img_files = ['obj_o{0:03d}.fits'.format(ii) for ii in fnum]
    reduce_fli.clean_images(img_files, rebin=1, sky_frame=sky_dir + 'pleiades_sky_204.fits')
    
    return
    

def reduce_pleiades_binned_closed():
    sky_dir = '/Users/jlu/data/imaka/2017_01_10/fli/Pleiades/'
    data_dir = '/Users/jlu/data/imaka/2017_01_10/fli/Pleiades/'
    os.chdir(data_dir)

    fnum = [8, 9, 12, 13, 16, 17, 22, 23, 26, 27, 32, 33, 36, 37]
    img_files = ['obj{0:03d}.fits'.format(ii) for ii in fnum]
    reduce_fli.clean_images(img_files, rebin=1, sky_frame=sky_dir + 'pleiades_sky_020.fits')

    fnum = [42, 43, 48, 49, 52, 53, 58, 59, 62, 63]
    img_files = ['obj{0:03d}.fits'.format(ii) for ii in fnum]
    reduce_fli.clean_images(img_files, rebin=1, sky_frame=sky_dir + 'pleiades_sky_046.fits')

    fnum = [66, 67, 72, 73, 76, 77, 80, 81, 86, 87, 90, 91]
    img_files = ['obj{0:03d}.fits'.format(ii) for ii in fnum]
    reduce_fli.clean_images(img_files, rebin=1, sky_frame=sky_dir + 'pleiades_sky_070.fits')

    fnum_1 = [94, 95, 100, 101, 104, 105]
    fnum_2 = [108, 109, 114, 115, 118, 119]
    img_files_1 = ['obj{0:03d}.fits'.format(ii) for ii in fnum_1]
    img_files_2 = ['obj_c{0:03d}.fits'.format(ii) for ii in fnum_2]
    img_files = img_files_1 + img_files_2
    reduce_fli.clean_images(img_files, rebin=1, sky_frame=sky_dir + 'pleiades_sky_098.fits')

    fnum = [122, 123, 128, 129, 132, 133, 136, 137]
    img_files = ['obj_c{0:03d}.fits'.format(ii) for ii in fnum]
    reduce_fli.clean_images(img_files, rebin=1, sky_frame=sky_dir + 'pleiades_sky_126.fits')

    fnum = [142, 143, 146, 147, 150, 151, 156, 157, 160, 161]
    img_files = ['obj_c{0:03d}.fits'.format(ii) for ii in fnum]
    reduce_fli.clean_images(img_files, rebin=1, sky_frame=sky_dir + 'pleiades_sky_140.fits')

    fnum = [164, 165, 170, 171, 174, 175, 178, 179, 182, 183, 188, 189, 192, 193,
                196, 197, 200, 201, 206, 207, 210, 211, 214, 215, 218, 219, 224, 225,
                228, 229, 232, 233, 236, 237, 242, 243, 246, 247, 250, 251]
    img_files = ['obj_c{0:03d}.fits'.format(ii) for ii in fnum]
    reduce_fli.clean_images(img_files, rebin=1, sky_frame=sky_dir + 'pleiades_sky_168.fits')

    return
    
def find_stars_pleiades_binned_open():
    data_dir = '/Users/jlu/data/imaka/2017_01_10/fli/Pleiades/'
    os.chdir(data_dir)
    
    fnum_1 = [10, 11, 14, 15, 18, 19, 24, 25, 28, 29, 34, 35]
    fnum_1 += [38, 39, 40, 41, 44, 45, 50, 51, 60, 61]
    fnum_1 += [64, 65, 68, 69, 74, 75, 78, 79, 82, 83, 88, 89]
    fnum_1 += [92, 93, 96, 97, 102, 103]
    img_files_1 = ['obj{0:03d}_bin_nobkg.fits'.format(ii) for ii in fnum_1]
    
    fnum_2 = [106, 107, 110, 111, 116, 117]
    fnum_2 = [120, 121, 124, 125, 130, 131]
    fnum_2 = [134, 135, 138, 139, 144, 145, 148, 149, 152, 153, 158, 159]
    fnum_2 = [162, 163, 166, 167, 172, 173, 176, 177, 180, 181, 184, 185, 190, 191, 194, 195]
    fnum_2 = [198, 199, 202, 203, 208, 209, 212, 213, 216, 217, 220,
                221, 226, 227, 230, 231, 234, 235, 238, 239, 244, 245, 248, 249, 252, 253]
    img_files_2 = ['obj_o{0:03d}_bin_nobkg.fits'.format(ii) for ii in fnum_2]

    img_files = img_files_1 + img_files_2
    
    reduce_fli.find_stars_bin(img_files, fwhm=2, threshold=6)

    return
    
def find_stars_pleiades_binned_closed():
    data_dir = '/Users/jlu/data/imaka/2017_01_10/fli/Pleiades/'
    os.chdir(data_dir)
    
    fnum = [8, 9, 12, 13, 16, 17, 22, 23, 26, 27, 32, 33, 36, 37]
    fnum = [42, 43, 48, 49, 52, 53, 58, 59, 62, 63]
    fnum = [66, 67, 72, 73, 76, 77, 80, 81, 86, 87, 90, 91]
    fnum_1 = [94, 95, 100, 101, 104, 105]
    img_files_1 = ['obj{0:03d}_bin_nobkg.fits'.format(ii) for ii in fnum_1]

    fnum_2 = [108, 109, 114, 115, 118, 119]
    fnum = [122, 123, 128, 129, 132, 133, 136, 137]
    fnum = [142, 143, 146, 147, 150, 151, 156, 157, 160, 161]
    fnum = [164, 165, 170, 171, 174, 175, 178, 179, 182, 183, 188, 189, 192, 193,
                196, 197, 200, 201, 206, 207, 210, 211, 214, 215, 218, 219, 224, 225,
                228, 229, 232, 233, 236, 237, 242, 243, 246, 247, 250, 251]
    img_files_2 = ['obj_c{0:03d}_bin_nobkg.fits'.format(ii) for ii in fnum_2]

    img_files = img_files_1 + img_files_2
    reduce_fli.find_stars_bin(img_files, fwhm=2, threshold=6)

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
    
