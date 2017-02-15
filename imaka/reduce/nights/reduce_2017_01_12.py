import pylab as plt
import numpy as np
from astropy.io import fits
from astropy import table
from astropy import units
import scipy
import glob
from imaka.reduce import reduce_fli
from imaka.reduce import calib
from imaka.reduce import util
import pdb
import os
#from flystar import match

root_dir = '/Volumes/g/lu/data/imaka/2017_01_12/fli/'

def make_sky():

    sky_dir = '/Volumes/g/lu/data/imaka/2017_01_12/fli/Pleiades/'
    out_dir = '/Volumes/g/lu/data/imaka/2017_01_12/fli/reduce/sky/'

    sky_num = [41, 42, 43, 71, 72, 73, 116, 117, 118, 119, 147, 148, 149, 177, 178, 179, 207, 208, 209, 243, 244, 245, 273, 274, 275, 303, 304, 305]
    sky_frames = ['{0:s}sky{1:03d}.fits'.format(sky_dir, ss) for ss in sky_num]
    calib.makedark(sky_frames, 'pleiades_sky.fits')

    return

def make_flat():
    # Didn't take flats, so just copy over the one from Friday night.
    old_flat = imaka_dir + '2017_01_13/fli/calib/flat.fits'
    new_flat = root_dir + 'calib/flat.fits'
    shutil.copyfile(old_flat, new_flat)
    return

def reduce_pleiades():
    sky_dir = root_dir + 'reduce/sky/'
    data_dir = root_dir + 'Pleiades/'
    flat_dir = root_dir + 'reduce/calib/'
    out_dir = root_dir + 'reduce/pleiades/'
    
    util.mkdir(out_dir)

    # Open Loop
    fnum1 = [56, 57, 58, 65, 66, 67, 77, 78, 79, 92, 93, 94, 98, 99, 100, 104, 105, 106]
    fnum2 = [113, 114, 115, 126, 127, 128, 135, 136, 137, 144, 145, 146, 156, 157, 158]
    fnum3 = [165, 166, 167, 174, 175, 176, 186, 187, 188, 195, 196, 197, 204, 205, 206]
    fnum4 = [216, 217, 218, 231, 232, 233, 240, 241, 242, 252, 253, 254, 261, 262, 263]
    fnum5 = [270, 271, 272, 280, 281, 286, 287, 292, 293, 299, 300, 301, 302, 312, 313]
    fnum6 = [314, 321, 322, 323]
    fnum = fnum1 + fnum2 + fnum3 + fnum4 + fnum5 + fnum6
    img_files = [data_dir + 'obj_o{0:03d}.fits'.format(ii) for ii in fnum]
    reduce_fli.clean_images(img_files, out_dir, rebin=1, sky_frame=sky_dir + 'pleiades_sky.fits', flat_frame=flat_dir + 'flat.fits')

    # TTF Closed Loop
    fnum1 = [32, 33, 34, 38, 39, 40, 50, 51, 52, 59, 60, 61, 68, 69, 70, 80, 81, 82, 86]
    fnum2 = [87, 88, 95, 96, 97, 107, 108, 109, 120, 121, 122, 129, 130, 131, 138, 139]
    fnum3 = [140, 150, 151, 152, 159, 160, 161, 168, 169, 170, 180, 181, 182, 189, 190]
    fnum4 = [191, 198, 199, 200, 210, 211, 212, 219, 220, 221, 225, 226, 227, 234, 235]
    fnum5 = [236, 246, 247, 248, 255, 256, 257, 264, 265, 266, 276, 277, 282, 283, 288]
    fnum6 = [289, 294, 295, 306, 307, 308, 315, 316, 317]
    fnum =  fnum1 + fnum2 + fnum3 + fnum4 + fnum5 + fnum6
    img_files = [data_dir + 'obj_ttf{0:03d}.fits'.format(ii) for ii in fnum]
    reduce_fli.clean_images(img_files, out_dir, rebin=1, sky_frame=sky_dir + 'pleiades_sky.fits', flat_frame=flat_dir + 'flat.fits')

    # Closed Loop
    fnum1 = [27, 28, 29, 30, 31, 35, 36, 37, 44, 45, 46, 47, 48, 49, 53, 54, 55, 62, 63, 64]
    fnum2 = [74, 75, 76, 89, 90, 91, 101, 102, 103, 110, 111, 112, 123, 124, 125, 132, 133]
    fnum3 = [134, 141, 142, 143, 153, 154, 155, 162, 163, 164, 171, 172, 173, 183, 184, 185]
    fnum4 = [192, 193, 194, 201, 202, 203, 213, 214, 215, 222, 223, 224, 228, 229, 230, 237]
    fnum5 = [238, 239, 249, 250, 251, 258, 259, 260, 267, 268, 269, 278, 279, 284, 285, 290]
    fnum6 = [291, 296, 297, 309, 310, 311, 318, 319, 320]
    fnum =  fnum1 + fnum2 + fnum3 + fnum4 + fnum5 + fnum6
    img_files = [data_dir + 'obj_c{0:03d}.fits'.format(ii) for ii in fnum]
    reduce_fli.clean_images(img_files, out_dir, rebin=1, sky_frame=sky_dir + 'pleiades_sky.fits', flat_frame = flat_dir + 'flat.fits')
    
    return
    
    
def find_stars_pleiades_open():
    reduce_dir = root_dir + 'reduce/pleiades/'
    
    fnum1 = [56, 57, 58, 65, 66, 67, 77, 78, 79, 92, 93, 94, 98, 99, 100, 104, 105, 106]
    fnum2 = [113, 114, 115, 126, 127, 128, 135, 136, 137, 144, 145, 146, 156, 157, 158]
    fnum3 = [165, 166, 167, 174, 175, 176, 186, 187, 188, 195, 196, 197, 204, 205, 206]
    fnum4 = [216, 217, 218, 231, 232, 233, 240, 241, 242, 252, 253, 254, 261, 262, 263]
    fnum5 = [270, 271, 272, 280, 281, 286, 287, 292, 293, 299, 300, 301, 302, 312, 313]
    fnum6 = [314, 321, 322, 323]
    fnum = fnum1 + fnum2 + fnum3 + fnum4 + fnum5 + fnum6
    img_files = [reduce_dir + 'obj_o{0:03d}_clean.fits'.format(ii) for ii in fnum]
    reduce_fli.find_stars(img_files, fwhm=5, threshold=6)

    return

def find_stars_pleiades_ttf():
    reduce_dir = root_dir + 'reduce/pleiades/'
    
    fnum1 = [32, 33, 34, 38, 39, 40, 50, 51, 52, 59, 60, 61, 68, 69, 70, 80, 81, 82, 86]
    fnum2 = [87, 88, 95, 96, 97, 107, 108, 109, 120, 121, 122, 129, 130, 131, 138, 139]
    fnum3 = [140, 150, 151, 152, 159, 160, 161, 168, 169, 170, 180, 181, 182, 189, 190]
    fnum4 = [191, 198, 199, 200, 210, 211, 212, 219, 220, 221, 225, 226, 227, 234, 235]
    fnum5 = [236, 246, 247, 248, 255, 256, 257, 264, 265, 266, 276, 277, 282, 283, 288]
    fnum6 = [289, 294, 295, 306, 307, 308, 315, 316, 317]
    fnum =  fnum1 + fnum2 + fnum3 + fnum4 + fnum5 + fnum6   
    img_files = [reduce_dir + 'obj_ttf{0:03d}_clean.fits'.format(ii) for ii in fnum]
    reduce_fli.find_stars(img_files, fwhm=5, threshold=6)

    return

def find_stars_pleiades_closed():
    reduce_dir = root_dir + 'reduce/pleiades/'
    
    fnum1 = [27, 28, 29, 30, 31, 35, 36, 37, 44, 45, 46, 47, 48, 49, 53, 54, 55, 62, 63, 64]
    fnum2 = [74, 75, 76, 89, 90, 91, 101, 102, 103, 110, 111, 112, 123, 124, 125, 132, 133]
    fnum3 = [134, 141, 142, 143, 153, 154, 155, 162, 163, 164, 171, 172, 173, 183, 184, 185]
    fnum4 = [192, 193, 194, 201, 202, 203, 213, 214, 215, 222, 223, 224, 228, 229, 230, 237]
    fnum5 = [238, 239, 249, 250, 251, 258, 259, 260, 267, 268, 269, 278, 279, 284, 285, 290]
    fnum6 = [291, 296, 297, 309, 310, 311, 318, 319, 320]
    fnum =  fnum1 + fnum2 + fnum3 + fnum4 + fnum5 + fnum6 
    img_files = [reduce_dir + 'obj_c{0:03d}_clean.fits'.format(ii) for ii in fnum]
    reduce_fli.find_stars(img_files, fwhm=3, threshold=6)

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
    

def stack_open():
    """
    Stack the open loop images.
    """
    fnum = [56, 57, 58]
    # fnum += [77, 78, 79, 92, 93, 94, 98, 99, 100, 104, 105, 106, 113, 114]
    # fnum += [115, 126, 127, 128, 135]
    img_files = ['obj_o{0:03d}_bin_nobkg.fits'.format(ii) for ii in fnum]
    starlists = ['obj_o{0:03d}_bin_nobkg_stars.txt'.format(ii) for ii in fnum]

    reduce_fli.shift_and_add(img_files, starlists, 'stack_open.fits', method='meanclip')
    
    return

def stack_ttf():
    """
    Stack the TTF closed loop images.
    """
    return

def stack_closed():
    """
    Stack the closed loop images.
    """
    return
