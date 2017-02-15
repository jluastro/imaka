import pylab as plt
import numpy as np
from astropy.io import fits
from astropy import table
from astropy import units
import glob
from imaka.reduce import reduce_fli
from imaka.reduce import calib
from imaka.reduce import util
import pdb
import os, shutil
#from flystar import match

root_dir = '/Volumes/g/lu/data/imaka/2017_01_11/fli/'

def make_sky():
    sky_dir = '/Volumes/g/lu/data/imaka/2017_01_11/fli/Pleiades/'

    sky_num = [23, 42, 60] #45 sec integration time
    sky_frames = ['{0:s}sky{1:03d}.fits'.format(sky_dir, ss) for ss in sky_num]
    calib.makedark(sky_frames, 'pleiades_sky_023.fits')
    
    sky_num = [25, 43, 61, 78, 79] #30 sec integration time 
    sky_frames = ['{0:s}sky{1:03d}.fits'.format(sky_dir, ss) for ss in sky_num]
    calib.makedark(sky_frames, 'pleiades_sky_025.fits')

    sky_num = [78, 79, 114, 115, 132, 133] #30 sec integration time
    sky_frames = ['{0:s}sky{1:03d}.fits'.format(sky_dir, ss) for ss in sky_num]
    calib.makedark(sky_frames, 'pleiades_sky_078.fits')
    
    sky_num = [195, 196, 197] #30 sec integration time
    sky_frames = ['{0:s}sky{1:03d}.fits'.format(sky_dir, ss) for ss in sky_num]
    calib.makedark(sky_frames, 'pleiades_sky_195.fits')

    return

def make_flat():
    # Didn't take flats, so just copy over the one from Friday night.
    old_flat = imaka_dir + '2017_01_13/fli/calib/flat.fits'
    new_flat = root_dir + 'calib/flat.fits'
    shutil.copyfile(old_flat, new_flat)
    return


def reduce_pleiades_binned_open():
    sky_dir = root_dir + 'reduce/sky/'
    data_dir = root_dir + 'Pleiades/'
    flat_dir = root_dir + 'reduce/calib/'
    out_dir = root_dir + 'reduce/pleiades/'

    util.mkdir(out_dir)
    os.chdir(data_dir)

    fnum = [9, 10, 13, 14, 17, 18, 21, 22, 28, 29, 32, 33, 36, 37, 40, 41, 46, 47, 50, 51, 54, 55, 58, 59]
    img_files = ['obj_o{0:03d}.fits'.format(ii) for ii in fnum]
    reduce_fli.clean_images(img_files, out_dir, rebin=1, sky_frame=sky_dir + 'pleiades_sky_023.fits', flat_frame=flat_dir + 'flat.fits')

    fnum = [64, 65, 68, 69, 72, 73, 76, 77, 82, 83, 86, 87, 90, 91, 94, 95, 100, 101, 104, 105, 108, 109, 112, \
            113, 118, 119, 122, 123, 126, 127, 130, 131]
    img_files = ['obj_o{0:03d}.fits'.format(ii) for ii in fnum]
    reduce_fli.clean_images(img_files, out_dir, rebin=1, sky_frame=sky_dir + 'pleiades_sky_025.fits', flat_frame=flat_dir + 'flat.fits')

    fnum = [141, 142, 149, 150]
    img_files = ['obj_o{0:03d}.fits'.format(ii) for ii in fnum]
    reduce_fli.clean_images(img_files, out_dir, rebin=1, sky_frame=sky_dir + 'pleiades_sky_078.fits', flat_frame=flat_dir + 'flat.fits')
    
    fnum = [179, 180, 185, 186, 193, 194, 200, 201, 206, 207, 212, 213, 218, 219]
    img_files = ['obj_o{0:03d}.fits'.format(ii) for ii in fnum]
    reduce_fli.clean_images(img_files, out_dir, rebin=1, sky_frame=sky_dir + 'pleiades_sky_195.fits', flat_frame=flat_dir + 'flat.fits')
     
    return
    

def reduce_pleiades_binned_closed():
    sky_dir = root_dir + 'reduce/sky/'
    data_dir = root_dir + 'Pleiades/'
    flat_dir = root_dir + 'reduce/flat/'
    out_dir = root_dir + 'reduce/pleiades/'

    os.chdir(data_dir)

    fnum1 = [7, 8, 12, 15, 16, 19, 20, 26, 27, 30, 31, 34, 35, 38, 39, 44, 45, 48, 49, 52]
    fnum2 = [53, 56, 57, 62, 63]
    fnum = fnum1 + fnum2
    img_files = ['obj_c{0:03d}.fits'.format(ii) for ii in fnum]
    reduce_fli.clean_images(img_files, out_dir, rebin=1, sky_frame=sky_dir + 'pleiades_sky_025.fits', flat_frame=flat_dir + 'flat.fits')
    
    fnum1 = [66, 67, 70, 71, 74, 75, 80, 81, 84, 85, 88, 89, 92, 93, 98, 99, 102, 103, 106]
    fnum2 = [107, 110, 111, 116, 117, 120, 121, 124, 125, 128, 129]
    fnum = fnum1 + fnum2
    img_files = ['obj_c{0:03d}.fits'.format(ii) for ii in fnum]
    reduce_fli.clean_images(img_files, out_dir, rebin=1, sky_frame=sky_dir + 'pleiades_sky_78.fits', flat_frame=flat_dir + 'flat.fits')

    fnum1 = [155, 156, 163, 164, 167, 168, 171, 172, 177, 178, 183, 184, 191, 192, 198, 199]
    fnum2 = [204, 205 ,210, 211, 216, 217, 222, 223, 226, 227, 230, 231]
    fnum = fnum1 + fnum2
    img_files = ['obj_c{0:03d}.fits'.format(ii) for ii in fnum]
    reduce_fli.clean_images(img_files, out_dir, rebin=1, sky_frame=sky_dir + 'pleiades_sky_195.fits', flat_frame=flat_dir + 'flat.fits')

    return


def reduce_pleiades_binned_tt():
    sky_dir = root_dir + 'reduce/sky/'
    data_dir = root_dir + 'Pleiades/'
    flat_dir = root_dir + 'reduce/flat/'
    out_dir = root_dir + 'reduce/pleiades/'
    os.chdir(data_dir)
    
    fnum = [134, 135, 136, 137, 138, 139, 140]
    img_files = ['obj_tt{0:03d}.fits'.format(ii) for ii in fnum]
    reduce_fli.clean_images(img_files, out_dir, rebin=1, sky_frame=sky_dir + 'pleiades_sky_078.fits', flat_frame=flat_dir + 'flat.fits')

    return


def reduce_pleiades_binned_ttf():
    sky_dir = root_dir + 'reduce/sky/'
    data_dir = root_dir + 'Pleiades/'
    flat_dir = root_dir + 'reduce/flat/'
    out_dir = root_dir + 'reduce/pleiades/'
    os.chdir(data_dir)
    
    fnum1 = [143, 144, 145, 146, 147, 148, 151, 152, 153, 154, 157, 157, 158, 161, 162, 165, 166]
    fnum2 = [169, 170, 173, 174, 175, 176, 181, 182, 187, 188, 189, 190, 202, 203, 208, 209, 214]
    fnum3 = [215, 220, 221, 224, 225, 228, 229]
    fnum = fnum1 + fnum2 + fnum3
    img_files = ['obj_ttf{0:03d}.fits'.format(ii) for ii in fnum]
    reduce_fli.clean_images(img_files, out_dir, rebin=1, sky_frame=sky_dir + 'pleiades_sky_195.fits', flat_frame=flat_dir + 'flat.fits')
    
    return

    
def find_stars_pleiades_binned_open():
    data_dir = root_dir + 'reduce/'
    os.chdir(data_dir)
      
    fnum1 = [9, 10, 13, 14, 17, 18, 21, 22, 28, 29, 32, 33, 36, 37, 40, 41, 46, 47, 50, 51]
    fnum2 = [54, 55, 58, 59, 64, 65, 68, 69, 72, 73, 76, 77, 82, 83, 86, 87, 90, 91, 94, 95]
    fnum3 = [100, 101, 104, 105, 108, 109, 112, 113, 118, 119, 122, 123, 126, 127, 130, 131]
    fnum4 = [141, 142, 149, 150, 179, 180, 185, 186, 193, 194, 200, 201, 206, 207, 212, 213, 218, 219]
    fnum = fnum1 + fnum2 + fnum3 + fnum4
    img_files = ['obj_o{0:03d}_clean.fits'.format(ii) for ii in fnum]
    reduce_fli.find_stars(img_files, fwhm=5, threshold=6)

    return
    
def find_stars_pleiades_closed():
    data_dir = root_dir + 'reduce/'
    os.chdir(data_dir)
    
    fnum1 = [7, 8, 12, 15, 16, 19, 20, 26, 27, 30, 31, 34, 35, 38, 39, 44, 45, 48, 49, 52]
    fnum2 = [53, 56, 57, 62, 63, 66, 67, 70, 71, 74, 75, 80, 81, 84, 85, 88, 89, 92, 93]
    fnum3 = [98, 99, 102, 103, 106, 107, 110, 111, 116, 117, 120, 121, 124, 125, 128, 129]
    fnum4 = [155, 156, 163, 164, 167, 168, 171, 172, 177, 178, 183, 184, 191, 192, 198, 199]
    fnum5 = [204, 205 ,210, 211, 216, 217, 222, 223, 226, 227, 230, 231]
    fnum = fnum1 #+ fnum2 + fnum3 + fnum4 + fnum5
    img_files = ['obj_c{0:03d}_clean.fits'.format(ii) for ii in fnum]
    reduce_fli.find_stars(img_files, fwhm=5, threshold=6)

    return

def find_stars_pleiades_tt():
    data_dir = root_dir + 'reduce/'
    os.chdir(data_dir)  
    
    fnum = [134, 135, 136, 137, 138, 139, 140]
    img_files = ['obj_tt{0:03d}_clean.fits'.format(ii) for ii in fnum]
    reduce_fli.find_stars(img_files, fwhm=5, threshold=6)
    
    return


def find_stars_pleiades_ttf():
    data_dir = root_dir + 'reduce/'
    os.chdir(data_dir)  
    #[143, 144, 145, 146, 147, 148, 151, 
    fnum1 = [152, 153, 154, 157, 157, 158, 161, 162, 165, 166]
    fnum2 = [169, 170, 173, 174, 175, 176, 181, 182, 187, 188, 189, 190, 202, 203, 208, 209, 214]
    fnum3 = [215, 220, 221, 224, 225, 228, 229]
    fnum = fnum1 + fnum2 + fnum3
    img_files = ['obj_ttf{0:03d}_clean.fits'.format(ii) for ii in fnum]
    reduce_fli.find_stars(img_files, fwhm=5, threshold=6)
    
    return

    
def calc_star_stats_open():
    reduce_dir = root_dir + 'reduce/'
    stats_dir = root_dir + 'stats/'
    os.chdir(data_dir)

    fnum1 = [9, 10, 13, 14, 17, 18, 21, 22, 28, 29, 32, 33, 36, 37, 40, 41, 46, 47, 50, 51]
    fnum2 = [54, 55, 58, 59, 64, 65, 68, 69, 72, 73, 76, 77, 82, 83, 86, 87, 90, 91, 94, 95]
    fnum3 = [100, 101, 104, 105, 108, 109, 112, 113, 118, 119, 122, 123, 126, 127, 130, 131]
    fnum4 = [141, 142, 149, 150, 179, 180, 185, 186, 193, 194, 200, 201, 206, 207, 212, 213, 218, 219]
    fnum = fnum1 + fnum2 + fnum3 + fum4
    img_files = ['obj_o{0:03d}_clean.fits'.format(ii) for ii in fnum]
    reduce_fli.calc_star_stats(img_files, output_stats=stats_dir + 'stats_open.fits')
    
    return

def calc_star_stats_closed():
    reduce_dir = '/Volumes/g/lu/data/imaka/2017_01_11/fli/reduce/'
    stats_dir = '/Volumes/g/lu/data/imaka/2017_01_11/fli/reduce/stats/'
    os.chdir(data_dir)
    
    fnum1 = [7, 8, 12, 15, 16, 19, 20, 26, 27, 30, 31, 34, 35, 38, 39, 44, 45, 48, 49, 52]
    fnum2 = [53, 56, 57, 62, 63, 66, 67, 70, 71, 74, 75, 80, 81, 84, 85, 88, 89, 92, 93]
    fnum3 = [98, 99, 102, 103, 106, 107, 110, 111, 116, 117, 120, 121, 124, 125, 128, 129]
    fnum4 = [155, 156, 163, 164, 167, 168, 171, 172, 177, 178, 183, 184, 191, 192, 198, 199]
    fnum5 = [204, 205 ,210, 211, 216, 217, 222, 223, 226, 227, 230, 231]
    fnum = fnum1 + fnum2 + fnum3 + fnum4 + fnum5
    img_files = ['obj_c{0:03d}_clean.fits'.format(ii) for ii in fnum]    
    reduce_fli.calc_star_stats(img_files, output_stats=stats_dir + 'stats_closed.fits')

    return
