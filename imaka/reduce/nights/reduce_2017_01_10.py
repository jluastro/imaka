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

root_dir = '/Volumes/g/lu/data/imaka/2017_01_10/fli/'

def make_sky():

    sky_raw_dir = root_dir + 'Pleiades/'
    sky_out_dir = root_dir + 'reduce/sky/'

    util.mkdir(sky_out_dir)
    
    # 30 second integration time
    sky_num = [20, 21, 30, 31, 46, 47, 56, 57, 70, 71, 84, 85, 98, 99, 112, 113, \
               126, 127, 140, 141, 154, 155, 168, 169, 186, 187]
    sky_frames = ['{0:s}sky{1:03d}.fits'.format(sky_dir, ss) for ss in sky_num]
    calib.makedark(sky_frames, 'pleiades_sky_30s.fits')

    # 45 second integration time
    sky_num = [204, 205, 222, 240, 241]
    sky_frames = ['{0:s}sky{1:03d}.fits'.format(sky_dir, ss) for ss in sky_num]
    calib.makedark(sky_frames, 'pleiades_sky_45s.fits')

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

    # Open Loop, old naming scheme, 30 s
    fnum1 = [10, 11, 14, 15, 18, 19, 24, 25, 28, 29, 34, 35, 38, 39, 40, 41, 44, 45, 50, 51, 55]
    fnum2 = [60, 61, 64, 65, 68, 69, 74, 75, 78, 79, 82, 83, 88, 89, 92, 93, 96, 97, 102, 103]
    fnum = fnum1 + fnum2 
    img_files = [data_dir + 'obj{0:03d}.fits'.format(ii) for ii in fnum]
    reduce_fli.clean_images(img_files, out_dir, rebin=1, sky_frame=sky_dir + 'pleiades_sky_30s.fits', flat_frame=flat_dir + 'flat.fits')

    # Open loop, new naming scheme, 30 s
    fnum1 = [106, 107, 110, 111, 116, 117, 120, 121, 124, 125, 130, 131, 134, 135, 138, 139, 144]
    fnum2 = [145, 148, 149, 152, 153, 158, 159, 162, 163, 166, 167, 172, 173, 176, 177, 180, 181]
    fnum = fnum1 + fnum2
    img_files = [data_dir + 'obj_o{0:03d}.fits'.format(ii) for ii in fnum]
    reduce_fli.clean_images(img_files, out_dir, rebin=1, sky_frame=sky_dir + 'pleiades_sky_30s.fits', flat_frame=flat_dir + 'flat.fits')
    
    # Open Loop, new naming scheme, 45 s
    fnum1 = [184, 185, 190, 191, 194, 195, 198, 199, 202, 203, 208, 209, 212, 213, 216, 217]
    fnum2 = [220, 221, 226, 227, 230, 231, 235, 238, 239, 244, 245, 248, 249, 252, 253]
    fnum = fnum1 + fnum2
    img_files = [data_dir + 'obj_o{0:03d}.fits'.format(ii) for ii in fnum]
    reduce_fli.clean_images(img_files, out_dir, rebin=1, sky_frame=sky_dir + 'pleiades_sky_45s.fits', flat_frame=flat_dir + 'flat.fits')

    # Closed Loop, old naming scheme, 30 s
    fnum1 = [8, 9, 12, 13, 16, 17, 22, 23, 26, 27, 32, 33, 36, 37, 42, 43, 48, 49, 52, 53, 58, 59]
    fnum2 = [62, 63, 66, 67, 72, 73, 76, 77, 80, 81, 86, 87, 90, 91, 94, 95, 100, 101, 104, 105]
    fnum =  fnum1 + fnum2
    img_files = [data_dir + 'obj{0:03d}.fits'.format(ii) for ii in fnum]
    reduce_fli.clean_images(img_files, out_dir, rebin=1, sky_frame=sky_dir + 'pleiades_sky_30s.fits', flat_frame = flat_dir + 'flat.fits')
    
    # Closed Loop, new naming scheme, 30 s
    fnum1 = [108, 109, 114, 115, 118, 119, 122, 123, 128, 129, 132, 133, 136, 137, 142, 143, 146]
    fnum2 = [147, 150, 151, 156, 157, 160, 161, 164, 165, 170, 171, 174, 175, 178, 179, 182, 183]
    fnum3 = [188, 189, 192, 193, 196, 197, 200, 201, 206, 207, 210, 211, 214, 215, 218, 219, 224]
    fnum4 = [225, 228, 229, 232, 233, 236, 237, 242, 243, 246, 247, 250, 251]
    fnum = fnum1 + fnum2 + fnum3 + fnum4
    img_files = [data_dir + 'obj_c{0:03d}.fits'.format(ii) for ii in fnum]
    reduce_fli.clean_images(img_files, out_dir, rebin=1, sky_frame=sky_dir + 'pleiades_sky_30s.fits', flat_frame = flat_dir + 'flat.fits')
    
    return
    
    
def find_stars_pleiades_open():
    reduce_dir = root_dir + 'reduce/pleiades/'
    
    # Old Naming Scheme
#     fnum1 = [10, 11, 14, 15, 18, 19, 24, 25, 28, 29, 34, 35, 38, 39, 40, 41, 44, 45, 50, 51, 55]
#     fnum2 = [60, 61, 64, 65, 68, 69, 74, 75, 78, 79, 82, 83, 88, 89, 92, 93, 96, 97, 102, 103]
#     fnum = fnum1 + fnum2 
#     img_files = [reduce_dir + 'obj{0:03d}_clean.fits'.format(ii) for ii in fnum]
#     reduce_fli.find_stars(img_files, fwhm=5, threshold=6)

    # New Naming Scheme
    #fnum1 = [106, 107, 110, 111, 116, 117, 120, 121, 124, 125, 130, 131, 134, 135, 138, 139, 144]
    #fnum2 = [145, 148, 149, 152, 153, 158, 159, 162, 163, 166, 167, 172, 173, 176, 177, 180, 181]
    fnum3 = [184, 185, 190, 191, 194, 195, 198, 199, 202, 203, 208, 209, 212, 213, 216, 217]
    fnum4 = [220, 221, 226, 227, 230, 231, 235, 238, 239, 244, 245, 248, 249, 252, 253]
    #fnum = fnum1 + fnum2 + fnum3 + fnum4
    fnum = fnum3 + fnum4
    img_files = [reduce_dir + 'obj_o{0:03d}_clean.fits'.format(ii) for ii in fnum]
    reduce_fli.find_stars(img_files, fwhm=5, threshold=6)

    return


def find_stars_pleiades_closed():
    reduce_dir = root_dir + 'reduce/pleiades/'
    
    # Old Naming Scheme
    fnum1 = [8, 9, 12, 13, 16, 17, 22, 23, 26, 27, 32, 33, 36, 37, 42, 43, 48, 49, 52, 53, 58, 59]
    fnum2 = [62, 63, 66, 67, 72, 73, 76, 77, 80, 81, 86, 87, 90, 91, 94, 95, 100, 101, 104, 105]
    fnum =  fnum1 + fnum2
    img_files = [reduce_dir + 'obj{0:03d}_clean.fits'.format(ii) for ii in fnum]
    reduce_fli.find_stars(img_files, fwhm=3, threshold=6)
    
    #New Naming Scheme
    fnum1 = [108, 109, 114, 115, 118, 119, 122, 123, 128, 129, 132, 133, 136, 137, 142, 143, 146]
    fnum2 = [147, 150, 151, 156, 157, 160, 161, 164, 165, 170, 171, 174, 175, 178, 179, 182, 183]
    fnum3 = [188, 189, 192, 193, 196, 197, 200, 201, 206, 207, 210, 211, 214, 215, 218, 219, 224]
    fnum4 = [225, 228, 229, 232, 233, 236, 237, 242, 243, 246, 247, 250, 251]
    fnum = fnum1 + fnum2 + fnum3 + fnum4
    img_files = [reduce_dir + 'obj_c{0:03d}_clean.fits'.format(ii) for ii in fnum]
    reduce_fli.find_stars(img_files, fwhm=3, threshold=6)
    return


def calc_star_stats_open():
    reduce_dir = '/Volumes/g/lu/data/imaka/2017_01_10/fli/reduce/pleiades/'
    stats_dir = '/Volumes/g/lu/data/imaka/2017_01_10/fli/reduce/stats/'

    # Old Naming Scheme
    fnum1 = [10, 11, 14, 15, 18, 19, 24, 25, 28, 29, 34, 35, 38, 39, 40, 41, 44, 45, 50, 51, 55]
    fnum2 = [60, 61, 64, 65, 68, 69, 74, 75, 78, 79, 82, 83, 88, 89, 92, 93, 96, 97, 102, 103]
    fnum = fnum1 + fnum2 
    img_files_old_name = [reduce_dir + 'obj{0:03d}_clean.fits'.format(ii) for ii in fnum]
    
    #New Naming Scheme
    fnum1 = [106, 107, 110, 111, 116, 117, 120, 121, 124, 125, 130, 131, 134, 135, 138, 139, 144]
    fnum2 = [145, 148, 149, 152, 153, 158, 159, 162, 163, 166, 167, 172, 173, 176, 177, 180, 181]
    fnum3 = [184, 185, 190, 191, 194, 195, 198, 199, 202, 203, 208, 209, 212, 213, 216, 217]
    fnum4 = [220, 221, 226, 227, 230, 231, 235, 238, 239, 244, 245, 248, 249, 252, 253]
    fnum = fnum1 + fnum2 + fnum3 + fnum4
    img_files_new_name = [reduce_dir + 'obj_o{0:03d}_clean.fits'.format(ii) for ii in fnum]
    
    img_files = img_files_old_name + img_files_new_name
    
    reduce_fli.calc_star_stats(img_files, output_stats=stats_dir + 'stats_open.fits')
    
    return


def calc_star_stats_closed():
    reduce_dir = '/Volumes/g/lu/data/imaka/2017_01_10/fli/reduce/pleiades/'
    stats_dir = '/Volumes/g/lu/data/imaka/2017_01_10/fli/reduce/stats/'

    # Old Naming Scheme
    fnum1 = [8, 9, 12, 13, 16, 17, 22, 23, 26, 27, 32, 33, 36, 37, 42, 43, 48, 49, 52, 53, 58, 59]
    fnum2 = [62, 63, 66, 67, 72, 73, 76, 77, 80, 81, 86, 87, 90, 91, 94, 95, 100, 101, 104, 105]
    fnum =  fnum1 + fnum2
    img_files_old_name = [reduce_dir + 'obj{0:03d}_clean.fits'.format(ii) for ii in fnum]
    
    #New Naming Scheme
    fnum1 = [108, 109, 114, 115, 118, 119, 122, 123, 128, 129, 132, 133, 136, 137, 142, 143, 146]
    fnum2 = [147, 150, 151, 156, 157, 160, 161, 164, 165, 170, 171, 174, 175, 178, 179, 182, 183]
    fnum3 = [188, 189, 192, 193, 196, 197, 200, 201, 206, 207, 210, 211, 214, 215, 218, 219, 224]
    fnum4 = [225, 228, 229, 232, 233, 236, 237, 242, 243, 246, 247, 250, 251]
    fnum = fnum1 + fnum2 + fnum3 + fnum4
    img_files_new_name = [reduce_dir + 'obj_c{0:03d}_clean.fits'.format(ii) for ii in fnum]
    
    img_files = img_files_old_name + img_files_new_name
    
    reduce_fli.calc_star_stats(img_files, output_stats=stats_dir + 'stats_closed.fits')
    
    return

