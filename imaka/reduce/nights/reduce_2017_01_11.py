import pylab as plt
import numpy as np
from astropy.io import fits
from astropy import table
from astropy import units
import glob
import reduce_fli
import calib
import pdb
import os

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

def reduce_pleiades_binned_open():
    sky_dir = '/Volumes/g/lu/data/imaka/2017_01_11/fli/Pleiades/'
    data_dir = '/Volumes/g/lu/data/imaka/2017_01_11/fli/Pleiades/'
    flat_dir = '/Volumes/g/lu/data/imaka/2017_01_13/fli/twilights/'
    out_dir = '/Volumes/g/lu/data/imaka/2017_01_11/fli/reduce/'
    os.chdir(data_dir)

    fnum = [9, 10, 13, 14, 17, 18, 21, 22, 28, 29, 32, 33, 36, 37, 40, 41, 46, 47, 50, 51, 54, 55, 58, 59]
    img_files = ['obj_o{0:03d}.fits'.format(ii) for ii in fnum]
    reduce_fli.flat_sky_reduction(img_files, output_path=out_dir, sky_frame=sky_dir + 'pleiades_sky_023.fits', flat_frame=flat_dir + 'flat.fits')

    fnum = [64, 65, 68, 69, 72, 73, 76, 77, 82, 83, 86, 87, 90, 91, 94, 95, 100, 101, 104, 105, 108, 109, 112, \
            113, 118, 119, 122, 123, 126, 127, 130, 131]
    img_files = ['obj_o{0:03d}.fits'.format(ii) for ii in fnum]
    reduce_fli.flat_sky_reduction(img_files, output_path=out_dir, sky_frame=sky_dir + 'pleiades_sky_025.fits', flat_frame=flat_dir + 'flat.fits')

    fnum = [141, 142, 149, 150]
    img_files = ['obj_o{0:03d}.fits'.format(ii) for ii in fnum]
    reduce_fli.flat_sky_reduction(img_files, output_path=out_dir, sky_frame=sky_dir + 'pleiades_sky_078.fits', flat_frame=flat_dir + 'flat.fits')
    
    fnum = [179, 180, 185, 186, 193, 194, 200, 201, 206, 207, 212, 213, 218, 219]
    img_files = ['obj_o{0:03d}.fits'.format(ii) for ii in fnum]
    reduce_fli.flat_sky_reduction(img_files, output_path=out_dir, sky_frame=sky_dir + 'pleiades_sky_195.fits', flat_frame=flat_dir + 'flat.fits')
     
    return
    

def reduce_pleiades_binned_closed():
    sky_dir = '/Volumes/g/lu/data/imaka/2017_01_11/fli/Pleiades/'
    data_dir = '/Volumes/g/lu/data/imaka/2017_01_11/fli/Pleiades/'
    flat_dir = '/Volumes/g/lu/data/imaka/2017_01_13/fli/twilights/'
    out_dir = '/Volumes/g/lu/data/imaka/2017_01_11/fli/reduce/'

    os.chdir(data_dir)

    fnum1 = [7, 8, 12, 15, 16, 19, 20, 26, 27, 30, 31, 34, 35, 38, 39, 44, 45, 48, 49, 52]
    fnum2 = [53, 56, 57, 62, 63]
    fnum = fnum1 + fnum2
    img_files = ['obj_c{0:03d}.fits'.format(ii) for ii in fnum]
    reduce_fli.flat_sky_reduction(img_files, output_path=out_dir, sky_frame=sky_dir + 'pleiades_sky_025.fits', flat_frame=flat_dir + 'flat.fits')
    
    fnum1 = [66, 67, 70, 71, 74, 75, 80, 81, 84, 85, 88, 89, 92, 93, 98, 99, 102, 103, 106]
    fnum2 = [107, 110, 111, 116, 117, 120, 121, 124, 125, 128, 129]
    fnum = fnum1 + fnum2
    img_files = ['obj_c{0:03d}.fits'.format(ii) for ii in fnum]
    reduce_fli.flat_sky_reduction(img_files, output_path=out_dir, sky_frame=sky_dir + 'pleiades_sky_78.fits', flat_frame=flat_dir + 'flat.fits')

    fnum1 = [155, 156, 163, 164, 167, 168, 171, 172, 177, 178, 183, 184, 191, 192, 198, 199]
    fnum2 = [204, 205 ,210, 211, 216, 217, 222, 223, 226, 227, 230, 231]
    fnum = fnum1 + fnum2
    img_files = ['obj_c{0:03d}.fits'.format(ii) for ii in fnum]
    reduce_fli.flat_sky_reduction(img_files, output_path=out_dir, sky_frame=sky_dir + 'pleiades_sky_195.fits', flat_frame=flat_dir + 'flat.fits')

    return


def reduce_pleiades_binned_tt():
    sky_dir = '/Volumes/g/lu/data/imaka/2017_01_11/fli/Pleiades/'
    data_dir = '/Volumes/g/lu/data/imaka/2017_01_11/fli/Pleiades/'
    flat_dir = '/Volumes/g/lu/data/imaka/2017_01_13/fli/twilights/'
    out_dir = '/Volumes/g/lu/data/imaka/2017_01_11/fli/reduce/'
    os.chdir(data_dir)
    
    fnum = [134, 135, 136, 137, 138, 139, 140]
    img_files = ['obj_tt{0:03d}.fits'.format(ii) for ii in fnum]
    reduce_fli.flat_sky_reduction(img_files, output_path=out_dir, sky_frame=sky_dir + 'pleiades_sky_078.fits', flat_frame=flat_dir + 'flat.fits')

    return


def reduce_pleiades_binned_ttf():
    sky_dir = '/Volumes/g/lu/data/imaka/2017_01_11/fli/Pleiades/'
    data_dir = '/Volumes/g/lu/data/imaka/2017_01_11/fli/Pleiades/'
    flat_dir = '/Volumes/g/lu/data/imaka/2017_01_13/fli/twilights/'
    out_dir = '/Volumes/g/lu/data/imaka/2017_01_11/fli/reduce/'
    os.chdir(data_dir)
    
    fnum1 = [143, 144, 145, 146, 147, 148, 151, 152, 153, 154, 157, 157, 158, 161, 162, 165, 166]
    fnum2 = [169, 170, 173, 174, 175, 176, 181, 182, 187, 188, 189, 190, 202, 203, 208, 209, 214]
    fnum3 = [215, 220, 221, 224, 225, 228, 229]
    fnum = fnum1 + fnum2 + fnum3
    img_files = ['obj_ttf{0:03d}.fits'.format(ii) for ii in fnum]
    reduce_fli.flat_sky_reduction(img_files, output_path=out_dir, sky_frame=sky_dir + 'pleiades_sky_195.fits', flat_frame=flat_dir + 'flat.fits')
    
    return

    
def find_stars_pleiades_binned_open():
    data_dir = '/Volumes/g/lu/data/imaka/2017_01_11/fli/reduce/'
    os.chdir(data_dir)
      
    fnum1 = [9, 10, 13, 14, 17, 18, 21, 22, 28, 29, 32, 33, 36, 37, 40, 41, 46, 47, 50, 51]
    fnum2 = [54, 55, 58, 59, 64, 65, 68, 69, 72, 73, 76, 77, 82, 83, 86, 87, 90, 91, 94, 95]
    fnum3 = [100, 101, 104, 105, 108, 109, 112, 113, 118, 119, 122, 123, 126, 127, 130, 131]
    fnum4 = [141, 142, 149, 150, 179, 180, 185, 186, 193, 194, 200, 201, 206, 207, 212, 213, 218, 219]
    fnum = fnum1 + fnum2 + fnum3 + fnum4
    img_files = ['obj_o{0:03d}_bin_nobkg.fits'.format(ii) for ii in fnum]
    reduce_fli.find_stars_bin(img_files, fwhm=8, threshold=6)

    return
    
def find_stars_pleiades_binned_closed():
    data_dir = '/Volumes/g/lu/data/imaka/2017_01_11/fli/reduce/'
    os.chdir(data_dir)
    
    fnum1 = [7, 8, 12, 15, 16, 19, 20, 26, 27, 30, 31, 34, 35, 38, 39, 44, 45, 48, 49, 52]
    fnum2 = [53, 56, 57, 62, 63, 66, 67, 70, 71, 74, 75, 80, 81, 84, 85, 88, 89, 92, 93]
    fnum3 = [98, 99, 102, 103, 106, 107, 110, 111, 116, 117, 120, 121, 124, 125, 128, 129]
    fnum4 = [155, 156, 163, 164, 167, 168, 171, 172, 177, 178, 183, 184, 191, 192, 198, 199]
    fnum5 = [204, 205 ,210, 211, 216, 217, 222, 223, 226, 227, 230, 231]
    fnum = fnum1 #+ fnum2 + fnum3 + fnum4 + fnum5
    img_files = ['obj_c{0:03d}_bin_nobkg.fits'.format(ii) for ii in fnum]
    reduce_fli.find_stars_bin(img_files, fwhm=4, threshold=6)

    return

def find_stars_pleiades_binned_tt():
    data_dir = '/Volumes/g/lu/data/imaka/2017_01_11/fli/reduce/'
    os.chdir(data_dir)  
    
    fnum = [134, 135, 136, 137, 138, 139, 140]
    img_files = ['obj_tt{0:03d}_bin_nobkg.fits'.format(ii) for ii in fnum]
    reduce_fli.find_stars_bin(img_files, fwhm=4, threshold=6)
    
    return


def find_stars_pleiades_binned_ttf():
    data_dir = '/Volumes/g/lu/data/imaka/2017_01_11/fli/reduce/'
    os.chdir(data_dir)  
    
    fnum1 = [143, 144, 145, 146, 147, 148, 151, 152, 153, 154, 157, 157, 158, 161, 162, 165, 166]
    fnum2 = [169, 170, 173, 174, 175, 176, 181, 182, 187, 188, 189, 190, 202, 203, 208, 209, 214]
    fnum3 = [215, 220, 221, 224, 225, 228, 229]
    fnum = fnum1 + fnum2 + fnum3
    img_files = ['obj_ttf{0:03d}_bin_nobkg.fits'.format(ii) for ii in fnum]
    reduce_fli.find_stars_bin(img_files, fwhm=4, threshold=6)
    
    return

    
def calc_star_stats_open():
    data_dir = '/Volumes/g/lu/data/imaka/2017_01_11/fli/reduce/'
    os.chdir(data_dir)

    fnum1 = [9, 10, 13, 14, 17, 18, 21, 22, 28, 29, 32, 33, 36, 37, 40, 41, 46, 47, 50, 51]
    fnum2 = [54, 55, 58, 59, 64, 65, 68, 69, 72, 73, 76, 77, 82, 83, 86, 87, 90, 91, 94, 95]
    fnum3 = [100, 101, 104, 105, 108, 109, 112, 113, 118, 119, 122, 123, 126, 127, 130, 131]
    fnum4 = [141, 142, 149, 150, 179, 180, 185, 186, 193, 194, 200, 201, 206, 207, 212, 213, 218, 219]
    fnum = fnum1 + fnum2 + fnum3 + fum4
    img_files = ['obj_o{0:03d}_bin_nobkg.fits'.format(ii) for ii in fnum]
    reduce_fli.calc_star_stats(img_files, output_stats='stats_open2.fits')
    
    return

def calc_star_stats_closed():
    data_dir = '/Volumes/g/lu/data/imaka/2017_01_11/fli/reduce/'
    os.chdir(data_dir)
    
    fnum1 = [7, 8, 12, 15, 16, 19, 20, 26, 27, 30, 31, 34, 35, 38, 39, 44, 45, 48, 49, 52]
    fnum2 = [53, 56, 57, 62, 63, 66, 67, 70, 71, 74, 75, 80, 81, 84, 85, 88, 89, 92, 93]
    fnum3 = [98, 99, 102, 103, 106, 107, 110, 111, 116, 117, 120, 121, 124, 125, 128, 129]
    fnum4 = [155, 156, 163, 164, 167, 168, 171, 172, 177, 178, 183, 184, 191, 192, 198, 199]
    fnum5 = [204, 205 ,210, 211, 216, 217, 222, 223, 226, 227, 230, 231]
    fnum = fnum1 + fnum2 + fnum3 + fnum4 + fnum5
    img_files = ['obj_c{0:03d}_bin_nobkg.fits'.format(ii) for ii in fnum]    
    reduce_fli.calc_star_stats(img_files, output_stats='stats_closed.fits')

    return


def calc_star_stats_tt():
    data_dir = '/Volumes/g/lu/data/imaka/2017_01_11/fli/reduce/'
    os.chdir(data_dir)
        
    fnum = [134, 135, 136, 137, 138, 139, 140]
    img_files = ['obj_tt{0:03d}_bin_nobkg.fits'.format(ii) for ii in fnum]    
    reduce_fli.calc_star_stats(img_files, output_stats='stats_tt.fits')
    
    return


def calc_star_stats_ttf():
    data_dir = '/Volumes/g/lu/data/imaka/2017_01_11/fli/reduce/'
    os.chdir(data_dir)
        
    fnum1 = [143, 144, 145, 146, 147, 148, 151, 152, 153, 154, 157, 157, 158, 161, 162, 165, 166]
    fnum2 = [169, 170, 173, 174, 175, 176, 181, 182, 187, 188, 189, 190, 202, 203, 208, 209, 214]
    fnum3 = [215, 220, 221, 224, 225, 228, 229]
    fnum = fnum1 + fnum2 + fnum3
    img_files = ['obj_ttf{0:03d}_bin_nobkg.fits'.format(ii) for ii in fnum]    
    reduce_fli.calc_star_stats(img_files, output_stats='stats_ttf.fits')
    
    return
