import pylab as plt
import numpy as np
from astropy.io import fits
from astropy import table
from astropy import units
import scipy
import glob
import reduce_fli
import calib
import util
import pdb
import os
from flystar import match

root_dir = '/Volumes/g/lu/data/imaka/2017_01_12/fli/'

def make_sky():

    sky_dir = '/Volumes/g/lu/data/imaka/2017_01_12/fli/Pleiades/'
    out_dir = '/Volumes/g/lu/data/imaka/2017_01_12/fli/reduce/sky/'

    sky_num = [41, 42, 43, 71, 72, 73, 116, 117, 118, 119, 147, 148, 149, 177, 178, 179, 207, 208, 209, 243, 244, 245, 273, 274, 275, 303, 304, 305]
    sky_frames = ['{0:s}sky{1:03d}.fits'.format(sky_dir, ss) for ss in sky_num]
    calib.makedark(sky_frames, 'pleiades_sky.fits')

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


def calc_star_stats():
    reduce_dir = root_dir + 'reduce/pleiades/'
    stats_dir = root_dir + 'reduce/stats/'

    # open loop (removed: 270)
    fnum1 = [56, 57, 58, 65, 66, 67, 77, 78, 79, 92, 93, 94, 98, 99, 100, 104, 105, 106]
    fnum2 = [113, 114, 115, 126, 127, 128, 135, 136, 137, 144, 145, 146, 156, 157, 158]
    fnum3 = [165, 166, 167, 174, 175, 176, 186, 187, 188, 195, 196, 197, 204, 205, 206]
    fnum4 = [216, 217, 218, 231, 232, 233, 240, 241, 242, 252, 253, 254, 261, 262, 263]
    fnum5 = [271, 272, 280, 281, 286, 287, 292, 293, 299, 300, 301, 302, 312, 313]
    fnum6 = [314, 321, 322, 323]
    fnum = fnum1 + fnum2 + fnum3 + fnum4 + fnum5 + fnum6
    img_files = ['{0:s}/obj_o{1:03d}_clean.fits'.format(reduce_dir, ii) for ii in fnum]
    reduce_fli.calc_star_stats(img_files, output_stats=stats_dir + 'stats_open.fits')

    # ttf (removed: 151, 159)
    fnum1 = [32, 33, 34, 38, 39, 40, 50, 51, 52, 59, 60, 61, 68, 69, 70, 80, 81, 82, 86]
    fnum2 = [87, 88, 95, 96, 97, 107, 108, 109, 120, 121, 122, 129, 130, 131, 138, 139]
    fnum3 = [140, 150, 152, 160, 161, 168, 169, 170, 180, 181, 182, 189, 190]
    fnum4 = [191, 198, 199, 200, 210, 211, 212, 219, 220, 221, 225, 226, 227, 234, 235]
    fnum5 = [236, 246, 247, 248, 255, 256, 257, 264, 265, 266, 276, 277, 282, 283, 288]
    fnum6 = [289, 294, 295, 306, 307, 308, 315, 316, 317]
    fnum =  fnum3 + fnum4 + fnum5 +fnum6 #fnum1 + fnum2 + fnum3 + fnum4 + fnum5 + fnum6   
    img_files = ['{0:s}/obj_ttf{1:03d}_clean.fits'.format(reduce_dir, ii) for ii in fnum]
    reduce_fli.calc_star_stats(img_files, output_stats=stats_dir + 'stats_ttf.fits')

    # closed loop (took out: 29)
    fnum1 = [27, 28, 30, 31, 35, 36, 37, 44, 45, 46, 47, 48, 49, 53, 54, 55, 62, 63, 64]
    fnum2 = [74, 75, 76, 89, 90, 91, 101, 102, 103, 110, 111, 112, 123, 124, 125, 132, 133]
    fnum3 = [134, 141, 142, 143, 153, 154, 155, 162, 163, 164, 171, 172, 173, 183, 184, 185]
    fnum4 = [192, 193, 194, 201, 202, 203, 213, 214, 215, 222, 223, 224, 228, 229, 230, 237]
    fnum5 = [238, 239, 249, 250, 251, 258, 259, 260, 267, 268, 269, 278, 279, 284, 285, 290]
    fnum6 = [291, 296, 297, 309, 310, 311, 318, 319, 320]
    fnum = fnum1 + fnum2 + fnum3 + fnum4 + fnum5 + fnum6 
    img_files = ['{0:s}/obj_c{1:03d}_clean.fits'.format(reduce_dir, ii) for ii in fnum]
    reduce_fli.calc_star_stats(img_files, output_stats=stats_dir + 'stats_closed.fits')
    
    return
