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

root_dir = '/Users/fatimaabdurrahman/Desktop/20170216/FLI/'


# def make_sky():
#     data_dir = root_dir + 'Pleiades/'
#     sky_dir = root_dir + 'reduce/sky/'

#     #util.mkdir(out_dir)

#     sky_num = np.arange(89, 110)
#     sky_frames = ['{0:s}sky_{1:04d}.fits'.format(data_dir, ss) for ss in sky_num]
#     calib.makedark(sky_frames, sky_dir+'pleiades_sky.fits')
    
#     return
        
    
def reduce_pleiades():
    sky_dir = root_dir + 'reduce/sky/'
    data_dir = root_dir + 'Pleiades/'
    flat_dir = root_dir + 'reduce/calib/'
    out_dir = root_dir + 'reduce/pleiades/'
    
    #util.mkdir(out_dir)

#     # Open Loop
#     fnum = [59, 62, 65, 68, 70, 76, 82, 88]#, 94]
#     img_files = [data_dir + 'obj{0:04d}_o.fits'.format(ii) for ii in fnum]
#     reduce_fli.clean_images(img_files, out_dir, rebin=1, sky_frame=sky_dir + 'pleiades_sky.fits', flat_frame=flat_dir + 'flat.fits')

#     # Closed Loop - a
#     fnum = [58, 61, 64, 75, 81, 87]#, 93]
#     img_files = [data_dir + 'obj{0:04d}_ca.fits'.format(ii) for ii in fnum]
#     reduce_fli.clean_images(img_files, out_dir, rebin=1, sky_frame=sky_dir + 'pleiades_sky.fits', flat_frame = flat_dir + 'flat.fits')

#     # Closed Loop - b
#     fnum = [71, 77, 83, 89]
#     img_files = [data_dir + "obj{0:04d}_cb.fits".format(ii) for ii in fnum]
#     reduce_fli.clean_images(img_files, out_dir, rebin=1, sky_frame=sky_dir + 'pleiades_sky.fits', flat_frame = flat_dir + 'flat.fits')

#     # Closed Loop - c
#     fnum = [72, 78, 84, 90]
#     img_files = [data_dir + "obj{0:04d}_cc.fits".format(ii) for ii in fnum]
#     reduce_fli.clean_images(img_files, out_dir, rebin=1, sky_frame=sky_dir + 'pleiades_sky.fits', flat_frame = flat_dir + 'flat.fits')

#     # Closed Loop - d
#     fnum = [73, 79, 85]#, 91]
#     img_files = [data_dir + "obj{0:04d}_cd.fits".format(ii) for ii in fnum]
#     reduce_fli.clean_images(img_files, out_dir, rebin=1, sky_frame=sky_dir + 'pleiades_sky.fits', flat_frame = flat_dir + 'flat.fits')

#     # Closed Loop - e
#     fnum = [74, 80, 86]#, 92]
#     img_files = [data_dir + "obj{0:04d}_ce.fits".format(ii) for ii in fnum]
#     reduce_fli.clean_images(img_files, out_dir, rebin=1, sky_frame=sky_dir + 'pleiades_sky.fits', flat_frame = flat_dir + 'flat.fits')
    
    # TT
    fnum = [60, 63, 69]
    img_files = [data_dir + "obj{0:04d}_tt.fits".format(ii) for ii in fnum]
    reduce_fli.clean_images(img_files, out_dir, rebin=1, sky_frame=sky_dir + 'pleiades_sky.fits', flat_frame = flat_dir + 'flat.fits')


    return
    
    
def find_stars_pleiades():
    data_dir = root_dir + 'reduce/pleiades/'

    # Open Loop
    fnum = [59, 62, 65, 68, 70, 76, 82, 88]#, 94]
    img_files = [data_dir + 'obj{0:04d}_o_clean.fits'.format(ii) for ii in fnum]
    reduce_fli.find_stars(img_files, fwhm=5, threshold=6)

    # Closed Loop - a
    fnum = [58, 61, 64, 75, 81, 87]#, 93]
    img_files = [data_dir + 'obj{0:04d}_ca_clean.fits'.format(ii) for ii in fnum]
    reduce_fli.find_stars(img_files, fwhm=5, threshold=6)

    # Closed Loop - b
    fnum = [71, 77, 83, 89]
    img_files = [data_dir + "obj{0:04d}_cb_clean.fits".format(ii) for ii in fnum]
    reduce_fli.find_stars(img_files, fwhm=5, threshold=6)

    # Closed Loop - c
    fnum = [72, 78, 84, 90]
    img_files = [data_dir + "obj{0:04d}_cc_clean.fits".format(ii) for ii in fnum]
    reduce_fli.find_stars(img_files, fwhm=5, threshold=6)

    # Closed Loop - d
    fnum = [73, 79, 85]#, 91]
    img_files = [data_dir + "obj{0:04d}_cd_clean.fits".format(ii) for ii in fnum]
    reduce_fli.find_stars(img_files, fwhm=5, threshold=6)

    # Closed Loop - e
    fnum = [74, 80, 86]#, 92]
    img_files = [data_dir + "obj{0:04d}_ce_clean.fits".format(ii) for ii in fnum]
    reduce_fli.find_stars(img_files, fwhm=5, threshold=6)
    
    # TT
    fnum = [60, 63, 69]
    img_files = [data_dir + "obj{0:04d}_tt_clean.fits".format(ii) for ii in fnum]
    reduce_fli.find_stars(img_files, fwhm=5, threshold=6)
    

    return


def calc_star_stats():
    data_dir = root_dir + 'reduce/pleiades/'
    stats_dir = root_dir +'reduce/stats/'

    
    # Open Loop
    fnum = [59, 62, 65, 68, 70, 76, 82, 88]#, 94]
    img_files = [data_dir + 'obj{0:04d}_o_clean.fits'.format(ii) for ii in fnum]
    reduce_fli.calc_star_stats(img_files, output_stats=stats_dir + 'stats_open.fits')

    # Closed Loop - a
    fnum = [58, 61, 64, 75, 81, 87]#, 93]
    img_files = [data_dir + 'obj{0:04d}_ca_clean.fits'.format(ii) for ii in fnum]
    reduce_fli.calc_star_stats(img_files, output_stats=stats_dir + 'stats_closeda.fits')

    # Closed Loop - b
    fnum = [71, 77, 83, 89]
    img_files = [data_dir + "obj{0:04d}_cb_clean.fits".format(ii) for ii in fnum]
    reduce_fli.calc_star_stats(img_files, output_stats=stats_dir + 'stats_closedb.fits')

    # Closed Loop - c
    fnum = [72, 78, 84, 90]
    img_files = [data_dir + "obj{0:04d}_cc_clean.fits".format(ii) for ii in fnum]
    reduce_fli.calc_star_stats(img_files, output_stats=stats_dir + 'stats_closedc.fits')

    # Closed Loop - d
    fnum = [73, 79, 85]#, 91]
    img_files = [data_dir + "obj{0:04d}_cd_clean.fits".format(ii) for ii in fnum]
    reduce_fli.calc_star_stats(img_files, output_stats=stats_dir + 'stats_closedd.fits')

    # Closed Loop - e
    fnum = [74, 80, 86]#, 92]
    img_files = [data_dir + "obj{0:04d}_ce_clean.fits".format(ii) for ii in fnum]
    reduce_fli.calc_star_stats(img_files, output_stats=stats_dir + 'stats_closede.fits')
    
    # TT
    fnum = [60, 63, 69]
    img_files = [data_dir + "obj{0:04d}_tt_clean.fits".format(ii) for ii in fnum]
    reduce_fli.calc_star_stats(img_files, output_stats=stats_dir + 'stats_tt.fits')

    
    return