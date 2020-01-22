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
import os, shutil
import pdb

root_dir = '/Users/fatimaabdurrahman/Desktop/Research/RUN5/20170521/FLI/'

#def make_flat():
     #use flat from 20170520, no flat taken this night
    

#def make_sky():
    #use sky from 20170521, no sky taken this night (actually check to make sure at least integration times match up)
    
    #return
 
def reduce_FLD2():
    sky_dir = root_dir + 'reduce/sky/' 
    data_dir = root_dir + 'FLD2_2/'
    flat_dir = root_dir + 'reduce/calib/'
    out_dir = root_dir + 'reduce/FLD2_2/'
    
    util.mkdir(out_dir)

    
    # Open Loop
    fnum = [51, 53, 56, 59, 62, 65, 68, 71, 74, 77, 80, 83, 86, 89, 92, 96, 99, 102]
    fnum += [105, 108, 111, 114, 117, 120, 123, 126, 129, 132, 135, 138, 141, 144]
    fnum += [147, 150, 153, 156, 159, 162, 165, 168, 171, 180, 183, 186, 189, 192]
    fnum += [195, 198, 201, 204, 207, 210, 213]
    img_files = [data_dir + 'obj{0:04d}_o.fits'.format(ii) for ii in fnum]
    reduce_fli.clean_images(img_files, out_dir, rebin=1, sky_frame=sky_dir + 'FLD2_2_sky.fits', flat_frame=flat_dir + 'flat.fits')


    # Closed
    fnum = [49, 52, 55, 58, 61, 64, 67, 70, 73, 76, 79, 82, 85, 88, 91, 95, 98, 101]
    fnum += [104, 107, 110, 113, 116, 119, 122, 125, 128, 131, 134, 137, 140, 143, 146]
    fnum += [149, 152, 155, 158, 161, 164, 167, 170, 173, 176, 179, 182, 185, 188, 191]
    fnum += [194, 197, 200, 203, 206, 209, 212]
    img_files = [data_dir + 'obj{0:04d}_c.fits'.format(ii) for ii in fnum]
    reduce_fli.clean_images(img_files, out_dir, rebin=1, sky_frame=sky_dir + 'FLD2_2_sky.fits', flat_frame = flat_dir + 'flat.fits')

    # Closed A
    fnum = [50, 54, 57, 60, 63, 66, 69, 72, 75, 78, 81, 84, 87, 90, 93, 97, 100, 103, 106]
    fnum += [109, 112, 115, 118, 121, 124, 127, 130, 133, 136, 139, 142, 145, 148, 151, 154]
    fnum += [157, 160, 163, 166, 169, 172, 175, 178, 181, 184, 187, 190, 193, 196, 199, 202]
    fnum += [205, 208, 211, 214]
    img_files = [data_dir + 'obj{0:04d}_cA.fits'.format(ii) for ii in fnum]
    reduce_fli.clean_images(img_files, out_dir, rebin=1, sky_frame=sky_dir + 'FLD2_2_sky.fits', flat_frame = flat_dir + 'flat.fits')

    return
    
def find_stars_FLD2():
    data_dir = root_dir + 'reduce/FLD2_2/'

    # Open 
    fnum = [51, 53, 56, 59, 62, 65, 68, 71, 74, 77, 80, 83, 86, 89, 92, 96, 99, 102]
    fnum += [105, 108, 111, 114, 117, 120, 123, 126, 129, 132, 135, 138, 141, 144]
    fnum += [147, 150, 153, 156, 159, 162, 165, 168, 171, 180, 183, 186, 189, 192]
    fnum += [195, 198, 201, 204, 207, 210, 213]
    img_files = [data_dir + 'obj{0:04d}_o_clean.fits'.format(ii) for ii in fnum]
    reduce_fli.find_stars(img_files, fwhm=8, threshold=10)
    
    #Closed 
    
    fnum = [49, 52, 55, 58, 61, 64, 67, 70, 73, 76, 79, 82, 85, 88, 91, 95, 98, 101]
    fnum += [104, 107, 110, 113, 116, 119, 122, 125, 128, 131, 134, 137, 140, 143, 146]
    fnum += [149, 152, 155, 158, 161, 164, 167, 170, 173, 176, 179, 182, 185, 188, 191]
    fnum += [194, 197, 200, 203, 206, 209, 212]
    img_files = [data_dir + 'obj{0:04d}_c_clean.fits'.format(ii) for ii in fnum]
    reduce_fli.find_stars(img_files, fwhm=6, threshold=6)    

    # Closed A
    fnum = [50, 54, 57, 60, 63, 66, 69, 72, 75, 78, 81, 84, 87, 90, 93, 97, 100, 103, 106]
    fnum += [109, 112, 115, 118, 121, 124, 127, 130, 133, 136, 139, 142, 145, 148, 151, 154]
    fnum += [157, 160, 163, 166, 169, 172, 175, 178, 181, 184, 187, 190, 193, 196, 199, 202]
    fnum += [205, 208, 211, 214]
    img_files = [data_dir + 'obj{0:04d}_cA_clean.fits'.format(ii) for ii in fnum]
    reduce_fli.find_stars(img_files, fwhm=6, threshold=6)    

    return


def calc_star_stats():
    data_dir = root_dir + 'reduce/FLD2_2/'
    stats_dir = root_dir +'reduce/stats/'

    
    # Open Loop
    fnum = [51, 53, 56, 59, 62, 65, 68, 71, 74, 77, 80, 83, 86, 89, 92, 96, 99, 102]
    fnum += [105, 108, 111, 114, 117, 120, 123, 126, 129, 132, 135, 138, 141, 144]
    fnum += [147, 150, 153, 156, 159, 162, 165, 168, 171, 180, 183, 186, 189, 192]
    fnum += [195, 198, 201, 204, 207, 210, 213]
    img_files = [data_dir + 'obj{0:04d}_o_clean.fits'.format(ii) for ii in fnum]
    reduce_fli.calc_star_stats(img_files, output_stats=stats_dir + 'stats_open.fits')

    # Closed 
    fnum = [49, 52, 55, 58, 61, 64, 67, 70, 73, 76, 79, 82, 85, 88, 91, 95, 98, 101]
    fnum += [104, 107, 110, 113, 116, 119, 122, 125, 128, 131, 134, 137, 140, 143, 146]
    fnum += [149, 152, 155, 158, 161, 164, 167, 170, 173, 176, 179, 182, 185, 188, 191]
    fnum += [194, 197, 200, 203, 206, 209, 212]
    img_files = [data_dir + "obj{0:04d}_c_clean.fits".format(ii) for ii in fnum]
    reduce_fli.calc_star_stats(img_files, output_stats=stats_dir + 'stats_closedA.fits')

    # Closed A
    fnum = [50, 54, 57, 60, 63, 66, 69, 72, 75, 78, 81, 84, 87, 90, 93, 97, 100, 103, 106]
    fnum += [109, 112, 115, 118, 121, 124, 127, 130, 133, 136, 139, 142, 145, 148, 151, 154]
    fnum += [157, 160, 163, 166, 169, 172, 175, 178, 181, 184, 187, 190, 193, 196, 199, 202]
    fnum += [205, 208, 211, 214]
    img_files = [data_dir + "obj{0:04d}_cA_clean.fits".format(ii) for ii in fnum]
    reduce_fli.calc_star_stats(img_files, output_stats=stats_dir + 'stats_closedA.fits')
    
    return 

def stack_FLD2():

    data_dir = root_dir + 'reduce/FLD2_2/'
    stats_dir = root_dir + 'reduce/stats/'
    stacks_dir = root_dir + 'reduce/stacks/'

    util.mkdir(stacks_dir)

    # Open Loop
    fnum = [51, 53, 56, 59, 62, 65, 68, 71, 74, 77, 80, 83, 86, 89, 92, 96, 99, 102]
    fnum += [105, 108, 111, 114, 117, 120, 123, 126, 129, 132, 135, 138, 141, 144]
    fnum += [147, 150, 153, 156, 159, 162, 165, 168, 171, 180, 183, 186, 189, 192]
    fnum += [195, 198, 201, 204, 207, 210, 213]
    open_images = [data_dir + 'obj{0:04d}_o_clean.fits'.format(ii) for ii in fnum]
    open_starlists = [data_dir + 'obj{0:04d}_o_clean_stars.txt'.format(ii) for ii in fnum]
    open_output_root = stacks_dir + 'FLD2_stack_open'
    reduce_fli.shift_and_add(open_images, open_starlists, open_output_root, method='mean')
    
    # Closed 
    fnum = [49, 52, 55, 58, 61, 64, 67, 70, 73, 76, 79, 82, 85, 88, 91, 95, 98, 101]
    fnum += [104, 107, 110, 113, 116, 119, 122, 125, 128, 131, 134, 137, 140, 143, 146]
    fnum += [149, 152, 155, 158, 161, 164, 167, 170, 173, 176, 179, 182, 185, 188, 191]
    fnum += [194, 197, 200, 203, 206, 209, 212]
    closed_images = [data_dir + 'obj{0:04d}_c_clean.fits'.format(ii) for ii in fnum]
    closed_starlists = [data_dir + 'obj{0:04d}_c_clean_stars.txt'.format(ii) for ii in fnum]
    closed_output_root = stacks_dir + 'FLD2_stack_closed'
    reduce_fli.shift_and_add(closed_images, closed_starlists, closed_output_root, method='mean')
    
    # Closed A
    fnum = [50, 54, 57, 60, 63, 66, 69, 72, 75, 78, 81, 84, 87, 90, 93, 97, 100, 103, 106]
    fnum += [109, 112, 115, 118, 121, 124, 127, 130, 133, 136, 139, 142, 145, 148, 151, 154]
    fnum += [157, 160, 163, 166, 169, 172, 175, 178, 181, 184, 187, 190, 193, 196, 199, 202]
    fnum += [205, 208, 211, 214]
    closed_images = [data_dir + 'obj{0:04d}_cA_clean.fits'.format(ii) for ii in fnum]
    closed_starlists = [data_dir + 'obj{0:04d}_cA_clean_stars.txt'.format(ii) for ii in fnum]
    closed_output_root = stacks_dir + 'FLD2_stack_closedA'
    reduce_fli.shift_and_add(closed_images, closed_starlists, closed_output_root, method='mean')

    return

def analyze_stacks():
    data_dir = root_dir + 'reduce/stacks/'
    stats_dir = root_dir + 'reduce/stats/'
    
    img_files = [data_dir + 'FLD2_stack_open.fits',
                 data_dir + 'FLD2_stack_closed.fits',
                 data_dir + 'FLD2_stack_closedA.fits']
    
    #Find stars in image
    reduce_fli.find_stars(img_files, fwhm=5, threshold=10, N_passes=2, plot_psf_compare=False)
    
    # Calc stats on all the stacked images
    reduce_fli.calc_star_stats(img_files, output_stats= stats_dir + 'stats_stacks.fits')

    return
    
