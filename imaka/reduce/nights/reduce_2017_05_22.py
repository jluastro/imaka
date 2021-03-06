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

root_dir = '/Users/fatimaabdurrahman/Desktop/Research/RUN5/20170522/FLI/'

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
    fnum = [124, 129, 135, 141, 147, 153, 159, 165, 171, 177, 183, 189, 195]
    fnum += [201, 207, 215, 221, 227, 233, 239, 245, 251, 257, 263, 269]
    img_files = [data_dir + 'obj{0:04d}_o.fits'.format(ii) for ii in fnum]
    reduce_fli.clean_images(img_files, out_dir, rebin=1, sky_frame=sky_dir + 'FLD2_2_sky.fits', flat_frame=flat_dir + 'flat.fits')


    # Closed
    fnum = [126, 128, 134, 140, 146, 152, 158, 164, 170, 176, 182, 188, 194]
    fnum += [200, 206, 212, 214, 220, 226, 232, 238, 244, 250, 256, 262, 268]
    img_files = [data_dir + 'obj{0:04d}_c.fits'.format(ii) for ii in fnum]
    reduce_fli.clean_images(img_files, out_dir, rebin=1, sky_frame=sky_dir + 'FLD2_2_sky.fits', flat_frame = flat_dir + 'flat.fits')

    # Closed A
    fnum = [125, 127, 130, 136, 142, 148, 154, 160, 166, 172, 178, 184, 190]
    fnum += [196, 202, 208, 216, 222, 228, 234, 240, 246, 252, 258, 264, 270]
    img_files = [data_dir + 'obj{0:04d}_cA.fits'.format(ii) for ii in fnum]
    reduce_fli.clean_images(img_files, out_dir, rebin=1, sky_frame=sky_dir + 'FLD2_2_sky.fits', flat_frame = flat_dir + 'flat.fits')

    # Closed B
    fnum = [131, 137, 143, 149, 155, 161, 167, 173, 179, 185, 191, 197, 203]
    fnum += [209, 217, 223, 229, 235, 241, 247, 253, 259, 265, 271]
    img_files = [data_dir + 'obj{0:04d}_cB.fits'.format(ii) for ii in fnum]
    reduce_fli.clean_images(img_files, out_dir, rebin=1, sky_frame=sky_dir + 'FLD2_2_sky.fits', flat_frame = flat_dir + 'flat.fits')

    # Closed C
    fnum = [132, 138, 144, 150, 156, 162, 168, 174, 180, 186, 192, 198, 204]
    fnum += [210, 218, 224, 230, 236, 242, 248, 254, 260, 266, 272]
    img_files = [data_dir + 'obj{0:04d}_cC.fits'.format(ii) for ii in fnum]
    reduce_fli.clean_images(img_files, out_dir, rebin=1, sky_frame=sky_dir + 'FLD2_2_sky.fits', flat_frame = flat_dir + 'flat.fits')

    # Closed D
    fnum = [133, 139, 145, 151, 157, 163, 169, 175, 181, 187, 193, 199, 205]
    fnum += [211, 219, 225, 231, 237, 243, 249, 255, 261, 267, 273]
    img_files = [data_dir + 'obj{0:04d}_cD.fits'.format(ii) for ii in fnum]
    reduce_fli.clean_images(img_files, out_dir, rebin=1, sky_frame=sky_dir + 'FLD2_2_sky.fits', flat_frame = flat_dir + 'flat.fits')
    
    
    return
    
def find_stars_FLD2():
    data_dir = root_dir + 'reduce/FLD2_2/'

    # Open Loop
    fnum = [124, 129, 135, 141, 147, 153, 159, 165, 171, 177, 183, 189, 195]
    fnum += [201, 207, 215, 221, 227, 233, 239, 245, 251, 257, 263, 269]
    img_files = [data_dir + 'obj{0:04d}_o_clean.fits'.format(ii) for ii in fnum]
    reduce_fli.find_stars(img_files, fwhm=8, threshold=10)
    
    #Closed 
    
    fnum = [126, 128, 134, 140, 146, 152, 158, 164, 170, 176, 182, 188, 194]
    fnum += [200, 206, 212, 214, 220, 226, 232, 238, 244, 250, 256, 262, 268]
    img_files = [data_dir + 'obj{0:04d}_c_clean.fits'.format(ii) for ii in fnum]
    reduce_fli.find_stars(img_files, fwhm=6, threshold=6)    

    # Closed A
    fnum = [125, 127, 130, 136, 142, 148, 154, 160, 166, 172, 178, 184, 190]
    fnum += [196, 202, 208, 216, 222, 228, 234, 240, 246, 252, 258, 264, 270]
    img_files = [data_dir + 'obj{0:04d}_cA_clean.fits'.format(ii) for ii in fnum]
    reduce_fli.find_stars(img_files, fwhm=6, threshold=6)    

    # Closed B
    fnum = [131, 137, 143, 149, 155, 161, 167, 173, 179, 185, 191, 197, 203]
    fnum += [209, 217, 223, 229, 235, 241, 247, 253, 259, 265, 271]
    img_files = [data_dir + 'obj{0:04d}_cB_clean.fits'.format(ii) for ii in fnum]
    reduce_fli.find_stars(img_files, fwhm=6, threshold=6)   

    # Closed C
    fnum = [132, 138, 144, 150, 156, 162, 168, 174, 180, 186, 192, 198, 204]
    fnum += [210, 218, 224, 230, 236, 242, 248, 254, 260, 266, 272]
    img_files = [data_dir + 'obj{0:04d}_C_clean.fits'.format(ii) for ii in fnum]
    reduce_fli.find_stars(img_files, fwhm=6, threshold=6)   

    # Closed D
    fnum = [133, 139, 145, 151, 157, 163, 169, 175, 181, 187, 193, 199, 205]
    fnum += [211, 219, 225, 231, 237, 243, 249, 255, 261, 267, 273]
    img_files = [data_dir + 'obj{0:04d}_D_clean.fits'.format(ii) for ii in fnum]
    reduce_fli.find_stars(img_files, fwhm=6, threshold=6)   
        
    return


def calc_star_stats():
    data_dir = root_dir + 'reduce/FLD2_2/'
    stats_dir = root_dir +'reduce/stats/'

    
    # Open Loop
    fnum = [124, 129, 135, 141, 147, 153, 159, 165, 171, 177, 183, 189, 195]
    fnum += [201, 207, 215, 221, 227, 233, 239, 245, 251, 257, 263, 269]
    img_files = [data_dir + 'obj{0:04d}_o_clean.fits'.format(ii) for ii in fnum]
    reduce_fli.calc_star_stats(img_files, output_stats=stats_dir + 'stats_open.fits')

    # Closed 
    fnum = [126, 128, 134, 140, 146, 152, 158, 164, 170, 176, 182, 188, 194]
    fnum += [200, 206, 212, 214, 220, 226, 232, 238, 244, 250, 256, 262, 268]
    img_files = [data_dir + "obj{0:04d}_c_clean.fits".format(ii) for ii in fnum]
    reduce_fli.calc_star_stats(img_files, output_stats=stats_dir + 'stats_closedA.fits')

    # Closed A
    fnum = [125, 127, 130, 136, 142, 148, 154, 160, 166, 172, 178, 184, 190]
    fnum += [196, 202, 208, 216, 222, 228, 234, 240, 246, 252, 258, 264, 270]
    img_files = [data_dir + "obj{0:04d}_cA_clean.fits".format(ii) for ii in fnum]
    reduce_fli.calc_star_stats(img_files, output_stats=stats_dir + 'stats_closedA.fits')

    # Closed B
    fnum = [131, 137, 143, 149, 155, 161, 167, 173, 179, 185, 191, 197, 203]
    fnum += [209, 217, 223, 229, 235, 241, 247, 253, 259, 265, 271]
    img_files = [data_dir + "obj{0:04d}_cB_clean.fits".format(ii) for ii in fnum]
    reduce_fli.calc_star_stats(img_files, output_stats=stats_dir + 'stats_closedB.fits')
    
    # Closed C
    fnum = [132, 138, 144, 150, 156, 162, 168, 174, 180, 186, 192, 198, 204]
    fnum += [210, 218, 224, 230, 236, 242, 248, 254, 260, 266, 272]
    img_files = [data_dir + "obj{0:04d}_cC_clean.fits".format(ii) for ii in fnum]
    reduce_fli.calc_star_stats(img_files, output_stats=stats_dir + 'stats_closedC.fits')

    # Closed D
    fnum = [133, 139, 145, 151, 157, 163, 169, 175, 181, 187, 193, 199, 205]
    fnum += [211, 219, 225, 231, 237, 243, 249, 255, 261, 267, 273]
    img_files = [data_dir + "obj{0:04d}_cD_clean.fits".format(ii) for ii in fnum]
    reduce_fli.calc_star_stats(img_files, output_stats=stats_dir + 'stats_closedD.fits')
    
    return 

def stack_FLD2():

    data_dir = root_dir + 'reduce/FLD2_2/'
    stats_dir = root_dir + 'reduce/stats/'
    stacks_dir = root_dir + 'reduce/stacks/'

    util.mkdir(stacks_dir)

    # Open Loop
    fnum = [124, 129, 135, 141, 147, 153, 159, 165, 171, 177, 183, 189, 195]
    fnum += [201, 207, 215, 221, 227, 233, 239, 245, 251, 257, 263, 269]
    open_images = [data_dir + 'obj{0:04d}_o_clean.fits'.format(ii) for ii in fnum]
    open_starlists = [data_dir + 'obj{0:04d}_o_clean_stars.txt'.format(ii) for ii in fnum]
    open_output_root = stacks_dir + 'FLD2_stack_open'
    reduce_fli.shift_and_add(open_images, open_starlists, open_output_root, method='mean')
    
    # Closed 
    fnum = [126, 128, 134, 140, 146, 152, 158, 164, 170, 176, 182, 188, 194]
    fnum += [200, 206, 212, 214, 220, 226, 232, 238, 244, 250, 256, 262, 268]
    closed_images = [data_dir + 'obj{0:04d}_c_clean.fits'.format(ii) for ii in fnum]
    closed_starlists = [data_dir + 'obj{0:04d}_c_clean_stars.txt'.format(ii) for ii in fnum]
    closed_output_root = stacks_dir + 'FLD2_stack_closed'
    reduce_fli.shift_and_add(closed_images, closed_starlists, closed_output_root, method='mean')
    
    # Closed A
    fnum = [125, 127, 130, 136, 142, 148, 154, 160, 166, 172, 178, 184, 190]
    fnum += [196, 202, 208, 216, 222, 228, 234, 240, 246, 252, 258, 264, 270]
    closed_images = [data_dir + 'obj{0:04d}_cA_clean.fits'.format(ii) for ii in fnum]
    closed_starlists = [data_dir + 'obj{0:04d}_cA_clean_stars.txt'.format(ii) for ii in fnum]
    closed_output_root = stacks_dir + 'FLD2_stack_closedA'
    reduce_fli.shift_and_add(closed_images, closed_starlists, closed_output_root, method='mean')

    # Closed B
    fnum = [131, 137, 143, 149, 155, 161, 167, 173, 179, 185, 191, 197, 203]
    fnum += [209, 217, 223, 229, 235, 241, 247, 253, 259, 265, 271]
    closed_images = [data_dir + 'obj{0:04d}_cB_clean.fits'.format(ii) for ii in fnum]
    closed_starlists = [data_dir + 'obj{0:04d}_cB_clean_stars.txt'.format(ii) for ii in fnum]
    closed_output_root = stacks_dir + 'FLD2_stack_closedB'
    reduce_fli.shift_and_add(closed_images, closed_starlists, closed_output_root, method='mean')

    # Closed C
    fnum = [132, 138, 144, 150, 156, 162, 168, 174, 180, 186, 192, 198, 204]
    fnum += [210, 218, 224, 230, 236, 242, 248, 254, 260, 266, 272]
    closed_images = [data_dir + 'obj{0:04d}_cC_clean.fits'.format(ii) for ii in fnum]
    closed_starlists = [data_dir + 'obj{0:04d}_cC_clean_stars.txt'.format(ii) for ii in fnum]
    closed_output_root = stacks_dir + 'FLD2_stack_closedC'
    reduce_fli.shift_and_add(closed_images, closed_starlists, closed_output_root, method='mean')

    # Closed D
    fnum = [133, 139, 145, 151, 157, 163, 169, 175, 181, 187, 193, 199, 205]
    fnum += [211, 219, 225, 231, 237, 243, 249, 255, 261, 267, 273]
    closed_images = [data_dir + 'obj{0:04d}_cD_clean.fits'.format(ii) for ii in fnum]
    closed_starlists = [data_dir + 'obj{0:04d}_cD_clean_stars.txt'.format(ii) for ii in fnum]
    closed_output_root = stacks_dir + 'FLD2_stack_closedD'
    reduce_fli.shift_and_add(closed_images, closed_starlists, closed_output_root, method='mean')
    
    return


def analyze_stacks():
    data_dir = root_dir + 'reduce/stacks/'
    stats_dir = root_dir + 'reduce/stats/'
    
    img_files = [data_dir + 'FLD2_stack_open.fits',
                 data_dir + 'FLD2_stack_closed.fits',
                 data_dir + 'FLD2_stack_closedA.fits',
                 data_dir + 'FLD2_stack_closedB.fits',
                 data_dir + 'FLD2_stack_closedC.fits',
                 data_dir + 'FLD2_stack_closedD.fits']
    
    #Find stars in image
    reduce_fli.find_stars(img_files, fwhm=5, threshold=10, N_passes=2, plot_psf_compare=False)
    
    # Calc stats on all the stacked images
    reduce_fli.calc_star_stats(img_files, output_stats= stats_dir + 'stats_stacks.fits')

    return
    
