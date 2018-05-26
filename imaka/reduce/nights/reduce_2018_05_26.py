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
from imaka.analysis import moffat
import os, shutil
import pdb

root_dir = '//Volumes/DATA5/imaka/20180526/FLI/'

sky_dir = root_dir + 'reduce/sky/' 
data_dir = root_dir + 'FLD2/'
flat_dir = root_dir + 'reduce/calib/'
out_dir = root_dir + 'reduce/FLD2/'
stats_dir = root_dir +'reduce/stats/'
stacks_dir = root_dir + 'reduce/stacks/'
    
fnum_o = [] #open loop img file numbers
fnum_c = [] #closed loop img file numbers


    
#def make_flat():    

#def make_sky():

 
def reduce_FLD2():
    
    util.mkdir(out_dir)

    # Open Loop
    img_files = [data_dir + 'obj{0:04d}_o.fits'.format(ii) for ii in fnum_o]
    reduce_fli.clean_images(img_files, out_dir, rebin=1, sky_frame=data_dir+'../twilight/sky03.fits', flat_frame=None)

    # Closed
    img_files = [data_dir + 'obj{0:04d}_fourWFS_c.fits'.format(ii) for ii in fnum_c]
    reduce_fli.clean_images(img_files, out_dir, rebin=1, sky_frame=data_dir + '../twilight/sky03.fits', flat_frame =None)

    return


def find_stars_FLD2():

    # Open Loop
    img_files = [out_dir + 'obj{0:04d}_o_clean.fits'.format(ii) for ii in fnum_o]
    reduce_fli.find_stars(img_files, fwhm=8, threshold=10)
    
    #Closed Loop
    img_files = [out_dir + 'obj{0:04d}_fourWFS_c_clean.fits'.format(ii) for ii in fnum_c]
    reduce_fli.find_stars(img_files, fwhm=6, threshold=6)    
        
    return


def calc_star_stats():
    
    # Open Loop
    img_files = [data_dir + 'obj{0:04d}_o_clean.fits'.format(ii) for ii in fnum_o]
    reduce_fli.calc_star_stats(img_files, output_stats=stats_dir + 'stats_open.fits')

    # Closed Loop
    img_files = [data_dir + "obj{0:04d}_c_clean.fits".format(ii) for ii in fnum_C]
    reduce_fli.calc_star_stats(img_files, output_stats=stats_dir + 'stats_closed.fits')

    return


def calc_mof_stats():

    # Open Loop
    stats_file = stats_dir + 'stats_open.fits'
    img_files = [data_dir + 'obj{0:04d}_o_clean.fits'.format(ii) for ii in fnum_o]
    moffat.fit_moffat(img_files, stats_file)

    # Closed Loop
    stats_file = stats_dir + 'stats_closed.fits'
    img_files = [data_dir + 'obj{0:04d}_c_clean.fits'.format(ii) for ii in fnum_c]
    moffat.fit_moffat(img_files, stats_file)


def stack_FLD2():

    util.mkdir(stacks_dir)

    # Open Loop
    open_images = [data_dir + 'obj{0:04d}_o_clean.fits'.format(ii) for ii in fnum_o]
    open_starlists = [data_dir + 'obj{0:04d}_o_clean_stars.txt'.format(ii) for ii in fnum_o]
    open_output_root = stacks_dir + 'Beehive-W_stack_open'
    reduce_fli.shift_and_add(open_images, open_starlists, open_output_root, method='mean')
    
    # Closed Loop
    closed_images = [data_dir + 'obj{0:04d}_c_clean.fits'.format(ii) for ii in fnum_c]
    closed_starlists = [data_dir + 'obj{0:04d}_c_clean_stars.txt'.format(ii) for ii in fnum_c]
    closed_output_root = stacks_dir + 'Beehive-w_stack_closed'
    reduce_fli.shift_and_add(closed_images, closed_starlists, closed_output_root, method='mean')
    
    return


#def analyze_stacks():

#    img_files = []
    
    #Find stars in image
#    reduce_fli.find_stars(img_files, fwhm=5, threshold=10, N_passes=2, plot_psf_compare=False)
    
    # Calc stats on all the stacked images
#    reduce_fli.calc_star_stats(img_files, output_stats= stats_dir + 'stats_stacks.fits')

    return
    
