
from astropy.io import fits
from os import listdir
import numpy as np
import glob 
from skimage.measure import block_reduce
from astropy.stats import sigma_clipped_stats



def rebin(a, bin_fac):
    
    #bins an array 'a' by a factor 'bin_fac'; sums pixels

    img_bin = block_reduce(a, (bin_fac, bin_fac), func=np.sum)
    return img_bin



def make_flat(flat_files, dark_files):

    ##Makes a flat from a list of twilights; specifically for the 2016-11-18 run's weird naming scheme
    
    flats = []
    
    for twilight in flat_files:
        segment = twilight.split("_")[1].split(".")[0] #cuts off file naming weirdness, leaves xxxx.fits (so I can search it in the darks)
    
        twi, twi_head = fits.getdata(twilight, header=True)
        med = np.median(twi, axis=None)
        if med < 20000: #checks for saturated images   
            for file in dark_files: #finds corresponding dark
                if segment in file:
                    dark, dark_head = fits.getdata(file, header=True)
                    if twi_head["EXPTIME"] == dark_head["EXPTIME"] and twi_head["BINFAC"] == dark_head["BINFAC"]:
                        twi_ds = twi - dark
                        mean, median, stdev = sigma_clipped_stats(twi_ds, sigma=3, axis=None)
                        twi_ds_norm = twi_ds / median #normalization; outlier rejection here later iterative mean with outlier rejection; find it in imaka analysis file)
                        flats.append(twi_ds_norm) 
    
    ave, master_flat, stdev = sigma_clipped_stats(flats, sigma=3, axis=0) #Combine all flats
    return master_flat.data



def make_dark(dark_files, bin_fac):

    #combine dark images of some bin factor (images are either binned 1 or 3)

    darks = []
    for file in dark_files:
        data, header = fits.getdata(file, header=True)
        if header['EXPTIME'] == 45 and header['BINFAC'] == bin_fac:
            darks.append(data)
    dark = np.median(darks, axis=0)
    return dark



def make_sky(sky_files, dark_array):
    
    #make sky from a list of sky files and a np.array dark

    skies = []
    for file in sky_files:
        data = fits.getdata(file)
        skies.append(data-dark_array) #subtracts appropriate dark

    ave, sky, stdev = sigma_clipped_stats(skies, sigma=3, iters=2, axis=0)
    return sky.data



def clean_image(object_file, dark_array, flat_array, sky_array):
    
    #combine image with dark, sky, and flat
    
    obj_data, obj_head = fits.getdata(object_file, header=True)
    
    if np.shape(obj_data) != np.shape(flat_array):
        flat = rebin(flat_array[:, 0:-1], np.shape(obj_data))
   
    final = (obj_data - dark_array - sky_array) / flat
    return final
