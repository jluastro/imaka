import os
import math
import pylab as plt
import numpy as np
from astropy.io import fits
from astropy import table
from astropy import units
import pdb
#from astroscrappy import detect_cosmics
import glob
import photutils
from photutils import psf
from photutils import morphology as morph
from photutils import DAOStarFinder
from photutils import CircularAperture
from photutils import CircularAnnulus
from photutils import aperture_photometry
from astropy.stats import gaussian_fwhm_to_sigma
from astropy.stats import sigma_clipped_stats
from astropy.modeling import models, fitting
from flystar import match
from flystar import align
from flystar import transforms
from imaka.reduce import calib
from imaka.reduce import util
#import ccdproc
from scipy.ndimage import interpolation
from scipy.ndimage import median_filter
import scipy.ndimage
from astropy.table import Table
from skimage.measure import block_reduce
from datetime import datetime
import pytz
import multiprocessing as mp
from itertools import repeat
import warnings
from astropy.utils.exceptions import AstropyWarning

np.seterr(all='ignore')

def rebin(a, bin_fac):
    """
    Bins an array 'a' by a factor 'bin_fac'; sums pixels
    """
    img_bin = block_reduce(a, (bin_fac, bin_fac), func=np.sum)
    
    return img_bin

def write_rebin(file, binfac):
    #Rewrites file binned by some factor 
    dat, hdr = fits.getdata(file, header=True)
    rebin_dat = rebin(dat, binfac)
    fits.writeto(file.split(".fits")[0]+"orig.fits", dat, hdr, overwrite=False)
    orig_binfac = hdr['BINFAC']
    new_binfac = orig_binfac * binfac
    hdr['BINFAC'] = new_binfac
    fits.writeto(file, rebin_dat, hdr, overwrite=True)
    return


def subtract_dark(image_file, dark_file):
    # Dark subtract an image.
    img, i_hdr = fits.getdata(image_file, header=True)
    drk, d_hdr = fits.getdata(dark_file, header=True)

    # Rescale the dark to match the integration time... this isn't perfect.
    drk_time = d_hdr['EXPTIME']
    img_time = i_hdr['EXPTIME']

    drk_scale = img_time / drk_time
    img_ds = img - (drk * drk_scale)

    return img_ds


def rescale_dark(dark_file, image_file, dark_save):
    # Rescales a dark to given images exposure time
    img, i_hdr = fits.getdata(image_file, header=True)
    drk, d_hdr = fits.getdata(dark_file, header=True)

    # Rescale the dark to match the integration time... this isn't perfect.
    drk_time = d_hdr['EXPTIME']
    img_time = i_hdr['EXPTIME']

    drk_scale = img_time / drk_time
    img_ds = drk * drk_scale
    
    d_hdr['EXPTIME'] = img_time
    
    fits.writeto(dark_save, img_ds, d_hdr)
    return
    
    


def fix_bad_pixels(img):
    crmask, cleanarr = detect_cosmics(img, sigclip=5, sigfrac=0.3, objlim=5.0, gain=1.0,
                                      readnoise=6.5, satlevel=65535.0, pssl=0, niter=4,
                                      sepmed=True, cleantype='meanmask', fsmode='median',
                                      psfmodel='gauss', psffwhm=2.5, psfsize=7,
                                      psfk=None, psfbeta=4.765, verbose=False)
                       
    return cleanarr


def clean_images(img_files, out_dir, rebin=1, sky_frame=None, flat_frame=None,
                     fix_bad_pixels=False, worry_about_edges=False):
    """
    Clean a stack of imaka images. Options include:
      - sky subtraction
      - flat fielding
      - bad pixel fixing
      - rebinning
    
    """
    print('\nREDUCE_FLI: clean_images()')

    # If we are doing sky subtraction with a separate sky, then load it up
    if sky_frame != None:
        sky = fits.getdata(sky_frame)

        if rebin != 1:
            sky = block_reduce(sky, (rebin, rebin), func=np.sum)
    else:
        sky = sky_frame

    # If we are doing flat fielding, then load it up
    if flat_frame != None:
        flat = fits.getdata(flat_frame)
    else:
        flat = flat_frame

    ####
    # Setup for parallel processing.
    ####
    cpu_count = _cpu_num()
    # Start the pool
    #pool = mp.Pool(cpu_count)
    pool = mp.Pool(1) ## DEBUG
    
    #####
    # Loop through the image stack
    #####
    # Run
    print(f'Cleaning images in parallel with {cpu_count} cores.')
    args = zip(img_files, repeat(sky), repeat(flat), repeat(out_dir),
                   repeat(rebin), repeat(fix_bad_pixels), repeat(worry_about_edges))
    pool.starmap(clean_and_save_single_image, args)

    pool.close()
    pool.join()

    return

def clean_and_save_single_image(img_file, sky, flat, out_dir, rebin, fix_bad_pixels, worry_about_edges):
    """
    Clean and save an image file. 
    """
    pid = mp.current_process().pid
    
    print(f'  p{pid} - Working on image: {img_file}')

    # Read image in
    img, hdr = fits.getdata(img_file, header=True)
    
    # Clean it
    print(f'  p{pid} - Cleaning image')
    img_clean, bad_pixels = clean_single_image(img, sky=sky, flat=flat, rebin=rebin,
                                                fix_bad_pixels=fix_bad_pixels,
                                                worry_about_edges=worry_about_edges)
    

    # Write it out
    print(f'  p{pid} - Writing image')
    file_dir, file_name = os.path.split(img_file)
    out_name = file_name.replace('.fits', '_clean.fits')
    fits.writeto(out_dir + out_name, img_clean, hdr, overwrite=True)

    if bad_pixels != None:
        mask_name = file_name.replace('.fits', '_mask.fits')
        fits.writeto(out_dir + mask_name, bad_pixels, hdr, overwrite=True)
    
    print(f'  p{pid} - Image stored: {out_dir}/{out_name}')

    return


def clean_single_image(img, sky=None, flat=None, rebin=1,
                           fix_bad_pixels=False, worry_about_edges=False, verbose=False):
    pid = mp.current_process().pid    

    # Bin by a factor of 10 (or as specified).
    if rebin != 1:
        print(f'  p{pid} Bin by a factor of {rebin}')
        img = block_reduce(img, (rebin, rebin), func=np.sum)

    # Subtract background.
    if sky is None:
        bkg_box_size = int(round(img.shape[0] / 10.))
        print(f'  p{pid} Subtract background with smoothing box size of {bkg_box_size}')
        blurred = median_filter(img, size=bkg_box_size)
        img = img - blurred.astype(float)
    else:
        img = img - sky

        # Calculate the mean background and make sure it isn't too negative.
        # This happens when the skies are taken with moonlight or twilight.
        mean, median, stddev = sigma_clipped_stats(img.flatten(), sigma_lower=5, sigma_upper=3, maxiters=5)

        if median < 0:
            img -= median
            
    # Flat field.
    if flat is not None:
        print(f'  p{pid} Dividing by flat')
        img = img / flat

    # Fix any inifinte or nan pixels... just set to 0
    idx = np.where(np.isfinite(img) == False)
    img[idx] = 0

    # Fix hot and dead pixels
    if fix_bad_pixels:
        print(f'  p{pid} Fixing bad pixels')
        bad_pixels, img = find_outlier_pixels(img, tolerance=5, worry_about_edges=worry_about_edges)
    else:
        bad_pixels = None
            
    return img, bad_pixels


def find_stars(img_files, fwhm=5, threshold=4, N_passes=2, plot_psf_compare=False, mask_file=None):
    """
    img_files - a list of image files.
    fwhm - First guess at the FWHM of the PSF for the first pass on the first image.
    threshold - the SNR threshold (mean + threshold*std) above which to search for sources.
    N_passes - how many times to find sources, recalc the PSF FWHM, and find sources again. 
    """
    print('\nREDUCE_FLI: find_stars()')
    
    # Prepare for parallel compute.
        
    cpu_count = _cpu_num()
    # Start the pool
    pool = mp.Pool(cpu_count)
    #pool = mp.Pool(1) ## DEBUG

    # Normalize each of the flats: in parallel
    print(f'  Finding stars in parallel with {cpu_count} cores.')
    args = zip(img_files, repeat(fwhm), repeat(threshold), repeat(N_passes), repeat(plot_psf_compare), repeat(mask_file))
    results = pool.starmap(find_stars_single, args)

    pool.close()
    pool.join()

    return

def find_stars_single(img_file, fwhm, threshold, N_passes, plot_psf_compare, mask_file):
    pid = mp.current_process().pid
    
    print(f'  p{pid} - Working on image: {img_file}')
    img, hdr = fits.getdata(img_file, header=True, ignore_missing_end=True)

    mask = fits.getdata(mask_file).astype('bool')
    
    img = np.ma.masked_array(img, mask=mask)
    
    fwhm_curr = fwhm
    
    # Calculate the bacgkround and noise (iteratively)
    print(f'  p{pid} - Calculating background')
    bkg_threshold_above = 1
    bkg_threshold_below = 3

    good_pix = np.where(np.isfinite(img))
    
    for nn in range(5):
        bkg_mean = img[good_pix].mean()
        bkg_std = img[good_pix].std()

        bad_hi = bkg_mean + (bkg_threshold_above * bkg_std)
        bad_lo = bkg_mean - (bkg_threshold_below * bkg_std)

        good_pix = np.where((img > bad_lo) & (img < bad_hi))
    
    bkg_mean = img[good_pix].mean()
    bkg_std = img[good_pix].std()
    img_threshold = threshold * bkg_std 
    print(f'    p{pid} - Bkg = {bkg_mean:.2f} +/- {bkg_std:.2f}')
    print(f'    p{pid} - Bkg Threshold = {img_threshold:.2f}')

    # Detect stars
    print(f'     p{pid} - Detecting Stars')

    # Each pass will have an updated fwhm for the PSF.
    for nn in range(N_passes):
        print(f'    p{pid} - Pass {nn:d} assuming FWHM = {fwhm_curr:.1f}')
        daofind = DAOStarFinder(fwhm=fwhm_curr, threshold = img_threshold, exclude_border=True)
        sources = daofind(img - bkg_mean, mask=mask)
        print(f'    p{pid} - {len(sources)} sources found, now fitting for FWHM.')

        # Calculate FWHM for each detected star.
        x_fwhm = np.zeros(len(sources), dtype=float)
        y_fwhm = np.zeros(len(sources), dtype=float)
        theta = np.zeros(len(sources), dtype=float)
    
        # We will actually be resampling the images for the Gaussian fits.
        resamp = 2
        
        cutout_half_size = int(round(fwhm_curr * 3.5))
        cutout_size = 2 * cutout_half_size + 1

        # Define variables to hold final averages PSFs
        final_psf_obs = np.zeros((cutout_size*resamp, cutout_size*resamp), dtype=float)
        final_psf_mod = np.zeros((cutout_size*resamp, cutout_size*resamp), dtype=float)
        final_psf_count = 0

        # Setup our gaussian fitter with some good initial guesses.
        sigma_init_guess = fwhm_curr * gaussian_fwhm_to_sigma
        g2d_model = models.Gaussian2D(1.0, cutout_half_size*resamp, cutout_half_size*resamp,
                                          sigma_init_guess*resamp, sigma_init_guess*resamp, theta=0,
                                          bounds={'x_stddev':[0, 20], 'y_stddev':[0, 20], 'amplitude':[0, 2]})
        c2d_model = models.Const2D(0.0)
        
        the_model = g2d_model + c2d_model
        the_fitter = fitting.LevMarLSQFitter()
        
        cut_y, cut_x = np.mgrid[:cutout_size, :cutout_size]
        
        for ss in range(len(sources)):
            x_lo = int(round(sources[ss]['xcentroid'] - cutout_half_size))
            x_hi = x_lo + cutout_size
            y_lo = int(round(sources[ss]['ycentroid'] - cutout_half_size))
            y_hi = y_lo + cutout_size

            cutout_tmp = img[y_lo:y_hi, x_lo:x_hi].astype(float)
            if ((cutout_tmp.shape[0] != cutout_size) | (cutout_tmp.shape[1] != cutout_size)):
                # Edge source... fitting is no good
                continue
        
            # Oversample the image
            cutout_resamp = scipy.ndimage.zoom(cutout_tmp, resamp, order = 1)
            cutout_resamp /= cutout_resamp.sum()
            cut_y_resamp, cut_x_resamp = np.mgrid[:cutout_size*resamp, :cutout_size*resamp]

            # Fit a 2D gaussian + constant
            with warnings.catch_warnings():
                # Suppress warnings... too many.
                warnings.simplefilter("ignore", category=UserWarning)
                warnings.simplefilter("ignore", category=AstropyWarning)
                g2d_params = the_fitter(the_model, cut_x_resamp, cut_y_resamp, cutout_resamp)
                
            g2d_image = g2d_params(cut_x_resamp, cut_y_resamp)

            # Catch bad fits and ignore. 
            if (np.isnan(g2d_params.x_mean_0.value) or
                (np.abs(g2d_params.x_mean_0.value) > (cutout_size * resamp)) or
                (np.abs(g2d_params.y_mean_0.value) > (cutout_size * resamp))):
                print(f'      p{pid} - Bad fit for {ss}')
                continue

            # Add to our average observed/model PSFs
            if sources['flux'][ss] > 1.9:
                final_psf_count += 1
                final_psf_obs += cutout_resamp
                final_psf_mod += g2d_image

            if (plot_psf_compare == True) and (x_lo > 200) and (y_lo > 200):
                plt.figure(4, figsize=(6, 4))
                plt.clf()
                plt.subplots_adjust(left=0.08)
                plt.imshow(cutout_resamp)
                plt.colorbar(fraction=0.25)
                plt.axis('equal')
                plt.title(f'Image (resamp={resamp:d})')
                plt.pause(0.05)
            
                plt.figure(5, figsize=(6, 4))
                plt.clf()
                plt.subplots_adjust(left=0.08)
                plt.imshow(g2d_image)
                plt.colorbar(fraction=0.25)
                plt.axis('equal')
                plt.title(f'Model (resamp={resamp:d})')
                plt.pause(0.05)
            
                plt.figure(6, figsize=(6, 4))
                plt.clf()
                plt.subplots_adjust(left=0.08)
                plt.imshow((cutout_resamp - g2d_image) / cutout_resamp)
                cbar = plt.colorbar(fraction=0.25)
                plt.axis('equal')
                cbar.set_label('% Residual')
                plt.pause(0.05)
                
                #pdb.set_trace()

            # Save the FWHM and angle.
            x_fwhm[ss] = g2d_params.x_stddev_0.value / gaussian_fwhm_to_sigma / resamp
            y_fwhm[ss] = g2d_params.y_stddev_0.value / gaussian_fwhm_to_sigma / resamp
            theta[ss] = g2d_params.theta_0.value
                
            # Some occasional display
            if (plot_psf_compare == True) and (ss % 250 == 0):
                plt.figure(2, figsize=(6, 4))
                plt.clf()
                plt.subplots_adjust(left=0.08)
                plt.imshow(final_psf_obs)
                plt.colorbar(fraction=0.25)
                plt.axis('equal')
                plt.title(f'Obs PSF (resamp = {resamp:d})')
                plt.pause(0.05)
                
                plt.figure(3, figsize=(6, 4))
                plt.clf()
                plt.subplots_adjust(left=0.08)
                plt.imshow(final_psf_mod)
                plt.colorbar(fraction=0.25)
                plt.axis('equal')
                plt.title(f'Mod PSF (resamp = {resamp:d})')
                plt.pause(0.05)

                print(f'    p{pid} - ss={ss} fwhm_x={x_fwhm[ss]:.1f} fwhm_y={y_fwhm[ss]:.1f}')

                

        sources['x_fwhm'] = x_fwhm
        sources['y_fwhm'] = y_fwhm
        sources['theta'] = theta

        # Save the average PSF (flux-weighted). Note we are making a slight mistake
        # here since each PSF has a different sub-pixel position... still same for both
        # obs and model.
        final_psf_obs /= final_psf_count
        final_psf_mod /= final_psf_count
        final_psf_obs /= final_psf_obs.sum()
        final_psf_mod /= final_psf_mod.sum()
        fits.writeto(img_file.replace('.fits', '_psf_obs.fits'), final_psf_obs, hdr, overwrite=True)
        fits.writeto(img_file.replace('.fits', '_psf_mod.fits'), final_psf_mod, hdr, overwrite=True)

        # Drop sources with flux (signifiance) that isn't good enough.
        # Empirically this is <1.2
        # Also drop sources that couldn't be fit.
        good = np.where((sources['flux'] > 1.9) & (sources['x_fwhm'] > 0) & (sources['y_fwhm'] > 0))[0]
        sources = sources[good]

        # Only use the brightest sources for calculating the mean. This is just for printing.
        idx = np.where(sources['flux'] > 5)[0]
        x_fwhm_med = np.median(sources['x_fwhm'][idx])
        y_fwhm_med = np.median(sources['y_fwhm'][idx])
        
        print(f'      p{pid} - Number of sources = {len(sources)}')
        print(f'      p{pid} - Median x_fwhm = {x_fwhm_med:.1f} +/- {sources["x_fwhm"].std():.1f}')
        print(f'      p{pid} - Median y_fwhm = {y_fwhm_med:.1f} +/- {sources["y_fwhm"].std():.1f}')

        fwhm_curr = np.mean([x_fwhm_med, y_fwhm_med])

        formats = {'xcentroid': '%8.3f', 'ycentroid': '%8.3f', 'sharpness': '%.2f',
                   'roundness1': '%.2f', 'roundness2': '%.2f', 'peak': '%10.1f',
                   'flux': '%10.6f', 'mag': '%6.2f', 'x_fwhm': '%5.2f', 'y_fwhm': '%5.2f',
                   'theta': '%6.3f'}
    
        sources.write(img_file.replace('.fits', '_stars.txt'), format='ascii.fixed_width',
                          delimiter=None, bookend=False, formats=formats, overwrite=True)

    
    return
    
def find_outlier_pixels(data, tolerance=3, worry_about_edges=True, median_filter_size=10):
    """Find the hot or dead pixels in a 2D dataset. 
    Parameters
    ----------
    data: numpy array (2D)
        image
    tolerance: float (default 3)
        number of standard deviations used to cutoff the hot pixels
    worry_about_edges: boolean (default True)
        If you want to ignore the edges and greatly speed up the code, then set
        worry_about_edges to False.
        
    Return
    ------
    The function returns a list of hot pixels and also an image with with hot pixels removed
    """
    from scipy.ndimage import median_filter
    blurred = median_filter(data, size=median_filter_size)
    difference = data - blurred.astype(float)
    threshold = tolerance * np.std(difference)

    # Find the hot pixels, but ignore the edges
    hot_pixels = np.nonzero((np.abs(difference[1:-1,1:-1]) > threshold) )
    hot_pixels = np.array(hot_pixels) + 1 # because we ignored the first row and first column

    fixed_image = np.copy(data) # This is the image with the hot pixels removed
    for y, x in zip(hot_pixels[0], hot_pixels[1]):
        fixed_image[y, x] = blurred[y, x]

    if worry_about_edges == True:
        height, width = np.shape(data)

        ###Now get the pixels on the edges (but not the corners)###

        #left and right sides
        for index in range(1,height-1):
            #left side:
            med  = np.median(data[index-1:index+2,0:2])
            diff = np.abs(data[index,0] - med)
            if diff>threshold: 
                hot_pixels = np.hstack(( hot_pixels, [[index],[0]]  ))
                fixed_image[index,0] = med

            #right side:
            med  = np.median(data[index-1:index+2,-2:])
            diff = np.abs(data[index,-1] - med)
            if diff>threshold: 
                hot_pixels = np.hstack(( hot_pixels, [[index],[width-1]]  ))
                fixed_image[index,-1] = med

        #Then the top and bottom
        for index in range(1,width-1):
            #bottom:
            med  = np.median(data[0:2,index-1:index+2])
            diff = np.abs(data[0,index] - med)
            if diff>threshold: 
                hot_pixels = np.hstack(( hot_pixels, [[0],[index]]  ))
                fixed_image[0,index] = med

            #top:
            med  = np.median(data[-2:,index-1:index+2])
            diff = np.abs(data[-1,index] - med)
            if diff>threshold: 
                hot_pixels = np.hstack(( hot_pixels, [[height-1],[index]]  ))
                fixed_image[-1,index] = med

        ###Then the corners###

        #bottom left
        med  = np.median(data[0:2,0:2])
        diff = np.abs(data[0,0] - med)
        if diff>threshold: 
            hot_pixels = np.hstack(( hot_pixels, [[0],[0]]  ))
            fixed_image[0,0] = med

        #bottom right
        med  = np.median(data[0:2,-2:])
        diff = np.abs(data[0,-1] - med)
        if diff>threshold: 
            hot_pixels = np.hstack(( hot_pixels, [[0],[width-1]]  ))
            fixed_image[0,-1] = med

        #top left
        med  = np.median(data[-2:,0:2])
        diff = np.abs(data[-1,0] - med)
        if diff>threshold: 
            hot_pixels = np.hstack(( hot_pixels, [[height-1],[0]]  ))
            fixed_image[-1,0] = med

        #top right
        med  = np.median(data[-2:,-2:])
        diff = np.abs(data[-1,-1] - med)
        if diff>threshold: 
            hot_pixels = np.hstack(( hot_pixels, [[height-1],[width-1]]  ))
            fixed_image[-1,-1] = med

    return hot_pixels, fixed_image


def calc_star_stats(img_files, output_stats='image_stats.fits', filt=None, fourfilt=False, starlists=None):
    """
    Calculate statistics for the Data Metrics table.
    """
    # Create arrays for all the final statistics.
    N_files = len(img_files)
    s_time_hst = np.zeros(N_files, dtype='S15')
    s_time_utc = np.zeros(N_files, dtype='S15')
    s_date_hst = np.zeros(N_files, dtype='S15')
    s_date_utc = np.zeros(N_files, dtype='S15')
    s_band = np.zeros(N_files, dtype='S15')
    s_binfac = np.zeros(N_files, dtype=float)
    s_ee25 = np.zeros(N_files, dtype=float)
    s_ee50 = np.zeros(N_files, dtype=float)
    s_ee80 = np.zeros(N_files, dtype=float)
    s_xfwhm = np.zeros(N_files, dtype=float)
    s_yfwhm = np.zeros(N_files, dtype=float)
    s_theta = np.zeros(N_files, dtype=float)
    s_fwhm = np.zeros(N_files, dtype=float)
    s_fwhm_std = np.zeros(N_files, dtype=float)
    s_NEA = np.zeros(N_files, dtype=float)
    s_NEA2 = np.zeros(N_files, dtype=float)
    s_emp_fwhm = np.zeros(N_files, dtype=float)
    s_emp_fwhm_std = np.zeros(N_files, dtype=float)
    # Make a column to designate filter position in four filter images
    s_quadrant = np.zeros(N_files, dtype='S15')

    ####
    # Setup for parallel processing.
    ####
    cpu_count = _cpu_num()
    
    # Start the pool
    pool = mp.Pool(cpu_count)
    
    results_async = []
    print(f'calc_stats in parallel with {cpu_count} cores.')
    
    # making file names
    for ii in range(N_files):
        # Select the image file and starlist to work on
        img_file = img_files[ii]

        if starlists == None:
            if filt==None:
                starlist = img_file.replace('.fits', '_stars.txt')
            else:
                starlist = img_file.replace('.fits', '_'+filt+'_stars.txt')
        elif starlists != None:
            starlist = starlists[ii]

    
        #####
        # Add calc for this starlist to the pool.
        #####
        results = pool.apply_async(calc_star_stats_single, (img_file, starlist, fourfilt))
        results_async.append(results)

    pool.close()
    pool.join()

    # Fetch the plate scale
    img, hdr = fits.getdata(img_files[0], header=True)
    plate_scale = util.get_plate_scale(img, hdr)
        

    for ii in range(N_files):
        #pdb.set_trace()
        results = results_async[ii].get()
        
        # Save results
        s_band[ii] = results['band']
        s_binfac[ii] = results['binfac']
        s_time_utc[ii] = results['time_utc']
        s_date_utc[ii] = results['date_utc']
        s_time_hst[ii] = results['time_hst']
        s_date_hst[ii] = results['date_hst']
        s_fwhm[ii] = results['fwhm']
        s_fwhm_std[ii] = results['fwhm_std']
        s_ee25[ii] = results['ee25']
        s_ee50[ii] = results['ee50']
        s_ee80[ii] = results['ee80']
        s_NEA[ii] = results['NEA']
        s_NEA2[ii] = results['NEA2']
        s_xfwhm[ii] = results['xfwhm']
        s_yfwhm[ii] = results['yfwhm']
        s_theta[ii] = results['theta']
        s_emp_fwhm[ii] = results['emp_fwhm']
        s_emp_fwhm_std[ii] = results['emp_fwhm_std']
        s_quadrant[ii] = results['quadrant']

    
    # FUTURE: Make a date array that holds UTC.

    stats = table.Table([img_files, s_band, s_binfac, s_date_utc, s_time_utc, s_date_hst, s_time_hst,
                             s_fwhm, s_fwhm_std, s_ee25, s_ee50, s_ee80,
                             s_NEA, s_NEA2, s_xfwhm, s_yfwhm, s_theta, s_emp_fwhm, s_emp_fwhm_std,
                             s_quadrant],
                             names=('Image', 'FILTER', 'BINFAC', 'DATE_UTC', 'TIME_UTC', 'DATE_HST', 'TIME_HST',
                                        'FWHM', 'FWHM_std', 'EE25', 'EE50', 'EE80',
                                        'NEA', 'NEA2', 'xFWHM', 'yFWHM', 'theta', 'emp_fwhm', 'emp_fwhm_std',
                                        'quadrant'),
                            meta={'name':'Stats Table', 'scale': plate_scale})
    
    ## I believe this is already covered in previous command
    #if fourfilt:
    #    quad_col = stats.Column(quadrant, name='quadrant')
    #    stats.add_column(quad_col)
    
    stats['FWHM'].format = '7.3f'
    stats['FWHM_std'].format = '7.3f'
    stats['EE25'].format = '7.3f'
    stats['EE50'].format = '7.3f'
    stats['EE80'].format = '7.3f'
    stats['NEA'].format = '7.3f'
    stats['NEA2'].format = '7.3f'
    stats['xFWHM'].format = '7.3f'
    stats['yFWHM'].format = '7.3f'
    stats['theta'].format = '7.3f'
    stats['emp_fwhm'].format = '7.3f'
    stats['emp_fwhm_std'].format = '7.3f'

    add_frame_number_column(stats)

    stats.write(output_stats, overwrite=True)
                        
    return

def calc_star_stats_single(img_file, starlist, is_four_filt):
    pid = mp.current_process().pid
    
    # radial bins for the EE curves
    max_radius = 3.0   # arcsec
    radii = np.arange(0.05, max_radius, 0.05)  # Arcsec

    # Load up the image to work on.
    print(f"p{pid} - calc_star_stats on image: {img_file}")
    img, hdr = fits.getdata(img_file, header=True)

    # Get times
    d_hst, t_hst, d_utc, t_utc = util.get_times(hdr)

    # Get the bin fraction from the header
    plate_scale = util.get_plate_scale(img, hdr)

        
    # Read in the appropriate starlist for this image. 
    stars = table.Table.read(starlist, format='ascii')
    N_stars = len(stars)

    # Figure out quadrant information (relevant of four-filter only)
    if is_four_filt:
        quad = util.get_four_filter_quadrant(starlist)
    else:
        quad = ''

        
    # Put the positions into an array for photutils work.
    coords = np.array([stars['xcentroid'], stars['ycentroid']])

    # Define the background annuli (typically from 2"-3"). Calculate mean background for each star.
    bkg_annuli = CircularAnnulus(coords.T, max_radius / plate_scale, (max_radius + 1) / plate_scale)
    bkg = aperture_photometry(img, bkg_annuli)
    bkg_mean = bkg['aperture_sum'] / bkg_annuli.area

    enc_energy = np.zeros((N_stars, len(radii)), dtype=float)
    int_psf2_all = np.zeros(N_stars, dtype=float)

    # Loop through radial bins and calculate EE
    for rr in range(len(radii)):
        radius_pixel = radii[rr] / plate_scale
        apertures = CircularAperture(coords.T, r=radius_pixel)
        phot_table = aperture_photometry(img, apertures)

        energy = phot_table['aperture_sum']

        bkg_sum = apertures.area * bkg_mean
        
        enc_energy[:, rr] = energy - bkg_sum

        # Calculate the sum((PSF - bkg))^2 -- have to do this for each star.
        # Only do this on the last radius measurement.
        if rr == (len(radii) - 1):
            for ss in range(N_stars):
                aperture_ss = CircularAperture(coords[:,ss], r=radius_pixel)
                
                phot2_table = aperture_photometry((img - bkg_mean[ss])**2, aperture_ss)

                int_psf2 = phot2_table['aperture_sum'][0]
                int_psf2 /= enc_energy[ss, rr]**2 # normalize

                int_psf2_all[ss] = int_psf2

    # Normalize all the curves so that the mean of the last 5 bins = 1
    enc_energy /= np.tile(enc_energy[:, -5:].mean(axis=1), (len(radii), 1)).T

    # Calculate the EE 25, 50, and 80% radii for each star separately for saving.
    i_ee25_rad = np.zeros(N_stars, dtype=float)
    i_ee50_rad = np.zeros(N_stars, dtype=float)
    i_ee80_rad = np.zeros(N_stars, dtype=float)
    for ss in range(N_stars):
        i_ee25_rad[ss] = radii[ np.where(enc_energy[ss] >= 0.25)[0][0] ]
        i_ee25_rad[ss] = radii[ np.where(enc_energy[ss] >= 0.50)[0][0] ]
        i_ee25_rad[ss] = radii[ np.where(enc_energy[ss] >= 0.80)[0][0] ]

    stars['ee25_rad'] = i_ee25_rad
    stars['ee50_rad'] = i_ee50_rad
    stars['ee80_rad'] = i_ee80_rad

    # Find the median EE curve. But first, trim to just the brightest stars.
    idx = np.where(stars['flux'] > 5)[0]
    if len(idx) == 0:
        # Didn't find any bright stars... use all of them.
        idx = np.arange(N_stars)
    enc_energy_final = np.median(enc_energy[idx], axis=0)

    # Plot and save the EE curve and data.
    img_dir_name, img_file_name = os.path.split(img_file)
    ee_dir = img_dir_name + '/ee/'
    util.mkdir(ee_dir)
    _ee_out = open(ee_dir + img_file_name.replace('.fits', '_ee.txt'), 'w')
    _ee_out.write('{0:10s}  {1:10s}\n'.format('#Radius', 'EE'))
    _ee_out.write('{0:10s}  {1:10s}\n'.format('#(arcsec)', '()'))
    for rr in range(len(radii)):
        _ee_out.write('{0:10.2f}  {1:10.4f}\n'.format(radii[rr], enc_energy_final[rr]))
    _ee_out.close()
    
    # Find the 50% and 80% EE values
    ee25_rad = radii[ np.where(enc_energy_final >= 0.25)[0][0]]
    ee50_rad = radii[ np.where(enc_energy_final >= 0.5)[0][0] ]
    ee80_rad = radii[ np.where(enc_energy_final >= 0.8)[0][0] ]

    # Find the median NEA in pixel^2.
    stars['nea2'] = (1.0 / int_psf2_all)       # inidividual stars
    nea2 = 1.0 / np.median(int_psf2_all[idx])  # combined

    # Calculate the NEA in a different way. (in pixel^2)
    stars['nea'] = 0.0  # make new column
    r_dr_2pi = 2.0 * math.pi * (radii[1:] / plate_scale)  * (np.diff(radii) / plate_scale)
    for ss in range(N_stars):
        stars['nea'][ss] = 1.0 / (np.diff(enc_energy[ss])**2 / r_dr_2pi).sum()
    nea = 1.0 / (np.diff(enc_energy_final)**2 / r_dr_2pi).sum()

    #plt.clf()
    #plt.plot(radii, enc_energy_final, 'k-')
    #plt.axvline(ee50_rad, color='r', linestyle='--', label='r(50% EE)')
    #plt.axvline(ee80_rad, color='g', linestyle='--', label='r(80% EE)')
    #plt.xlabel('Radius (arcsec)')
    #plt.ylabel('Encircled Energy')
    #plt.legend(loc='lower right')
    #plt.pause(0.05)
    #plt.savefig(ee_dir + img_file_name.replace('.fits', '_ee.png'))

    # Calculate the average FWHM.
    xfwhm = stars['x_fwhm'][idx].mean()
    yfwhm = stars['y_fwhm'][idx].mean()
    fwhm = np.mean([stars['x_fwhm'][idx], stars['y_fwhm'][idx]])
    fwhm_std = np.std([stars['x_fwhm'][idx], stars['y_fwhm'][idx]])
    theta = stars['theta'][idx].mean()
    
    # calculate emperical FWHM 
    emp_FWHM_list = np.zeros(N_stars, dtype=float)

    for jj in range(N_stars):
        # Make a 21x21 patch centered on centroid, oversample and interpolate
        x_cent = int(round(float(coords[0][jj])))
        y_cent = int(round(float(coords[1][jj])))
        
        if (y_cent-10 > 0 and x_cent-10 > 0 and
            y_cent+10 < np.shape(img)[0] and x_cent+10 < np.shape(img)[1]):
            
            one_star = img[y_cent-10 : y_cent+10+1, x_cent-10 : x_cent+10+1]  # Odd box, with center in middle pixel.
            over_samp_5 = scipy.ndimage.zoom(one_star, 5, order = 1)

            # # Make an array with the radius at each pixel.
            # y_1d = np.arange(over_samp_5.shape[0])
            # x_1d = np.arange(over_samp_5.shape[1])
            # y_2d, x_2d = np.meshgrid(y_1d, x_1d)
            # r_2d = np.hypot(x_2d, y_2d)

            # Find the pixels where the flux is a above half max value.
            max_flux = np.amax(over_samp_5) 
            half_max = max_flux / 2.0
            idx = np.where(over_samp_5 >= half_max)
        
            # Find the equivalent circle diameter for the area of pixels.
            #    Area = \pi * (FWHM / 2.0)**2
            area_count = len(idx[0]) / 5**2   # area in pix**2 -- note we went back to raw pixels (not oversampled)
            emp_FWHM = 2.0 * (area_count / np.pi)**0.5
            emp_FWHM_list[jj] = emp_FWHM

    stars['fwhm_emp'] = emp_FWHM_list
    
    # Find the median empirical FWHM of all stars. But first, trim to just the brightest stars.
    idx = np.where(stars['flux'] > 5)[0]
    if len(idx) == 0:
        # Didn't find any bright stars... use all of them.
        idx = np.arange(N_stars)
    med_emp_FWHM = np.median(emp_FWHM_list[idx])
    std_emp_FWHM = np.std(emp_FWHM_list[idx])

    band = util.get_filter(hdr)
    binfac = util.get_bin_factor(img, hdr)

    # Save the individual star stats. (and plate scale info)
    stars.meta['scale'] = plate_scale
    new_starlist = img_file.replace('.fits', '_stars_stats.fits')
    stars.write(new_starlist, overwrite=True)

    results = {}
    results['band'] = band
    results['binfac'] = binfac
    results['time_utc'] = t_utc
    results['date_utc'] = d_utc
    results['time_hst'] = t_hst
    results['date_hst'] = d_hst
    results['fwhm'] = fwhm
    results['fwhm_std'] = fwhm_std
    results['ee25'] = ee25_rad
    results['ee50'] = ee50_rad
    results['ee80'] = ee80_rad
    results['NEA'] = nea
    results['NEA2'] = nea2
    results['xfwhm'] = xfwhm
    results['yfwhm'] = yfwhm
    results['theta'] = theta
    results['emp_fwhm'] = med_emp_FWHM
    results['emp_fwhm_std'] = std_emp_FWHM
    results['quadrant'] = quad
    
    return results

def add_frame_number_column(stats_table):
    # Get the frame numbers for plotting.
    frame_num = np.zeros(len(stats_table), dtype=int)

    for ii in range(len(stats_table)):

        full_name = stats_table['Image'][ii]
        last_bit = full_name.split("/")[-1]
        number = '0'
        found_digits = False
        
        for character in last_bit:
            if character.isdigit() == True:
                found_digits = True
                digit = (character)
                number += digit
            else:
                if found_digits:
                    # Condition is already found digits,
                    # and found first letter after that.
                    # If so, then we should stop.
                    break

        frame_num[ii] = int(number)

    stats_table['Index'] = frame_num

    return


def shift_and_add(img_files, starlists, output_root, method='mean', clip_sigma=None):
    """
    Take a stack of images and starlists, calculate the 
    """
    N_files = len(img_files)

    # Loop through all the starlists and get the transformations.
    # We will align all frames to the first frame. We can repeat
    # this process. 
    shift_trans = get_transforms_from_starlists(starlists)

    # Make some arrays to store the stack of images... this is a big array.
    # But it allows us to do some clever sigma clipping or median combine.
    # We only do this if asked (becaues it is RAM intensive).
    if method != 'mean':
        img, hdr = fits.getdata(img_files[0], header=True)
        
        image_stack = np.zeros((N_files, img.shape[0], img.shape[1]), dtype=float)
        image_cover = np.zeros((N_files, img.shape[0], img.shape[1]), dtype=int)

        ccddata_arr = []

    # Now loop through all the images and stack them.
    for ii in range(N_files):
        # Load up the image to work on.
        print("Shifting image: ", img_files[ii])
        img, hdr = fits.getdata(img_files[ii], header=True)

        # Make a coverage map that will also get shifted.
        img_covers = np.ones(img.shape, dtype=int)

        # Pull out the shifts
        shiftx = shift_trans[ii].px.c0_0.value
        shifty = shift_trans[ii].py.c0_0.value

        # Make an array to hold our final image.
        if method == 'mean':
            if ii == 0:
                final_image = np.zeros(img.shape, dtype=float)
                final_count = np.zeros(img.shape, dtype=int)

            # Shift the image and the coverage image and add to total
            final_image += interpolation.shift(img, (shifty, shiftx), order=1, cval=0)
            final_count += interpolation.shift(img_covers, (shifty, shiftx), order=1, cval=0)
        else:
            # image_stack[ii, :, :] = interpolation.shift(img, (shifty, shiftx), order=1, cval=0)
            # image_cover[ii, :, :] = interpolation.shift(img_covers, (shifty, shiftx), order=1, cval=0)
            img_tmp = interpolation.shift(img, (shifty, shiftx), order=1, cval=0)
            img_cnt_tmp = interpolation.shift(img_covers, (shifty, shiftx), order=1, cval=0)

            ccddata_cur = ccdproc.CCDData(img_tmp, unit=units.adu)

            ccddata_arr.append(ccddata_cur)

            
    # Done with all the images... do clipping, combining, or averaging.
    # Depends on method.
    if method == 'mean':
        final_image /= final_count
        final_image[final_count == 0] = 0

    if method == 'median':
        final_image = np.median(image_stack, axis=0)

    if method == 'meanclip':
        # avg = image_stack.mean(axis=0)
        # std = image_stack.std(axis=0)

        combiner = ccdproc.Combiner(ccddata_arr)
        combiner.sigma_clipping()
        final_image = combiner.average_combine()

    # Save file (note we are just using the last hdr... not necessarily the best)
    fits.writeto(output_root + '.fits', final_image, hdr, overwrite=True)
    
    return
    

def get_transforms_from_starlists(starlists):
    N_files = len(starlists)
    N_brite = 50
    N_passes = 2

    # Read in the first image and use this as the initial
    stars_ref = read_starlist(starlists[0])

    for nn in range(N_passes):
        # These will be the reference positions in the next pass... calculate them in this pass.
        x_ref_avg = np.zeros(len(stars_ref), dtype=float)
        y_ref_avg = np.zeros(len(stars_ref), dtype=float)
        n_ref_avg = np.zeros(len(stars_ref), dtype=int)

        trans = []

        for ii in range(len(starlists)):
            # Load up the corresponding starlist.
            stars = read_starlist(starlists[ii])

            t = align.initial_align(stars, stars_ref, briteN=N_brite, transformModel=transforms.Shift, req_match=3)

            idx1, idx2 = align.transform_and_match(stars, stars_ref, t, dr_tol=5, dm_tol=None)

            t, N_trans = align.find_transform(stars[idx1], None, stars_ref[idx2], transModel=transforms.Shift)

            trans.append(t)

            # Do a final matching to update reference positions.
            idx1, idx2 = align.transform_and_match(stars, stars_ref, t, dr_tol=5, dm_tol=None)

            xnew, ynew = t.evaluate(stars[idx1]['x'], stars[idx1]['y'])

            x_ref_avg[idx2] += xnew
            y_ref_avg[idx2] += ynew
            n_ref_avg[idx2] += 1

        if nn < (N_passes - 1):
            # Find the average position of any matched stars in the reference epoch.
            # This will improve our reference list later on. Trim down to those stars
            # in more than one frame.
            x_ref_avg /= n_ref_avg
            y_ref_avg /= n_ref_avg
        
            idx = np.where(n_ref_avg > 1)[0]
            stars_ref = stars_ref[idx]
        
            stars_ref['x'] = x_ref_avg[idx]
            stars_ref['y'] = y_ref_avg[idx]

    return trans


def read_starlist(starlist):
    """
    Read in a starlist and change the column names to be useful
    with flystar.
    """
    stars = table.Table.read(starlist, format='ascii.fixed_width'')

    stars.rename_column('xcentroid', 'x')
    stars.rename_column('ycentroid', 'y')
    stars.rename_column('mag', 'm')

    return stars


def _cpu_num():
    # How many cores to use for parallel functions
    # returns a number
    # Use N_cpu - 2 so we leave 2 for normal
    # operations. 
    cpu_count = mp.cpu_count()
    if (cpu_count > 2):
        cpu_count -= 2
        cpu_count = cpu_count if cpu_count <= 8 else 8
        
    # reurn 1 #(DEBUG)
    return cpu_count
