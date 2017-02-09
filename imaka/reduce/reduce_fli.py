import math
import pylab as plt
from PIL import Image
import numpy as np
from astropy.io import fits
from astropy import table
from astropy import units
import glob
import photutils
from astropy.convolution import Gaussian2DKernel
from astropy.stats import gaussian_fwhm_to_sigma
import poppy
import pdb
from astroscrappy import detect_cosmics
import pylab as plt
from PIL import Image
import numpy as np
from astropy.io import fits
from astropy import table
from astropy import units
import glob
import photutils
from photutils import psf
from photutils import morphology as morph
from photutils import DAOStarFinder
from photutils import CircularAperture
from photutils import CircularAnnulus
from photutils import aperture_photometry
from astropy.convolution import Gaussian2DKernel
from astropy.stats import gaussian_fwhm_to_sigma
from astropy.stats import sigma_clipped_stats
from astropy.modeling import models, fitting
import poppy
import pdb
from flystar import match
from flystar import align
from flystar import transforms
from imaka.reduce import calib
import ccdproc
from scipy.ndimage import interpolation
import scipy.ndimage
from astropy.table import Table


scale = 40.0 # mas/pixel

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

def fix_bad_pixels(img):
    crmask, cleanarr = detect_cosmics(img, sigclip=5, sigfrac=0.3, objlim=5.0, gain=1.0,
                                      readnoise=6.5, satlevel=65535.0, pssl=0, niter=4,
                                      sepmed=True, cleantype='meanmask', fsmode='median',
                                      psfmodel='gauss', psffwhm=2.5, psfsize=7,
                                      psfk=None, psfbeta=4.765, verbose=False)
                       
    return cleanarr


def clean_images(img_files, rebin=10, sky_frame=None):
    from scipy.ndimage import median_filter

    sky_bin = None   # Indiciates we haven't loaded the sky frame yet.

    for ii in range(len(img_files)):
        print('Working on image: ', img_files[ii])
        img, hdr = fits.getdata(img_files[ii], header=True)

        # # Fix hot and dead pixels
        # print('\t Fix bad pixels')
        # bad_pixels, img_cr = find_outlier_pixels(img, tolerance=5, worry_about_edges=False)
        img_cr = img
 
        # Bin by a factor of 10 (or as specified).
        print('\t Bin by a factor of ', rebin)
        from skimage.measure import block_reduce
        img_bin = block_reduce(img_cr, (rebin, rebin), func=np.sum)
        fits.writeto(img_files[ii].replace('.fits', '_bin.fits'), img_bin, hdr, clobber=True)

        # Subtract background.
        if sky_frame == None:
            bkg_box_size = int(round(img_bin.shape[0] / 10.))
            print('\t Subtract background with smoothing box size of ', bkg_box_size)
            blurred = median_filter(img_bin, size=bkg_box_size)
            img_final = img_bin - blurred.astype(float)
        else:
            if sky_bin == None:
                sky = fits.getdata(sky_frame)
                sky_bin = block_reduce(sky, (rebin, rebin), func=np.sum)
                
            img_final = img_bin - sky_bin

        # Save the final image.
        fits.writeto(img_files[ii].replace('.fits', '_bin_nobkg.fits'), img_final, hdr, clobber=True)

    return

def makeflat(flat_files, dark_files):
    
    ##Dark subtracts and combines flats

    # Read in the dark images
    print('\nReading in dark files...')
    darks = np.array([fits.getdata(dark) for dark in dark_files])

    # Read in the flat images
    print('\nReading in flat files...')
    flats = np.array([fits.getdata(flat) for flat in flat_files])

    #Dark subtraction
    print('\nSubtracting darks from flats...')
    flat_ds = flats - darks 

    #Normalize each of the flats:
    print('\nNormalizing each flat...')
    norm_flats = []
    for i in range(len(flat_ds)):
        mean, median, stdev = sigma_clipped_stats(flat_ds[i], sigma=3, iters=2, axis=None)
        norm = flat_ds[i] / median
        norm_flats.append(norm)

    # Sigma clip  (about the median)
    print('\nSigma clipping (3 sigma, 2 iterations) about the median...')
    mean, median, stdev = sigma_clipped_stats(norm_flats, sigma=3, iters=2, axis=0)

    # Median combine the masked images.
    print('\nMedian combining all images...')
    flat = median.data
    
    return flat

def rebin(a, bin_fac):
    
    #bins an array 'a' by a factor 'bin_fac'; sums pixels

    img_bin = block_reduce(a, (bin_fac, bin_fac), func=np.sum)
    return img_bin


def flat_sky_reduction(img_files, output_path, sky_frame, flat_frame):
    
    ##Subtracts sky and applies flat
    
    for ii in range(len(img_files)):
        print('Working on image: ', img_files[ii])
        img, hdr = fits.getdata(img_files[ii], header=True)
        
        sky = fits.getdata(sky_frame)
        flat = fits.getdata(flat_frame)
        
        img_final = (img - sky) / flat
        
        
        fits.writeto(output_path + img_files[ii].replace('.fits', '_bin_nobkg.fits').split('/')[-1], img_final, hdr, clobber=True)
        
    return


def find_stars_bin(img_files, fwhm=5, threshold=4, N_passes=2, plot_psf_compare=False):
    """
    img_files - a list of image files.
    fwhm - First guess at the FWHM of the PSF for the first pass on the first image.
    threshold - the SNR threshold (mean + threshold*std) above which to search for sources.
    N_passes - how many times to find sources, recalc the PSF FWHM, and find sources again. 
    """
    for ii in range(len(img_files)):
        print("Working on image: ", img_files[ii])
        img, hdr = fits.getdata(img_files[ii], header=True)

        fwhm_curr = fwhm

        # Calculate the bacgkround and noise (iteratively)
        print("\t Calculating background")
        bkg_threshold = 3
        for nn in range(5):
            if nn == 0:
                bkg_mean = img.mean()
                bkg_std = img.std()
            else:
                bkg_mean = img[good_pix].mean()
                bkg_std = img[good_pix].std()

            bad_hi = bkg_mean + (bkg_threshold * bkg_std)
            bad_lo = bkg_mean - (bkg_threshold * bkg_std)

            good_pix = np.where((img < bad_hi) & (img > bad_lo))
        
        bkg_mean = img[good_pix].mean()
        bkg_std = img[good_pix].std()
        img_threshold = threshold * bkg_std
        print('     Bkg = {0:.2f} +/- {1:.2f}'.format(bkg_mean, bkg_std))
        print('     Bkg Threshold = {0:.2f}'.format(img_threshold))

        # Detect stars
        print('     Detecting Stars')

        # Each pass will have an updated fwhm for the PSF.
        for nn in range(N_passes):
            print('     Pass {0:d} assuming FWHM = {1:.1f}'.format(nn, fwhm_curr))
            daofind = DAOStarFinder(fwhm=fwhm_curr, threshold = img_threshold, exclude_border=True)
            sources = daofind(img - bkg_mean)

            # Calculate FWHM for each detected star.
            x_fwhm = np.zeros(len(sources), dtype=float)
            y_fwhm = np.zeros(len(sources), dtype=float)
            theta = np.zeros(len(sources), dtype=float)
        
            cutout_half_size = int(round(fwhm_curr * 3))
            cutout_size = 2 * cutout_half_size

            final_psf_obs = np.zeros((cutout_size, cutout_size), dtype=float)
            final_psf_mod = np.zeros((cutout_size, cutout_size), dtype=float)
            final_psf_count = 0
        
            cutouts = np.zeros((len(sources), cutout_size, cutout_size), dtype=float)
            sigma_init_guess = fwhm_curr * gaussian_fwhm_to_sigma
            g2d_model = models.Gaussian2D(1.0, cutout_half_size, cutout_half_size,
                                              sigma_init_guess, sigma_init_guess, theta=0,
                                              bounds={'x_stddev':[0, 20], 'y_stddev':[0, 20]})
            g2d_fitter = fitting.LevMarLSQFitter()
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
            
                cutouts[ss] = cutout_tmp
                cutouts[ss] /= cutouts[ss].sum()

                # Fit an elliptical gaussian to the cutout image.
                g2d_params = g2d_fitter(g2d_model, cut_x, cut_y, cutouts[ss])
                g2d_image = g2d_params(cut_x, cut_y)
                
                final_psf_count += 1
                final_psf_obs += cutout_tmp
                final_psf_mod += g2d_image

                if plot_psf_compare == True:
                    plt.figure(4)
                    plt.clf()
                    plt.imshow(cutouts[ss])
                    plt.pause(0.05)
                
                    plt.figure(5)
                    plt.clf()
                    plt.imshow(g2d_image)
                    plt.pause(0.05)
                
                    plt.figure(6)
                    plt.clf()
                    plt.imshow(cutouts[ss] - g2d_image)
                    plt.pause(0.05)
                
                    pdb.set_trace()

                x_fwhm[ss] = g2d_params.x_stddev.value / gaussian_fwhm_to_sigma
                y_fwhm[ss] = g2d_params.y_stddev.value / gaussian_fwhm_to_sigma
                theta[ss] = g2d_params.theta.value

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
            fits.writeto(img_files[ii].replace('.fits', '_psf_obs.fits'), final_psf_obs, hdr, clobber=True)
            fits.writeto(img_files[ii].replace('.fits', '_psf_mod.fits'), final_psf_mod, hdr, clobber=True)

            # Drop sources with flux (signifiance) that isn't good enough.
            # Empirically this is <1.2
            good = np.where(sources['flux'] > 1.2)[0]
            sources = sources[good]

            # Only use the brightest sources for calculating the mean. This is just for printing.
            idx = np.where(sources['flux'] > 5)[0]
            x_fwhm_med = np.median(sources['x_fwhm'][idx])
            y_fwhm_med = np.median(sources['y_fwhm'][idx])
            
            print('        Number of sources = ', len(sources))
            print('        Median x_fwhm = {0:.1f} +/- {1:.1f}'.format(x_fwhm_med,
                                                                     sources['x_fwhm'].std()))
            print('        Median y_fwhm = {0:.1f} +/- {1:.1f}'.format(y_fwhm_med,
                                                                     sources['y_fwhm'].std()))

            fwhm_curr = np.mean([x_fwhm_med, y_fwhm_med])


            formats = {'xcentroid': '%8.3f', 'ycentroid': '%8.3f', 'sharpness': '%.2f',
                       'roundness1': '%.2f', 'roundness2': '%.2f', 'peak': '%10.1f',
                       'flux': '%10.6f', 'mag': '%6.2f', 'x_fwhm': '%5.2f', 'y_fwhm': '%5.2f',
                       'theta': '%6.3f'}
        
            sources.write(img_files[ii].replace('.fits', '_stars.txt'), format='ascii.fixed_width',
                          delimiter=None, bookend=False, formats=formats)

        
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
    pdb.set_trace()

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

    
        
def calc_star_stats(img_files, output_stats='image_stats.fits'):
    """
    Calculate statistics for the Data Metrics table.
    """
    plate_scale_orig = 0.04 # " / pixel

    # radial bins for the EE curves
    max_radius = 3.0
    radii = np.arange(0.05, max_radius, 0.05)

    # Create arrays for all the final statistics.
    N_files = len(img_files)
    s_time = np.zeros(N_files, dtype='S15')
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
    
    for ii in range(N_files):
        # Load up the image to work on.
        print("Working on image: ", img_files[ii])
        img, hdr = fits.getdata(img_files[ii], header=True)

        s_time[ii] = hdr['TIMEOBS']

        # Get the bin fraction from the header
        bin_factor = hdr['BINFAC']
        plate_scale = plate_scale_orig * bin_factor

        # Load up the corresponding starlist.
        starlist = img_files[ii].replace('.fits', '_stars.txt')
        stars = table.Table.read(starlist, format='ascii')
        N_stars = len(stars)

        # Put the positions into an array for photutils work.
        coords = np.array([stars['xcentroid'], stars['ycentroid']])

        # Define the background annuli (typically from 2"-3"). Calculate mean background for each star.
        bkg_annuli = CircularAnnulus(coords, max_radius / plate_scale, (max_radius + 1) / plate_scale)
        bkg = aperture_photometry(img, bkg_annuli)
        bkg_mean = bkg['aperture_sum'] / bkg_annuli.area()

        enc_energy = np.zeros((N_stars, len(radii)), dtype=float)
        int_psf2_all = np.zeros(N_stars, dtype=float)

        # Loop through radial bins and calculate EE
        for rr in range(len(radii)):
            radius_pixel = radii[rr] / plate_scale
            apertures = CircularAperture(coords, r=radius_pixel)
            phot_table = aperture_photometry(img, apertures)

            energy = phot_table['aperture_sum']

            bkg_sum = apertures.area() * bkg_mean
            
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

        # Find the median EE curve. But first, trim to just the brightest stars.
        idx = np.where(stars['flux'] > 5)[0]
        if len(idx) == 0:
            # Didn't find any bright stars... use all of them.
            idx = np.arange(N_stars)
        enc_energy_final = np.median(enc_energy[idx], axis=0)

        # Plot and save the EE curve and data.
        _ee_out = open('ee/' + img_files[ii].replace('.fits', '_ee.txt'), 'w')
        _ee_out.write('{0:10s}  {1:10s}\n'.format('#Radius', 'EE'))
        _ee_out.write('{0:10s}  {1:10s}\n'.format('#(arcsec)', '()'))
        for rr in range(len(radii)):
            _ee_out.write('{0:10.2f}  {1:10.4f}\n'.format(radii[rr], enc_energy_final[rr]))
        _ee_out.close()

        # Find the 50% and 80% EE values
        ee50_rad = radii[ np.where(enc_energy_final >= 0.5)[0][0] ]
        ee80_rad = radii[ np.where(enc_energy_final >= 0.8)[0][0] ]

        # Find the median NEA. Convert into arcsec^2
        nea2 = 1.0 / np.median(int_psf2_all[idx])
        nea2 *= plate_scale**2

        # Calculate the NEA in a different way.
        nea = 1.0 / (np.diff(enc_energy_final)**2 / (2.0 * math.pi * radii[1:] * np.diff(radii))).sum()

        
        plt.clf()
        plt.plot(radii, enc_energy_final, 'k-')
        plt.axvline(ee50_rad, color='r', linestyle='--', label='r(50% EE)')
        plt.axvline(ee80_rad, color='g', linestyle='--', label='r(80% EE)')
        plt.xlabel('Radius (arcsec)')
        plt.ylabel('Encircled Energy')
        plt.legend(loc='lower right')
        plt.pause(0.05)
        plt.savefig('ee/' +  img_files[ii].replace('.fits', '_ee.png'))

        # Calculate the average FWHM.
        xfwhm = stars['x_fwhm'][idx].mean()
        yfwhm = stars['y_fwhm'][idx].mean()
        fwhm = np.mean([stars['x_fwhm'][idx], stars['y_fwhm'][idx]])
        fwhm_std = np.std([stars['x_fwhm'][idx], stars['y_fwhm'][idx]])
        theta = stars['theta'][idx].mean()
        
        # calculate emperical FWHM 

        emp_FWHM_list = []
       
        for jj in range(len(stars)):
           
            # Make a 20x20 patch centered on centroid, oversample and interpolate
            x_cent = float(coords[0][jj]); y_cent = float(coords[1][jj])
            one_star = img[y_cent-10 : y_cent+10, x_cent-10 : x_cent+10]
            over_samp_5 = scipy.ndimage.zoom(one_star, 5, order = 1)

            # Find area of stars above half max and calculate equivalnent circle diameter
            max_flux = np.amax(over_samp_5) 
            half_max = max_flux / 2.0
            idx = np.where(over_samp_5 >= half_max)
            area_count = len(idx[0])
            emp_FWHM = (2.0* ((area_count / np.pi) ** 0.5)) / 5.0
            emp_FWHM_list.append(emp_FWHM)
 
        # Take median of all stars in image
        med_emp_FWHM = np.median(emp_FWHM_list)
        std_emp_FWHM = np.std(emp_FWHM_list)

        s_ee50[ii] = ee50_rad
        s_ee80[ii] = ee80_rad
        s_xfwhm[ii] = xfwhm
        s_yfwhm[ii] = yfwhm
        s_theta[ii] = theta
        s_fwhm[ii] = fwhm
        s_fwhm_std[ii] = fwhm_std
        s_NEA[ii] = nea
        s_NEA2[ii] = nea2
        s_emp_fwhm[ii] = med_emp_FWHM
        s_emp_fwhm_std[ii] = std_emp_FWHM

    stats = table.Table([img_files, s_time, s_fwhm, s_fwhm_std, s_ee50, s_ee80, s_NEA, s_NEA2, s_xfwhm, s_yfwhm, s_theta, s_emp_fwhm, s_emp_fwhm_std],
                                names=('Image', 'TIME', 'FWHM', 'FWHM_std', 'EE50', 'EE80', 'NEA', 'NEA2', 'xFWHM', 'yFWHM', 'theta', 'emp_fwhm', 'emp_fwhm_std'),
                            meta={'name':'Stats Table'})
    
    stats['FWHM'].format = '7.3f'
    stats['FWHM_std'].format = '7.3f'
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
    stats.write(output_stats.replace('.fits', '.csv'), format='csv') # Auto overwrites
                        
    return


def add_frame_number_column(stats_table):
    # Get the frame numbers for plotting.
    frame_num = np.zeros(len(stats_table), dtype=int)
    for ii in range(len(stats_table)):
        foo = stats_table['Image'][ii].index('bin')

        frame_num[ii] = int(stats_table['Image'][ii][foo - 4:foo - 1])
        frame_num_col = table.Column(frame_num, name='Index')

    stats_table.add_column(frame_num_col, index=1)

    return


def shift_and_add(img_files, starlists, output_root, method='mean', clip_sigma=None):
    """
    Take a stack of images and starlists, calculate the 
    """
    plate_scale_orig = 0.04 # " / pixel
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
        shiftx = shift_trans[ii].px[0]
        shifty = shift_trans[ii].py[0]

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

    fits.writeto(output_root + '.fits', final_image.data)
    

    return
    

def get_transforms_from_starlists(starlists):
    plate_scale_orig = 0.04 # " / pixel
    N_files = len(starlists)
    N_brite = 7
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
    stars = table.Table.read(starlist, format='ascii')

    stars.rename_column('xcentroid', 'x')
    stars.rename_column('ycentroid', 'y')
    stars.rename_column('mag', 'm')

    return stars


def calc_empirical_FWHM(img_files, star_lists):
    
    #Input: image files, obj_x000_bin_nobkg_stars.txt file !!must be run after find_stars!!
    #output: writes new file obj_x_000_bin_nobkg_stars_FWHM.txt

    #read in image data
    images = np.array([fits.getdata(image) for image in img_files])

    for ii in range(len(images)):
        print("Working on image: ", img_files[ii].split('/')[-1])

        #read in centroid coordinates from _bin_nobkg_stars.txt
        with open(star_lists[ii], "r") as f:
            xcents = []; ycents = []; ids    = []
            for line in f:
                pieces = line.split('   ')
                x_cent = pieces[1]; y_cent = pieces[2]; iden = pieces[0]
                xcents.append(x_cent); ycents.append(y_cent); ids.append(iden)
        xcents.pop(0); ycents.pop(0), ids.pop(0)

        FWHM_list = []
        peaks = []
        
        for i in range(len(xcents)):
           
            # Make a 20x20 patch centered on centroid, oversample and interpolate
            x_cent = float(xcents[i]); y_cent = float(ycents[i])
            one_star = images[ii][y_cent-10 : y_cent+10, x_cent-10 : x_cent+10]
            over_samp_5 = scipy.ndimage.zoom(one_star, 5, order = 1)

            
            # Find area of stars above half max and calculate equivalnent circle diameter
            max_flux = np.amax(over_samp_5) 
            peaks.append(max_flux)
            half_max = max_flux / 2.0
            idx = np.where(over_samp_5 >= half_max)
            area_count = len(idx[0])
            FWHM = (2.0* ((area_count / np.pi) ** 0.5)) / 5.0
            FWHM_list.append(FWHM)
        
        print("     Sources calculated:", len(FWHM_list))
        print("     Writing text file...")
        
        # Write data to file
        t = Table([ids, xcents, ycents, FWHM_list, peaks],\
                  names=('id', 'xcentroid', 'ycentroid', 'FWHM', 'peak'), meta={'name': 'empirical FWHM'})
                            
        t.write(img_files[ii].replace('.fits', '_stars_FWHM.txt'), format='ascii.fixed_width', delimiter=None, bookend=False)

    return
