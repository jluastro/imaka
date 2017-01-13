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
from imaka.reduce import reduce_fli

from imaka.reduce import calib


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
                                              sigma_init_guess, sigma_init_guess, theta=0)
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

            x_fwhm_med = np.median(sources['x_fwhm'])
            y_fwhm_med = np.median(sources['y_fwhm'])
            
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
    s_ee50 = np.zeros(len(img_files), dtype=float)
    s_ee80 = np.zeros(len(img_files), dtype=float)
    s_xfwhm = np.zeros(len(img_files), dtype=float)
    s_yfwhm = np.zeros(len(img_files), dtype=float)
    s_theta = np.zeros(len(img_files), dtype=float)
    s_fwhm = np.zeros(len(img_files), dtype=float)
    s_fwhm_std = np.zeros(len(img_files), dtype=float)
    s_NEA = np.zeros(len(img_files), dtype=float)
    
    for ii in range(len(img_files)):
        # Load up the image to work on.
        print("Working on image: ", img_files[ii])
        img, hdr = fits.getdata(img_files[ii], header=True)

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
            # if rr == (len(radii) - 1):
            if radii[rr] == 1.0:
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
        nea = 1.0 / int_psf2_all[idx].mean()
        nea *= plate_scale**2

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

        s_ee50[ii] = ee50_rad
        s_ee80[ii] = ee80_rad
        s_xfwhm[ii] = xfwhm
        s_yfwhm[ii] = yfwhm
        s_theta[ii] = theta
        s_fwhm[ii] = fwhm
        s_fwhm_std[ii] = fwhm_std
        s_NEA[ii] = nea

    stats = table.Table([img_files, s_fwhm, s_fwhm_std, s_ee50, s_ee80, s_NEA, s_xfwhm, s_yfwhm, s_theta],
                                names=('Image', 'FWHM', 'FWHM_std', 'EE50', 'EE80', 'NEA', 'xFWHM', 'yFWHM', 'theta'),
                            meta={'name':'Stats Table'})

    stats.write(output_stats, overwrite=True)
                        
    return
