import numpy as np 
from astropy import units as u
from astropy.nddata import CCDData
from astropy.io import fits
from astropy.modeling import models
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
#from flystar import match
#from flystar import align
#from flystar import transforms
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
from imaka.reduce import reduce_fli
from astropy.table import Column


def treat_overscan(files):
    """
    Performs bias/dark subtraction from overscan images on STA camera
    Assumes STA 8x2 'cell' structure, each with its own overscan region
    Works independent of image binning
    Input: an array of files from STA camera
    Output: writes new images with '_scan' suffix in name
    """
    for file in files:
        
        # Read in data and divide into 'cells'
        print('Working on file: ', file)
        dat, hdr = fits.getdata(file, header=True)
        rows = np.vsplit(dat, 2)
        one_row = np.hstack([rows[0], rows[1]])
        cells = np.hsplit(one_row, 16)

        # Define the 'image' and 'overscan' region 
        # Done with fractions, so independent of binning
        imgs = []
        scans = []
        clean_cells = []

        orig_wid = np.shape(dat)[1]        # Width of original image
        cell_wid = orig_wid / 8            # Width of one cell (1/8 of image)
        scan_wid = int(cell_wid * (3/25))  # Width of overscan region (3/25 of cell)
        data_wid = int(cell_wid * (22/25)) # Width of data region (22/25 of cell)

        # For each cell, take the median of each overscan row and 
        # subtract that from the corresponding row in the data
        for ii in range(16):
            scan, img = np.split(cells[ii], [scan_wid], axis=1)
            scans.append(scan)
            imgs.append(img)
            med_col = np.median(scan, axis=1)
            med_cell = np.array([med_col,]*data_wid).transpose()
            clean_cell = img - med_cell 
            clean_cells.append(clean_cell)

        # Rebuild image with clean cells
        clean_one_row = np.hstack(clean_cells)
        clean_rows = np.hsplit(clean_one_row, 2)
        clean_image = np.vstack(clean_rows)

        # Write new image
        new_filename = file[:-5] + '_scan.fits'
        fits.writeto(new_filename, clean_image, hdr, overwrite=True)
    
    return

def clean_image(image_file):
    """
    Cleans a single image.  Cleaning includes:
        -Gain application
        -Cosmic ray removal
    
    Returns: astropy.nddata.ccddata.CCDData object
    
    Header needs: 'GAINI' (optional, but good); 'EXPTIME', 'BINFAC'
    """

    # Read in image and header
    img, hdr = fits.getdata(image_file, header=True)
    if 'GAINI' in hdr:
        gain = hdr['GAINI'] * u.electron/u.adu
    else:
        print("Gain not found in header: assuming 1 e/adu")
        gain = 1 *u.electron/u.adu

    readnoise = 5 * u.electron #From STA1600LN data sheet

    # Define data and uncertainty
    data = CCDData(img, unit=u.adu)
    data_with_deviation = ccdproc.create_deviation(data, gain=gain, readnoise=readnoise)
    data_with_deviation.header['exposure'] = float(hdr['EXPTIME'])

    # Apply gain to data
    gain_corrected = ccdproc.gain_correct(data_with_deviation, gain)

    # Clean cosmic rays
    cr_cleaned = ccdproc.cosmicray_lacosmic(gain_corrected, sigclip=5)

    return cr_cleaned


def create_bias(bias_files):
    """
    Cleans and median combines an input list of bias files with sigma clipping 3-σ about median
    
    Input: a list of bias image fits files 

    Returns: Saves master bias frame as 'master_bias.fits' in same directory as bias frames
    """
    
    # Read in bias images
    combiner = []
    for file in bias_files:
        print("Working on file: ", file.split("/")[-1])
        ccd = clean_image(file)
        combiner.append(ccd)

    # Sigma clip 3-σ about median
    print("Sigma clipping")
    combiner.sigma_clipping(func=np.ma.median)

    # Combine frames
    print("Combining frames")
    combined_median = combiner.median_combine()

    # Write to fits file with header info of last file
    print("Writing file")
    bias_dir = file.split('bias_')[0]
    combined_median.write(bias_dir+'master_bias.fits')

    return


def calc_star_stats(img_files, output_stats='image_stats.fits', filt=None, fourfilt=False, starlists=None):
    """
    Calculate statistics for the Data Metrics table.
    """
    plate_scale_orig = 0.016 # " / pixel

    # radial bins for the EE curves
    max_radius = 3.0
    radii = np.arange(0.05, max_radius, 0.05)

    # Create arrays for all the final statistics.
    N_files = len(img_files)
    s_time_hst = np.zeros(N_files, dtype='S15')
    s_time_utc = np.zeros(N_files, dtype='S15')
    s_date_hst = np.zeros(N_files, dtype='S15')
    s_date_utc = np.zeros(N_files, dtype='S15')
    band = np.zeros(N_files, dtype='S15')
    binfac = np.zeros(N_files, dtype=float)
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

    if fourfilt==True:
        # Make a column to designate filter position in four filter images
        quadrant = np.zeros(N_files, dtype='S15')
    
    for ii in range(N_files):
        # Load up the image to work on.
        print("Working on image: ", img_files[ii])
        img, hdr = fits.getdata(img_files[ii], header=True)

        # Make dates and times in UT and HST
        s_time_hst[ii] = hdr['SYS-TIME']
        s_date_hst[ii] = hdr['SYS-DATE']
        s_date_utc[ii] = hdr['DATE-OBS']

        time = hdr['UT']
        h, m, s = time.split(':')
        new_time = h+":"+m+":"+s[:2]
        s_time_utc[ii] = new_time
        

        # Get the bin fraction from the header
        bin_factor = hdr['CCDBIN1']
        plate_scale_orig = 0.016
        plate_scale = plate_scale_orig * bin_factor

        if starlists == None:
            # Load up the corresponding starlist.
            if filt==None:
                starlist = img_files[ii].replace('.fits', '_stars.txt')
            else:
                starlist = img_files[ii].replace('.fits', '_'+filt+'_stars.txt')
        elif starlists != None:
            starlist = starlists[ii]
        stars = table.Table.read(starlist, format='ascii')
        N_stars = len(stars)

        if fourfilt==True:
            # Get filter position in case of four filter data
            name_strings = starlist.split("_")
            filt_name    = name_strings[-3]
            filt_order   = name_strings[-2]

            if filt_name == filt_order[0]:
                quad = 'NW'
            elif filt_name == filt_order[1]:
                quad = 'NE'
            elif filt_name == filt_order[2]:
                quad = 'SE'
            elif filt_name == filt_order[3]:
                quad = 'SW'

            quadrant[ii] = quad
            
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
        img_dir_name, img_file_name = os.path.split(img_files[ii])
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

        # Find the median NEA. Convert into arcsec^2
        nea2 = 1.0 / np.median(int_psf2_all[idx])
        nea2 *= plate_scale**2

        # Calculate the NEA in a different way.
        nea = 1.0 / (np.diff(enc_energy_final)**2 / (2.0 * math.pi * radii[1:] * np.diff(radii))).sum()

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
            if y_cent-10 > 0 and x_cent-10 > 0 and y_cent+10<np.shape(img)[0] and x_cent+10<np.shape(img)[1]:
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
 
        # Find the median empirical FWHM of all stars. But first, trim to just the brightest stars.
        idx = np.where(stars['flux'] > 5)[0]
        if len(idx) == 0:
            # Didn't find any bright stars... use all of them.
            idx = np.arange(N_stars)
        med_emp_FWHM = np.median(emp_FWHM_list[idx])
        std_emp_FWHM = np.std(emp_FWHM_list[idx])

        band[ii] = "I"#hdr['FILTER']
        binfac[ii] = hdr['CCDBIN1']
        s_ee25[ii] = ee25_rad
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

    
    # Make a date array that holds UTC.

    stats = table.Table([img_files, band, binfac, s_date_utc, s_time_utc, s_date_hst, s_time_hst,
                             s_fwhm, s_fwhm_std, s_ee25, s_ee50, s_ee80,
                             s_NEA, s_NEA2, s_xfwhm, s_yfwhm, s_theta, s_emp_fwhm, s_emp_fwhm_std],
                             names=('Image', 'FILTER', 'BINFAC', 'DATE_UTC', 'TIME_UTC', 'DATE_HST', 'TIME_HST',
                                        'FWHM', 'FWHM_std', 'EE25', 'EE50', 'EE80',
                                        'NEA', 'NEA2', 'xFWHM', 'yFWHM', 'theta', 'emp_fwhm', 'emp_fwhm_std'),
                            meta={'name':'Stats Table'})

    if fourfilt==True:
        quad_col = Column(quadrant, name='quadrant')
        stats.add_column(quad_col)
    
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

    reduce_fli.add_frame_number_column(stats)
    
    stats.write(output_stats, overwrite=True)
    #stats.write(output_stats.replace('.fits', '.csv'), format='csv') # Auto overwrites
                        
    return


def add_focus(stats_files):
    """
    Retrieves focus value from image headers for
    each frame in a stats_file, and creates a 
    new column with these values.  For STA camera.
    """
    for file in stats_files:
        tab = Table.read(file)
        N_frames = len(tab)
        focci = []
        for ii in range(N_frames):
            frame_file = tab['Image'][ii]
            hdr = fits.getheader(frame_file)
            focus = hdr['FOCUS']
            focci.append(focus)
        col_focus = Column(name='Focus', data=focci)
        tab.add_column(col_focus)
        tab.write(file, overwrite=True)

    return 


def four_filt_split(starlists, filt_order):
    
    """
    For STA camera with four filters.
    Splits starlist into four starlists,
    each corresponding to one quadrant/filter
    Inputs: starlists: array of STA camera starlists;
            filt_order: a four character string
            giving filter labels from in order:
            NW, NE, SE, SW, eg: 'BVIR'
    Output: Four new starlists for each input
            starlist given with ending:
            R_BVIR_stars.txt' for R filter and
            filter order 'BVIR'
    """

    filt_1 = filt_order[0]
    filt_2 = filt_order[1]
    filt_3 = filt_order[2]
    filt_4 = filt_order[3]

    N = len(starlists)
    for ii in range(N):
        print("Splitting file", ii+1, "of", N)
        print("\t", starlists[ii])

        # Read in image and starlists
        stars = table.Table.read(starlists[ii], format='ascii')
        x = stars['xcentroid']
        y = stars['ycentroid']

        # Define quadrants with indicies
        img_half = 5280 / 2

        ind_1 = np.where((x>img_half) & (y<img_half))
        ind_2 = np.where((x<img_half) & (y<img_half))
        ind_3 = np.where((x<img_half) & (y>img_half))
        ind_4 = np.where((x>img_half) & (y>img_half))
    
        # Write new files
        list_root = starlists[ii].split('stars.txt')[0]
        stars[ind_1].write(list_root + filt_1 + '_' + filt_order + '_stars.txt', format='ascii', overwrite=True)
        stars[ind_2].write(list_root + filt_2 + '_' + filt_order + '_stars.txt', format='ascii', overwrite=True)
        stars[ind_3].write(list_root + filt_3 + '_' + filt_order + '_stars.txt', format='ascii', overwrite=True)
        stars[ind_4].write(list_root + filt_4 + '_' + filt_order + '_stars.txt', format='ascii', overwrite=True)
    
    return
