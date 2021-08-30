### moffat.py - Takes reduced images and existing stats files and conducts moffat fitting, adds data to stats tables and makes model PSFs

import os
import numpy as np
from astropy.io import fits
from astropy.table import Table, Column
from astropy.modeling import models, fitting
from astropy.modeling.models import custom_model
from imaka.analysis import plot_stats
from imaka.reduce import util
import multiprocessing as mp
from itertools import repeat
import matplotlib.pyplot as plt
import scipy
@custom_model



def Elliptical_Moffat2D(x, y, \
                        N_sky = 0., amplitude = 1., phi=0., power = 1.,\
                        x_0 = 0., y_0 = 0., width_x = 1., width_y = 1.):
    """
    A custom astropy model for a two dimensional elliptical moffat function.  
    N_sky: a constant background value
    Amplitude: A
    phi: rotation angle (in radians)
    power: slope of model (beta)
    x_0, y_0: star's centroid coordinates
    width_x, width_y: core widths (alpha)
    """

    c = np.cos(phi)
    s = np.sin(phi)
    A = (c / width_x) ** 2 + (s / width_y)**2
    B = (s / width_x) ** 2 + (c/ width_y)**2
    C = 2 * s * c * (1/width_x**2 - 1/width_y**2)
    denom = (1 + A * (x-x_0)**2 + B * (y-y_0)**2 + C*(x-x_0)*(y-y_0))**power

    return N_sky + amplitude / denom 


def fit_moffat(img_files, stats_file, x_guess=5, y_guess=5, flux_percent=0.9, starlists=None):
    """
    Conduct Moffat fit on data and add outputs to stats tables
    """

    # Create arrays for all the final statistics.
    N_files   = len(img_files)
    mof_stars = np.zeros(N_files, dtype=float) #Number of stars used in medians
    N_sky     = np.zeros(N_files, dtype=float)
    amplitude = np.zeros(N_files, dtype=float)
    phi       = np.zeros(N_files, dtype=float)
    power     = np.zeros(N_files, dtype=float)
    x_0       = np.zeros(N_files, dtype=float)
    y_0       = np.zeros(N_files, dtype=float)
    width_x   = np.zeros(N_files, dtype=float)
    width_y   = np.zeros(N_files, dtype=float)

    N_sky_std     = np.zeros(N_files, dtype=float)
    amplitude_std = np.zeros(N_files, dtype=float)
    phi_std       = np.zeros(N_files, dtype=float)
    power_std     = np.zeros(N_files, dtype=float)
    x_0_std       = np.zeros(N_files, dtype=float)
    y_0_std       = np.zeros(N_files, dtype=float)
    width_x_std   = np.zeros(N_files, dtype=float)
    width_y_std   = np.zeros(N_files, dtype=float)

    ####
    # Setup for parallel processing.
    ####
    # Use N_cpu - 2 so we leave 2 for normal operations. 
    cpu_count = mp.cpu_count()
    if (cpu_count > 2):
        cpu_count -= 2

    # Start the pool
    pool = mp.Pool(cpu_count)
    results_async = []
    print(f'fit_moffat in parallel with {cpu_count} cores.')
    
    
    for ii in range(N_files):
        img_file = img_files[ii]
        
        # Load up the corresponding starlist.
        if starlists==None:
            starlist = img_file.replace('.fits', '_stars_stats.fits')
        else:
            starlist = starlists[ii]

        #####
        # Add calc for this starlist to the pool.
        #####
        results = pool.apply_async(fit_moffat_single, (img_file, starlist, flux_percent))
        results_async.append(results)

    pool.close()
    pool.join()

    for ii in range(N_files):
        results = results_async[ii].get()
        
        N_sky[ii]     = results['N_sky']
        amplitude[ii] = results['amplitude']
        phi[ii]       = results['phi']
        power[ii]     = results['power']
        x_0[ii]       = results['x_0']
        y_0[ii]       = results['y_0']
        width_x[ii]   = results['width_x']
        width_y[ii]   = results['width_y']

        N_sky_std[ii]     = results['N_sky_std']
        amplitude_std[ii] = results['amplitude_std']
        phi_std[ii]       = results['phi_std']
        power_std[ii]     = results['power_std']
        x_0_std[ii]       = results['x_0_std']
        y_0_std[ii]       = results['y_0_std']
        width_x_std[ii]   = results['width_x_std']
        width_y_std[ii]   = results['width_y_std']

        mof_stars[ii] = results['mof_stars']
        

    # Read in existing stats table and append new data
    stats = Table.read(stats_file)
    
    stats['N Stars'] = mof_stars
    stats['N Sky'] = N_sky
    stats['N Sky std'] = N_sky_std
    stats['Amplitude'] = amplitude
    stats['Amplitude std'] = amplitude_std
    stats['Phi'] = phi
    stats['Phi std'] = phi_std
    stats['Beta'] = power
    stats['Beta std'] = power_std
    stats['Minor Alpha'] = width_x
    stats['Minor Alpha std'] = width_x_std
    stats['Major Alpha'] = width_y
    stats['Major Alpha std'] = width_y_std

    stats_file_root, stats_file_ext = os.path.splitext(stats_file)
    stats.write(stats_file_root.split("_mdp")[0] + stats_file_ext, overwrite=True)
    
    return

def fit_moffat_single(img_file, starlist, flux_percent):
    pid = mp.current_process().pid
    
    print(f"p{pid} - Fitting moffat for {img_file} and {starlist}")
    
    # Load up the image to work on.
    img, hdr = fits.getdata(img_file, header=True)

    stars = Table.read(starlist, format='ascii')
    N_stars = len(stars)

    # Put the positions into an array 
    coords = np.array([stars['xcentroid'], stars['ycentroid']])

    # Make empty arrays for trimmed star sample
    x_cents = []
    y_cents = []
    N_good_stars = 0

    # Lets figure out which stars to actually fit.
    flux_sorted = np.sort(stars['flux'])                  # array of increasing flux
    sdx_flux_cut = int(len(flux_sorted) * flux_percent)   # set the cut at brightest YY% of all stars.
    flux_cut = flux_sorted[sdx_flux_cut]                  # this yields the flux lower limit.

    # Trim star sample for edge sources, saturation, and low flux
    cut_size = int(40 / util.get_bin_factor(img, hdr))
    dcut = cut_size / 2.0
    xlo = (stars['xcentroid'] - dcut).astype('int')
    xhi = xlo + cut_size + 1
    ylo = (stars['ycentroid'] - dcut).astype('int')
    yhi = ylo + cut_size + 1
    
    gdx = np.where((xlo > 0) & (xhi < img.shape[1]) &
                   (ylo > 0) & (yhi < img.shape[0]) &
                   (stars['flux'] > flux_cut) &
                   (stars['peak'] < 20000))[0]
    
    if len(gdx) < 2:
        gdx = np.where((xlo > 0) & (xhi < img.shape[1]) &
                       (ylo > 0) & (yhi < img.shape[0]) &
                       (stars['peak'] < 20000))[0]

    N_good_stars = len(gdx)
    stars = stars[gdx]
    x_cents = stars['xcentroid']
    y_cents = stars['ycentroid']
    xlo = xlo[gdx]
    xhi = xhi[gdx]
    ylo = ylo[gdx]
    yhi = yhi[gdx]
    print(f'p{pid} -   Moffat fitting for {N_good_stars}')
    
    # Create lists for each image's stars' parameters
    N_sky_list     = np.zeros(N_good_stars, dtype=float)
    amplitude_list = np.zeros(N_good_stars, dtype=float)
    phi_list       = np.zeros(N_good_stars, dtype=float)
    power_list     = np.zeros(N_good_stars, dtype=float)
    x_0_list       = np.zeros(N_good_stars, dtype=float)
    y_0_list       = np.zeros(N_good_stars, dtype=float)
    width_x_list   = np.zeros(N_good_stars, dtype=float)
    width_y_list   = np.zeros(N_good_stars, dtype=float)

    # We will fit on oversampled images.  This is the over-sampling factor.
    over_samp = 2
    
    final_psf_mof = np.zeros(((cut_size + 1) * over_samp,
                              (cut_size + 1) * over_samp), dtype=float)
    
    for jj in range(N_good_stars):
        # Make image cut for good star
        x_cent = x_cents[jj]
        y_cent = y_cents[jj]
        image_cut = img[ylo[jj]:yhi[jj], xlo[jj]:xhi[jj]]

        over_samp_cut = scipy.ndimage.zoom(image_cut, over_samp, order=1)

        # Conduct moffat fit
        y, x = np.mgrid[:over_samp_cut.shape[0], :over_samp_cut.shape[1]]
        z = over_samp_cut
        m_init = Elliptical_Moffat2D(N_sky = 0,
                                         amplitude = np.amax(z),
                                         x_0 = over_samp_cut.shape[0] / 2,
                                         y_0 = over_samp_cut.shape[1] / 2,
                                         width_x = 4.55 * over_samp,
                                         width_y = 4.17 * over_samp)
        fit_m = fitting.LevMarLSQFitter()
        m = fit_m(m_init, x, y, z)
             
        N_sky_list[jj]     = m.N_sky.value
        amplitude_list[jj] = m.amplitude.value
        power_list[jj]     = m.power.value
        x_0_list[jj]       = m.x_0.value / over_samp
        y_0_list[jj]       = m.y_0.value / over_samp
        if m.width_x.value < m.width_y.value:
            width_x_list[jj]    = m.width_x.value / over_samp
            width_y_list[jj]    = m.width_y.value / over_samp
            phi_list[jj]       = m.phi.value 
        else:
            width_x_list[jj]    = m.width_y.value / over_samp
            width_y_list[jj]    = m.width_x.value / over_samp
            phi_list[jj]       = m.phi.value + (np.pi/2)

        final_psf_mof += m(x, y)

    # Save the average PSF (flux-weighted). Note we are making a slight 
    # mistake here since each PSF has a different sub-pixel position
    final_psf_mof /= N_good_stars
    fits.writeto(img_file.replace('.fits', '_psf_mof_oversamp{0:d}.fits'.format(over_samp)),
                 final_psf_mof, hdr, overwrite=True)

    # Save and updated list of stars with all their moffat fits.
    stars['N Sky'] = N_sky_list
    stars['Amplitude'] = amplitude_list
    stars['Phi'] = phi_list
    stars['Beta'] = power_list
    stars['Minor Alpha'] = width_x_list
    stars['Major Alpha'] = width_y_list

    stars_file_root, stars_file_ext = os.path.splitext(starlist)
    stars.write(starlist.replace('.fits', '_mdp.fits'), overwrite=True)

    
    # Take median values of parameters and put in intial lists
    results = {}
    results['N_sky']     = np.median(N_sky_list)
    results['amplitude'] = np.median(amplitude_list)
    results['phi']       = np.median(phi_list)
    results['power']     = np.median(abs(power_list))
    results['x_0']       = np.median(x_0_list)
    results['y_0']       = np.median(y_0_list)
    results['width_x']   = np.median(abs(width_x_list))
    results['width_y']   = np.median(abs(width_y_list))

    results['N_sky_std']     = np.std(N_sky_list)
    results['amplitude_std'] = np.std(amplitude_list)
    results['phi_std']       = np.std(phi_list)
    results['power_std']     = np.std(abs(power_list))
    results['x_0_std']       = np.std(x_0_list)
    results['y_0_std']       = np.std(y_0_list)
    results['width_x_std']   = np.std(abs(width_x_list))
    results['width_y_std']   = np.std(abs(width_y_list))

    results['mof_stars'] = int(N_good_stars)

    return results


def calc_mof_fwhm(stats_file, filt=False, plate_scale=0.016):

    """
    Takes stats_<type>.fits file and outputs four arrays:
        minor FWHM
        minor FWHM uncertainty
        major FWHM
        major FWHM uncertainty
    
    all in arcsec.
    
    If filt=True, data is scaled with filter data to 500 nm

    Input
    ----------
    stats_file : str
        Name of the stats FITS table to read.

    filt : boolean
        Rescale the values to 500 nm

    plate_scale : float
        The plate scale used to convert from pixels to arcsec (units = ''/pix).
        Note that this will be modified by the BINFAC keyword at the top of
        the stats table. 
    """
    data = Table.read(stats_file)

    filters = np.array(data['FILTER'])
    bin_fac = np.array(data['BINFAC'])
    N_stars = np.array(data['N Stars'])
    beta = np.array(data['Beta'])
    beta_std = np.array(data['Beta std'])
    alpha_min = np.array(data['Minor Alpha'])
    alpha_min_std = np.array(data['Minor Alpha std'])
    alpha_maj = np.array(data['Major Alpha'])
    alpha_maj_std = np.array(data['Major Alpha std'])

    # Calculate calibration factors
    calib = plate_scale * bin_fac

    if filt==True:
        wvs = plot_stats.filter2wv(filters)
        calib *= ((wvs/filt)**(1/5))

    #TEMPORARILY REMOVED BIN AN PLATE SCALE CALIB
    FWHM_min = 2 * alpha_min * np.sqrt((2**(1/beta))-1) * calib
    FWHM_maj = 2 * alpha_maj * np.sqrt((2**(1/beta))-1) * calib
    
    # Calculate uncertainties of median parameters
    sig_alpha_min = (alpha_min_std / np.sqrt(N_stars)) * np.sqrt((np.pi * N_stars) / (2 * (N_stars-1)))
    sig_alpha_maj = (alpha_maj_std / np.sqrt(N_stars)) * np.sqrt((np.pi * N_stars) / (2 * (N_stars-1)))
    sig_beta = (beta_std / np.sqrt(N_stars)) * np.sqrt((np.pi * N_stars) / (2 * (N_stars-1)))

    # Calculate partial derivatives for error propogation
    del_alpha = 2 * np.sqrt((2**(1/beta))-1)
    del_beta_min = alpha_min * (((2**(1/beta))-1)**(-0.5)) * (-np.log(2)*(2**(1/beta))/(beta**2))
    del_beta_maj = alpha_maj * (((2**(1/beta))-1)**(-0.5)) * (-np.log(2)*(2**(1/beta))/(beta**2))
    
    # Calculate uncertainties of FWHMs
    sig_FWHM_min = np.sqrt(((del_alpha*sig_alpha_min)**2)+((del_beta_min*sig_beta)**2)) * calib
    sig_FWHM_maj = np.sqrt(((del_alpha*sig_alpha_maj)**2)+((del_beta_maj*sig_beta)**2)) * calib
    
    return FWHM_min, sig_FWHM_min, FWHM_maj, sig_FWHM_maj


def rm_mof(stats_files):
    """
    Removes moffat fit columns from stats files
    Input: List of stats files
    Output: overwrites same file in same location w/o moffat cols
    """
    for file in stats_files:
        dat = Table.read(file)
        cols = dat.colnames
        if 'N Stars' in cols:
            dat.remove_columns(['N Stars', 'N Sky', 'N Sky std', 'Amplitude', 'Amplitude std', 'Phi', 'Phi std', \
                            'Beta', 'Beta std', 'Minor Alpha', 'Minor Alpha std', 'Major Alpha', 'Major Alpha std'])
            dat.write(file, overwrite=True)
        else:
            pass
    returnn


def combine_table(stats_file, filt=False, plate_scale=0.016):
    '''
    given stats_file with moffat parameters, calculates
    and adds minor and major moffat FWHM and errors
    '''
    
    dat = Table.read(stats_file)
    FWHM_min, sig_FWHM_min, FWHM_maj, sig_FWHM_maj = calc_mof_fwhm(stats_file, filt=filt, plate_scale=plate_scale)
    dat.add_columns([Column(FWHM_min), Column(sig_FWHM_min), Column(FWHM_maj), Column(sig_FWHM_maj)], \
                   names=['FWHM_min', 'sig_FWHM_min', 'FWHM_maj', 'sig_FWHM_maj'])
    return dat
