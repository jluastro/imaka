### moffat.py - Takes reduced images and existing stats files and conducts moffat fitting, adds data to stats tables and makes model PSFs

import os
import numpy as np
from astropy.io import fits
from astropy.table import Table, Column
from astropy.modeling import models, fitting
from astropy.modeling.models import custom_model

import matplotlib.pyplot as plt
@custom_model





def Elliptical_Moffat2D(x, y, \
                        N_sky = 0., amplitude = 1., phi=0., power = 1.,\
                        x_0 = 0., y_0 = 0., width_x = 1., width_y = 1.):
    """
    A custom astropy model for a two dimensional elliptical moffat function.  
    N_sky: a constant background value
    Amplitude: A
    phi: rotation angle (in radians?)
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


def fit_moffat(img_files, stats_file, x_guess=5, y_guess=5, ):
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
    
    for ii in range(N_files):
        # Load up the image to work on.
        print("Working on image: ", img_files[ii])
        img, hdr = fits.getdata(img_files[ii], header=True)

        # Load up the corresponding starlist.
        starlist = img_files[ii].replace('.fits', '_stars.txt')
        stars = Table.read(starlist, format='ascii')
        N_stars = len(stars)

        # Put the positions into an array 
        coords = np.array([stars['xcentroid'], stars['ycentroid']])

        # Make empty arrays for trimmed star sample
        x_cents = []
        y_cents = []
        N_good_stars = 0
        
        for jj in range(N_stars):
            # Trim star sample for edge sources, saturation, and low flux
            cut_size = 40
            x_cent = int(round(float(coords[0][jj])))
            y_cent = int(round(float(coords[1][jj])))
            flux   = np.array(stars['flux'])[jj]
            peak   = np.array(stars['peak'])[jj]
            if y_cent - cut_size/2 > 0 and x_cent - cut_size/2 > 0 and y_cent + cut_size/2 < np.shape(img)[0] and x_cent + cut_size/2<np.shape(img)[1] \
              and flux > 30 and peak < 20000:
                N_good_stars += 1
                x_cents.append(x_cent)
                y_cents.append(y_cent)

        print("Using", N_good_stars, "of", N_stars, "stars in calculation")
        
        # Create lists for each image's stars' parameters
        N_sky_list     = np.zeros(N_good_stars, dtype=float)
        amplitude_list = np.zeros(N_good_stars, dtype=float)
        phi_list       = np.zeros(N_good_stars, dtype=float)
        power_list     = np.zeros(N_good_stars, dtype=float)
        x_0_list       = np.zeros(N_good_stars, dtype=float)
        y_0_list       = np.zeros(N_good_stars, dtype=float)
        width_x_list   = np.zeros(N_good_stars, dtype=float)
        width_y_list   = np.zeros(N_good_stars, dtype=float)

        final_psf_mof = np.zeros((cut_size, cut_size), dtype=float)
        
        for jj in range(N_good_stars):
            # Make image cut for good star
            x_cent = x_cents[jj]
            y_cent = y_cents[jj]
            image_cut = img[int(y_cent-cut_size*0.5):int(y_cent+cut_size*0.5), int(x_cent-cut_size*0.5):int(x_cent+cut_size*0.5)]

            # Conduct moffat fit
            y, x = np.mgrid[:cut_size, :cut_size]
            z = image_cut
            m_init = Elliptical_Moffat2D(N_sky = 0, amplitude=np.amax(z),  x_0=cut_size/2, y_0=cut_size/2, width_x = 4.55, width_y=4.17)
            fit_m = fitting.LevMarLSQFitter()
            m = fit_m(m_init, x, y, z)
                 
            N_sky_list[jj]     = m.N_sky.value
            amplitude_list[jj] = m.amplitude.value
            phi_list[jj]       = m.phi.value
            power_list[jj]     = m.power.value
            x_0_list[jj]       = m.x_0.value
            y_0_list[jj]       = m.y_0.value
            if m.width_x.value < m.width_y.value:
                width_x_list[jj]    = m.width_x.value
                width_y_list[jj]    = m.width_y.value
            else:
                width_x_list[jj]    = m.width_y.value
                width_y_list[jj]    = m.width_x.value

            final_psf_mof += m(x, y)

        
        # Take median values of parameters and put in intial lists
        N_sky[ii]     = np.median(N_sky_list)
        amplitude[ii] = np.median(amplitude_list)
        phi[ii]       = np.median(phi_list)
        power[ii]     = np.median(power_list)
        x_0[ii]       = np.median(x_0_list)
        y_0[ii]       = np.median(y_0_list)
        width_x[ii]   = np.median(width_x_list)
        width_y[ii]   = np.median(width_y_list)

        N_sky_std[ii]     = np.std(N_sky_list)
        amplitude_std[ii] = np.std(amplitude_list)
        phi_std[ii]       = np.std(phi_list)
        power_std[ii]     = np.std(power_list)
        x_0_std[ii]       = np.std(x_0_list)
        y_0_std[ii]       = np.std(y_0_list)
        width_x_std[ii]   = np.std(width_x_list)
        width_y_std[ii]   = np.std(width_y_list)

        mof_stars[ii] = int(N_good_stars)
        # Save the average PSF (flux-weighted). Note we are making a slight 
        # mistake here since each PSF has a different sub-pixel position

        final_psf_mof /= N_good_stars
        fits.writeto(img_files[ii].replace('.fits', '_psf_mof.fits'), final_psf_mof, hdr, clobber=True)

    # Read in existing stats table and append new data
    stats = Table.read(stats_file)
    
    col_mof_stars = Column(name='N Stars', data = mof_stars)
    col_N_sky = Column(name='N Sky', data=N_sky)
    col_N_sky_std = Column(name='N Sky std', data=N_sky_std)
    col_amplitude = Column(name='Amplitude', data=amplitude)
    col_amplitude_std = Column(name='Amplitude std', data=amplitude_std)
    col_phi = Column(name='Phi', data=phi)
    col_phi_std = Column(name='Phi std', data=phi_std)
    col_beta = Column(name='Beta', data = power)
    col_beta_std = Column(name='Beta std', data=power_std)
    col_alpha_min = Column(name='Minor Alpha', data = width_x)
    col_alpha_min_std = Column(name='Minor Alpha std', data = width_x_std)
    col_alpha_maj = Column(name='Major Alpha', data = width_y)
    col_alpha_maj_std = Column(name='Major Alpha std', data = width_y_std)

    stats.add_columns([col_mof_stars, col_N_sky, col_N_sky_std, col_amplitude,\
                           col_amplitude_std, col_phi, col_phi_std, col_beta, \
                           col_beta_std, col_alpha_min, col_alpha_min_std, \
                           col_alpha_maj, col_alpha_maj_std])

    stats_file_root, stats_file_ext = os.path.splitext(stats_file)
    stats.write(stats_file_root.split("_mdp")[0] + stats_file_ext, overwrite=True)
    
    return
