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
import os
import numpy as np
from astropy.io import fits
from astropy.table import Table, Column
from astropy.modeling import models, fitting, Fittable2DModel, Parameter
from astropy.modeling.models import custom_model
from imaka.analysis import plot_stats
import matplotlib.pyplot as plt


FLOAT_EPSILON = float(np.finfo(np.float32).tiny)

# @custom_model
# def Elliptical_Moffat2D(x, y, \
#                         N_sky = 0., amplitude = 1., phi=0., power = 1.,\
#                         x_0 = 0., y_0 = 0., width_x = 1., width_y = 1.):
#     """
#     A custom astropy model for a two dimensional elliptical moffat function.  
#     N_sky: a constant background value
#     Amplitude: A
#     phi: rotation angle (in radians?)
#     power: slope of model (beta)
#     x_0, y_0: star's centroid coordinates
#     width_x, width_y: core widths (alpha)
#     """

#     width_x = np.abs(width_x)
#     width_y = np.abs(width_y)
    
#     c = np.cos(phi)
#     s = np.sin(phi)
#     A = (c / width_x)** 2 + (s / width_y)**2
#     B = (s / width_x)** 2 + (c / width_y)**2
#     C = 2 * s * c * (1/width_x**2 - 1/width_y**2)
#     denom = (1 + A * (x-x_0)**2 + B * (y-y_0)**2 + C*(x-x_0)*(y-y_0))**power

#     return N_sky + amplitude / denom


class Elliptical_Moffat2D(Fittable2DModel):
    r"""
    A custom astropy model for a two dimensional elliptical moffat function.  

    Parameters
    -----------
    N_sky : float
        a constant background value
    Amplitude : float
        surface brightness
    phi : float
        rotation angle (in radians)
    power : float
        slope of model (beta)
    x_0 : float
        x centroid of star
    y_0 : float
        y centroid of star
    width_x : float
        core width (alpha)
    width_y : float
        core width (alpha)

    See Also
    --------
    Gaussian2D, Moffat2D

    """

    N_sky = Parameter(default=0)
    amplitude = Parameter(default=1)
    phi = Parameter(default=0)
    power = Parameter(default=1)
    x_0 = Parameter(default=0)
    y_0 = Parameter(default=0)
    width_x = Parameter(default=1)
    width_y = Parameter(default=1)

    def __init__(self, N_sky=N_sky.default, amplitude=amplitude.default,
                     phi=phi.default, power=power.default, x_0=x_0.default, y_0=y_0.default,
                     width_x=width_x.default, width_y=width_y.default, **kwargs):

        # Ensure widths makes sense if its bounds are not explicitly set.
        # widths must be non-zero and positive.
        # TODO: Investigate why setting this in Parameter above causes
        #       convolution tests to hang.
        kwargs.setdefault('bounds', {})
        kwargs['bounds'].setdefault('width_x', (0, None))
        kwargs['bounds'].setdefault('width_y', (0, None))
        #kwargs['bounds'].setdefault('phi', (-2*3.14, 2.*3.14))

        super().__init__(
            N_sky=N_sky, amplitude=amplitude, phi=phi, power=power,
            x_0=x_0, y_0=y_0, width_x=width_x, width_y=width_y, **kwargs)
        
    
    @classmethod
    def evaluate(cls, x, y, N_sky, amplitude, phi, power, x_0, y_0, width_x, width_y):
        """Two dimensional Elliptical Moffat profile function at an 
        arbitrary angle."""
    
        cosp = np.cos(phi)
        sinp = np.sin(phi)
        A = (cosp / width_x)** 2 + (sinp / width_y)**2
        B = (sinp / width_x)** 2 + (cosp / width_y)**2
        C = 2 * sinp * cosp * (1/width_x**2 - 1/width_y**2)
        denom = (1 + A * (x-x_0)**2 + B * (y-y_0)**2 + C*(x-x_0)*(y-y_0))**power

        return N_sky + (amplitude / denom)

    @property
    def input_units(self):
        if self.x_0.unit is None:
            return None
        else:
            return {'x': self.x_0.unit,
                    'y': self.y_0.unit}

    def _parameter_units_for_data_units(self, inputs_unit, outputs_unit):
        # Note that here we need to make sure that x and y are in the same
        # units otherwise this can lead to issues since rotation is not well
        # defined.
        if inputs_unit['x'] != inputs_unit['y']:
            raise UnitsError("Units of 'x' and 'y' inputs should match")
        return OrderedDict([('x_0', inputs_unit['x']),
                            ('y_0', inputs_unit['x']),
                            ('width_x', inputs_unit['x']),
                            ('width_y', inputs_unit['y']),
                            ('phi', u.rad),
                            ('amplitude', outputs_unit['z'])])

    

def find_stars(img_file, fwhm=5, threshold=4, N_passes=2, plot=False):
    """
    img_file - an image file.
    fwhm - First guess at the FWHM of the PSF for the first pass on the first image.
    threshold - the SNR threshold (mean + threshold*std) above which to search for sources.
    N_passes - how many times to find sources, recalc the PSF FWHM, and find sources again. 
    """
    print('\nREDUCE_FLI: find_stars()')

    # Unbinned plate scale
    ps = 0.04 # arcsecs / pix

    img, hdr = fits.getdata(img_file, header=True)
    binfac = hdr['BINFAC']
    ps = ps * binfac

    fwhm_curr = fwhm
    N = 0
    
    # Calculate the bacgkround and noise (iteratively)
    print("    Calculating background")
    bkg_threshold_above = 1
    bkg_threshold_below = 3
    for nn in range(5):
        if nn == 0:
            bkg_mean = img.mean()
            bkg_std = img.std()
        else:
            bkg_mean = img[good_pix].mean()
            bkg_std = img[good_pix].std()

        bad_hi = bkg_mean + (bkg_threshold_above * bkg_std)
        bad_lo = bkg_mean - (bkg_threshold_below * bkg_std)

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
        print('     Pass {0:d} assuming FWHM = {1:.1f} pix'.format(nn, fwhm_curr))
        daofind = DAOStarFinder(fwhm=fwhm_curr, threshold = img_threshold, exclude_border=True)
        sources = daofind(img - bkg_mean)

        # Calculate FWHM for each detected star.
        min_fwhm = np.zeros(len(sources), dtype=float)
        maj_fwhm = np.zeros(len(sources), dtype=float)
        phi = np.zeros(len(sources), dtype=float)
        elon = np.zeros(len(sources), dtype=float)

        cutout_half_size = int(round(fwhm_curr * 3))
        cutout_size = 2 * cutout_half_size
        cutouts = np.zeros((len(sources), cutout_size, cutout_size), dtype=float)

        alpha_init_guess = fwhm_curr / 0.87 #fwhm with beta=4
        m_init = Elliptical_Moffat2D(N_sky = 0, amplitude=1.,
                                         x_0=cutout_half_size, y_0=cutout_half_size, 
                                         width_x = alpha_init_guess, width_y=alpha_init_guess)
        
        fit_m = fitting.LevMarLSQFitter()
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

            cutouts[ss] = cutout_tmp # - bkg_mean ###Might be a source of issue in reduced images
            cutouts[ss] /= cutouts[ss].sum()

            # Fit an elliptical moffat to the cutout image.
            mof_params = fit_m(m_init, cut_x, cut_y, cutouts[ss])
    
            phi[ss] = mof_params.phi.value % (2*np.pi)
            
            if mof_params.width_x.value < mof_params.width_y.value:
                x = mof_params.width_x.value
                y = mof_params.width_y.value
            else:
                x = mof_params.width_y.value
                y = mof_params.width_x.value
                phi[ss] += np.pi/2
            beta = mof_params.power.value
            
            min_fwhm[ss] = 2*x*np.sqrt((2**(1/beta))-1)
            maj_fwhm[ss] = 2*y*np.sqrt((2**(1/beta))-1)
            elon[ss] = maj_fwhm[ss] / min_fwhm[ss]
            N += 1

        # Drop sources with flux (signifiance) that isn't good enough.
        # Empirically this is <1.2
        #good = np.where(sources['flux'] > 1.9)[0]
        #sources = sources[good]

        # Only use the brightest sources for calculating the mean. This is just for printing.
        min_fwhm_mean, min_fwhm_med, min_fwhm_std = sigma_clipped_stats(min_fwhm)
        maj_fwhm_mean, maj_fwhm_med, maj_fwhm_std = sigma_clipped_stats(maj_fwhm)
        elon_mean, elon_med, elon_std = sigma_clipped_stats(elon)
        phi_mean, phi_med, phi_std = sigma_clipped_stats(phi)

        print('        Number of sources = ', len(sources))
        print('        Minor fwhm = {0:.2f} +/- {1:.2f} arcsec'.format(min_fwhm_med*ps,
                                                                 min_fwhm_std*ps))
        print('        Major fwhm = {0:.2f} +/- {1:.2f} arcsec'.format(maj_fwhm_med*ps,
                                                                 maj_fwhm_std*ps))
        print('        Elongation = {0:.2f} +/- {1:.2f}'.format(elon_med,
                                                                 elon_std))
        print('        Pos. Angle = {0:.2f} +/- {1:.2f} rad'.format(phi_med,
                                                                 phi_std))

        fwhm_curr = np.mean([min_fwhm_med, maj_fwhm_med])

        
    if plot == True:
        plt.figure(figsize=(8,6))
        plt.clf()
        plt.subplots_adjust(left=0.13)
        vmin = (min_fwhm_med - 3*min_fwhm_std) * ps
        vmax = (min_fwhm_med + 3*min_fwhm_std) * ps
        print(vmin, vmax)

        # Drop outliers:
        idx = np.where((min_fwhm*ps >= vmin) & (min_fwhm*ps <= vmax))
        plt.scatter(sources['xcentroid'][idx], sources['ycentroid'][idx],
                        c=min_fwhm[idx]*ps, s=50, edgecolor='none',
                        vmin=vmin, vmax=vmax)
        plt.xlabel('X (pixels)')
        plt.ylabel('Y (pixels)')
        plt.title(img_file)
        plt.colorbar(label='FWHM (")')
        plt.axis('equal')
        

    return min_fwhm_med*ps, min_fwhm_std*ps, maj_fwhm_med*ps, maj_fwhm_std*ps, elon_med*ps, elon_std*ps
    
