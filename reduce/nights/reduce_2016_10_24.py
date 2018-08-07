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
from astropy.convolution import Gaussian2DKernel
from astropy.stats import gaussian_fwhm_to_sigma
from photutils import Background, detect_threshold, detect_sources
from photutils import source_properties, properties_table
import poppy
import pdb

scale = 23.0 # mas/pixel


def load_image(image_file):
    if image_file.endswith('fits'):
        arr = fits.getdata(image_file)
    else:
        img = Image.open(image_file)
        arr = np.array(img, dtype=np.int64)

    return arr

def convert_to_fits(image_file, clobber=True):
    """Convert from a PGM file to a FITS file.
    """
    img = load_image(image_file)

    fits.writeto(image_file.replace('pgm', 'fits'), img, clobber=clobber)

    return

def convert_all_to_fits(directory, clobber=True):
    images = glob.glob(directory + "*.pgm")

    for img in images:
        convert_to_fits(img, clobber=clobber)

    return

def calc_fwhm(image_file):
    stddev = 1.0 # pixels
    round_sigma = (stddev * stddev)**0.5
    psf_guess = psf.IntegratedGaussianPRF(flux=1, sigma=round_sigma)
    psf_guess.flux.fixed = psf_guess.x_0.fixed = psf_guess.y_0.fixed = False
    psf_guess.x_0.sigma = True

    img = load_image(image_file)

    outtabi = psf.psf_photometry(img, intab, psf_guess, fitshape)

    pdb.set_trace()

    return

def calc_fwhm_on_bright_star(image_file, print=True, fwhm_init=2.0):
    """Calculate the FWHM on a single bright star (either open or closed loops).

    image_file -- either FITS or PGM
    fwhm_init -- (def=2) pixels for FWHM initial guess

    """
    
    img = load_image(image_file)
    
    # Calculate the bacgkround
    bkg = photutils.Background(img, img.shape, filter_shape=(1,1), method='median')

    threshold = bkg.background + (30.0 * bkg.background_rms)

    sigma = 2.0 * gaussian_fwhm_to_sigma    # FWHM = 2. pixels
    kernel = Gaussian2DKernel(sigma, x_size=3, y_size=3)
    kernel.normalize()
    segm = detect_sources(img, threshold, npixels=5, filter_kernel=kernel)    

    props = source_properties(img, segm)
    tbl = properties_table(props)

    # Check for junk stars (cosmic rays)
    idx = np.where((tbl['semimajor_axis_sigma'] > 1) & (tbl['semiminor_axis_sigma'] > 1))[0]
    tbl = tbl[idx]
    
    tbl['image_name'] = image_file

    if print == True:
        reformat_source_table(tbl)
        print_source_table(tbl)
        
    return tbl


    
def calc_fwhm_on_bright_star_stack(directory, wildcard="*.pgm"):
    
    images = glob.glob(directory + wildcard)

    tbl_final = None

    for ii in range(len(images)):
        image_file = images[ii]

        tbl = calc_fwhm_on_bright_star(image_file, print=False)
        
        if ii == 0:
            tbl_final = tbl
        else:
            tbl_final = table.vstack([tbl_final, tbl])

    reformat_source_table(tbl_final)
    print_source_table(tbl_final)

    return


def reformat_source_table(tbl_final):
    # Add a few columns
    tbl_final['fwhm_major'] = 2.355 * tbl_final['semimajor_axis_sigma']
    tbl_final['fwhm_minor'] = 2.355 * tbl_final['semiminor_axis_sigma']
    tbl_final['fwhm_major_mas'] = 2.355 * tbl_final['semimajor_axis_sigma'] * scale
    tbl_final['fwhm_minor_mas'] = 2.355 * tbl_final['semiminor_axis_sigma'] * scale
        
    # Set some column formatting:
    tbl_final['xcentroid'].format = '%7.3f'
    tbl_final['ycentroid'].format = '%7.3f'
    tbl_final['semimajor_axis_sigma'].format = '%6.3f'
    tbl_final['semiminor_axis_sigma'].format = '%6.3f'
    tbl_final['fwhm_major'].format = '%6.3f'
    tbl_final['fwhm_minor'].format = '%6.3f'
    tbl_final['fwhm_major_mas'].format = '%6.0f'
    tbl_final['fwhm_minor_mas'].format = '%6.0f'
    tbl_final['eccentricity'].format = '%6.3f'
    tbl_final['orientation'].format = '%6.3f'
    tbl_final['ellipticity'].format = '%6.3f'

    # Convert orientation to degrees.
    tbl_final['orientation'] = np.degrees(tbl_final['orientation'])
    tbl_final['orientation'].unit = units.degree
    tbl_final['fwhm_major_mas'].unit = 'mas'
    tbl_final['fwhm_minor_mas'].unit = 'mas'

    return tbl_final


def print_source_table(tbl_final):
    print(tbl_final['id','xcentroid','ycentroid','source_sum',
                        'fwhm_major', 'fwhm_minor',
                        'fwhm_major_mas', 'fwhm_minor_mas',
                        'orientation','ellipticity','image_name'])

    if (len(tbl_final) > 1):
        print('Average +/- Standard Deviation over stack:')
        hdr = '{0:3s}  {1:15s}   {2:15s}  '
        hdr += '{3:9s}  {4:9s}  '
        hdr += '{5:14s}  {6:14s}'
        print(hdr.format('ID', 'xcentroid', 'ycentroid', 'FWHM_major', ' FWHM_minor',
                        'Orietnation', 'Ellipticity'))
        fmt = '{0:3.0f}  {1:6.3f} +- {2:5.3f}   {3:6.3f} +- {4:5.3f}  '
        fmt += '{5:3.0f} +/- {6:2.0f}  {7:3.0f} +/- {8:2.0f}  '
        fmt += '{9:6.3f} +/- {10:4.3f}  {11:6.3f} +/- {12:4.3f}'
        print(fmt.format(tbl_final['id'].mean(),
                        tbl_final['xcentroid'].mean(), tbl_final['xcentroid'].std(), 
                        tbl_final['ycentroid'].mean(), tbl_final['ycentroid'].std(),
                        tbl_final['fwhm_major_mas'].mean(), tbl_final['fwhm_major_mas'].std(),
                        tbl_final['fwhm_minor_mas'].mean(), tbl_final['fwhm_minor_mas'].std(),
                        tbl_final['orientation'].mean(), tbl_final['orientation'].std(),
                        tbl_final['ellipticity'].mean(), tbl_final['ellipticity'].std()))
    
    return tbl_final

def stack_images(directory, wildcard="*.pgm", clobber=True):
    images = glob.glob(directory + wildcard)

    N_images = len(images)

    # Calculate the weight per image. 
    weight = 1.0 / N_images

    for ii in range(N_images):
        image_file = images[ii]
        img = load_image(image_file)

        if ii == 0:
            final_image = img * weight
        else:
            final_image += img * weight


    tmp = images[0].split('-')
    filename = '_'.join(tmp[0:-2]) + '.fits'

    # Chop the outer 2 pixels around the image... bad.
    final_image = final_image[2:-2, 2:-2]
    fits.writeto(filename, final_image, clobber=clobber)
    fits.setval(filename, 'PIXELSCL', value=1.0)
    
def compare_curve_of_growth_open_closed(open_loop_image, closed_loop_image):
    otab = calc_fwhm_on_bright_star(open_loop_image, print=False, fwhm_init=20)
    ctab = calc_fwhm_on_bright_star(closed_loop_image, print=False, fwhm_init=2)

    xpix_o = int(np.round(otab['xcentroid']))
    ypix_o = int(np.round(otab['ycentroid']))
    xpix_c = int(np.round(ctab['xcentroid']))
    ypix_c = int(np.round(ctab['ycentroid']))

    # Read in the images and background subtract
    hdu_o = fits.open(open_loop_image)
    hdu_c = fits.open(closed_loop_image)
    img_o = hdu_o[0].data.copy()
    img_c = hdu_c[0].data.copy()

    bkg_o = photutils.Background(img_o, img_o.shape, filter_shape=(1,1), method='median')
    bkg_c = photutils.Background(img_c, img_c.shape, filter_shape=(1,1), method='median')

    hdu_o[0].data -= bkg_o.background
    hdu_c[0].data -= bkg_c.background

    # Trim down to the central
    ihalf = 75
    xlo_o = xpix_o - ihalf
    xhi_o = xpix_o + ihalf
    ylo_o = ypix_o - ihalf
    yhi_o = ypix_o + ihalf
    hdu_o[0].data = hdu_o[0].data[ylo_o:yhi_o, xlo_o:xhi_o]
    xpix_o -= xlo_o
    ypix_o -= ylo_o

    ihalf = 75
    xlo_c = xpix_c - ihalf
    xhi_c = xpix_c + ihalf
    ylo_c = ypix_c - ihalf
    yhi_c = ypix_c + ihalf
    hdu_c[0].data = hdu_c[0].data[ylo_c:yhi_c, xlo_c:xhi_c]
    xpix_c -= xlo_c
    ypix_c -= ylo_c
    plt.figure(3)
    plt.imshow(hdu_o[0].data)
    
    
    rr_o, prof_o, ee_o = poppy.utils.radial_profile(hdu_o, EE=True,
                                                        center=[xpix_o, ypix_o],
                                                        normalize='total')
    rr_c, prof_c, ee_c = poppy.utils.radial_profile(hdu_c, EE=True,
                                                        center=[xpix_c, ypix_c],
                                                        normalize='total')
    # Plot the radial profiles.
    plt.figure(1)
    plt.clf()
    plt.plot(rr_o, prof_o, label='Open')
    plt.plot(rr_c, prof_c, label='Closed')
    plt.legend()
    plt.xlim(0, 75)
    plt.show()

    plt.figure(2)
    plt.clf()
    plt.plot(rr_o, ee_o, label='Open')
    plt.plot(rr_c, ee_c, label='Closed')
    plt.legend(loc='upper left')
    plt.xlim(0, 75)
    plt.ylim(0, 1)
    plt.show()
    
    
    
