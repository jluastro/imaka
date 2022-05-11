import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from astropy.modeling import models, fitting
from scipy.ndimage import interpolation, rotate
from math import radians, degrees
from astropy.modeling.models import custom_model
from astropy.stats import gaussian_fwhm_to_sigma
from matplotlib.patches import Rectangle
import pdb

@custom_model
def Elliptical_Moffat2D(x, y, \
                        amplitude = 1., phi=0., power = 1.,\
                        x_0 = 0., y_0 = 0., width_x = 1., width_y = 1.):
    c = np.cos(phi)
    s = np.sin(phi)
    A = (c / width_x) ** 2 + (s / width_y)**2
    B = (s / width_x) ** 2 + (c/ width_y)**2
    C = 2 * s * c * (1/width_x**2 - 1/width_y**2)
    denom = (1 + A * (x-x_0)**2 + B * (y-y_0)**2 + C*(x-x_0)*(y-y_0))**power
    return amplitude / denom



def plot_samples(open_file, x_o, y_o, closed_file, x_c, y_c, cutout_size):
    
    open_img = fits.getdata(open_file); closed_img = fits.getdata(closed_file)

    # Make open cutout
    x_cent = x_o; y_cent = y_o
    x_lo = int(round(x_cent - cutout_size/2)); x_hi = x_lo + cutout_size
    y_lo = int(round(y_cent - cutout_size/2)); y_hi = y_lo + cutout_size
    cutout_tmp = open_img[y_lo:y_hi, x_lo:x_hi].astype(float)
    open_cut = cutout_tmp
    mid_cut = open_cut.shape/2
    open_cut /= open_cut.max()
    #open_cut /= open_cut[mid_cut[0], mid_cut[1]]
    
    #Make closed cutout
    x_cent = x_c; y_cent = y_c
    x_lo = int(round(x_cent - cutout_size/2)); x_hi = x_lo + cutout_size
    y_lo = int(round(y_cent - cutout_size/2)); y_hi = y_lo + cutout_size
    cutout_tmp = closed_img[y_lo:y_hi, x_lo:x_hi].astype(float)
    closed_cut = cutout_tmp
    mid_cut = closed_cut.shape/2
    closed_cut /= closed_cut.max()
    #closed_cut /= open_cut[mid_cut[0], mid_cut[1]]
    
    # Plot sources
    plt.figure(1, figsize=(8, 4)); plt.suptitle("Original Sample Sources", fontsize=20)
    plt.subplot(1,2,1); plt.imshow(open_cut, origin='lower', interpolation='nearest'); plt.title("Open Loop Sample Source")
    plt.subplot(1,2,2); plt.imshow(closed_cut, origin='lower', interpolation='nearest'); plt.title("Closed Loop Sample Source")
    
    return open_cut, closed_cut

  
def model_fit(model, image, amp, x_wid=5, y_wid=5, angle=0, fignumber=1, title=None):

    cutout_size = np.shape(image)[0]
    # Open Fit
    y_o, x_o = np.mgrid[:cutout_size, :cutout_size]
    z_o = image
    
    if model == "Gaussian":
        p_init_o = models.Gaussian2D(np.amax(z_o), cutout_size/2, cutout_size/2, x_stddev = x_wid, y_stddev=y_wid, theta=angle)
    elif model == "Moffat":
        p_init_o = Elliptical_Moffat2D(amplitude=np.amax(z_o), x_0=cutout_size/2, y_0=cutout_size/2, width_x = x_wid, width_y=y_wid)
    fit_p_o = fitting.LevMarLSQFitter()
    p_o = fit_p_o(p_init_o, x_o, y_o, z_o)
    residual_o = np.sum((z_o - p_o(x_o, y_o))**2)
    PSF_mean = np.mean(z_o)
    diff = z_o - PSF_mean
    FUV = residual_o / np.sum((diff)**2)

    return p_o



def model_plot(model, image, amp, x_wid=5, y_wid=5, angle=0, fignumber=1, title=None):

    cutout_size = np.shape(image)[0]
    # Open Fit
    y_o, x_o = np.mgrid[:cutout_size, :cutout_size]
    z_o = image
    
    if model == "Gaussian":
        g2d_model = models.Gaussian2D(np.amax(z_o), cutout_size/2, cutout_size/2,
                                         x_stddev = x_wid, y_stddev=y_wid, theta=angle)
        c2d_model = models.Const2D(0.0)
        p_init_o = g2d_model + c2d_model
    elif model == "Moffat":
        g2d_model = Elliptical_Moffat2D(amplitude=np.amax(z_o), x_0=cutout_size/2,
                                           y_0=cutout_size/2, width_x = x_wid, width_y=y_wid)
        c2d_model = models.Const2D(0.0)
        p_init_o = g2d_model + c2d_model

    fit_p_o = fitting.LevMarLSQFitter()
    # weights = 1.0 / z_o**0.5
    # weights = z_o**0.5
    # weights = np.ones(z_o.shape, dtype=float)
    p_o = fit_p_o(p_init_o, x_o, y_o, z_o,
                      epsilon=1e-12, acc=1e-12, maxiter=300,
                      weights=None)
    
    # FWHM values?
    if model == "Gaussian":
        x_fwhm = p_o.x_stddev_0.value / gaussian_fwhm_to_sigma
        y_fwhm = p_o.y_stddev_0.value / gaussian_fwhm_to_sigma 
    elif model == "Moffat":
        # BUG this isn't quite right for a moffat fit...
        x_fwhm =  p_o.width_x_0.value 
        y_fwhm =  p_o.width_y_0.value 
    
    #Images
    model_img = p_o(x_o, y_o)
    diff_img = z_o - model_img
    
    #Images cut:
    half_size = image.shape[0]/2
    x_lo = int(round(half_size - x_fwhm))
    x_hi = x_lo + int(round(x_fwhm*2))
    y_lo = int(round(half_size - y_fwhm))
    y_hi = y_lo + int(round(y_fwhm*2))
    #cut size for metrics
    z_o_cut = z_o[y_lo:y_hi, x_lo:x_hi].astype(float)
    diff_img_cut = diff_img[y_lo:y_hi, x_lo:x_hi].astype(float)
    
    # Fit metrics
    residual_o = np.sum(diff_img_cut**2)
    PSF_mean = np.mean(z_o_cut)
    fresid = np.median(np.abs(diff_img_cut / z_o_cut)) # mean fractional residual
    FVU = np.sum(diff_img**2) / np.sum((z_o_cut - PSF_mean)**2)  # fraction of variance unexplained
    
    # Plotting
    plt.figure(fignumber, figsize=(12, 3)); plt.suptitle(title)

    vmin = z_o.min()
    vmax = z_o.max()

    # left panel
    plt.subplot(1,4,1)
    plt.imshow(z_o, origin='lower',
                   vmin=vmin, vmax=vmax)
    plt.gca().add_patch(Rectangle((x_lo, y_lo),x_hi-x_lo,y_hi-y_lo,
                    edgecolor='red',
                    facecolor='none',
                    lw=2))
    plt.colorbar(fraction=0.046, pad=0.05)
    plt.title("Data")
    
    # left middle panel
    plt.subplot(1,4,2)
    plt.imshow(model_img, origin='lower',
                   vmin=vmin, vmax=vmax)
    plt.gca().add_patch(Rectangle((x_lo, y_lo),x_hi-x_lo,y_hi-y_lo,
                    edgecolor='red',
                    facecolor='none',
                    lw=2))
    plt.colorbar(fraction=0.046, pad=0.04); plt.title("Model")
    
    # right middle panel
    plt.subplot(1,4,3)
    plt.imshow(z_o - model_img, origin='lower',
                   vmin=-vmax/6, vmax=vmax/6)
    plt.gca().add_patch(Rectangle((x_lo, y_lo),x_hi-x_lo,y_hi-y_lo,
                    edgecolor='red',
                    facecolor='none',
                    lw=1))
    plt.colorbar(fraction=0.046, pad=0.04)
    plt.title("Data-Model")

    # right panel
    plt.subplot(1,4,4)
    plt.imshow((z_o - model_img) / z_o,
                   origin='lower',
                   vmax=1, vmin=-1)
    plt.gca().add_patch(Rectangle((x_lo, y_lo),x_hi-x_lo,y_hi-y_lo,
                    edgecolor='red',
                    facecolor='none',
                    lw=1))
    plt.colorbar(fraction=0.046, pad=0.04)
    plt.title("Fract. Resid.")
    
    plt.tight_layout()

    print("Residuals Squared Summed               - ", title, ": ",'{:.2e}'.format(residual_o))
    print("Fraction of Variance Unexplained (FVU) - ", title, ": ",'{:.2e}'.format(FVU))
    print("Median Fractional Residual Per Pix     - ", title, ": ",'{:.2e}'.format(fresid))
    
    if model == "Moffat":
        print("beta: ", p_o.power_0.value)
    print('\nFit Info:')
    print('nfev = ', fit_p_o.fit_info['nfev'])
    print('ierr = ', fit_p_o.fit_info['ierr'])
    print('mesg = ', fit_p_o.fit_info['message'])
    print('\nParams:')

    for k, v in p_o._parameters_.items():
        val = v.value
        unit = 'pix'
        
        if 'phi' in k:
            val = val % 360.0
            unit = 'deg'
        if 'theta' in k:
            val = np.rad2deg(val) % 360.0
            unit = 'deg'
        if 'amp' in k:
            unit = 'flux'
            
        print(f'{k:15s} = {val:10.6f} {unit:5s}')

    print('\n')
        

    return p_o

def model_plot_double(model, image, amp, x_wid_0=5, y_wid_0=5, x_wid_1=5, y_wid_1=5,
                          angle=0, fignumber=1, title=None):

    cutout_size = np.shape(image)[0]
    # Open Fit
    y_o, x_o = np.mgrid[:cutout_size, :cutout_size]
    z_o = image
    
    if model == "Gaussian":
        gaus1_o  = models.Gaussian2D(np.amax(z_o), cutout_size/2, cutout_size/2,
                                         x_stddev=x_wid_0, y_stddev=y_wid_0, theta=angle)
        gaus2_o  = models.Gaussian2D(np.amax(z_o), cutout_size/2, cutout_size/2,
                                         x_stddev=x_wid_1, y_stddev=y_wid_1, theta=angle)
        p_init_o = gaus1_o + gaus2_o
    elif model == "Moffat":
        mof1_o  = Elliptical_Moffat2D(amplitude=np.amax(z_o), x_0=cutout_size/2,
                                          y_0=cutout_size/2, width_x = x_wid_0, width_y=y_wid_0)
        mof2_o  = Elliptical_Moffat2D(amplitude=np.amax(z_o), x_0=cutout_size/2,
                                          y_0=cutout_size/2, width_x = x_wid_1, width_y=y_wid_1)
        p_init_o = mof1_o + mof2_o
    fit_p_o = fitting.LevMarLSQFitter()
    p_o = fit_p_o(p_init_o, x_o, y_o, z_o)
    residual_o = np.sum((z_o - p_o(x_o, y_o))**2)

    PSF_mean = np.mean(z_o)
    diff = z_o - PSF_mean
    FUV = residual_o / np.sum((diff)**2)

    # Plotting
    plt.figure(fignumber, figsize=(12, 2.5)); plt.suptitle(title)
    plt.subplot(1,3,1)
    plt.imshow(z_o, origin='lower', interpolation='nearest')
    plt.colorbar( size="5%")
    plt.title("Data")
    
    plt.subplot(1,3,2)
    plt.imshow(p_o(x_o, y_o), origin='lower', interpolation='nearest')
    plt.colorbar(); plt.title("Model")
    plt.subplot(1,3,3)
    plt.imshow((z_o - p_o(x_o, y_o)) /z_o, origin='lower', interpolation='nearest')
    plt.colorbar( size="5%"); plt.title("Residual")
    
    plt.tight_layout()
    print("Least Squares Sum - ", title, ": ",'{:.2e}'.format(residual_o))
    print("FUV - ", title, ": ",'{:.2e}'.format(FUV))

    return p_o

def model_plot_double_test(model, image, amp, x_wid_0=5, y_wid_0=5, x_wid_1=5, y_wid_1=5,
                          angle=0, fignumber=1, title=None):

    cutout_size = np.shape(image)[0]
    # Open Fit
    y_o, x_o = np.mgrid[:cutout_size, :cutout_size]
    z_o = image
    
    if model == "Gaussian":
        gaus1_o  = models.Gaussian2D(np.amax(z_o), cutout_size/2, cutout_size/2,
                                         x_stddev=x_wid_0, y_stddev=y_wid_0, theta=angle)
        gaus2_o  = models.Gaussian2D(np.amax(z_o), cutout_size/2, cutout_size/2,
                                         x_stddev=x_wid_1, y_stddev=y_wid_1, theta=angle)
        p_init_o = gaus1_o + gaus2_o
    elif model == "Moffat":
        mof1_o  = Elliptical_Moffat2D(amplitude=np.amax(z_o), x_0=cutout_size/2,
                                          y_0=cutout_size/2, width_x = x_wid_0, width_y=y_wid_0)
        mof2_o  = Elliptical_Moffat2D(amplitude=np.amax(z_o), x_0=cutout_size/2,
                                          y_0=cutout_size/2, width_x = x_wid_1, width_y=y_wid_1)
        p_init_o = mof1_o + mof2_o
    fit_p_o = fitting.LevMarLSQFitter()
    p_o = fit_p_o(p_init_o, x_o, y_o, z_o)
    residual_o = np.sum((z_o - p_o(x_o, y_o))**2)

    PSF_mean = np.mean(z_o)
    diff = z_o - PSF_mean
    FUV = residual_o / np.sum((diff)**2)

    # Plotting
    resamp=2
    
    plt.figure(4, figsize=(12,3))
    plt.clf()
    plt.subplot(1,3,1)
    plt.subplots_adjust(left=0.08)
    plt.imshow(p_o(x_o, y_o), origin='lower', interpolation='nearest')
    plt.colorbar(fraction=0.25)
    #plt.axis('equal')
    plt.title(f'Image (resamp={resamp:d})')
    #plt.pause(0.05)

    plt.subplot(1,3,2)
    plt.subplots_adjust(left=0.08)
    plt.imshow(p_o(x_o, y_o), origin='lower', interpolation='nearest')
    plt.colorbar(fraction=0.25)
    #plt.axis('equal')
    plt.title(f'Model (resamp={resamp:d})')
    #plt.pause(0.05)

    plt.subplot(1,3,3)
    plt.subplots_adjust(left=0.08)
    plt.imshow(z_o - p_o(x_o, y_o), origin='lower', interpolation='nearest')# take out outliers?
    cbar = plt.colorbar(fraction=0.25)
    #plt.axis('equal')
    cbar.set_label('% Residual')
    plt.suptitle(f"Source fit")
    plt.tight_layout()
    #plt.pause(0.05)
    
    #plt.figure(fignumber, figsize=(12, 2.5)); plt.suptitle(title)
    #plt.subplot(1,3,1)
    #plt.imshow(z_o, origin='lower', interpolation='nearest')
    #plt.colorbar()
    #plt.title("Data")
    
    #plt.subplot(1,3,2)
    #plt.imshow(p_o(x_o, y_o), origin='lower', interpolation='nearest')
    #plt.colorbar(); plt.title("Model")
    #plt.subplot(1,3,3)
    #plt.imshow(z_o - p_o(x_o, y_o), origin='lower', interpolation='nearest')
    #plt.colorbar(); plt.title("Residual")
    
    #plt.tight_layout()
    #print("Least Squares Sum - ", title, ": ",'{:.2e}'.format(residual_o))
    #print("FUV - ", title, ": ",'{:.2e}'.format(FUV))

    return p_o

def plot_xy_psf(image, angle):
    #enter angle in radians
    test = rotate(image, degrees(angle), axes=(1, 0), reshape=True, output=None, order=3, mode='constant', cval=0.0, prefilter=True)
    col_sums = []
    for i in range(np.shape(test)[0]):
        col_sums.append(np.sum(test[i,:]))
    max_col_ind = np.argmax(col_sums)
    row_sums = []
    for i in range(np.shape(test)[0]):
        row_sums.append(np.sum(test[:,i]))
    max_row_ind = np.argmax(row_sums)
    plt.figure(1, figsize=(12, 4))
    plt.subplot(1,3,1), plt.title("Image"); plt.imshow(test); plt.plot([0, 64], [max_col_ind, max_col_ind], 'r-'); plt.plot([max_row_ind, max_row_ind], [0, 64] , 'm-')
    plt.subplot(1,3,2), plt.title("x slice"); plt.plot(np.arange(0, np.shape(test)[0]), test[max_col_ind,:], 'r-')
    plt.subplot(1,3,3), plt.title("y slice"); plt.plot(np.arange(0, np.shape(test)[0]), test[:,max_row_ind], 'm-')
    
    return
