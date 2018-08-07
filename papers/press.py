import pdb
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors
from astropy.io import fits
from astropy.table import Table
from imaka.reduce import reduce_fli
from flystar import match
from datetime import datetime
from matplotlib import dates as mp_dates
from matplotlib import ticker
import subprocess

def run_find_stars():
    work_dir = '/Users/jlu/work/imaka/pleiades/press_release/'
    
    img_files = ['west_stack_open.fits', 'west_stack_ttf.fits', 'west_stack_closed.fits']
    for ii in range(len(img_files)):
        img_files[ii] = work_dir + img_files[ii]
    
    reduce_fli.find_stars(img_files, fwhm=5, threshold=4, N_passes=2)

    # Calc stats on all the stacked images
    reduce_fli.calc_star_stats(img_files, output_stats=work_dir + 'stats_stacks.fits')
    
    return

def first_light_images():
    """
    Plot up the first light images for the press release.
    """
    work_dir = '/Users/jlu/work/imaka/pleiades/press_release/'

    img_op = fits.getdata(work_dir + 'west_stack_open.fits')
    img_tf = fits.getdata(work_dir + 'west_stack_ttf.fits')
    img_cl = fits.getdata(work_dir + 'west_stack_closed.fits')

    # Use the centroids in each of the images to shift the images.
    t_op = Table.read(work_dir + 'west_stack_open_stars.txt', format='ascii')
    t_tf = Table.read(work_dir + 'west_stack_ttf_stars.txt', format='ascii')
    t_cl = Table.read(work_dir + 'west_stack_closed_stars.txt', format='ascii')

    # Trim out any stars with FWHM > fwhm_max and flux < flux_min
    fwhm_max = 20
    flux_min = 5
    idx_op = np.where((t_op['x_fwhm'] < fwhm_max) & (t_op['flux'] > flux_min))[0]
    idx_tf = np.where((t_tf['x_fwhm'] < fwhm_max) & (t_tf['flux'] > flux_min))[0]
    idx_cl = np.where((t_cl['x_fwhm'] < fwhm_max) & (t_cl['flux'] > flux_min))[0]
    t_op = t_op[idx_op]
    t_tf = t_tf[idx_tf]
    t_cl = t_cl[idx_cl]

    # Add zeropoints (completely arbitrary)
    ZP = 10.0
    t_op['mag'] += ZP
    t_tf['mag'] += ZP
    t_cl['mag'] += ZP
    
    # Matched catalogs
    mdx_c, mdx_o, dr, dm = match.match(t_cl['xcentroid'], t_cl['ycentroid'], t_cl['mag'],
                                       t_op['xcentroid'], t_op['ycentroid'], t_op['mag'],
                                       dr_tol=10)
    t_op_match = t_op[mdx_o]
    t_cl_match = t_cl[mdx_c]

    x_sh = t_cl_match['xcentroid'] - t_op_match['xcentroid']
    y_sh = t_cl_match['ycentroid'] - t_op_match['ycentroid']

    x_sh_med = np.median(x_sh)
    y_sh_med = np.median(y_sh)

    fmt = 'Closed - Open Shifts: X = {0:5.2f} +/- {1:5.2f}, Y = {2:5.2f} +/- {3:5.2f}'
    print(fmt.format(x_sh_med, x_sh.std(), y_sh_med, y_sh.std()))

    # Define image extent
    img_height = float(img_cl.shape[0])
    img_width = float(img_cl.shape[1])

    img_extent_cl = [0, img_cl.shape[1], 0, img_cl.shape[0]]
    img_extent_op = [x_sh_med, img_cl.shape[1] + x_sh_med, y_sh_med, img_cl.shape[0] + y_sh_med]

    img_limits_cl = [0, img_cl.shape[1], 0, img_cl.shape[0]]
    img_limits_op = img_limits_cl
    
    plot_image_with_zoom_insets(img_cl, img_extent_cl, img_limits_cl,
                                    "Closed Loop", work_dir + 'full_closed.png')
    plot_image_with_zoom_insets(img_op, img_extent_op, img_limits_op,
                                    "Open Loop", work_dir + 'full_open.png')

    cmd1 = 'convert full_open.png full_closed.png -delay 100 -morph 10 full.gif'
    cmd2 = 'convert full.gif -coalesce -duplicate 1,-2-1 -quiet -layers OptimizePlus -loop 0 full_patrol.gif'
    subprocess.call(cmd1.split())
    subprocess.call(cmd2.split())

    return

def plot_image_with_zoom_insets(img, extent, limits, text_label, output_name):
    img_height = float(img.shape[0])
    img_width = float(img.shape[1])

    # norm = plt.Normalize(30, 2000)
    # norm = plt.Normalize(img_cl.min(), img.max())
    norm = colors.LogNorm(vmin=30, vmax=2000)

    star1 = [339, 1207]  # X, Y, Top Left
    star2 = [1570, 1676] # Top Right
    star3 = [2382, 958]  # Bottom Right
    star4 = [182, 475]   # Bottom Left
    stars = [star1, star2, star3, star4]

    axes_x = [-0.025, 0.825, 0.825, -0.025]
    axes_y = [0.8, 0.8, 0.0, 0.0]

    sub_norm_min = [1, 1, -1, -1]
    sub_norm_max = [900, 730, 4200, 210]
    
    # Make the large scale closed-loop image.
    fig = plt.figure(1)
    fig.set_size_inches(3*img_width/img_height, 3, forward=False)
    plt.clf()
    plt.subplots_adjust(top=1, bottom=0, right=1, left=0, hspace=0, wspace=0)
    ax0 = plt.axes([0, 0, 1, 1])
    ax0.set_axis_off()
    fig.add_axes(ax0)
    ax0.imshow(img+50, extent=extent, cmap=cm.Greys_r, norm=norm)
    ax0.set_xlim(limits[0], limits[1])
    ax0.set_ylim(limits[2], limits[3])

    # text_label0 = "'imaka Ground-Layer Adaptive Optics"
    text_label0 = "'imaka"
    ax0.text(0.5*img_width, 0.15*img_height, text_label0,
                 color='w', horizontalalignment='center', fontsize=14)
    ax0.text(0.5*img_width, 0.05*img_height, text_label,
                 color='w', horizontalalignment='center', fontsize=14)

    for ss in range(len(stars)):
        sub_xlo = stars[ss][0] - 15
        sub_xhi = sub_xlo + 31
        sub_ylo = stars[ss][1] - 15
        sub_yhi = sub_ylo + 31

        ext_xlo = extent[0] + sub_xlo
        ext_xhi = ext_xlo + 31
        ext_ylo = extent[2] + sub_ylo
        ext_yhi = ext_ylo + 31

        lim_xlo = limits[0] + sub_xlo + 2
        lim_xhi = lim_xlo + 27
        lim_ylo = limits[2] + sub_ylo + 2
        lim_yhi = lim_ylo + 27
        
        # Make an inset
        sub_img = img[sub_ylo:sub_yhi, sub_xlo:sub_xhi]
        sub_extent = [ext_xlo, ext_xhi, ext_ylo, ext_yhi]
        sub_norm = plt.Normalize(sub_norm_min[ss], sub_norm_max[ss])

        # print(sub_xlo, sub_xhi, sub_ylo, sub_yhi)
        # print(ext_xlo, ext_xhi, ext_ylo, ext_yhi)
        # print(lim_xlo, lim_xhi, lim_ylo, lim_yhi)

        ax1_center = [(axes_x[ss] + 0.1) * img_width, (axes_y[ss] + 0.1) * img_height]
        ax0.arrow(ax1_center[0], ax1_center[1],
                  (stars[ss][0] - ax1_center[0])*0.9, (stars[ss][1] - ax1_center[1])*0.9,
                  width=0.001, head_width=10, head_length=10, fc='w', ec='w')
        ax1 = plt.axes([axes_x[ss], axes_y[ss], 0.2, 0.2], frameon=False, aspect='equal')
        ax1.set_axis_off()
        ax1.imshow(sub_img, extent=sub_extent, norm=sub_norm)
        ax1.set_xlim(lim_xlo, lim_xhi)
        ax1.set_ylim(lim_ylo, lim_yhi)

    plt.savefig(output_name, pad_inches=0.05, dpi=300)
    
    return


def plot_stats():
    stats_dir = '/Users/jlu/work/imaka/pleiades/20170112/fli/reduce/stats/'
    work_dir  = '/Users/jlu/work/imaka/pleiades/press_release/'

    colors = ['b', 'g', 'r', 'c', 'm', 'k', 'yellow', 'purple', 'orange']
    scale = 0.12

    stats_o = Table.read(stats_dir + 'stats_open.fits')
    stats_c = Table.read(stats_dir + 'stats_closed.fits')

    hst_dt_o = []
    hst_dt_c = []

    for ii in range(len(stats_o)):
        time_tmp = stats_o['DATE_HST'][ii] + ' ' + stats_o['TIME_HST'][ii]
        hst_dt_o.append( datetime.strptime(time_tmp, '%Y-%m-%d %H:%M:%S') )
    for ii in range(len(stats_c)):
        time_tmp = stats_c['DATE_HST'][ii] + ' ' + stats_c['TIME_HST'][ii]
        hst_dt_c.append( datetime.strptime(time_tmp, '%Y-%m-%d %H:%M:%S') )

    # Read in Olivier's Stats
    ao_telem = fits.getdata(work_dir + 'profile_data_20170112.fits')
    ao_time_tmp = ao_telem[:, 0]
    ao_seeing_all = ao_telem[:, 1]
    ao_seeing_gl = ao_telem[:, 2]

    ao_time = []
    for ii in range(len(ao_telem)):
        hour = int( math.floor(ao_time_tmp[ii]) )
        minute = int( round((ao_time_tmp[ii] % 1) * 60.0) )

        if minute >= 60:
            minute -= 60
            hour += 1
            
        if ao_time_tmp[ii] < 24.0:
            time_str = '{0:s} {1:02d} {2:02d}'.format(stats_c['DATE_HST'][0],
                                                          hour, minute)
        else:
            hour -= 24
            time_str = '{0:s} {1:02d} {2:02d}'.format(stats_c['DATE_HST'][-1],
                                                          hour, minute)
            
        ao_time.append( datetime.strptime(time_str, '%Y-%m-%d %H %M') )

    time_fmt = mp_dates.DateFormatter('%H:%M')
    time_loc = mp_dates.MinuteLocator(byminute=[0])

    #####
    # EMPIRICAL FWHM PLOT
    #####
    plt.figure(1)
    plt.clf()
    plt.subplots_adjust(bottom=0.2)
    plt.plot(hst_dt_o, stats_o['emp_fwhm']*scale, color='r', marker='o', label='Open Loop',
                 alpha=1.0, linestyle='none')
    plt.plot(hst_dt_c, stats_c['emp_fwhm']*scale, color='g', marker='o', label='Closed Loop',
                 alpha=1.0, linestyle='none')

    # plt.plot(ao_time, ao_seeing_all, label='Seeing')
    # plt.plot(ao_time, ao_seeing_gl, label='GL Seeing')
    plt.gca().xaxis.set_major_formatter(time_fmt)
    plt.gca().xaxis.set_major_locator(time_loc)
    plt.legend(loc='upper right', numpoints=1)
    plt.xticks(rotation=35)
    plt.ylim(0, 1.4)
    xlo = datetime(hst_dt_c[0].year, hst_dt_c[0].month, hst_dt_c[0].day,
                            19, 0, 0)
    xhi = datetime(hst_dt_c[-1].year, hst_dt_c[-1].month, hst_dt_c[-1].day,
                            0, 30, 0)
    plt.xlim(xlo, xhi)
    plt.xlabel("Time (hours)", labelpad=10)
    plt.ylabel('Image Resolution (arcsec)', labelpad=10)
    plt.savefig(work_dir + 'fwhm_vs_time_2012_01_12.png')
    
    return
    
