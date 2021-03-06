import pylab as plt
import numpy as np
from astropy.io import fits
from astropy import table
from astropy import units
from astropy.stats import sigma_clipped_stats
import scipy
import glob
from imaka.reduce import reduce_fli
from imaka.reduce import calib
from imaka.reduce import util
import os

root_dir = '/Volumes/g/lu/data/imaka/2017_01_13/fli/'

def make_sky():
    sky_raw_dir = root_dir + 'Pleiades/'
    sky_out_dir = root_dir + 'reduce/sky/'

    util.mkdir(sky_out_dir)

    sky_num = [46, 47, 48, 72, 73, 74, 99, 100, 101]
    sky_frames = ['{0:s}sky{1:03d}.fits'.format(sky_raw_dir, ss) for ss in sky_num]
    calib.makedark(sky_frames, sky_out_dir + 'pleiades_sky_46.fits')
    
    sky_num = [598, 599, 600, 626, 627, 1038, 1039, 1040, 1077, 1078, 1079, \
              1109, 1110, 1111, 1112, 1113, 1114, 1115, 1116, 1117]
    sky_frames = ['{0:s}sky{1:03d}.fits'.format(sky_raw_dir, ss) for ss in sky_num]
    calib.makedark(sky_frames, sky_out_dir + 'pleiades_sky_598.fits')

    return

def make_flat():
    flat_raw_dir = root_dir + "twilights/"
    dark_raw_dir = root_dir + "twilights/"
    flat_out_dir = root_dir + "reduce/calib/"

    util.mkdir(flat_out_dir)
    
    flat_num = [0, 1, 2, 3, 4, 5, 6]
    flat_frames = ['{0:s}twi_{1:03d}.fits'.format(flat_raw_dir, ss) for ss in flat_num]
    dark_frames = ['{0:s}dark_{1:03d}.fits'.format(flat_raw_dir, ss) for ss in flat_num]
    calib.makeflat(flat_frames, dark_frames, flat_out_dir + 'flat_1.fits')

    flat_num = [8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, \
               26, 27, 28, 29]
    flat_frames = ['{0:s}twi_{1:03d}.fits'.format(flat_raw_dir, ss) for ss in flat_num]
    dark_frames = ['{0:s}dark_{1:03d}.fits'.format(flat_raw_dir, ss) for ss in flat_num]
    calib.makeflat(flat_frames, dark_frames, flat_out_dir + 'flat_2.fits')

    flat1 = fits.getdata(flat_out_dir + 'flat_1.fits')
    flat2 = fits.getdata(flat_out_dir + 'flat_2.fits')

    flat1_rebin3 = reduce_fli.rebin(flat1, 3)[:,:-1]/9   ## WHAT IS THIS REVERSING/DIVIDNG DOING? 
    mean, median, stdev = sigma_clipped_stats([flat1_rebin3, flat2], sigma=3, iters=2, axis=0)
    hdu = fits.PrimaryHDU(median.data)
    hdu.writeto(flat_out_dir + 'flat.fits')

    return
        
def reduce_pleiades():
    sky_dir = root_dir + 'reduce/sky/'
    data_dir = root_dir + 'Pleiades/'
    flat_dir = root_dir + 'reduce/calib/'
    out_dir = root_dir + 'reduce/pleiades/'

    util.mkdir(out_dir)

    #####
    # Open Loop First Half
    #####
    ###What to do with obj_o594-597??? look at images to see orientation
    ###obj_o640 has a 0.5s integration time...whut
    ###obj_o613.fits in log but not in data.... also 624
    fnum1 = [34, 35, 36, 43, 44, 45, 56, 57, 58, 59, 68, 69, 70, 71, 83, 84, 85, 86]
    fnum2 = [95, 96, 97, 98, 480, 582, 583, 584, 585]
    fnum = fnum1 + fnum2
    img_files = ['{0:s}/obj_o{1:03d}.fits'.format(data_dir, ii) for ii in fnum]
    reduce_fli.clean_images(img_files, out_dir, rebin=1,
                                sky_frame=sky_dir + 'pleiades_sky_46.fits',
                                flat_frame=flat_dir + 'flat.fits')

    #####
    # Open Loop Second Half
    #####
    fnum1 = [609, 610, 611, 612, 622, 623, 624, 636, 637, 638, 639, 1010, 1011]
    fnum2 = [1012, 1013, 1022, 1023, 1024, 1025, 1034, 1035, 1036, 1037, 1049, 1052, 1061]
    fnum3 = [1062, 1063, 1064, 1073, 1074, 1075, 1076, 1088, 1089, 1090, 1091, 1104, 1105]
    fnum4 = [1106, 1107]
    fnum = fnum1 + fnum2 + fnum3 + fnum4
    img_files = ['{0:s}/obj_o{1:03d}.fits'.format(data_dir, ii) for ii in fnum]
    reduce_fli.clean_images(img_files, out_dir, rebin=1,
                                sky_frame=sky_dir + 'pleiades_sky_598.fits',
                                flat_frame=flat_dir + 'flat.fits')

    
    #####
    # TTF Closed Loop First Half
    #####
    ###obj_ttf541.fits in log but not data
    fnum1 = [37, 38, 39, 49, 50, 51, 60, 61, 62, 63, 75, 76, 77, 78, 87, 88, 89, 90]
    fnum2 = [586, 587, 588, 589]
    fnum = fnum1 + fnum2
    img_files = ['{0:s}/obj_ttf{1:03d}.fits'.format(data_dir, ii) for ii in fnum]
    reduce_fli.clean_images(img_files, out_dir, rebin=1,
                                sky_frame=sky_dir + 'pleiades_sky_46.fits',
                                flat_frame=flat_dir + 'flat.fits')
    
    
    #####
    # TTF Closed Loop Second Half
    #####
    ###obj_ttf617 not in file but in log
    fnum1 = [601, 602, 603, 604, 614, 615, 616, 629, 630, 631, 1014, 1015, 1016, 1017]
    fnum2 = [1026, 1027, 1028, 1029, 1041, 1042, 1043, 1044, 1053, 1054, 1055, 1056, 1065]
    fnum3 = [1066, 1067, 1068, 1080, 1081, 1082, 1083, 1092, 1093, 1094, 1095]
    fnum = fnum1 + fnum2 + fnum3
    img_files = ['{0:s}/obj_ttf{1:03d}.fits'.format(data_dir, ii) for ii in fnum]
    reduce_fli.clean_images(img_files, out_dir, rebin=1,
                                sky_frame=sky_dir + 'pleiades_sky_598.fits',
                                flat_frame=flat_dir + 'flat.fits')

    
    #####
    # Closed Loop First Half
    #####
    #obj_c102-103, 407-416, 573-574 have 1s or 0.1s integration times... not reduced    
    fnum1 = [31, 32, 33, 40, 41, 42, 52, 53, 54, 55, 64, 65, 66, 67, 79, 80, 81, 82, 91, 92]
    fnum2 = [93, 94, 417, 448, 479, 577, 578, 579, 580, 581, 590, 591, 592, 593]
    fnum = fnum1 + fnum2
    img_files = ['{0:s}/obj_c{1:03d}.fits'.format(data_dir, ii) for ii in fnum]
    reduce_fli.clean_images(img_files, out_dir, rebin=1,
                                sky_frame=sky_dir + 'pleiades_sky_46.fits',
                                flat_frame=flat_dir + 'flat.fits')

    
    #####
    # Closed Loop Second Half
    #####
    ###obj_c641 has 0.5s integration time    
    ###obj_c621.fits doesn't exist in file, only in log
    fnum1 = [605, 606, 607, 608, 618, 619, 620, 632, 633, 634, 635, 1006, 1007, 1008]
    fnum2 = [1009, 1018, 1019, 1020, 1021, 1030, 1031, 1032, 1033, 1045, 1046, 1047, 1048]
    fnum3 = [1057, 1058, 1059, 1060, 1069, 1070, 1071, 1072, 1084, 1085, 1086, 1087, 1096]
    fnum4 = [1097, 1098, 1099]
    fnum = fnum1 + fnum2 + fnum3 + fnum4
    img_files = ['{0:s}/obj_c{1:03d}.fits'.format(data_dir, ii) for ii in fnum]
    reduce_fli.clean_images(img_files, out_dir, rebin=1,
                                sky_frame=sky_dir + 'pleiades_sky_598.fits',
                                flat_frame=flat_dir + 'flat.fits')
    
    
    return
    
    
def find_stars_pleiades_open():
    reduce_dir = root_dir + 'reduce/pleiades/'
    
    ##########
    # first half of night
    ##########
    fnum1 = [34, 35, 36, 43, 44, 45, 56, 57, 58, 59, 68, 69, 70, 71, 83, 84, 85, 86]
    fnum2 = [95, 96, 97, 98, 480, 582, 583, 584, 585]
    fnum = fnum1 + fnum2 
    img_files = ['{0:s}/obj_o{1:03d}_clean.fits'.format(reduce_dir, ii) for ii in fnum]
    
    reduce_fli.find_stars(img_files, fwhm=5, threshold=6)

    ##########
    # second half of night
    ##########
    fnum1 = [609, 610, 611, 612, 622, 623, 624, 636, 637, 638, 639, 1010, 1011]
    fnum2 = [1012, 1013, 1022, 1023, 1024, 1025, 1034, 1035, 1036, 1037, 1049, 1052, 1061]
    fnum3 = [1062, 1063, 1064, 1073, 1074, 1075, 1076, 1088, 1089, 1090, 1091, 1104, 1105]
    fnum4 = [1106, 1107]
    fnum = fnum1 + fnum2 + fnum3 + fnum4 
    img_files = ['{0:s}/obj_o{1:03d}_clean.fits'.format(reduce_dir, ii) for ii in fnum]
    
    reduce_fli.find_stars(img_files, fwhm=5, threshold=6)

    return


def find_stars_pleiades_ttf():
    reduce_dir = root_dir + 'reduce/pleiades/'
    
    ##########
    # first half of the night
    ##########
    fnum1 = [37, 38, 39, 49, 50, 51, 60, 61, 62, 63, 75, 76, 77, 78, 87, 88, 89, 90]
    fnum2 = [586, 587, 588, 589]
    fnum = fnum1 + fnum2
    img_files = ['{0:s}/bj_ttf{1:03d}_clean.fits'.format(reduce_dir, ii) for ii in fnum]
    
    reduce_fli.find_stars(img_files, fwhm=5, threshold=6)

    ##########
    # second half of the night
    ##########
    fnum1 = [601, 602, 603, 604, 614, 615, 616, 629, 630, 631, 1014, 1015, 1016, 1017]
    fnum2 = [1026, 1027, 1028, 1029, 1041, 1042, 1043, 1044, 1053, 1054, 1055, 1056, 1065]
    fnum3 = [1066, 1067, 1068, 1080, 1081, 1082, 1083, 1092, 1093, 1094, 1095]
    fnum = fnum1 + fnum2 + fnum3
    img_files = ['{0:s}/obj_ttf{1:03d}_clean.fits'.format(reduce_dir, ii) for ii in fnum]
    
    reduce_fli.find_stars(img_files, fwhm=5, threshold=6)    
    
    return


def find_stars_pleiades_closed():
    reduce_dir = root_dir + 'reduce/pleiades/'
    
    ##########
    # first half of night
    ##########
    fnum1 = [31, 32, 33, 40, 41, 42, 52, 53, 54, 55, 64, 65, 66, 67, 79, 80, 81, 82, 91, 92]
    fnum2 = [93, 94, 417, 448, 479, 577, 578, 579, 580, 581, 590, 591, 592, 593]
    fnum = fnum1 + fnum2
    img_files = ['{0:s}/obj_c{1:03d}_clean.fits'.format(reduce_dir, ii) for ii in fnum]
    
    reduce_fli.find_stars(img_files, fwhm=3, threshold=6)

    ##########
    # second half of night
    ##########
    fnum1 = [605, 606, 607, 608, 618, 619, 620, 632, 633, 634, 635, 1006, 1007, 1008]
    fnum2 = [1009, 1018, 1019, 1020, 1021, 1030, 1031, 1032, 1033, 1045, 1046, 1047, 1048]
    fnum3 = [1057, 1058, 1059, 1060, 1069, 1070, 1071, 1072, 1084, 1085, 1086, 1087, 1096]
    fnum4 = [1097, 1098, 1099]
    fnum = fnum1 + fnum2 + fnum3 + fnum4
    img_files = ['{0:s}/obj_c{1:03d}_clean.fits'.format(reduce_dir, ii) for ii in fnum]
    
    reduce_fli.find_stars(img_files, fwhm=3, threshold=6)
    
    return


def calc_star_stats():
    reduce_dir = root_dir + 'reduce/pleiades/'
    stats_dir = root_dir + 'reduce/stats/'

    # open first half of night
    fnum1 = [34, 35, 36, 43, 44, 45, 56, 57, 58, 59, 68, 69, 70, 71, 83, 84, 85, 86]
    fnum2 = [95, 96, 97, 98, 480, 582, 583, 584, 585] 
    fnum = fnum1 + fnum2 
    img_files = ['{0:s}/obj_o{1:03d}_clean.fits'.format(reduce_dir, ii) for ii in fnum]
    reduce_fli.calc_star_stats(img_files, output_stats=stats_dir + 'stats_open1.fits')

    # open second half of night
    fnum1 = [609, 610, 611, 612, 622, 623, 624, 636, 637, 638, 639, 1010, 1011]
    fnum2 = [1012, 1013, 1022, 1023, 1024, 1025, 1034, 1035, 1036, 1037, 1049, 1052, 1061]
    fnum3 = [1062, 1063, 1064, 1073, 1074, 1075, 1076, 1088, 1089, 1090, 1091, 1104, 1105]
    fnum4 = [1106, 1107]
    fnum = fnum1 + fnum2 + fnum3 + fnum4   
    img_files = ['{0:s}/obj_o{1:03d}_clean.fits'.format(reduce_dir, ii) for ii in fnum]
    reduce_fli.calc_star_stats(img_files, output_stats=stats_dir + 'stats_open2.fits')

    # # ttf first half of the night
    # fnum1 = [37, 38, 39, 49, 50, 51, 60, 61, 62, 63, 75, 76, 77, 78, 87, 88, 89, 90]
    # fnum2 = [586, 587, 588, 589]
    # fnum = fnum1 + fnum2
    # img_files = ['{0:s}/obj_ttf{1:03d}_clean.fits'.format(reduce_dir, ii) for ii in fnum]
    # reduce_fli.calc_star_stats(img_files, output_stats=stats_dir + 'stats_ttf1.fits')

    # # ttf second half of the night
    # fnum1 = [601, 602, 603, 604, 614, 615, 616, 629, 630, 631, 1014, 1015, 1016, 1017]
    # fnum2 = [1026, 1027, 1028, 1029, 1041, 1042, 1043, 1044, 1053, 1054, 1055, 1056, 1065]
    # fnum3 = [1066, 1067, 1068, 1080, 1081, 1082, 1083, 1092, 1093, 1094, 1095]
    # fnum = fnum1 + fnum2 + fnum3
    # img_files = ['{0:s}/obj_ttf{1:03d}_clean.fits'.format(reduce_dir, ii) for ii in fnum]
    # reduce_fli.calc_star_stats(img_files, output_stats=stats_dir + 'stats_ttf2.fits')
    
    # # closed first half of the night
    # fnum1 = [31, 32, 33, 40, 41, 42, 52, 53, 54, 55, 64, 65, 66, 67, 79, 80, 81, 82, 91, 92]
    # fnum2 = [93, 94, 417, 448, 479, 577, 578, 579, 580, 581, 590, 591, 592, 593]
    # fnum = fnum1 + fnum2
    # img_files = ['{0:s}/obj_c{1:03d}_clean.fits'.format(reduce_dir, ii) for ii in fnum]
    # reduce_fli.calc_star_stats(img_files, output_stats=stats_dir + 'stats_closed1.fits')
        
    # # closed second half of night
    # fnum1 = [605, 606, 607, 608, 618, 619, 620, 632, 633, 634, 635, 1006, 1007, 1008]
    # fnum2 = [1009, 1018, 1019, 1020, 1021, 1030, 1031, 1032, 1033, 1045, 1046, 1047, 1048]
    # fnum3 = [1057, 1058, 1059, 1060, 1069, 1070, 1071, 1072, 1084, 1085, 1086, 1087, 1096]
    # fnum4 = [1097, 1098, 1099]
    # fnum = fnum1 + fnum2 + fnum3 + fnum4
    # img_files = ['{0:s}/obj_c{1:03d}_clean.fits'.format(reduce_dir, ii) for ii in fnum]
    # reduce_fli.calc_star_stats(img_files, output_stats=stats_dir + 'stats_closed2.fits')

    # plot_stats()
    
    return
    

def compare_fwhm_list():
    data_dir = '/Users/jlu/data/imaka/2017_01_13/fli/Pleiades/'
    os.chdir(data_dir)
    
    # o_list = np.arange(163, 173) # Open loop star lists
    # c_list = np.arange(153, 163) # Closed loop star lists
    # o_list = np.arange(163, 173) # Open loop star lists
    # c_list = np.arange(153, 163) # Closed loop star lists
    o_list = [10, 11, 14, 15] # open
    c_list = [8, 9, 12, 13]   # closed
    
    plt.ion()

    for ii in range(len(o_list)):
        open_list = '{0:s}/obj{1:03d}_clean_stars.txt'.format(reduce_dir, o_list[ii])
        closed_list = '{0:s}/obj{1:03d}_clean_stars.txt'.format(reduce_dir, c_list[ii])

        compare_fwhm(open_list, closed_list, scale=3*0.04, flux_min=1)
        pdb.set_trace()

def compare_fwhm(open_list, closed_list, scale=0.04, flux_min=2.0):
    topen = table.Table.read(open_list, format='ascii')
    tclose = table.Table.read(closed_list, format='ascii')

    print("N_open = {0:3d}  N_close = {1:3d} in original lists".format(len(topen), len(tclose)))

    # Trim out any stars with FWHM > 20
    idx_o = np.where(topen['x_fwhm'] < 20)[0]
    idx_c = np.where(tclose['x_fwhm'] < 20)[0]

    topen = topen[idx_o]
    tclose = tclose[idx_c]
    
    print("N_open = {:3d}  N_close = {:3d} after trimming fwhm<20".format(len(topen), len(tclose)))

    # Trim out any stars with flux < flux_min
    idx_o = np.where(topen['flux'] > flux_min)[0]
    idx_c = np.where(tclose['flux'] > flux_min)[0]

    topen = topen[idx_o]
    tclose = tclose[idx_c]
    
    print("N_open = {:3d}  N_close = {:3d} after trimming low flux sources".format(len(topen), len(tclose)))
    
    m_c = np.ones(len(tclose))
    m_o = np.ones(len(topen))

    idx_c, idx_o, dr, dm = match.match(tclose['xcentroid'], tclose['ycentroid'], m_c,
                                       topen['xcentroid'], topen['ycentroid'], m_o,
                                       dr_tol=10)
    # Matched catalogs
    to_match = topen[idx_o]
    tc_match = tclose[idx_c]
    
    print("N_open = {:3d}  N_close = {:3d} after matching".format(len(to_match), len(tc_match)))

    # Plot
    plt.figure(1)
    plt.clf()
    plt.plot(to_match['x_fwhm']*scale, tc_match['x_fwhm']*scale, 'r.', label='X')
    plt.plot(to_match['y_fwhm']*scale, tc_match['y_fwhm']*scale, 'b.', label='Y')
    plt.plot([0, 10], [0, 10], 'k--')
    plt.xlabel('FWHM in Open Loop (")')
    plt.ylabel('FWHM in Closed Loop (")')
    plt.axis('equal')

    max_fwhm = 1.5*np.mean([to_match['x_fwhm'], to_match['y_fwhm'], tc_match['x_fwhm'], tc_match['y_fwhm']])
    plt.ylim(0, max_fwhm*scale)
    plt.xlim(0, max_fwhm*scale)
    plt.legend(numpoints=1, loc='upper left')
    plt.pause(0.05)

    plt.figure(2)
    plt.clf()
    x_ratio = tc_match['x_fwhm'] / to_match['x_fwhm']
    y_ratio = tc_match['y_fwhm'] / to_match['y_fwhm']
    
    plt.plot(to_match['x_fwhm']*scale, x_ratio, 'r.', label='X')
    plt.plot(to_match['y_fwhm']*scale, y_ratio, 'b.', label='Y')
    plt.axhline(1, linestyle='--', color='black')
    plt.xlabel('FWHM in Open Loop (")')
    plt.ylabel('Closed / Open FWHM')

    max_fwhm = 1.5*np.mean([to_match['x_fwhm'], to_match['y_fwhm'], tc_match['x_fwhm'], tc_match['y_fwhm']])
    plt.xlim(0, max_fwhm*scale)
    plt.ylim(0, 1.5)
    plt.legend(numpoints=1, loc='upper left')
    plt.pause(0.05)

    x_fwhm_o_med = np.median(to_match['x_fwhm'])    
    y_fwhm_o_med = np.median(to_match['y_fwhm'])    
    x_fwhm_c_med = np.median(tc_match['x_fwhm'])    
    y_fwhm_c_med = np.median(tc_match['y_fwhm'])
    
    print('Open Loop Stats:')
    print('\t Median x_fwhm = {0:.1f} +/- {1:.1f}'.format(x_fwhm_o_med * scale,
                                                          to_match['x_fwhm'].std() * scale))
    print('\t Median y_fwhm = {0:.1f} +/- {1:.1f}'.format(y_fwhm_o_med * scale,
                                                          to_match['y_fwhm'].std() * scale))

    print('Closed Loop Stats:')
    print('\t Median x_fwhm = {0:.1f} +/- {1:.1f}'.format(x_fwhm_c_med * scale,
                                                          tc_match['x_fwhm'].std() * scale))
    print('\t Median y_fwhm = {0:.1f} +/- {1:.1f}'.format(y_fwhm_c_med * scale,
                                                          tc_match['y_fwhm'].std() * scale))

    print('Fractional Improvement')
    print('\t x_fwhm closed / open = {0:.2f}'.format(x_fwhm_c_med / x_fwhm_o_med))
    print('\t y_fwhm closed / open = {0:.2f}'.format(y_fwhm_c_med / y_fwhm_o_med))
    


def stack_pleiades_west_ttf():
    """
    Stack all the images from the night for the Open loop, TTF closed loop, and GLAO closed loop
    images. Calculate statstics on all three.
    """
    red_dir = root_dir + 'reduce/pleiades/'
    
    # open second half of night (on west field
    open_img_nums = [609, 610, 611, 612, 622, 623, 624, 636, 637, 638, 639, 1010, 1011]
    open_img_nums += [1012, 1013, 1022, 1023, 1024, 1025, 1034, 1035, 1036, 1037, 1049, 1052, 1061]
    open_img_nums += [1062, 1063, 1064, 1073, 1074, 1075, 1076, 1088, 1089, 1090, 1091, 1104, 1105]
    open_img_nums += [1106, 1107]
    open_images = ['{0:s}/obj_o{1:03d}_clean.fits'.format(red_dir, ii) for ii in open_img_nums]
    open_starlists = ['{0:s}/obj_o{1:03d}_clean.txt'.format(red_dir, ii) for ii in open_img_nums]
    open_output_root = 'west_stack_open'
    reduce_fli.shift_and_add(open_images, open_starlists, open_output_root, method='mean')

    # TTF closed second half of night (on west field)
    ttf_img_nums = [601, 602, 603, 604, 614, 615, 616, 629, 630, 631, 1014, 1015, 1016, 1017]
    ttf_img_nums += [1026, 1027, 1028, 1029, 1041, 1042, 1043, 1044, 1053, 1054, 1055, 1056, 1065]
    ttf_img_nums += [1066, 1067, 1068, 1080, 1081, 1082, 1083, 1092, 1093, 1094, 1095]
    ttf_images = ['{0:s}/obj_ttf{1:03d}_clean.fits'.format(red_dir, ii) for ii in ttf_img_nums]
    ttf_starlists = ['{0:s}/obj_ttf{1:03d}_clean_stars.txt'.format(red_dir, ii) for ii in ttf_img_nums]
    ttf_output_root = 'west_stack_ttf'
    reduce_fli.shift_and_add(ttf_images, ttf_starlists, ttf_output_root, method='mean')

    # GLAO closed second half of night (on west field
    closed_img_nums = [605, 606, 607, 608, 618, 619, 620, 632, 633, 634, 635, 1006, 1007, 1008]
    closed_img_nums += [1009, 1018, 1019, 1020, 1021, 1030, 1032, 1033, 1045, 1046, 1047, 1048]
    closed_img_nums += [1057, 1058, 1059, 1060, 1069, 1070, 1071, 1072, 1084, 1085, 1086, 1087, 1096]
    closed_img_nums += [1097, 1098, 1099]
    closed_images = ['{0:s}/obj_c{1:03d}_clean.fits'.format(red_dir, ii) for ii in closed_img_nums]
    closed_starlists = ['{0:s}/obj_c{1:03d}_clean_stars.txt'.format(red_dir, ii) for ii in closed_img_nums]
    closed_output_root = 'west_stack_closed'
    reduce_fli.shift_and_add(closed_images, closed_starlists, closed_output_root, method='mean')

    return


def analyze_stacks():
    img_files = ['west_stack_open.fits', 'west_stack_ttf.fits', 'west_stack_closed.fits']
    
    find_stars(img_files, fwhm=5, threshold=4, N_passes=2, plot_psf_compare=True)
    
    # Calc stats on all the stacked images
    reduce_fli.calc_star_stats(img_files, output_stats='stats_stacks.fits')

    return
