from datetime import datetime
import numpy as np
from imaka.analysis import plot_stats
from astropy.table import Table, Column, vstack, hstack
from astropy.io import fits
import os
from astropy import units as u
from astropy.coordinates import SkyCoord, AltAz
from astropy.coordinates import EarthLocation
from astropy.time import Time
import pdb



def append_alt_data(stats_file, alt_file, date):  

    stats = Table.read(stats_file)
    new_table = fits.getdata(alt_file)

    time = new_table[:,0]

    alt_date = np.zeros(len(stats), dtype='S15')
    alt_time = np.zeros(len(stats), dtype='S15')
    int_seeing = np.zeros(len(stats), dtype=float)
    alt_0 = np.zeros(len(stats), dtype=float)
    alt_60 = np.zeros(len(stats), dtype=float)
    alt_120 = np.zeros(len(stats), dtype=float)
    alt_180 = np.zeros(len(stats), dtype=float)
    alt_240 = np.zeros(len(stats), dtype=float)
    alt_300 = np.zeros(len(stats), dtype=float)
    alt_360 = np.zeros(len(stats), dtype=float)
    alt_420 = np.zeros(len(stats), dtype=float)
    alt_480 = np.zeros(len(stats), dtype=float)
    alt_560 = np.zeros(len(stats), dtype=float)
    alt_4000 = np.zeros(len(stats), dtype=float)

    for ii in range(len(stats)):
        hst_str = stats[ii]['DATE_UTC']+","+stats[ii]['TIME_UTC']
        time_obj = datetime.datetime.strptime(hst_str, "%Y-%m-%d,%H:%M:%S")

        time_diffs = []
        for jj in range(len(new_table)):

            hour = int(time[jj]//1)
            minute_dec = (time[jj]-hour)*60
            minute = int(minute_dec//1)
            second = int(round((minute_dec-minute)*60, 0))
            if hour >= 24:
                hour = hour-24
            if minute == 60:
                hour = hour + 1
                minute = 0
            if second == 60:
                minute = minute + 1
                second = 0
            HST_date_str = date[0:4]+"-"+date[4:6]+"-"+date[6:8]
            HST_time_str = str(hour)+":"+str(minute)+":"+str(second)
            comp_time = datetime.datetime.strptime(HST_date_str+","+HST_time_str, '%Y-%m-%d,%H:%M:%S')

            time_diff = ((time_obj-comp_time)).total_seconds()
            time_diffs.append(abs(time_diff))

        min_index = np.argmin(time_diffs)

        alt_date[ii] = HST_date_str
        alt_time[ii] = HST_time_str
        int_seeing[ii] = new_table[:,1][min_index]
        alt_0[ii] = new_table[:,2][min_index]
        alt_60[ii] = new_table[:,3][min_index]
        alt_120[ii] = new_table[:,4][min_index]
        alt_180[ii] = new_table[:,5][min_index]
        alt_240[ii] = new_table[:,6][min_index]
        alt_300[ii] = new_table[:,7][min_index]
        alt_360[ii] = new_table[:,8][min_index]
        alt_420[ii] = new_table[:,9][min_index]
        alt_480[ii] = new_table[:,10][min_index]
        alt_560[ii] = new_table[:,11][min_index]
        alt_4000[ii] = new_table[:,12][min_index]
        
    ##Make columns to append to table
    col_date = Column(name='alt_HST_date', data=alt_date)
    col_time = Column(name='alt_HST_time', data=alt_time)
    col_int_seeing =Column(name='int_seeing', data=int_seeing)
    col_alt_0 =Column(name='alt_0', data=alt_0)
    col_alt_60 =Column(name='alt_60', data=alt_60)
    col_alt_120 =Column(name='alt_120', data=alt_120)
    col_alt_180 =Column(name='alt_180', data=alt_180)
    col_alt_240 =Column(name='alt_240', data=alt_240)
    col_alt_300 =Column(name='alt_300', data=alt_300)
    col_alt_360 =Column(name='alt_360', data=alt_360)
    col_alt_420 =Column(name='alt_420', data=alt_420)
    col_alt_480 =Column(name='alt_480', data=alt_480)
    col_alt_560 =Column(name='alt_560', data=alt_560)
    col_alt_4000 =Column(name='alt_4000', data=alt_4000)

    ## make new stats file with appended data

    stats.add_columns([col_date, col_time, col_int_seeing, col_alt_0, col_alt_60, col_alt_120, \
                      col_alt_180, col_alt_240, col_alt_300, col_alt_360, \
                      col_alt_420, col_alt_480, col_alt_560, col_alt_4000])

    stats_file_root, stats_file_ext = os.path.splitext(stats_file)
    stats.write(stats_file_root + '_alt' + stats_file_ext, overwrite=True)
    
    return



def append_altaz(stats_file, img_file_dir):
    
    mauna_kea = EarthLocation(lat=19.8206*u.deg, lon=-155.4681*u.deg, height=4207*u.m)
    
    # Read in stats table and image files in that table
    stats_table = Table.read(stats_file)
    img_files = stats_table['Image']
    
    # Make some empty columns for new data
    az = np.zeros(len(img_files), dtype=float)
    alt = np.zeros(len(img_files), dtype=float)
    HA = []
    # Calculate new data from each header
    for ii in range(len(img_files)):
        file = img_file_dir + img_files[ii].split("/")[-1]
        
        # Read in header data
        data, hdr = fits.getdata(file, header=True)
        DEC = hdr['DEC']
        RA = hdr['RA']
        UTIM = hdr['UTIME']
        date = hdr['DATEOBS']
        
        # Define observation parameters
        time = Time(date[6:]+"-"+date[0:2]+"-"+date[3:5]+ " " +UTIM, scale='utc', location=mauna_kea)
        mauna_kea = EarthLocation(lat=19.8206*u.deg, lon=-155.4681*u.deg, height=4207*u.m)
        c = SkyCoord(RA*u.hour, DEC*u.deg)
        
        # Calculate Alt/Az coordinates (both output in degrees)
        alt_az = c.transform_to(AltAz(obstime=time, location=mauna_kea))
        string = alt_az.to_string('decimal').split(" ")
        az[ii] = string[0]
        alt[ii] = string[1]
        
        # Calculate local sidereal time and hour angle:
        lst = time.sidereal_time('apparent')
        hour_angle = lst - c.ra
        HA.append(hour_angle.to_string("hour"))
    
    HA = np.array(HA)
    
    # Make columns to append to table
    col_az = Column(name='AZ', data=az)
    col_alt = Column(name='ALT', data=alt)
    col_HA = Column(name='HA', data=HA)
    
    stats_table.add_columns([col_az, col_alt, col_HA])
    stats_file_root, stats_file_ext = os.path.splitext(stats_file)
    stats_table.write(stats_file_root + '_altaz' + stats_file_ext, overwrite=True)
    
    return





def match_cols(open_file, closed_file, comp_col):
    """
    time matches colums of data from different stats files, e.g. matches FWHM column from 
    a given open stats file to those in a closed file from the same night.
    
    inputs:
        open_file: the stats file that you are matching to
        closed_file: the comparison file to match
        comp_col: the column you want (e.g., NEA1, emp_FWHM, etc)
        
    outputs:
        UTC_TIME: a time average of each two matched data points
        UTC_DATE: corresponding date to above time
        data_1: the column of whatever comp_col described from the base file
        data_2: the column of whatever comp_col described from the comparison file 
    """
            
    stats1 = Table.read(open_file)
    stats2 = Table.read(closed_file)
    
    a = 0    
    if len(stats1) > len(stats2):
        stats1 = Table.read(closed_file)
        stats2 = Table.read(open_file)
        a+=1

    UTC_TIME = np.zeros(len(stats1), dtype='S15')
    UTC_DATE = np.zeros(len(stats1), dtype='S15')
    data_1 = np.zeros(len(stats1),dtype=float)
    data_2 = np.zeros(len(stats1), dtype=float)
    data_1_err = np.zeros(len(stats1), dtype=float)
    data_2_err = np.zeros(len(stats1), dtype=float)

    for ii in range(len(stats1)):
        data_1[ii] = stats1[comp_col][ii]
        data_1_err[ii] = stats1['emp_fwhm_std'][ii]

        dt_str = stats1['DATE_UTC'][ii] + ' ' + stats1['TIME_UTC'][ii]
        dt_utc = datetime.datetime.strptime(dt_str, '%Y-%m-%d %H:%M:%S')

        times2 = []
        for jj in range(len(stats2)):
            dt_str2 = stats2['DATE_UTC'][jj] + ' ' + stats2['TIME_UTC'][jj]
            dt_utc2 = datetime.datetime.strptime(dt_str2, '%Y-%m-%d %H:%M:%S')
            diff = abs(dt_utc-dt_utc2)
            times2.append(diff)

        min_id = np.argmin(times2)
        data_2[ii] = stats2[comp_col][min_id]
        data_2_err[ii] = stats2['emp_fwhm_std'][min_id]

        time2_str = stats2['DATE_UTC'][min_id] + ' ' + stats2['TIME_UTC'][min_id]
        time2_utc = datetime.datetime.strptime(time2_str, '%Y-%m-%d %H:%M:%S')
        diff_time = time2_utc - dt_utc
        ave_time = dt_utc + (diff_time/2)


        date_time = ((ave_time.isoformat().split('T')))
        UTC_TIME[ii] = str(date_time[1])[:8] #cuts of .5 second if present
        UTC_DATE[ii] = str(date_time[0])
	
    if a == 1:
        data_open = data_2
        data_open_err = data_2_err
        data_closed = data_1
        data_closed_err = data_1_err
    else:
        data_open = data_1
        data_open_err = data_1_err
        data_closed = data_2
        data_closed_err = data_2_err

    return  UTC_TIME, UTC_DATE, data_open, data_closed, data_open_err, data_closed_err


def add_filter(file, filter_name):
    
    ## Edits fits file header and adds a 'FILTER' card with the given value 'filter_name'.
    ## Meant to use on runs before run 5 that didn't have filter info in header so updated reduction pipeline can 
    ## pull info from header rather than have it hardcoded
    
    hdulist = fits.open(file, mode='update')
    hdulist[0].header.set('FILTER', filter_name)
    hdulist.flush()
    hdulist.close()
    
    return



def week_table(data_dir_root, stats_dir_end, labels):

    #creates a table for a week of observing showing nightly averages:
    # UT Date | average open fwhm | average closed fwhm | average DIMM seeing | average MASS seeing
    #inputs:
         #data_dir_root: directory containing nightly data files, e.g.:  "/Users/fatimaabdurrahman/Desktop/Research/RUN4/"
         #stats_dir_root: the rest of the path to the stats directory after date, e.g.: "/FLI/reduce/stats/"
         #labels: a list of two-element entries, which each have a date string and a closed stat file name, e.g.: [['20170214', 'closed'], ['20170215', 'closed'], ['20170216', 'closeda'],['20170217', 'closeda'], ['20170218', 'closeda']]

    # function writes a file "nightly_aves.fits" in the data_dir_root directory of output table
         
    date = []
    open_ave = np.zeros(len(labels), dtype=float)
    closed_ave = np.zeros(len(labels), dtype=float)
    DIMM_ave = np.zeros(len(labels), dtype=float)
    MASS_ave = np.zeros(len(labels), dtype=float)

    for i in range(len(labels)):
        open_file = data_dir_root+labels[i][0]+stats_dir_end+'stats_open.fits'
        print("Open File: ", open_file)
        closed_file = data_dir_root+labels[i][0]+stats_dir_end+'stats_'+labels[i][1]+'.fits'
        print("Closed File: ", closed_file)
        open_data = Table.read(open_file)
        closed_data = Table.read(closed_file)
        print(open_data)
        print(closed_data)

        scale = 0.04

        o_data = np.array(open_data['emp_fwhm'])
        o_binfac = np.array(open_data['BINFAC'])
        o_filt =  plot_stats.filter2wv(np.array(open_data['FILTER']))
        open_ave[i] = scale * np.mean(o_data * o_binfac * (500/o_filt)**(-1/5))
        
        c_data = np.array(closed_data['emp_fwhm'])
        c_binfac = np.array(closed_data['BINFAC'])
        c_filt =  plot_stats.filter2wv(np.array(closed_data['FILTER']))
        closed_ave[i] = scale * np.mean(c_data * c_binfac * (500/c_filt)**(-1/5))
        
        MASS_ave[i] = np.mean(np.concatenate((np.array(closed_data['MASS']), np.array(open_data['MASS'])), axis=0))
        DIMM_ave[i] = np.mean(np.concatenate((np.array(closed_data['DIMM']), np.array(open_data['DIMM'])), axis=0))
        
        string = labels[i][0][0:4]+"-"+labels[i][0][4:6]+"-"+labels[i][0][6:8]
        date.append(string)

    t = Table([date, open_ave, closed_ave, DIMM_ave, MASS_ave], names=('HST_Date', 'ave_open_fwhm', 'ave_closed_fwhm', 'ave_DIMM', 'ave_MASS'))
    t.write(data_dir_root+"nightly_aves.fits", overwrite=True)

    
    return



def read_starlist(starlist_file):
    data=[]
    f = open(starlist_file, 'r')
    for line in f:
        data.append(line)
    f.close()

    new_data =[]
    for line in data:
        new_line = []
        for element in line.split(" "):
            if element != '':
                new_line.append(element)
        new_data.append(new_line)

    array = np.array(new_data)

    x_cents_str = array[:,1][1:]
    y_cents_str = array[:,2][1:]
    x_fwhm_str = array[:,11][1:]
    y_fwhm_str = array[:,12][1:]
    roundness_str = array[:,4][1:]

    x_cents=np.zeros(len(x_cents_str), dtype=float)
    y_cents=np.zeros(len(x_cents_str), dtype=float)
    x_fwhm=np.zeros(len(x_cents_str), dtype=float)
    y_fwhm=np.zeros(len(x_cents_str), dtype=float)
    roundness=np.zeros(len(x_cents_str), dtype=float)

    for i in range(len(x_cents_str)):
        if 0<float(x_fwhm_str[i])<45 and 0<float(y_fwhm_str[i])<45:
            x_cents[i]=float(x_cents_str[i])
            y_cents[i]=float(y_cents_str[i])
            x_fwhm[i]=float(x_fwhm_str[i])
            y_fwhm[i]=float(y_fwhm_str[i])
            roundness[i]=float(roundness_str[i])



    fwhm = (x_fwhm + y_fwhm)/2

    #define coordinate center of image and star distances from center
    y_cent = 3066/2
    x_cent = 4088/2
    dist=(((x_cents-x_cent)**2)+((y_cents-y_cent)**2))**0.5
    
    return x_cents, y_cents, fwhm, x_fwhm, y_fwhm, roundness, dist

def combine_all_stats(data_root='/Users/jlu/data/imaka/'):
    stats_files = [
                   'RUN3/20170110/FLI/reduce/stats/stats_closed_mdp_alt.fits',
                   'RUN3/20170110/FLI/reduce/stats/stats_open_mdp_alt.fits',
                   'RUN3/20170111/FLI/reduce/stats/stats_closed_mdp_alt.fits',
                   'RUN3/20170111/FLI/reduce/stats/stats_open_mdp_alt.fits',
                   'RUN3/20170113/FLI/reduce/stats/stats_closed1_mdp_alt.fits',
                   'RUN3/20170113/FLI/reduce/stats/stats_open1_mdp_alt.fits',
                   'RUN3/20170113/FLI/reduce/stats/stats_closed2_mdp_alt.fits',
                   'RUN3/20170113/FLI/reduce/stats/stats_open2_mdp_alt.fits',
                   'RUN4/20170214/FLI/reduce/stats/stats_open_mdp_alt.fits',
                   'RUN4/20170214/FLI/reduce/stats/stats_closed_mdp_alt.fits',
                   # 'RUN4/20170215/FLI/reduce/stats/stats_open_mdp_alt.fits',
                   # 'RUN4/20170215/FLI/reduce/stats/stats_closed_mdp_alt.fits',
                   'RUN4/20170216/FLI/reduce/stats/stats_open_mdp_alt.fits',
                   'RUN4/20170216/FLI/reduce/stats/stats_closeda_mdp_alt.fits',
                   'RUN4/20170217/FLI/reduce/stats/stats_open_mdp_alt.fits',
                   'RUN4/20170217/FLI/reduce/stats/stats_closeda_mdp_alt.fits',
                   'RUN5/20170517/FLI/reduce/stats/stats_open_mdp_alt.fits',
                   'RUN5/20170518/FLI/reduce/stats/stats_closed_mdp_alt.fits',
                   'RUN5/20170518/FLI/reduce/stats/stats_open_mdp_alt.fits',
                   'RUN5/20170519/FLI/reduce/stats/stats_closed_mdp_alt.fits',
                   'RUN5/20170519/FLI/reduce/stats/stats_open_mdp_alt.fits',
                   'RUN5/20170520/FLI/reduce/stats/stats_closedA_mdp_alt.fits',
                   'RUN5/20170520/FLI/reduce/stats/stats_open_mdp_alt.fits',
                   'RUN5/20170521/FLI/reduce/stats/stats_closed_mdp_alt.fits',
                   'RUN5/20170521/FLI/reduce/stats/stats_open_mdp_alt.fits'
                   # 'RUN5/20170522/FLI/reduce/stats/stats_open_mdp_alt.fits',
                   # 'RUN5/20170522/FLI/reduce/stats/stats_closed_mdp_alt.fits',
                   ]

    for ii in range(len(stats_files)):
        print('Appending ', stats_files[ii])
        _stats = Table.read(data_root + stats_files[ii])

        foo1 = stats_files[ii].index('stats_')
        foo2 = stats_files[ii].index('_mdp')
        reconstruct = stats_files[ii][foo1+6:foo2]

        if 'closed' in stats_files[ii]:
            _stats['loop_status'] = 'closed'
            _stats['recon'] = reconstruct
        else:
            _stats['loop_status'] = 'open'
            _stats['recon'] = 'open'

        if ii == 0:
            final_table = _stats
        else:
            final_table = vstack([final_table, _stats])


    # Also make a matched table for open vs. closed plots.
    # While we are at it, calculate a few extra columns.

    # Calculate GL seeing from MASS/DIMM
    gl_see = (final_table['DIMM']**(5.0/3.0) - final_table['MASS']**(5.0/3.0))**(3.0/5.0)
    final_table['GL_MASSDIMM'] = gl_see

    # Fix BINFAC
    idx = np.where(final_table['BINFAC'].mask == True)[0]
    final_table['BINFAC'][idx] = 3.0

    # Fix Filter
    idx = np.where(final_table['FILTER'].mask == True)[0]
    print('replacing filter', len(idx))
    final_table['FILTER'][idx] = 'R'
    
    # Convert empirical FWHM measurements to " accounting for binning.
    scale = 0.04
    final_table['emp_fwhm'] *= scale * final_table['BINFAC']

    opcl = final_table['loop_status'] == 'open'
    i_op = np.where(opcl == True)
    i_cl = np.where(opcl == False)
    s_op = final_table[i_op]
    s_cl = final_table[i_cl]
    
    # Files to store the open-only and closed-only tables. 
    op_file = data_root + 'stats_aggregate_open_mdp_alt.fits'
    cl_file = data_root + 'stats_aggregate_closed_mdp_alt.fits'
    
    # Make an index array
    idx_all = np.arange(len(final_table))
    idx_op = idx_all[i_op]
    idx_cl = idx_all[i_cl]

    s_op['Index'] = idx_op
    s_cl['Index'] = idx_cl

    s_op.write(op_file, overwrite=True)
    s_cl.write(cl_file, overwrite=True)

    # Matching
    time, date, idx_op_match, idx_cl_match, err1, err2 = match_cols(op_file, cl_file, 'Index')
    idx_op_match = idx_op_match.astype('int')
    idx_cl_match = idx_cl_match.astype('int')

    # Horizontally stack
    s_op_match = final_table[idx_op_match]
    s_cl_match = final_table[idx_cl_match]

    s_opcl_match = hstack([s_op_match, s_cl_match], join_type='exact', table_names=['o', 'c'])

    final_table.write(data_root + 'stats_aggregate_mdp_alt.fits', overwrite=True)
    s_opcl_match.write(data_root + 'stats_aggregate_matched_mdp_alt.fits', overwrite=True)

    return
    
    

