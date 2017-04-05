from datetime import datetime
import numpy as np
from astropy.table import Table, Column
from astropy.io import fits
import datetime
import os




def append_alt_data(stats_file, alt_file, date):  

    #This function takes a stats file and a file of Olivier's altitude seeing data and matches
    #them, making a new stats file.  There aren't as many altitude data points as frames in stats
    #files, missing values are zeros
    
    #Inputs: 
        #stats_file: the original stats file
        #alt_file: altitide data file
        #date: a string of the UT date in format yyyymmdd, eg '20170112'
        
    #outputs: a new fits file with '_alt' added to the original file name
    
    #####
    
    ## Read in original stats table and altitude data

    stats = Table.read(stats_file)
    new_table = fits.getdata(alt_file) #altitude data

    time = new_table[:,0]

    ##Make some empty columns for new data

    alt_date = np.zeros(len(stats), dtype=str)
    alt_time = np.zeros(len(stats), dtype=str)
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


    # Convert times in alititude data list to ut datetime objects to compare to
    # original stats table times

    for ii in range(len(new_table)):
        hour = int(time[ii]//1)
        minute_dec = (time[ii]-hour)*60
        minute = int(minute_dec//1)
        second = int(round((minute_dec-minute)*60, 0))
        if hour >= 24:
            hour = hour-24
        HST_date_str = date[0:4]+"-"+date[4:6]+"-"+date[6:8]
        HST_time_str = str(hour)+":"+str(minute)+":"+str(second)
        comp_time = datetime.datetime.strptime(HST_date_str+","+HST_time_str, '%Y-%m-%d,%H:%M:%S')
        noon = datetime.time(12, 0, 0)
        del_day = datetime.timedelta(days=1)
        if comp_time.time() < noon: #comparison time
            comp_time += del_day

        ##For one altitude time, look through original stats table and find closest    
        time_diffs = []
        for jj in range(len(stats)):
            hst_str = stats[jj]['DATE_HST']+","+stats[jj]['TIME_HST']
            time_obj = datetime.datetime.strptime(hst_str, "%Y-%m-%d,%H:%M:%S")
            if time_obj > comp_time:
                time_diff = time_obj - comp_time
            else:
                time_diff = comp_time - time_obj
            time_diffs.append(time_diff)

        min_index = np.argmin(time_diffs)

        ## Add matching data to new table
        alt_date[min_index] = HST_date_str
        alt_time[min_index] = HST_time_str
        int_seeing[min_index] = new_table[:,1][ii]
        alt_0[min_index] = new_table[:,2][ii]
        alt_60[min_index] = new_table[:,3][ii]
        alt_120[min_index] = new_table[:,4][ii]
        alt_180[min_index] = new_table[:,5][ii]
        alt_240[min_index] = new_table[:,6][ii]
        alt_300[min_index] = new_table[:,7][ii]
        alt_360[min_index] = new_table[:,8][ii]
        alt_420[min_index] = new_table[:,9][ii]
        alt_480[min_index] = new_table[:,10][ii]
        alt_560[min_index] = new_table[:,11][ii]
        alt_4000[min_index] = new_table[:,12][ii]

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






def match_cols(base_file, comp_file, comp_col):
    
    #time matches colums of data from different stats files, e.g. matches FWHM column from 
    #a given open stats file to those in a closed file from the same night.
    
    #inputs:
        #base_file: the stats file that you are matching to
        #comp_file: the comparison file to match
        #comp_col: the column you want (e.g., NEA1, emp_FWHM, etc)
        
    #outputs:
        #UTC_TIME: a time average of each two matched data points
        #UTC_DATE: corresponding date to above time
        #data_1: the column of whatever comp_col described from the base file
        #data_2: the column of whatever comp_col described from the comparison file 
            
    stats1 = Table.read(base_file)
    stats2 = Table.read(comp_file)
        
    if len(stats1) > len(stats2):
        stats1 = Table.read(comp_file)
        stats2 = Table.read(base_file)

    UTC_TIME = np.zeros(len(stats1), dtype='S15')
    UTC_DATE = np.zeros(len(stats2), dtype='S15')
    data_open = np.zeros(len(stats1),dtype=float)
    data_closed = np.zeros(len(stats1), dtype=float)

    for ii in range(len(stats1)):
        data_open[ii] = stats1[comp_col][ii]

        dt_str = stats1['DATE_UTC'][ii] + ' ' + stats1['TIME_UTC'][ii]
        dt_utc = datetime.datetime.strptime(dt_str, '%Y-%m-%d %H:%M:%S')

        times2 = []
        for jj in range(len(stats2)):
            dt_str2 = stats2['DATE_UTC'][jj] + ' ' + stats2['TIME_UTC'][jj]
            dt_utc2 = datetime.datetime.strptime(dt_str2, '%Y-%m-%d %H:%M:%S')
            diff = abs(dt_utc-dt_utc2)
            times2.append(diff)

        min_id = np.argmin(times2)
        data_closed[ii] = stats2[comp_col][min_id]

        time2_str = stats2['DATE_UTC'][min_id] + ' ' + stats2['TIME_UTC'][min_id]
        time2_utc = datetime.datetime.strptime(time2_str, '%Y-%m-%d %H:%M:%S')
        diff_time = time2_utc - dt_utc
        ave_time = dt_utc + (diff_time/2)


        date_time = ((ave_time.isoformat().split('T')))
        UTC_TIME[ii] = str(date_time[1])[:8] #cuts of .5 second if present
        UTC_DATE[ii] = str(date_time[0])

    return  UTC_TIME, UTC_DATE, data_open, data_closed





