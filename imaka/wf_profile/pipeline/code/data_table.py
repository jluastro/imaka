# Eden McEwen
# July 16, 2020
# made to store functions used in generate data table

#imports
import os
import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.table import Table
from astropy.time import Time
from datetime import datetime, date, time, timezone
import pytz

## Functions for Timezones

def hst_to_utc(date):
    """
    Returns the Universal Coordinated Time calculated from Hawaii Standard Time
    """
    tz = pytz.timezone('US/Hawaii')
    return tz.normalize(tz.localize(date)).astimezone(pytz.utc)

def str_to_utc(d_str):
    d_dt = str_to_datetime(d_str)
    return d_dt.astimezone(pytz.utc)

def str_to_datetime(string):
    #assuming this is coming in in the format of the data fits
    string = string.replace("-", "")
    l_str = [string[:4]] + [string[i:i+2] for i in range(4, len(string), 2)]
    l_int = [int(i) for i in l_str]
    return datetime(l_int[0], l_int[1], l_int[2], l_int[3], l_int[4], l_int[5])

######################################################
# Pulls from a list of fits and returns a datatable
######################################################

#returns dataframe for all run files given
def df_gen_main(run_files):
    #generate empty dataframe
    df_main = pd.DataFrame()
    for run in run_files:
        dates = read_file(run)[1:]
        df_run = pd.DataFrame()

        # finda all fits associated with this run:
        rname = run.replace(run_input, "").replace(".txt", "")
        r = rname.replace("RUN", "")
        #run_fits = []
        for d in dates:
            fits_out_p = "/home/emcewen/out/"+d+"/fits/"
            d_fits = []
            if os.path.isdir(fits_out_p):
                fs = os.listdir(fits_out_p)
                fits_out_files = [fits_out_p + fn for fn in fs if fnmatch.fnmatch(fn, '*.fits')]
                d_fits.extend(fits_out_files)
            df_date = data_out_fits(d_fits)
            df_date["runfile"] = run
            df_date["run"] = r
            df_date.to_csv("csv/{}.csv".format(d))
            df_run = df_run.append(df_date, ignore_index = True)
        df_run.to_csv("csv/{}.csv".format(rname))
        df_main=df_main.append(df_run, ignore_index = True)
    # returns dataframe
    return df_main






#pulls from one fits file
def data_out_fits(out_fits):
    # input a list of all the corr data products
    dfObj = pd.DataFrame(columns=['dataname', 'OBSDATE', 'DATETIME', 
                                  'infits', 'outfits', "outdir", 'targetfile', 'nwfs',
                                  'tmax', 'ssub', 'ttsub',
                                  'wfs', 'TSAMPLE', 'FSAMPLE', 'TEXP', 'EMGAIN',
                                 'ra', 'dec', 'mag', 'datavalid'])
    # this is one data file
    for f in out_fits:
        if not os.path.exists(f):
            continue
        hdulist = fits.open(f)
        hdr = hdulist[0].header
        
        #files and identifiers 
        dataname = hdr.get('DATANAME')
        obsdate = hdr.get('OBSDATE')
        datetime = hdr.get('DATETIME')
        
        #file info
        infits = hdr.get('DATAFILE')
        outfits = f
        outdir = hdr.get('OUTPATH')
        targetfile = hdr.get('TARGET')
        
        #data product information
        tmax = hdr.get('TMAX')
        s_sub = hdr.get('SSUB')
        tt_sub = hdr.get('TTSUB')
        
        # Retrieve info from all wfs
        # granted, this isn't actually accurate :///
        ccor = hdr.get("CCOR") 
        acor = hdr.get("ACOR")
        
        nwfs = hdr.get('NWFS')
        #all info common among wfs
        data_dic = {'dataname':dataname, 'OBSDATE':obsdate, 'DATETIME': datetime, 
                                  'infits': infits, 'outfits': outfits, "outdir":outdir, 'targetfile':targetfile, 
                                  'nwfs':nwfs, 'tmax': tmax, 'ssub': s_sub, 'ttsub':tt_sub}
        
        
        for wfs in range(nwfs):
            hdr_t = hdulist[wfs+1].header
            # camera info
            TSAMPLE = hdr_t.get('TSAMPLE') 
            FSAMPLE = hdr_t.get('FSAMPLE') 
            TEXP = hdr_t.get('TEXP')
            EMGAIN = hdr_t.get('EMGAIN')
            # pointing info
            ra = hdr_t.get('ra')
            dec = hdr_t.get('dec')
            mag = hdr_t.get('mag')
            datavalid = hdr_t.get('WFSVALID')
            
            wfs_dic = {'wfs':wfs, 'TSAMPLE':TSAMPLE, 'FSAMPLE':FSAMPLE, 'TEXP':TEXP, 'EMGAIN':EMGAIN,
                                 'ra':ra, 'dec':dec, 'mag':mag, 'datavalid':datavalid}
            
            # adding in a column to the data frame 
            tmp = data_dic.copy()
            tmp.update(wfs_dic)
            dfObj = dfObj.append(tmp, ignore_index=True)
            
    return dfObj














##############################################
################### Pulling from CSVs


def data_wind_cfht_slow(datetimes, wind_CFHT_p):
    #generate a list of wind_cfht dat for supplied UT date and times
    df_cft = pd.DataFrame(columns=['cfht_dt', 'cft_wspd', 'cft_wdir', 
                                  'cft_avg_wspd', 'cft_avg_wdir'])
    for dt in datetimes:
        print("Date data: ", dt)
        year = dt.year
        wind_CFHT = wind_CFHT_p + "cfht-wx.{}.csv".format(year)
        df = pd.read_csv(wind_CFHT, memory_map=True)
        
        names, data = zip(*pd.to_datetime(df["datetime"]).iteritems())
        c_datetime = min(data, key=lambda x:abs(x-dt))
        print("min: ", c_datetime)
        
        #df["datetime"] = pd.to_datetime(df["datetime"])
        #date_idx = np.abs(df['datetime']-dt).idxmin()
        #print("min 2: ", date_idx)
        
        is_date = pd.to_datetime(df["datetime"])==c_datetime
        df_date = df[is_date]
        
        dt_dic = {'cfht_dt':df_date, 'cft_wspd': df.loc[df_date,'wind_speed(kts)'], 
                  'cft_wdir': df.loc[df_date,'wind_direction(dec)'],
                  'cft_avg_wspd': df['wind_speed(kts)'].iloc[df_date-avg_dt : df_date+avg_dt].mean(),
                  'cft_avg_wdir': df['wind_direction(dec)'].iloc[df_date-avg_dt : df_date+avg_dt].mean(),
        }
        
        df_cft = df_cft.append(dt_dic, ignore_index=True)
    return df_cft

# importing the CFHT ground data
# this means converting times, I think
# need to find file, convert UT time to HST, pull from the line that has the closest time 
# them pull 5 and 6 (wind_speed(kts) wind_direction(dec))
        
def data_wind_cfht(datetimes, wind_CFHT_p):
    df_cft = pd.DataFrame(columns=['cfht_dt', 'cft_wspd', 'cft_wdir', 
                                  'cft_avg_wspd', 'cft_avg_wdir'])
    avg_dt = 10
    year = 0  
    for dt in datetimes:
        if dt.year != year:
            year = dt.year
            wind_CFHT = wind_CFHT_p + "cfht-wx.{}.csv".format(year)
            df = pd.read_csv(wind_CFHT, memory_map=True, index_col ='datetime')
        df.index = pd.to_datetime(df.index)
        
        df = df.groupby(level=0).first()
               
        #getting the specific row
        spot = df.index.get_loc(dt, method='nearest')
        #getting the rown index
        idx = df.index[df.index.get_loc(dt, method='nearest')]
        
        dt_dic = {'cfht_dt':idx, 'cft_wspd': df.loc[idx,'wind_speed(kts)'], 
                  'cft_wdir': df.loc[idx,'wind_direction(dec)'],
                  'cft_avg_wspd': df['wind_speed(kts)'].iloc[spot-avg_dt : spot+avg_dt].mean(),
                  'cft_avg_wdir': df['wind_direction(dec)'].iloc[spot-avg_dt : spot+avg_dt].mean(),
        }
        df_cft = df_cft.append(dt_dic, ignore_index=True)
    return df_cft


################## Plotting code 







