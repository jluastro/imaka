import os, sys
from . import util
from astropy.io import fits
from astropy.stats import sigma_clip, sigma_clipped_stats
import numpy as np
import pdb

module_dir = os.path.dirname(__file__)

def makedark(dark_files, output):
    """
    Make dark image for data. Makes a calib/ directory
    and stores all output there. All output and temporary files
    will be created in a darks/ subdirectory. 

    files: integer list of the files. Does not require padded zeros.
    output: output file name. Include the .fits extension.
    """
    # Catch if there are no files sent in.
    if len(dark_files) == 0:
        raise RuntimeError('No files passed into makedark')
        
    # Output directory is the same as the input directory.
    dark_dir = os.path.dirname(dark_files[0]) + '/'
    
    _out = dark_dir + output
    _outlis = dark_dir + 'dark.lis'
    util.rmall([_out, _outlis])

    # Print out useful info.
    print('Making dark frame: ' + output)
    for ff in range(len(dark_files)):
        filedir, filename = os.path.split(dark_files[ff])
        print('\t' + filename)

    # Read in the dark images
    print('\nReading in files')
    darks = np.array([fits.getdata(dark_file) for dark_file in dark_files])

    # Sigma clip  (about the median)
    print('Sigma clipping and median combining (3 sigma lower, 2 sigma upper, 5 iterations) about the median.')
    dark_mean, dark_median, dark_std = sigma_clipped_stats(darks, sigma_lower=3, sigma_upper=2, iters=5, axis=0)
    dark = dark_median
    
    # Save the output files.
    # Get the header from the first image.
    hdr = fits.getheader(dark_files[0])
    fits.writeto(_out, dark.data, header=hdr)

    outlis_file = open(_outlis, 'w')
    for ff in range(len(dark_files)):
        outlis_file.write(dark_files[ff] + '\n')
    outlis_file.close()
    
    return

def makeflat(flat_files, dark_files, output):
    """
    Make flats image for data. 

    flat_files: array fits file names
    dark_files: array fits file names matched, 1:1, to the flats
    """
    # Catch if there are no files sent in.
    if (len(flat_files) == 0) or (len(dark_files) == 0):
        raise RuntimeError('No files passed into makeflat')
        
    # Output directory is the same as the input directory.
    dark_dir = os.path.dirname(dark_files[0]) + '/'
    flat_dir = os.path.dirname(flat_files[0]) + '/'
    
    _out = flat_dir + output
    _outlis = flat_dir + 'flat.lis'
    util.rmall([_out, _outlis])

    # Print out useful info.
    print('Making flat frame: ' + output)
    for ff in range(len(dark_files)):
        filedir1, filename1 = os.path.split(dark_files[ff])
        filedir2, filename2 = os.path.split(flat_files[ff])
        print('\t Flat: ' + filename2 + '  Dark: ' + filename1)

    # Read in the dark images
    print('\nReading in dark files')
    darks = np.array([fits.getdata(dark_file) for dark_file in dark_files])
    
    # Read in the flat images
    print('\nReading in flat files')
    flats = np.array([fits.getdata(flat_file) for flat_file in flat_files])

    flats_ds = flats - darks

    # Normalize each of the flats:
    norms = flats_ds.sum(axis=(1,2))

    pdb.set_trace()
    
    for ff in range(len(flat_files)):
        flats_ds[ff]
    # Sigma clip  (about the median)
    print('Sigma clipping (3 sigma, 2 iterations) about the median.')
    darks_masked = sigma_clip(darks, sigma=3, iters=2, axis=0)
    print('\t Masked {0:d} total pixels'.format(darks_masked.mask.sum()))

    # Median combine the masked images.
    print('Median combining')
    dark = np.ma.median(darks_masked, axis=0)
    
    # Save the output files.
    # Get the header from the first image.
    hdr = fits.getheader(dark_files[0])
    fits.writeto(_out, dark.data, header=hdr)

    outlis_file = open(_outlis, 'w')
    for ff in range(len(dark_files)):
        outlis_file.write(dark_files[ff] + '\n')
    outlis_file.close()
    
    return



# def makeflat(onFiles, offFiles, output, normalizeFirst=False):
#     """
#     Make flat field image for NIRC2 data. Makes a calib/ directory
#     and stores all output there. All output and temporary files
#     will be created in a flats/ subdirectory. 
    
#     onFiles: integer list of lamps ON files. Does not require padded zeros.
#     offFiles: integer list of lamps OFF files. Does not require padded zeros.

#     If only twilight flats were taken (as in 05jullgs), use these flats as
#     the onFiles, and use 0,0 for offFiles. So the reduce.py file should look
#     something like this: onFiles = range(22, 26+1) and offFiles = range(0,0)
#     The flat will then be made by doing a median combine using just the
#     twilight flats.
#     output: output file name. Include the .fits extension.
#     """
#     redDir = os.getcwd() + '/'
#     curDir = redDir + 'calib/'
#     flatDir = util.trimdir(curDir + 'flats/')
#     rawDir = util.trimdir(os.path.abspath(redDir + '../raw') + '/')

#     util.mkdir(curDir)
#     util.mkdir(flatDir)

#     _on = flatDir + 'lampsOn.fits'
#     _off = flatDir + 'lampsOff.fits'
#     _norm = flatDir + 'flatNotNorm.fits'
#     _out = flatDir + output
#     _onlis = flatDir + 'on.lis'
#     _offlis = flatDir + 'off.lis'
#     _onNormLis = flatDir + 'onNorm.lis'

#     util.rmall([_on, _off, _norm, _out, _onlis, _offlis, _onNormLis])

#     lampson = [rawDir + 'n' + str(i).zfill(4) + '.fits' for i in onFiles]
#     lampsoff = [rawDir + 'n' + str(i).zfill(4) + '.fits' for i in offFiles]
#     lampsonNorm = [flatDir + 'norm' + str(i).zfill(4) + '.fits' for i in onFiles]
#     util.rmall(lampsonNorm)

#     if (len(offFiles) != 0):
#         f_on = open(_onlis, 'w')
#         f_on.write('\n'.join(lampson) + '\n')
#         f_on.close()
#         f_on = open(_offlis, 'w')
#         f_on.write('\n'.join(lampsoff) + '\n')
#         f_on.close()
#         f_onn = open(_onNormLis, 'w')
#         f_onn.write('\n'.join(lampsonNorm) + '\n')
#         f_onn.close()
    
#         # Combine to make a lamps on and lamps off
#         ir.unlearn('imcombine')
#         ir.imcombine.combine = 'median'
#         ir.imcombine.reject = 'sigclip'
#         ir.imcombine.nlow = 1
#         ir.imcombine.nhigh = 1
#         ir.imcombine('@' + _offlis, _off)

#         # Check if we should normalize individual flats first
#         # such as in the case of twilight flats.
#         if normalizeFirst:
#             f_on = open(_offlis, 'w')
#             f_on.write('\n'.join(lampsoff) + '\n')
#             f_on.close()
            
#             # Subtract "off" from individual frames
#             ir.imarith('@'+_onlis, '-', _off, '@'+_onNormLis)
            
#             # Scale them and combine
#             ir.imcombine.scale = 'median'
#             ir.imcombine('@' + _onNormLis, _norm)
#         else:
#             # Combine all "on" frames
#             ir.imcombine('@' + _onlis, _on)

#             # Now do lampsOn - lampsOff
#             ir.imarith(_on, '-', _off, _norm)


#         # Normalize the final flat
#         ir.module.load('noao', doprint=0, hush=1)
#         ir.module.load('imred', doprint=0, hush=1)
#         ir.module.load('generic', doprint=0, hush=1)
#         orig_img = fits.getdata(_norm)
#         orig_size = (orig_img.shape)[0]
#         if (orig_size >= 1024):
#             flatRegion = '[100:900,513:950]'
#         else:
#             flatRegion = ''
#         ir.normflat(_norm, _out, sample=flatRegion)

#     else:
#         f_on = open(_onlis, 'w')
#         f_on.write('\n'.join(lampson) + '\n')
#         f_on.close()

#         # Combine twilight flats
#         ir.unlearn('imcombine')
#         ir.imcombine.combine = 'median'
#         ir.imcombine.reject = 'sigclip'
#         ir.imcombine.nlow = 1
#         ir.imcombine.nhigh = 1
#         if normalizeFirst:
#             # Scale them
#             ir.imcombine.scale = 'median'
#         ir.imcombine('@' + _onlis, _norm)

#         # Normalize the flat
#         ir.module.load('noao', doprint=0, hush=1)
#         ir.module.load('imred', doprint=0, hush=1)
#         ir.module.load('generic', doprint=0, hush=1)
#         flatRegion = '[100:900,513:950]'
#         ir.normflat(_norm, _out, sample=flatRegion)

# def makemask(dark, flat, output):
#     """Make bad pixel mask for NIRC2 data. Makes a calib/ directory
#     and stores all output there. All output and temporary files
#     will be created in a masks/ subdirectory. 
    
#     @param dark: The full relative path to a dark file. This is used to
#         construct a hot pixel mask. Use a long (t>20sec) exposure dark.
#     @type dark: str
#     @param flat: The full relative path to a flat file. This is used to 
#         construct a dead pixel mask. The flat should be normalized.
#     @type flat: str
#     @param output: output file name. This will be created in the masks/
#         subdirectory.
#     @type output: str
#     """
#     redDir = os.getcwd() + '/'
#     calDir = redDir + 'calib/'
#     maskDir = util.trimdir(calDir + 'masks/')
#     flatDir = util.trimdir(calDir + 'flats/')
#     darkDir = util.trimdir(calDir + 'darks/')
#     rawDir = util.trimdir(os.path.abspath(redDir + '../raw') + '/')
#     dataDir = util.trimdir(os.path.abspath(redDir + '../..') + '/')

#     util.mkdir(calDir)
#     util.mkdir(maskDir)

#     _out = maskDir + output
#     _dark = darkDir + dark
#     _flat = flatDir + flat
#     _nirc2mask = module_dir + '/masks/nirc2mask.fits'

#     util.rmall([_out])

#     # Make hot pixel mask
#     whatDir = redDir + dark
#     print(whatDir)

#     text_output = ir.imstatistics(_dark, fields="mean,stddev", 
# 				  nclip=10, format=0, Stdout=1)
#     print(text_output)
#     values = text_output[0].split()
#     hi = float(values[0]) + (10.0 * float(values[1]))

#     img_dk = fits.getdata(_dark)
#     hot = img_dk > hi

#     # Make dead pixel mask
#     text_output = ir.imstatistics(_flat, fields="mean,stddev", 
# 				  nclip=10, format=0, Stdout=1)
#     values = text_output[0].split()
#     #lo = float(values[0]) - (15.0 * float(values[1]))
#     # If flat is normalized, then lo should be set to 0.5
#     lo = 0.5
#     hi = float(values[0]) + (15.0 * float(values[1]))

#     img_fl = fits.getdata(_flat)
#     dead = np.logical_or(img_fl > hi, img_fl < lo)
    
#     # We also need the original NIRC2 mask (with cracks and such)
#     nirc2mask = fits.getdata(_nirc2mask)

#     # Combine into a final supermask. Use the flat file just as a template
#     # to get the header from.
#     ofile = fits.open(_flat)
    
#     if ((hot.shape)[0] == (nirc2mask.shape)[0]):
#         mask = hot + dead + nirc2mask
#     else:
#         mask = hot + dead
#     mask = (mask != 0)
#     unmask = (mask == 0)
#     ofile[0].data[unmask] = 0
#     ofile[0].data[mask] = 1
#     ofile[0].writeto(_out, output_verify='silentfix')
    

