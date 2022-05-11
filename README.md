imaka
====================

INSTALL (from git)

## Dependencies:

 * python (preferably via Anaconda, as it includes necessary pacakges like astropy).
 * astropy
 * photutils
 * ccdproc
 * flystar (<https://github.com/jluastro/imaka>)


## Summary of Basic reduction and Analysis

### Part 1: Make sky and flat images and use them to reduce science images

**Skies and Darks:** These are from calibration data taken at the beginning and sometimes end of the night.
- `make_flat()`: takes in Darks and creates a flat from them
  - **Requires:** input of flat numbers and locations
  - **Outputs:** 
    - `sta/reduce/calib/flats.fits` 
    - `sta/reduce/calib/flats.list`
  - Will scan flat files if needed (in flats original location)
  - Creates a mask from this flat to cut noisy edges: sta/reduce/calib/mask.fits
- `make_sky()`: Takes in sky files and combines in a master sky
  - **Requires:** input of sky numbers and locations
  - **Outputs**: 
    - `sta/reduce/sky/fld2_sky.fits`
    - `sta/reduce/sky/fld2_sky.list`
  - scans sky fits if needed (in fits original location)

**Reducing**: Reduction takes into account both skies and darks for the final science image
- `reduce_fld2()`
  - **Requires:** flats and sky fits
  - Will treat overscan regions for science images: `sta/Fld2/ *_scan.fits`
  - Will clean scanned files: `sta/reduce/Fld2/*_scan_clean.fits`

### Part 2: Find stars in clean images and calculate stats

**Star Finding:** Uses the DAO Starfinder
- `find_stars_fld2()`
  - **Requires**: mask, cleaned image files
  - **Outputs**:
    - `sta/reduce/Fld2/*_scan_clean_stars.txt` : a list of all sources centroids, peak brightness, and fwhm, etc.
    - `sta/reduce/Fld2/*_scan_clean_pfs_mod.fits` : average model psf
    - `sta/reduce/Fld2/*_scan_clean_pfs_obs.fits` : average observed psf

**Star Stats:** 
This is run from the reduce file on a list of clean files. Each starlist found will give 
- `calc_star_stats()`
  - **Requires**: stars.txt files, cleaned images
  - **Calculates**: From starlists
    - emperical FWHM
    - encircled energy (EE) at 25, 50 and 80 
    - NEA (? encircled energy and plate scale radius…)
  - **Outputs**:
    - `sta/reduce/Fld2/*_stars_stats.fits` : saved files of parameters
    - `sta/reduce/Fld2/stats/stats<KEY>.FITS` : summary of all files passed in
    - `sta/reduce/Fld2/ee/ee*` : encircled energy profile saved
- `fit_moffat()` : A PSF fit that is more extended than a gaussian.
  - **Requires**: start lists, star stats
  - **Outputs**:
    - `sta/reduce/Fld2/*_stars_stats_mdp.fits` : Save and updated list of stars with all their moffat fits.
    - `sta/reduce/Fld2/*_psf_mof_oversamp#.fits:` Save the average PSF (flux-weighted). oversampling set at 2

### Part 3: Make a stack of images for each mode and analyze stacks

- `append_massdimm()`
  - **Requires**: stats files in /stats/ folder
  - Pulls MASS DIMM data based on time
  - **Outputs**:
    - reduce/stats/stats*_mdp.fits
- `stack()`
  - **Requires**: stars.txt, cleaned images
  - Loop through all the starlists to get transformations, then shift images to those starlists
  - **Outputs**:
    - `reduce/stacks/fld2_stack_<key>.fits` : an image formed through stacking medians
- `analyze_stacks()`
  - **Requires**: stacked images, mask
  - Runs starfinding, star stats, and moffat fitting on stacked images  
  - **Outputs**:
    - `reduce/stacks/fld2_stack_<key>_<file>.fits` : general psf_mod, psf_obs, psf_mof_oversample, stars.txt, star_stats, star_stats_mdp

 
## Summary of Four Filter Reduction and Analysis

**General goal:** Know filter rotation for each file, split individual starlists by filter, combine STATS on starlists

- `split_filters()`
  - **Requires**: stars.txt
  - split star lists by color, saving each in own txt, noting orientation
  - iters by rotation key
  - **Outputs**
    - `reduce/Fld/*<filt>_<order>_stars.txt`

- `calc_fourfilt_stats()`
  - **Requires**: stats, clean images, starlists
  - for each suffix, for each color, collect color’s star stats
  - **Outputs**
    - `reduce/stats/stats_<key>_<color>.fits`

