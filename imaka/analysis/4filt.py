from astropy.table import Table
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.table import Column
import matplotlib.patches as patches
from imaka.analysis import moffat

def plot_stability(files_o, files_c):
    plt.figure(figsize=(8,10))
    for i in range(4):

        # Make divisions by focus
        dat = Table.read(files_o[i])
        ind2 = np.where(dat['Focus']==str(-2))
        ind3 = np.where(dat['Focus']==str(-3.5))
        if i == 0:
            ind = ind3
        else:
            ind = ind2
            
    list_o = []
    list_c = []

    for i in range(4):

        # Make divisions by focus
        dat = Table.read(files_o[i])
        ind2 = np.where(dat['Focus']==str(-2))
        ind3 = np.where(dat['Focus']==str(-3.5))
        if i == 0:
            ind = ind3
        else:
            ind = ind2

        FWHM_min_o, sig_FWHM_min, FWHM_maj, sig_FWHM_maj = moffat.calc_mof_fwhm(files_o[i], filt=False)
        FWHM_min_c, sig_FWHM_min, FWHM_maj, sig_FWHM_maj = moffat.calc_mof_fwhm(files_c[i], filt=False)
        list_o.append(FWHM_min_o[ind]*0.12)
        list_c.append(FWHM_min_c[ind]*0.12)

    fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(7, 5))

    # rectangular box plot
    bplot1 = axes.boxplot(list_o, sym='', vert=True, patch_artist=True, labels=['Open Loop'] * 4)   # fill with color
    for patch in bplot1['boxes']:
         patch.set_alpha(0.3);
    for element in ['boxes', 'whiskers', 'fliers', 'means', 'medians', 'caps']:
        plt.setp(bplot1[element], color='blue')

    bplot2 = axes.boxplot(list_c, sym='', vert=True, patch_artist=True, labels=['Closed Loop'] * 4)
    for patch in bplot2['boxes']:
        patch.set_alpha(0.3)
    for element in ['boxes', 'whiskers', 'fliers', 'means', 'medians', 'caps']:
        plt.setp(bplot2[element], color='red')

    # Formatting Axes and Labels
    axes.yaxis.grid(True)
    axes.set_xticks([y+1 for y in range(4)], )
    axes.set_xlabel('Observation Band', fontsize=18)
    axes.set_ylabel('Minor FWHM (arcsec)', fontsize=18)

    plt.setp(axes, xticks=[y+1 for y in range(4)],
             xticklabels=labels)
    axes.tick_params(axis='both', which='major', labelsize=16)
    plt.tight_layout()

    bplot1['boxes'][0].set_label('Open Loop')
    bplot2['boxes'][0].set_label('Closed Loop')
    axes.legend(loc=2)

    axes.text(3.9,0.93,'AO-Off', color='k', fontsize=16)
    axes.text(3.2,0.93,'AO-On', color='k', fontsize=16)
    axes.add_patch(patches.Rectangle((3.85, 0.92), 0.6, 0.05, color='blue', alpha=0.3))
    axes.add_patch(patches.Rectangle((3.13, 0.92), 0.6, 0.05, color='red', alpha=0.3))

    plt.savefig('Four_Filter_Stability.png')
    plt.show()

    return
