from astropy.table import Table
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.table import Column
import matplotlib.patches as patches
from imaka.analysis import moffat
from astropy.modeling import models, fitting


def plot_stability(files_o, files_c):

    labels = ['B', 'V', 'R', 'I']

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
    bplot1 = axes.boxplot(list_o, sym='', vert=True, patch_artist=True)   # fill with color
    for patch in bplot1['boxes']:
         patch.set_alpha(0.3);
    for element in ['boxes', 'whiskers', 'fliers', 'means', 'medians', 'caps']:
        plt.setp(bplot1[element], color='blue')

    bplot2 = axes.boxplot(list_c, sym='', vert=True, patch_artist=True)
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

    axes.text(3.9,0.93,'AO-Off', color='k', fontsize=16)
    axes.text(3.2,0.93,'AO-On', color='k', fontsize=16)
    axes.add_patch(patches.Rectangle((3.85, 0.92), 0.6, 0.05, color='blue', alpha=0.3))
    axes.add_patch(patches.Rectangle((3.13, 0.92), 0.6, 0.05, color='red', alpha=0.3));

    plt.savefig('Four_Filter_Stability.png');
    plt.show();

    return


def focus_comp(files_o, files_c):

    plt.figure(figsize=(15,6))
    colors = ['b', 'g', 'r', 'purple']
    waves = [445, 551, 658, 806]
    ps = 0.12

    ax0 = plt.subplot2grid((2,7), (0,0), colspan=3, rowspan=2)
    ax1 = plt.subplot2grid((2,7), (0,3), colspan=2)
    ax2 = plt.subplot2grid((2,7), (0,5), colspan=2)
    ax3 = plt.subplot2grid((2,7), (1,3), colspan=2)
    ax4 = plt.subplot2grid((2,7), (1,5), colspan=2)

    for i in range(4):

        FWHM_min_o, sig_FWHM_min_o, FWHM_maj_o, sig_FWHM_maj_o = moffat.calc_mof_fwhm(files_o[i], filt=False)
        FWHM_min_c, sig_FWHM_min_c, FWHM_maj_c, sig_FWHM_maj_c = moffat.calc_mof_fwhm(files_c[i], filt=False)    

        dat = Table.read(files_o[i])
        ind2 = np.where(dat['Focus']==str(-2))
        ind3 = np.where(dat['Focus']==str(-3.5))
        FWHM_2 = np.mean(FWHM_min_o[ind2])*0.12
        FWHM_3 = np.mean(FWHM_min_o[ind3])*0.12
        FWHM_0 = np.mean(FWHM_min_o)*0.12
        err_2 = (np.std(FWHM_min_o[ind2])/np.sqrt(len(FWHM_min_o[ind2])))*0.12
        err_3 = (np.std(FWHM_min_o[ind3])/np.sqrt(len(FWHM_min_o[ind3])))*0.12
        err_0 = (np.std(FWHM_min_o)/np.sqrt(len(FWHM_min_o)))*0.12
        ax0.errorbar(waves[i]+10, FWHM_2, yerr=err_2, fmt='x', color=colors[i])
        ax0.errorbar(waves[i]-10, FWHM_3, yerr=err_3, fmt='+', color=colors[i])
        ax0.errorbar(waves[i], FWHM_0, yerr=err_0, fmt='*', color=colors[i])


        dat = Table.read(files_c[i])
        ind2 = np.where(dat['Focus']==str(-2))
        ind3 = np.where(dat['Focus']==str(-3.5))
        FWHM_2 = np.mean(FWHM_min_c[ind2])*0.12
        FWHM_3 = np.mean(FWHM_min_c[ind3])*0.12
        FWHM_0 = np.mean(FWHM_min_c)*0.12
        err_2 = (np.std(FWHM_min_c[ind2])/np.sqrt(len(FWHM_min_c[ind2])))*0.12
        err_3 = (np.std(FWHM_min_c[ind3])/np.sqrt(len(FWHM_min_c[ind3])))*0.12
        err_0 = (np.std(FWHM_min_c)/np.sqrt(len(FWHM_min_c)))*0.12
        ax0.errorbar(waves[i]+10, FWHM_2, yerr=err_2, fmt='x', color=colors[i])
        ax0.errorbar(waves[i]-10, FWHM_3, yerr=err_3, fmt='+', color=colors[i])
        ax0.errorbar(waves[i], FWHM_0, yerr=err_0, fmt='*', color=colors[i])

    ax0.plot(0,0, 'kx', label='Focus=-2.0')
    ax0.plot(0,0, 'k+', label='Focus=-3.5')
    ax0.plot(0,0, 'k*', label='Combined')
    ax0.axhline(0.578, color='k', linestyle='--')
    ax0.legend(loc=3, fontsize=14)
    ax0.axis([410, 840, 0.45, 0.7])
    ax0.text(600, 0.583, 'AO-Off', fontsize=14)
    ax0.text(600, 0.565, 'AO-On', fontsize=14)

    ax0.xaxis.set_tick_params(labelsize=14)
    ax0.yaxis.set_tick_params(labelsize=14)
    ax0.set_xlabel('Observation Wavelength (nm)', fontsize=16)
    ax0.set_ylabel('Minor FWHM (arcsec)', fontsize=16)

    #####

    bins = np.arange(3*ps, 10*ps, 0.001)

    for i in range(4):
        dat = Table.read(files_o[i])
        ind2 = np.where(dat['Focus']==str(-2))
        FWHM_2 = dat['emp_fwhm'][ind2]*ps
        ax1.hist(FWHM_2, bins=bins, color=colors[i], histtype='step', normed=True, cumulative=True)
        ax1.axvline(np.median(FWHM_2), color=colors[i], alpha=0.5)
    ax1.set_ylabel('Focus=-2.0', fontsize=16)
    ax1.set_title('AO-Off', fontsize=18)
    ax1.set_xlim(4*ps, 9*ps)

    for i in range(4):
        dat = Table.read(files_o[i])
        ind2 = np.where(dat['Focus']==str(-3.5))
        FWHM_2 = dat['emp_fwhm'][ind2]*ps
        ax3.hist(FWHM_2, bins=bins, color=colors[i], histtype='step', normed=True, cumulative=True)
        ax3.axvline(np.median(FWHM_2), color=colors[i], alpha=0.5)
    ax3.set_ylabel('Focus=-3.5', fontsize=16)
    ax3.set_xlabel('Minor FWHM (arcsec)', fontsize=16)
    ax3.set_xlim(4*ps, 9*ps)

    for i in range(4):
        dat = Table.read(files_c[i])
        ind2 = np.where(dat['Focus']==str(-2))
        FWHM_2 = dat['emp_fwhm'][ind2]*ps
        ax2.hist(FWHM_2, bins=bins, color=colors[i], histtype='step', normed=True, cumulative=True)
        ax2.axvline(np.median(FWHM_2), color=colors[i], alpha=0.5)
    ax2.set_title('AO-On', fontsize=18)
    ax2.set_xlim(3*ps, 6.5*ps)
    ax2.set_ylabel('Focus=-3.5', fontsize=16)


    for i in range(4):
        dat = Table.read(files_c[i])
        ind2 = np.where(dat['Focus']==str(-3.5))
        FWHM_2 = dat['emp_fwhm'][ind2]*ps
        ax4.hist(FWHM_2, bins=bins, color=colors[i], histtype='step', normed=True, cumulative=True)
        ax4.axvline(np.median(FWHM_2), color=colors[i], alpha=0.5)
    ax4.set_xlabel('Minor FWHM (arcsec)', fontsize=16)
    ax4.set_xlim(3*ps, 6.5*ps)
    ax4.set_ylabel('Focus=-2.0', fontsize=16)

    
    plt.subplots_adjust(hspace=0, wspace=1);
    plt.savefig('filt_comp.png');

    return


def minmaj_FWHM(files_o, files_c):
    
    labels = ['B', 'V', 'R', 'I']
    ps=0.12
    plt.figure(figsize=(9,8))
    for i in range(4):

        # Make divisions by focus
        dat = Table.read(files_o[i])
        ind2 = np.where(dat['Focus']==str(-2))
        ind3 = np.where(dat['Focus']==str(-3.5))
        if i == 0:
            ind = ind3
        else:
            ind = ind2

        # Read in Moffat fit data
        FWHM_min_o, sig_FWHM_min_o, FWHM_maj_o, sig_FWHM_maj_o = moffat.calc_mof_fwhm(files_o[i], filt=False)
        FWHM_min_c, sig_FWHM_min_c, FWHM_maj_c, sig_FWHM_maj_c = moffat.calc_mof_fwhm(files_c[i], filt=False)    

        # Plot FWHM distributions
        plt.subplot(4,2,(i*2)+1)
        bins = np.arange(2*ps,10*ps,0.001*ps)

        plt.hist(FWHM_min_o[ind]*ps, bins=bins, histtype='step', cumulative=True, normed=True, color='b', label='AO Off, min')
        plt.hist(FWHM_maj_o[ind]*ps, bins=bins, histtype='step', cumulative=True, normed=True, color='b', linestyle='--', label='AO Off, Maj')
        plt.hist(FWHM_min_c[ind]*ps, bins=bins, histtype='step', cumulative=True, normed=True, color='r', label='AO On, min')
        plt.hist(FWHM_maj_c[ind]*ps, bins=bins, histtype='step', cumulative=True, normed=True, color='r', linestyle='--', label='AO On, maj')
        plt.axvline(np.median(FWHM_min_o[ind])*ps, color='b', alpha=0.5)
        plt.axvline(np.median(FWHM_maj_o[ind])*ps, color='b', alpha=0.5, linestyle='--')
        plt.axvline(np.median(FWHM_min_c[ind])*ps, color='r', alpha=0.5)
        plt.axvline(np.median(FWHM_maj_c[ind])*ps, color='r', alpha=0.5, linestyle='--')

        plt.yticks(fontsize=14)
        plt.ylabel(labels[i], fontsize=16)
        plt.xlim(2.9*ps, 9.5*ps)
        plt.yticks([0, 0.2, 0.4, 0.6, 0.8]) 
        if i == 0:
            plt.title('Minor and Major FWHM', fontsize=18)
        if i == 3:
            plt.xlabel('FWHM (arcsec)', fontsize=16)
            plt.xticks(fontsize=14)
            plt.legend(fontsize=9, loc=4)
        else:
            plt.xticks([], [])

        plt.subplot(4,2,(i*2)+2)
        bins = np.arange(1,2,0.05)
        elon_o = FWHM_maj_o / FWHM_min_o
        elon_c = FWHM_maj_c / FWHM_min_c
        plt.hist(elon_o[ind], bins=bins, color='b', alpha=0.5, normed=True, label='AO-Off')
        plt.hist(elon_c[ind], bins=bins, color='r', alpha=0.5, normed=True, label='AO-On')

        plt.xlim(1, 1.8)
        plt.yticks(fontsize=14)
        plt.axvline(np.median(elon_o), color='Navy')
        plt.axvline(np.median(elon_c), color='Firebrick')
        if i == 0:
            plt.title('Elongation', fontsize=18)
            plt.legend(fontsize=14)
        if i == 3:
            plt.xlabel('Elongation Factor', fontsize=16)
            plt.xticks(fontsize=14)
        else:
            plt.xticks([], [])


    plt.subplots_adjust(hspace=0);
    plt.savefig('4filt_fwhm.png');
    return


def plot_gain(files_o, files_c):

    plt.figure(1, figsize=(10,4))
    waves = [445, 551, 658, 806]
    ps=0.12

    for i in range(4):
        plt.subplot(121)
        FWHM_min_o, sig_FWHM_min_o, FWHM_maj_o, sig_FWHM_maj_o = moffat.calc_mof_fwhm(files_o[i], filt=False)
        FWHM_min_c, sig_FWHM_min_c, FWHM_maj_c, sig_FWHM_maj_c = moffat.calc_mof_fwhm(files_c[i], filt=False)    
        if i==0:
            plt.errorbar(waves[i], np.median(FWHM_min_o)*ps, yerr=np.std(FWHM_min_o)*ps/np.sqrt(len(FWHM_min_o)), fmt='bo', label='AO Off, Min')
            plt.errorbar(waves[i], np.median(FWHM_min_c)*ps, yerr=np.std(FWHM_min_c)*ps/np.sqrt(len(FWHM_min_c)), fmt='ro', label='AO On, Min')
            plt.errorbar(waves[i], np.median(FWHM_maj_o)*ps, yerr=np.std(FWHM_maj_o)*ps/np.sqrt(len(FWHM_maj_o)), fmt='bs', label='AO Off, Maj')
            plt.errorbar(waves[i], np.median(FWHM_maj_c)*ps, yerr=np.std(FWHM_maj_c)*ps/np.sqrt(len(FWHM_maj_c)), fmt='rs', label='AO On, Maj')

        else:
            plt.errorbar(waves[i], np.median(FWHM_min_o)*ps, yerr=np.std(FWHM_min_o)*ps/np.sqrt(len(FWHM_min_o)), fmt='bo')
            plt.errorbar(waves[i], np.median(FWHM_min_c)*ps, yerr=np.std(FWHM_min_c)*ps/np.sqrt(len(FWHM_min_c)), fmt='ro')
            plt.errorbar(waves[i], np.median(FWHM_maj_o)*ps, yerr=np.std(FWHM_maj_o)*ps/np.sqrt(len(FWHM_maj_o)), fmt='bs')
            plt.errorbar(waves[i], np.median(FWHM_maj_c)*ps, yerr=np.std(FWHM_maj_c)*ps/np.sqrt(len(FWHM_maj_c)), fmt='rs')

        plt.subplot(122)
        gain_min = np.median(FWHM_min_o) / np.median(FWHM_min_c)
        gain_maj = np.median(FWHM_maj_o) / np.median(FWHM_maj_c)
        if i == 0:
            plt.plot(waves[i], gain_min, 'ko', label='Minor')
            plt.plot(waves[i], gain_maj, 'ks', label='Major')
        else:
            plt.plot(waves[i], gain_min, 'ko')
            plt.plot(waves[i], gain_maj, 'ks')


    plt.subplot(121)
    plt.xlabel('Observation Wavelength (nm)', fontsize=16)
    plt.ylabel('Minor FWHM (arcsec)', fontsize=16)
    plt.legend(fontsize=14, loc='best')
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.legend(fontsize=14)

    plt.subplot(122)
    plt.xlabel('Observation Wavelength (nm)', fontsize=16)
    plt.ylabel('Correction Factor', fontsize=16)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.legend(fontsize=14)

    plt.savefig('colors_arcsec.png')
    plt.tight_layout();
    
    return


def power_model(files_o, files_c):

    waves = [445, 551, 658, 806]
    ps = 0.12

    fwhm_o = []
    fwhm_c = []
    fwhm_o_err = []
    fwhm_c_err = []
    for i in range(4):
        FWHM_min_o, sig_FWHM_min_o, FWHM_maj_o, sig_FWHM_maj_o = moffat.calc_mof_fwhm(files_o[i], filt=False)
        FWHM_min_c, sig_FWHM_min_c, FWHM_maj_c, sig_FWHM_maj_c = moffat.calc_mof_fwhm(files_c[i], filt=False)   
        fwhm_o.append(np.median(FWHM_min_o)*ps)
        fwhm_o_err.append(np.std(FWHM_min_o)*ps/np.sqrt(len(FWHM_min_o)))
        fwhm_c.append(np.median(FWHM_min_c)*ps)
        fwhm_c_err.append(np.std(FWHM_min_c)*ps/np.sqrt(len(FWHM_min_c)))

    plt.errorbar(waves, fwhm_o, yerr=fwhm_o_err, fmt='bo', label='AO-Off Data')
    plt.errorbar(waves, fwhm_c, yerr=fwhm_c_err, fmt='ro', label='AO-On Data')

    init = models.PowerLaw1D(amplitude=1., x_0=1., alpha=1.)
    fit = fitting.LevMarLSQFitter()
    p_o = fit(init, waves, fwhm_o, weights=1.0/np.array(fwhm_o_err))
    p_c = fit(init, waves, fwhm_c, weights=1.0/np.array(fwhm_c_err))

    x = np.linspace(445, 806, 100)
    plt.plot(x, p_o(x), 'b-', label='AO-Off Model')
    plt.plot(x, p_c(x), 'r-', label='AO-On Model');

    χ2_o = np.sum(((p_o(waves)-fwhm_o)/fwhm_o_err)**2) 
    χ2_c = np.sum(((p_c(waves)-fwhm_c)/fwhm_c_err)**2)

    α_o = p_o.alpha.value
    α_c = p_c.alpha.value

    plt.text(450, 0.61, 'χ$^2$='+str(np.round(χ2_o,2)), color='b', fontsize=14)
    plt.text(450, 0.5, 'χ$^2$='+str(np.round(χ2_c,2)), color='r', fontsize=14)
    plt.text(450, 0.63, 'α='+str(np.round(α_o,2)), color='b', fontsize=14)
    plt.text(450, 0.52, 'α='+str(np.round(α_c,2)), color='r', fontsize=14)

    plt.xlabel('Observation Wavelength (nm)', fontsize=16)
    plt.ylabel('Minor FWHM (arcsec)', fontsize=16)
    plt.title('Wavelength Dependence Model', fontsize=18)
    plt.legend()
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14);
    
    return