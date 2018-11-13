#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 20 13:23:34 2018

@author: X-phile

EXPL_PLOTS() creates explanatory plots of parameters pulled from Spitzer, MEarth, and LCO transit data.
Parameters include:
    ser - transit numbers
    rp - radius of planet in terms of stellar radius
    a - semi-major axis
    t0 - time of mid-transit
    secosw - sqrt(eccentricity)*cos(omega)
    sesinw - sqrt(eccentricity)*sin(omega)
    q1 - limb darkening parameter a
    q2 - limb darkening parameter b
    slope - slope of light curve
    offset - offset of light curve
    depth - depth of transit
    duration - duration of transit
    ecc - eccentricity of orbit
    w - argument of periapsis
    max - list of best values for each parameter, can be used to make violin plots

"""
import numpy as np
import glob
import os
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
from brokenaxes import brokenaxes
import pickle
import seaborn as sb

global dataset


def make_plot(dataset, ser, x_label, y_label, x_data, y_data, y_err, spitzer_avg, mearth_avg, LCO_avg, x_err=None, ptitle=None, ftitle=None):
    colors = cm.rainbow(np.linspace(0, 1, len(x_data)))
    colors_1 = ['#a50026', '#d73027', '#f46d43', '#fdae61', '#fee090', '#ffffbf', '#e0f3f8', '#abd9e9', '#74add1', '#4575b4', '#313695', '#000066']
    colors_2 = ['#40004b', '#762a83', '#9970ab', '#c2a5cf', '#e7d4e8', '#f7f7f7', '#d9f0d3', '#a6dba0', '#5aae61', '#1b7837', '#00441b', '#003300']
#    colors picked using ColorBrewer2.0, http://colorbrewer2.org/


    ## Parameter A v. transit number

    if ser is x_data:
        print 'path b1: no z values, param vs. time'

        ##Split x axis
        plt.figure()
        if 'S' in dataset:
            bax = brokenaxes.brokenaxes(xlims=((185., 200.), (235., 237.)), hspace=.05)
        elif 'M' in dataset:
            bax = brokenaxes.brokenaxes(xlims=((85., 90.), (170., 195.)), hspace=.05)
        elif 'L' in dataset:
            bax = brokenaxes.brokenaxes(xlims=((155., 160.), (175., 192.), (280., 290.)), hspace=.05)

        bax.scatter(x_data, y_data, color=colors, marker='o', label=ser, linewidth=.2, edgecolors='k')

        if x_err is not None:
            bax.errorbar(x_data, y_data, np.array(x_err).T, np.array(y_err).T,  capsize=3, fmt='.', c='k', alpha=.4, zorder=0)
        else:
            bax.errorbar(x_data, y_data, np.array(y_err).T, capsize=3, fmt='.', c='k', alpha=.4, zorder=0)

        bax.set_ylabel(y_label, fontsize=14, labelpad = 45)
        bax.set_xlabel('Transit', fontsize=14, labelpad = 30)
#        bax.autoscale()

        ## Build legend
        recs = []
        for i in range(0,len(colors)):
            recs.append(mpatches.Rectangle((0,0),1,1,fc=colors[i]))
        lgd = bax.legend(handles = recs, labels = np.ndarray.tolist(ser), loc=6, title='Transit Number', shadow=False, bbox_to_anchor=(1, 0.5))

        ## Draw refline from average values
        if 'S' in str(dataset):
            for i in spitzer_avg:
                if np.mean(y_data)==i:
                    bax.axhline(y=i, xmin=-10, xmax=25, linestyle='dashed', color='r', alpha=.6)
                    print 'drawing spitzer avg (red)'
        elif 'M' in str(dataset):
            for i in mearth_avg:
                if np.mean(y_data)==i:
                    bax.axhline(y=i, xmin=-10, xmax=25, linestyle='dashed', color='b', alpha=.6)
                    print 'drawing mearth avg (blue)'
        elif 'L' in dataset:
            for i in LCO_avg:
                if np.mean(y_data)==i:
                    bax.axhline(y=i, xmin=-10, xmax=25, linestyle='dashed', color='m', alpha=.6)
                    print 'drawing LCO avg (mauve)'

        if ptitle is None:
            plt.title(y_label + ' v. ' + x_label, fontsize=16)
        else:
            plt.title(ptitle, fontsize=16)

        if ftitle is None:
            print 'no title given'
            x_label = x_label.replace(' ', '-')
            y_label = y_label.replace(' ', '-')
            ftitle = y_label + '-v-' + x_label
        print ftitle
        if 'S' in str(dataset):
            plt.savefig('../plots/'+'S_' + ftitle + '_plot.pdf', bbox_extra_artists=(lgd,), bbox_inches='tight')
        elif 'M' in str(dataset):
            plt.savefig('../plots/'+'M_' + ftitle + '_plot.pdf', bbox_extra_artists=(lgd,), bbox_inches='tight')
        elif 'L' in dataset:
            plt.savefig('../plots/'+'L_' + ftitle + '_plot.pdf', bbox_extra_artists=(lgd,), bbox_inches='tight')

    #==========================================================#
    ## Parameter A v. Parameter B

    else:
        print 'path b2: no z values, param vs. param'
        plt.figure()
        plt.scatter(x_data, y_data, color=colors, linewidth=.3, edgecolors='k')
        if x_err is not None:
            plt.errorbar(x_data, y_data, np.array(x_err).T, np.array(y_err).T,  capsize=3, fmt='.', c='k', alpha=.4, zorder=0)
            print 'including x error'
        else:
            plt.errorbar(x_data, y_data, np.array(y_err).T,  fmt='.', c='k', alpha=.4, zorder=0)
        plt.xlabel(x_label, fontsize = 14)
        plt.ylabel(y_label, fontsize = 14)
        plt.ylim(np.mean(y_data)-.35*np.mean(y_data), np.mean(y_data)+.35*np.mean(y_data))

        ## Build legend
        recs = []
        for i in range(0,len(colors)):
            recs.append(mpatches.Rectangle((0,0),1,1,fc=colors[i]))
            plt.legend(recs, np.ndarray.tolist(ser), loc=6, title='Transit Number', shadow=False, bbox_to_anchor=(1, 0.5))

        if ptitle is None:
            plt.title(y_label + ' v. ' + x_label, fontsize = 16)
        else:
            plt.title(ptitle, fontsize = 16)

        ## Draw refline from average values
        if 'S' in str(dataset):
            for i in spitzer_avg:
                if np.mean(y_data)==i:
                    plt.axhline(y=i, xmin=-10, xmax=25, linestyle='dashed', color='b', alpha=.6)
                    print 'drawing spitzer avg (red)'
        elif 'M' in str(dataset):
            for i in mearth_avg:
                if np.mean(y_data)==i:
                    plt.axhline(y=i, xmin=-10, xmax=25, linestyle='dashed', color='r', alpha=.6)
                    print 'drawing mearth avg (blue)'
        elif 'L' in dataset:
            for i in LCO_avg:
                if np.mean(y_data)==i:
                    bax.axhline(y=i, xmin=-10, xmax=25, linestyle='dashed', color='m', alpha=.6)
                    print 'drawing LCO avg (mauve)'

        if ftitle is None:
            print 'no title given'
            x_label = x_label.replace(' ', '-')
            y_label = y_label.replace(' ', '-')
            ftitle = y_label + '-v-' + x_label
        print ftitle
        if 'S' in str(dataset):
            plt.savefig('../plots/'+'S_' + ftitle + '_plot.pdf', bbox_inches='tight')
        elif 'M' in str(dataset):
            plt.savefig('../plots/'+'M_' + ftitle + '_plot.pdf', bbox_inches='tight')
        elif 'L' in dataset:
            plt.savefig('../plots/'+'L_' + ftitle + '_plot.pdf', bbox_inches='tight')


def violin(varyparams, vplot, d=''):
    ## varyparams = all parameters that emcee is varying
    ## vplot = all params to be plotted here
    if 'K' in d:
        pickles = glob.glob('/Users/X-phile/Dropbox/k225_transits/code/pickle*K2_k2sc')
    else:
        pickles = glob.glob('/Users/X-phile/Dropbox/k225_transits/code/pickle*_')
    print pickles
    ser = []
    dset = [] ##doesn't work for S and M >> put more if statements (ugh)
    reduc = []
    pdata = [None]*len(pickles)
    i = 0

    ## Unpickle and organize
    for i, p in enumerate(pickles):
        pdata[i] = pickle.load(open(p))
        h, t = os.path.split(pickles[i])
        ser.append(t.split('_')[1])
        dset.append(t.split('_')[2])
        if 'K2' in t.split('_')[2]:
            reduc.append(t.split('_')[3])

    ## Assign serials for GJ1132, Megastack data
    for i in range(0,len(ser)):
        if 'GJ' in ser[i]:
            ser[i] = '1132.'
            dset[i] = 'GJ'
        if 'Megastack' in ser[i]:
            ser[i] = '-1.'

    ## Convert serials from string to int
    ser = [int(float(ser[i])) for i in range(0,len(pdata))]

    ## Sort everything by increasing transit number
    bigthing = zip(ser, pdata, dset, reduc)
    bigthing.sort(key=lambda x: x[0])
#    print bigthing

    ser, pdata, dset, reduc = zip(*bigthing)
    ser = list(ser)
    pdata = list(pdata)
    dset = list(dset)
    reduc = list(reduc)

#    print dset

    ## Make plot for each parameter
    for v in vplot:
        print v
        ind = varyparams.index(v)
        pplot = [pdata[i][:,ind] for i in range(0,len(ser))]

        colors = [None]*len(dset)
        for i in range(len(dset)):
            if dset[i]=='spitzer':
                colors[i] = '#579BC3'
            elif dset[i]=='mearth':
                colors[i] = '#DACF33'
            elif dset[i]=='LCO':
                colors[i] = '#8184F6'
            elif dset[i]=='K2':
                colors[i] = '#f9ae1f'
            else:
                colors[i] = '#70BF81'

        ## O-C plot: observed - calculated, data - t0*numtransits
        if 't0' in v:
            print 'correcting t0s'

            ## List of myt0 values zipped to serial
            if 'K2' in dset:
                ktwo_data = open('K2_ktwo_slope_file.txt', 'r')
                k2sc_data = open('K2_ktwo_slope_file.txt', 'r')

                data = [ktwo_data, k2sc_data]
                tel = ['spitz', 'mearth', 'ktwo', 'k2sc']

                ## Parse data files for each dataset
                nums, t0_list, tel_list = [], [], []
                for i in range(0, len(data)):
                    dset = data[i]
                    lines = dset.readlines()
                    for bit in lines:
                        nums.append(bit.split(';')[0])
                        t0_list.append(bit.split(';')[3])
                        tel_list.append(tel[i])
                    dset.close() ##

                nums = list(nums)
                t0s = list(t0s)
                reduc = list(reduc)


            else:
                s_data = open('spitzer_slope_file.txt', 'r')
                m_data = open('mearth_slope_file.txt', 'r')

                data = [s_data, m_data, ktwo_data, k2sc_data]
                tel = ['spitz', 'mearth', 'ktwo', 'k2sc']


                ## Parse data files for each dataset
                ser, t0_list, tel_list = [], [], []
                for i in range(0, len(data)):
                    dset = data[i]
                    lines = dset.readlines()
                    for bit in lines:
                        ser.append(bit.split(';')[0])
                        t0_list.append(bit.split(';')[3])
                        tel_list.append(tel[i])
                    dset.close() ##

                ## Remove Megastack things from everywhere -- hardcode option, change later?
                if 'S' in dset or 'M' in dset:
                    ser = list(ser[2:])
                    pdata = list(pdata[2:])
                    pplot = list(pplot[2:])
                    dset = list(dset[2:])
                    colors = list(colors[2:])

                nums, t0s = zip(*timeslist)
                nums = list(nums)
                t0s = list(t0s)

            ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
            ## Center time values by t0
            if ser==nums:
                pplot = list(np.array(pplot[i]) - np.array(t0s[i]) for i in range(0,len(pplot)))
            ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##

        ## Change colors based on data reduction
        for i in range(0, len(colors)):
            if 'k2sc' in reduc[i]:
                colors[i] = '#70BF81'

        ## Make plots
        if len(ser)>40:
            fig, (ax1, ax2, ax3) = plt.subplots(3,1, figsize=(15,15), sharey=True)
            sb.violinplot(data=pplot[:14], inner = 'quartiles', cut = 0, palette = colors, ax = ax1)
            sb.violinplot(data=pplot[14:29], inner = 'quartiles', cut = 0, palette = colors, ax = ax2)
            sb.violinplot(data=pplot[29:], inner = 'quartiles', cut = 0, palette = colors, ax = ax3)
            ax1.set(xticklabels=ser[:14])
            ax2.set(xticklabels=ser[14:29])
            ax3.set(xticklabels=ser[29:])
            ax3.set_xlabel('Transit Number')
            ax2.set_ylabel(str(v))
            ax1.set_ylim(-.01, .01)
            fig.suptitle(str(v) + ' Distribution')
        elif len(ser)>11:  ##2 vertically stacked panels for large number of transits
            fig, (ax1, ax2) = plt.subplots(2,1, sharey=True, figsize=(15,10))
            sb.violinplot(data=pplot[:12], inner = 'quartiles', cut = 0, palette = colors, ax = ax1)
            sb.violinplot(data=pplot[12:], inner = 'quartiles', cut = 0, palette = colors, ax = ax2)
            ax1.set(xticklabels=ser[:12])
            ax2.set(xticklabels=ser[12:])
            ax2.set_xlabel('Transit Number')
            ax1.set_ylabel(str(v))
            fig.suptitle(str(v) + ' Distribution')
        else:
            fig, ax1 = plt.subplots(figsize=(15,5))
            sb.violinplot(data=pplot, inner = 'quartiles', cut = 0, palette = colors, ax = ax1)
            ax1.set(xticklabels=ser)
            ax1.set_xlabel('Transit Number')
            ax1.set_ylabel(str(v))
            fig.suptitle(str(v) + ' Distribution')
            plt.tight_layout()

        #Build legend
        if 'K2' in dset:
            patch1 = mpatches.Rectangle((0,0),1,1,fc='#70BF81')
            patch2 = mpatches.Rectangle((0,0),1,1,fc='#f9ae1f')
            ax1.legend(loc=1, handles=(patch1, patch2), labels=('K2SC', 'KTWO'))
        else:
            patch1 = mpatches.Rectangle((0,0),1,1,fc='#579BC3')
            patch2 = mpatches.Rectangle((0,0),1,1,fc='#DACF33')
            ax1.legend(loc=1, handles=(patch1, patch2), labels=('Spitzer', 'MEarth'))

        # For specific vplots, set limits
        if 'rp' in v:
            ax1.set_ylim(0.09, 0.12)
#        if 't0' in v:
#            ax1.set_ylim(-0.004, 0.005)

        plt.show()
        fig.savefig('../plots/' + str(v) + '_violin.png')


#8184f6 lilac
#a7efa7 light green
#70bf81 soft green
#579bc3 soft blue
#dacf33 gold
#bea823 darker gold