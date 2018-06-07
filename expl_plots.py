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


def violin(varyparams, vplot):
    pickles = glob.glob('/Users/X-phile/Dropbox/k225_transits/code/pickle*')
    ser = []
    dset = []
    pdata = [None]*len(pickles)
    i = 0

    ## Unpickle and organize
    for i, p in enumerate(pickles):
        pdata[i] = pickle.load(open(p))
        h, t = os.path.split(pickles[i])
        ser.append(t.split('_')[1])
        dset.append(t.split('_')[2])

    ## Fiddle with GJ, Megastack things
    for i in range(0,len(ser)):
        if 'GJ' in ser[i]:
            ser[i] = '1132.'
            dset[i] = 'GJ'
        if 'Megastack' in ser[i]:
            ser[i] = '0000.'

    ## Convert serials from string to int
    ser = [int(float(ser[i])) for i in range(0,len(pdata))]

    ## Sort everything by increasing transit number
    bigthing = zip(ser, pdata, dset)
    bigthing.sort(key=lambda x: x[0])
    ser, pdata, dset = zip(*bigthing)
    ser = list(ser)
    pdata = list(pdata)
    dset = list(dset)

    print dset

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
            else:
                colors[i] = '#70BF81'

        ## O-C plot: observed - calculated, data - t0*numtransits
        if 't0' in v:
            print 'correcting t0s'

            ## List of myt0 values zipped to serial
            a = [2457721.16224, 2457728.13136, 2457731.61593, 2457735.10049, 2457742.06962, 2457745.55418, \
                2457749.03875, 2457752.52331, 2457756.00788, 2457884.93674, 2457668.89378,  \
                2457675.86290, 2457682.83203, 2457689.80116, 2457696.77029, 2457703.73942,  \
                2457710.70854, 2457717.67767, 2457724.6468, 2457731.61593, 2457365.73671, 2457372.70584]
            b = [189, 191, 192, 193, 195, 196, 197, 198, 199, 236, 174, 176, 178, 180, 182, 184, 186, 188, 190, 192, 87, 89]
            timeslist = zip(b,a)
            timeslist.sort()

            ## Remove Megastack things from everywhere
            ser = list(ser[2:])
            pdata = list(pdata[2:])
            pplot = list(pplot[2:])
            dset = list(dset[2:])
            colors = list(colors[2:])

            ## Convert t0 to O-C
            nums, t0s = zip(*timeslist)
            nums = list(nums)
            t0s = list(t0s)

            if ser==nums:
                print 'doing the thing'
                pplot = list(np.array(pplot[i]) - np.array(t0s[i]) for i in range(0,len(pplot)))


        ## Make plots
        fig, ax = plt.subplots(figsize=(15,5))
        sb.violinplot(data=pplot, inner = 'quartiles', cut = 0, palette = colors, ax = ax)
        ax.set(xticklabels=ser)
        ax.set_xlabel('Transit Number')
        ax.set_ylabel(str(v))
        ax.set_title(str(v) + ' Distribution')

        if 'rp' in v:
            ax.set_ylim(0.08, 0.14)
        if 't0' in v:
            ax.set_ylim(-0.005, 0.005)

        plt.show()
        ax.figure.savefig('../plots/' + str(v) + '_violin.png')


#8184f6 lilac
#a7efa7 light green
#70bf81 soft green
#579bc3 soft blue
#dacf33 gold
#bea823 darker gold