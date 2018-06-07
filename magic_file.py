#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 15 09:59:43 2017

@author: X-phile
"""
from __future__ import division, unicode_literals
import numpy as np
import glob
import functionsCollection as fc
#import magic_file as m
import expl_plots as ex
from astropy.io import ascii
import gc
import matplotlib.pyplot as plt

import plotwalkers as pw
import os
computer = os.getenv('COMPUTERNAME')


'''
<Description: add later>

'''

def call_things(x, y, yerr, serial, conv_flux=True):
    time, magnitude, magError = fc.prelim(x, y, yerr, conv_flux=conv_flux)
    print 'finished prelim'
    rp_rstar, a_rstar, myt0, params, time, magnitude, magError, thismodel, transit_num = fc.initialModel(time, magnitude, magError, serial)
    print 'finished initial model'
    samples, sampler, variables = fc.emceeModel(rp_rstar, a_rstar, myt0, params, time, magnitude, magError, thismodel, serial, transit_num)
    print 'finished emcee model'
    fc.cornerPlot(samples, serial, variables)
    print 'finshed corner plot'
    pw.walkers(sampler, serial, variables)
    print 'finished plotting walkers'



#==============================================================================
# FUNCTION 1: FOR INDIVIDUAL FILES
#==============================================================================

def ind_analysis(list_of_files, dataset):
    ask_ind = input('Would you like to analyze individual datasets? Type ' + '"' + 'Y' + '"' + ' for yes, ' + '"' + 'N' + '"' + ' for no: ')
    if 'N' in ask_ind:
        print 'Ok, no idividual analysis.'
        return 0

    for file in list_of_files:
        print file

        ## Grab data from each file
        if 'M' in dataset:
            x, y, yerr = np.genfromtxt(file, delimiter=',', dtype=None, unpack=True, usecols=(0,1,2))

            t0 = 2457062.57964
            per = 3.484564
            numtransits = round((x[0]-t0)/per)
            myt0 = t0 + per*numtransits
            print myt0
            serial = str(numtransits) + '_mearth'   #example: 89_mearth
            print "serial:", serial
            conv_flux = True

        elif 'S' in dataset:
            data = ascii.read(file, delimiter = '\s', header_start = 0, data_start = 1, guess = False)
            x = np.array(data['BMJD_OBS'])+2400000. + .5
            y = data['Normalized_Flux']

            t0 = 2457062.57964
            per = 3.484564
            numtransits = round((x[0]-t0)/per)
            myt0 = t0 + per*numtransits
            print myt0

            serial = str(numtransits) + '_spitzer'   #example: 189_spitzer
            print serial

            ## Trim data to window around transit
            use = np.abs(x-myt0)>.0208
            xsample = x[use]
            ysample = y[use]

            use_cut = (ysample<1.02) & (ysample>.98)
            ysample = ysample[use_cut]
            xsample = xsample[use_cut]

            ## Set error to std
            std = np.std(ysample)
            yerr = np.full(len(x), std)
            conv_flux = False
            print 'no need to convert to flux'

        elif 'L' in dataset:
            data = ascii.read(file, delimiter = '\s', header_start = 0, data_start = 1, guess = False)
            x = np.array(data['Time'])
            y = np.array(data['Fl'])
            yerr = np.array(data['Fl_err'])

            t0 = 2457062.57964
            per = 3.484564
            numtransits = round((x[0]-t0)/per)
            myt0 = t0 + per*numtransits
            print myt0
            serial = str(numtransits) + '_LCO'   #example: 89_mearth
            print "serial:", serial
            conv_flux = False ##need to normalize?

        elif 'GJ' in dataset:
            data = ascii.read(file, delimiter = '\s', guess = False, format = 'cds')
            x = np.array(data['BJD']) + 2400000. + 0.5
            y = np.array(data['NFlux'])
            yerr = np.array(data['e_NFlux'])
            inst = np.array(data['Inst'])

            t0 = 2457184.55804
            per = 1.6289246

            spitz = inst=='Spitzer'
            mrth = inst=='MEarth'

            ## MEarth transits
            mx = x[mrth]
            my = y[mrth] - 1
            myerr = yerr[mrth]

            mx = mx%per #phase folding

            serial = 'GJ1132_mearth'
            conv_flux = True

            dataset = 'M'
            fc.dataset = dataset
            ex.dataset = dataset

            call_things(mx, my, myerr, serial, conv_flux=conv_flux)

            ## Spitzer transits
            sx = x[spitz]
            sy = y[spitz]
            syerr = yerr[spitz]

            sx = sx%per

            serial = 'GJ1132_spitzer'
            conv_flux = False

            dataset = 'S'
            fc.dataset = dataset
            ex.dataset = dataset

            call_things(sx, sy, syerr, serial, conv_flux=conv_flux)

            break

        ## Analyze data
        call_things(x, y, yerr, serial, conv_flux=conv_flux)

        gc.collect()
        print "garbage collected"
    print 'finished!'


#==============================================================================
# FUNCTION 2: FOR MEGASTACK
#==============================================================================

def megastack_analysis(list_of_files, dataset):
    ask_megastack = input('Would you like to make a megastack of your data? Type ' + '"' + 'Y' + '"' + ' for yes, ' + '"' + 'N' + '"' + ' for no: ')
    if 'N' in ask_megastack:
        print 'Ok, no megastack.'
        return 0
    if 'Y' in ask_megastack:
        if 'L' in dataset:
            print 'Don\'t do that'
            return 0
        else:
            print 'Megashark!'

    ## Correct slope and offset of individual datasets; stack corrected data
    x, y, yerr = [], [], []
    conv_flux=False
    if 'M' in dataset:
        f=open('mearth_slope_file.txt',"r")
        lines=f.readlines()
#        print lines
        ser, slope_list, offset_list = [], [], []
        for bit in lines:
            ser.append(bit.split(';')[0])
            slope_list.append(bit.split(';')[8])
            offset_list.append(bit.split(';')[9])
        f.close()
#        print 'ser: ', ser
        print 'slope file parsed'

        ## Determine time of mid transit for first datafile
        t = np.genfromtxt(list_of_files[0], delimiter=',', dtype=None, unpack=True, usecols=(0))
        t0 = 2457062.57964
        per = 3.484564
        numtransit = round((t[0]-t0)/per)
        first_t0 = t0 + per*numtransit
        count = 0

        for file in list_of_files:
            t, m, mE = np.genfromtxt(file, delimiter=',', dtype=None, unpack=True, usecols=(0,1,2))

            t0 = 2457062.57964
            per = 3.484564
            numtransit = round((t[0]-t0)/per)
            myt0 = t0 + per*numtransit

            ## Grab best slope and offset values
            for num in ser:
                if int(numtransit)==int(float(num)):
                    i = ser.index(str(numtransit))
                    best_slope = float((slope_list[i][1:-2]).split(',')[0])
                    best_offset = float((offset_list[i][1:-2]).split(',')[0])

                    print 'best values: ', best_slope, best_offset

            ## Convert to flux so linearization works
            m = np.power(10, m/-2.5)
            mE = m*mE/-2.5*np.log(10)

            ## Linearize data
            linmodel = best_slope*(t-float(myt0)) + best_offset
            print 'linmodel: ', linmodel
            print 'm before linearization: ', m

            m = m/linmodel
            print 'm post linearization: ', m

            ## Center all datasets besides first at first_t0
            if count!=1:
                t = t-(myt0-first_t0)
#                print 'myt0: ', myt0
#                print 'first t0: ', first_t0

            x.extend(t)
            y.extend(m)
            yerr.extend(mE)

            count+=1

    elif 'S' in dataset:

        ## Determine time of mid transit for first datafile
        data = ascii.read(list_of_files[0], delimiter = '\s', header_start = 0, data_start = 1, guess = False)
        t = np.array(data['BMJD_OBS'])+2400000. + .5

        t0 = 2457062.57964
        per = 3.484564
        numtransit = round((t[0]-t0)/per)
        first_t0 = t0 + per*numtransit
        count = 0

        ## Iterate through all datafiles to grab and stack data
        for file in list_of_files:
            count = count+1
            data = ascii.read(file, delimiter = '\s', header_start = 0, data_start = 1, guess = False)
            t = np.array(data['BMJD_OBS'])+2400000. + .5
            m = data['Normalized_Flux']

            t0 = 2457062.57964
            per = 3.484564
            numtransit = round((t[0]-t0)/per)
            myt0 = t0 + per*numtransit

            ## Center all datasets besides first at first_t0
            if count!=1:
                t = t-(myt0-first_t0)

            ## Set yerr as stdev of data points in window around transit
            use = np.abs(t-myt0)>.0208
            xsample = t[use]
            ysample = m[use]
            use_cut = (ysample<1.02) & (ysample>.98)
            ysample = ysample[use_cut]
            xsample = xsample[use_cut]

            std = np.std(t)
            mE = np.full(len(t), std)

            #no need to linearize, less interference with Spitzer data
            x.extend(t)
            y.extend(m)
            yerr.extend(mE)

    print 'values gotten', len(x), len(y), len(yerr)
    print 'finished removing slope and offset from individual stack'

    x = np.array(x)
    y = np.array(y)
    yerr = np.array(yerr)

    print 'all done! \n'

    if 'M' in dataset:
        serial = 'Megastack_mearth'
    elif 'S' in dataset:
        serial = 'Megastack_spitzer'
    elif 'L' in dataset:
        serial = 'Megastack_LCO'

    ## Analyze data
    call_things(x, y, yerr, serial, conv_flux=conv_flux)

    gc.collect()
    print "garbage collected"
    print 'finished!!!'

#==============================================================================
# EXPLANATORY PLOTS
#==============================================================================
#   To create explanatory plots, you first need to individually analyze all data
#   files in your dataset.

def visual_analysis(dataset):

    ask_plot = input('Would you like to make explanatory plots? Type \'Y\' for yes, \'N\' for no: ')
    if 'N' in ask_plot:
        print 'OK then, no helpful plots today…'
        return 0

    ## Read in parameters from slope file
    if 'M' in dataset:
        f=open('mearth_slope_file.txt',"r")
    elif 'S' in dataset:
        f=open('spitzer_slope_file.txt',"r")
    lines=f.readlines()
    ser, rp_list, a_list, t0_list, ecw_list, esw_list, limb1_list, limb2_list, \
                    slope_list, offset_list, dur_list, depth_list, ecc_list, w_list, max_list \
                    = [], [], [], [], [], [], [], [], [], [], [], [], [], [], []
    for bit in lines:
        ser.append(bit.split(';')[0])
        rp_list.append(bit.split(';')[1])
        a_list.append(bit.split(';')[2])
        t0_list.append(bit.split(';')[3])
        ecw_list.append(bit.split(';')[4])
        esw_list.append(bit.split(';')[5])
        limb1_list.append(bit.split(';')[6])
        limb2_list.append(bit.split(';')[7])
        slope_list.append(bit.split(';')[8])
        offset_list.append(bit.split(';')[9])
        dur_list.append(bit.split(';')[10])
        depth_list.append(bit.split(';')[11])
        ecc_list.append(bit.split(';')[12])
        w_list.append(bit.split(';')[13])
        max_list.append(bit.split(';')[14])
    f.close()

    ## Change serial from '89_mearth' to 89.0 (float)
    i = 0
    for num in ser:
        num = float(num.split('_')[0])
        ser[i] = num
        i+=1

    rp_err, a_err, t0_err, ecw_err, esw_err, limb1_err, limb2_err, slope_err, \
                            offset_err, dur_err, dep_err, ecc_err, w_err = \
                            [], [], [], [], [], [], [], [], [], [], [], [], []

    ## Grab best values
    i = 0
    while i<len(ser):
        rp_list[i] = np.array(rp_list[i].strip('[]').split(','), dtype=float)
        a_list[i] = np.array(a_list[i].strip('[]').split(','), dtype='float')
        t0_list[i] = np.array(t0_list[i].strip('[]').split(','), dtype='float')
        ecw_list[i] = np.array(ecw_list[i].strip('[]').split(','), dtype='float')
        esw_list[i] = np.array(esw_list[i].strip('[]').split(','), dtype='float')
        limb1_list[i] = np.array(limb1_list[i].strip('[]').split(','), dtype='float')
        limb2_list[i] = np.array(limb2_list[i].strip('[]').split(','), dtype='float')
        slope_list[i] = np.array(slope_list[i].strip('[]').split(','), dtype='float')
        offset_list[i] = np.array(offset_list[i].strip('[]').split(','), dtype='float')
        dur_list[i] = np.array(dur_list[i].strip('[]').split(','), dtype='float')
        depth_list[i] = np.array(depth_list[i].strip('[]').split(','), dtype='float')
        ecc_list[i] = np.array(ecc_list[i].strip('[]').split(','), dtype='float')
        w_list[i] = np.array(w_list[i].strip('[]').split(','), dtype='float')
        max_list[i] = np.array(max_list[i].strip('[]').split(','), dtype='float')
        i+=1

    ## Grab error values
    i = 0
    while i<len(ser):
        rp_err.append([rp_list[i][2], rp_list[i][1]])
        rp_list[i] = rp_list[i][0]
        a_err.append([a_list[i][2], a_list[i][1]])
        a_list[i] = a_list[i][0]
        t0_err.append([t0_list[i][2], t0_list[i][1]])
        t0_list[i] = t0_list[i][0]
        ecw_err.append([ecw_list[i][2], ecw_list[i][1]])
        ecw_list[i] = ecw_list[i][0]
        esw_err.append([esw_list[i][2], esw_list[i][1]])
        esw_list[i] = esw_list[i][0]
        limb1_err.append([limb1_list[i][2], limb1_list[i][1]])
        limb1_list[i] = limb1_list[i][0]
        limb2_err.append([limb2_list[i][2], limb2_list[i][1]])
        limb2_list[i] = limb2_list[i][0]
        slope_err.append([slope_list[i][2], slope_list[i][1]])
        slope_list[i] = slope_list[i][0]
        offset_err.append([offset_list[i][2], offset_list[i][1]])
        offset_list[i] = offset_list[i][0]
        dur_err.append([dur_list[i][2], dur_list[i][1]])
        dur_list[i] = dur_list[i][0]
        dep_err.append([depth_list[i][2], depth_list[i][1]])
        depth_list[i] = depth_list[i][0]
        ecc_err.append([ecc_list[i][2], ecc_list[i][1]])
        ecc_list[i] = ecc_list[i][0]
        w_err.append([w_list[i][2], w_list[i][1]])
        w_list[i] = w_list[i][0]
        i+=1

    ## Order parameters by transit number
    index = np.argsort(ser)

    ser = np.array(ser)[index]
    rp_list = np.array(rp_list)[index]
    a_list = np.array(a_list)[index]
    t0_list = np.array(t0_list)[index]
    ecw_list = np.array(ecw_list)[index]
    esw_list = np.array(esw_list)[index]
    limb1_list = np.array(limb1_list)[index]
    limb2_list = np.array(limb2_list)[index]
    slope_list = np.array(slope_list)[index]
    offset_list = np.array(offset_list)[index]
    dur_list = np.array(dur_list)[index]
    depth_list = np.array(depth_list)[index]
    ecc_list = np.array(ecc_list)[index]
    w_list = np.array(w_list)[index]
    max_list = np.array(max_list)[index]

    ## Calculate average values for each parameter
    spitzer_avg = [np.mean(rp_list), np.mean(a_list), np.mean(t0_list), np.mean(ecw_list), \
                   np.mean(esw_list), np.mean(limb1_list), np.mean(limb2_list), np.mean(slope_list), \
                   np.mean(offset_list), np.mean(dur_list), np.mean(depth_list)]
    mearth_avg = [np.mean(rp_list), np.mean(a_list), np.mean(t0_list), np.mean(ecw_list), \
                  np.mean(esw_list), np.mean(limb1_list), np.mean(limb2_list), np.mean(slope_list), \
                  np.mean(offset_list), np.mean(dur_list), np.mean(depth_list)]
    LCO_avg = [np.mean(rp_list), np.mean(a_list), np.mean(t0_list), np.mean(ecw_list), \
               np.mean(esw_list), np.mean(limb1_list), np.mean(limb2_list), np.mean(slope_list), \
               np.mean(offset_list), np.mean(dur_list), np.mean(depth_list)]


    ## Pick which plots you would like to make here.
    ## Helpful guide for order of parameters:
    #===================================================================================================================================================#
    # SET  SER  x_label  y_label  x_data  y_data  y_err  spitzer_avg  mearth_avg  LCO_avg(x_err)  (ptitle)  (ftitle)
    #===================================================================================================================================================#
    ## where ftitle is the title of the file, and ptitle is the title of the plot as shown in the pdf.

    print "transit numbers", ser

    ##B1 -- param v. transit: PLANET RADIUS
    ex.make_plot(dataset, ser, 'Transit', 'Planet radius', ser, rp_list, rp_err, spitzer_avg, mearth_avg, LCO_avg, ftitle='rp')
    ex.make_plot(dataset, ser, 'Transit', 'Eccentricity', ser, ecc_list, ecc_err, spitzer_avg, mearth_avg, LCO_avg, ftitle='ecc')
    ex.make_plot(dataset, ser, 'Transit', 'Time of Mid-Transit', ser, t0_list, t0_err, spitzer_avg, mearth_avg, LCO_avg, ftitle='t0')
    ex.make_plot(dataset, ser, 'Transit', 'Transit Depth', ser, depth_list, dep_err, spitzer_avg, mearth_avg, LCO_avg, ftitle='depth')
    ex.make_plot(dataset, ser, 'Transit', 'Transit Duration', ser, dur_list, dur_err, spitzer_avg, mearth_avg, LCO_avg, ftitle='duration')

    ##B2 -- param v. param: RADIUS V DURATION
#    ex.make_plot(dataset, ser, 'Planet radius', 'Transit duration', rp_list, dur_list, dur_err, spitzer_avg, LCO_avg, mearth_avg, x_err=rp_err, ftitle='dur-v-rp')

#==============================================================================
# ACTION ZONE
#==============================================================================

#separate as main file?
dataset = input('Which dataset are you using? Type \"S\" for Spitzer, \"M\" for MEarth, \"L\" for LCO: ')
fc.dataset = dataset
ex.dataset = dataset

if computer == 'Juno':
    if 'M' in dataset:
        print 'Analyzing MEarth data'
        list_of_files = glob.glob('/Users/enewton/Dropbox (MIT)/k225_transits/data/MEarth_*.dat')
        print list_of_files
    elif 'S' in dataset:
        print 'Analyzing Spitzer data'
        list_of_files = glob.glob('/Users/enewton/Dropbox (MIT)/k225_transits/data/K2-25_*.ascii')
        print list_of_files

elif computer == 'Vesta':
    if 'M' in dataset:
        print 'Analyzing MEarth data'
        list_of_files = glob.glob('/Users/ellie/Dropbox (MIT)/k225_transits/data/MEarth_s*.dat')
        print list_of_files
    elif 'S' in dataset:
        print 'Analyzing Spitzer data'
        list_of_files = glob.glob('/Users/ellie/Dropbox (MIT)/k225_transits/data/K2-25_*.ascii')
        print list_of_files

else:
    if 'M' in dataset:
        print 'Analyzing MEarth data'
        list_of_files = glob.glob('/Users/X-phile/Dropbox/k225_transits/data/MEarth_*.dat')
        print list_of_files
    elif 'S' in dataset:
        print 'Analyzing Spitzer data'
        list_of_files = glob.glob('/Users/X-phile/Dropbox/k225_transits/data/K2-25_*.ascii')
        print list_of_files
    elif 'L' in dataset:
        print 'Analyzing LCO data'
        list_of_files = glob.glob('/Users/X-phile/Dropbox/k225_transits/data/LCO_K225/K2-25b_*.ascii')
        print list_of_files
    elif 'GJ' in dataset:
        list_of_files = glob.glob('/Users/X-phile/Dropbox/k225_transits/data/GJ1132_data.txt')
        print list_of_files

onefile = input('Would you like to only run one file? Enter either \"filename\" or \"N\"')
if 'N' not in onefile:
    list_of_files = [onefile]

#///\\\///\\\///\\\///\\\///\\\///\\\#

ind_analysis(list_of_files, dataset)
megastack_analysis(list_of_files, dataset)
#visual_analysis(dataset)

#///\\\///\\\///\\\///\\\///\\\///\\\#
