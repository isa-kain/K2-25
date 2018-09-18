from __future__ import division, unicode_literals
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
import glob
import os

import functionsCollection as fc
import expl_plots as ex


'''
<Description: add later>

'''

#==============================================================================
# SET PREFERENCES
#==============================================================================

analyze_individual = True   #analyze individual data files
analyze_stacked = False     #stack multiple data files and analyze the stack
dataset = 'S'               #'S' for Spitzer, 'M' for MEarth

t0 = 2457062.57964 
per = 3.484564
conv_flux = False  #if True, will convert from magnitude to flux
                   #Example files: Spitzer given in flux, MEarth in magnitude.

#==============================================================================
# READ IN DATA
#==============================================================================

if 'M' in dataset:
    print 'Analyzing MEarth data'
    list_of_files = glob.glob(os.getcwd()+'/example_MEarth.dat')

elif 'S' in dataset:
    print 'Analyzing Spitzer data'
    list_of_files = glob.glob(os.getcwd()+'/example_Spitzer.ascii)


#==============================================================================
# ANALYZE INDIVIDUAL FILES
#==============================================================================
if analyze_individual:
    for f in list_of_files:
        ## Grab data from each file
        if 'M' in dataset:
            x, y, yerr = np.genfromtxt(f, delimiter=',', dtype=None, unpack=True, usecols=(0,1,2))

            numtransits = round((x[0]-t0)/per)
            midtransit = t0 + per*numtransits

            serial = str(numtransits) + '_mearth'   #example: 89_mearth

        elif 'S' in dataset:
            data = ascii.read(file, delimiter = '\s', header_start = 0, data_start = 1, guess = False)
            x = np.array(data['BMJD_OBS'])+2400000. + .5
            y = data['Normalized_Flux']

            numtransits = round((x[0]-t0)/per)
            midtransit = t0 + per*numtransits

            serial = str(numtransits) + '_spitzer'   #example: 189_spitzer

            ## Trim data to window around transit -- Spitzer example file has lots of out-of-transit data
            use = np.abs(x-midtransit)>.0208
            xsample = x[use]
            ysample = y[use]

            use_cut = (ysample<1.02) & (ysample>.98)
            ysample = ysample[use_cut]
            xsample = xsample[use_cut]

            ## Set error to standard deviation
            std = np.std(ysample)
            yerr = np.full(len(x), std)

        ## Analyze data
        fc.call_things(x, y, yerr, serial, conv_flux=conv_flux)

#==============================================================================
# ANALYZE STACKED FILES
#==============================================================================
if analyze_stacked:
    x, y, yerr = [], [], []
    
    if 'M' in dataset:
        f=open('mearth_slope_file.txt',"r")
        lines=f.readlines()
        ser, slope_list, offset_list = [], [], []
        for bit in lines:
            ser.append(bit.split(';')[0])
            slope_list.append(bit.split(';')[8])
            offset_list.append(bit.split(';')[9])
        f.close()

        ## Determine time of mid transit for first datafile
        t = np.genfromtxt(list_of_files[0], delimiter=',', dtype=None, unpack=True, usecols=(0))
        
        numtransit = round((t[0]-t0)/per)
        first_t0 = t0 + per*numtransit
        count = 0

        ## Iterate through list of files
        for f in list_of_files:
            t, m, mE = np.genfromtxt(f, delimiter=',', dtype=None, unpack=True, usecols=(0,1,2))

            numtransit = round((t[0]-t0)/per)
            midtransit = t0 + per*numtransit

            ## Grab best slope and offset values
            for num in ser:
                if int(numtransit)==int(float(num)):
                    i = ser.index(str(numtransit))
                    best_slope = float((slope_list[i][1:-2]).split(',')[0])
                    best_offset = float((offset_list[i][1:-2]).split(',')[0])

                    print 'best values: ', best_slope, best_offset

            ## Convert to flux now so linearization works
            if conv_flux:
                m = np.power(10, m/-2.5)
                mE = m*mE/-2.5*np.log(10)

            ## Linearize data
            linmodel = best_slope*(t-float(midtransit)) + best_offset
            m = m/linmodel

            ## Center all datasets besides first at first_t0
            if count!=1:
                t = t-(midtransit-first_t0)

            ## Append newly corrected data to lists
            x.extend(t)
            y.extend(m)
            yerr.extend(mE)

            count+=1

            
    elif 'S' in dataset:

        ## Determine time of mid transit for first datafile
        data = ascii.read(list_of_files[0], delimiter = '\s', header_start = 0, data_start = 1, guess = False)
        t = np.array(data['BMJD_OBS'])+2400000. + .5

        numtransit = round((t[0]-t0)/per)
        first_t0 = t0 + per*numtransit
        count = 0

        ## Iterate through all datafiles to grab and stack data
        for file in list_of_files:
            count = count+1
            data = ascii.read(file, delimiter = '\s', header_start = 0, data_start = 1, guess = False)
            t = np.array(data['BMJD_OBS'])+2400000. + .5
            m = data['Normalized_Flux']

            numtransit = round((t[0]-t0)/per)
            midtransit = t0 + per*numtransit

            ## Center all datasets besides first at first_t0
            if count!=1:
                t = t-(midtransit-first_t0)

            ## Set error as stdev of out-of-transit data
            use = np.abs(t-midtransit)>.0208
            xsample = t[use]
            ysample = m[use]
            use_cut = (ysample<1.02) & (ysample>.98)
            ysample = ysample[use_cut]
            xsample = xsample[use_cut]

            std = np.std(t)
            mE = np.full(len(t), std)

            ## Slopes negligible with Spitzer data, don't linearize
            x.extend(t)
            y.extend(m)
            yerr.extend(mE)

    ## Make into nice arrays
    x = np.array(x)
    y = np.array(y)
    yerr = np.array(yerr)

    ## Set serials
    if 'M' in dataset:
        serial = 'Megastack_mearth'
    elif 'S' in dataset:
        serial = 'Megastack_spitzer'

    ## Analyze data
    fc.call_things(x, y, yerr, serial, conv_flux=False)


