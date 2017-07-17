#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 15 09:59:43 2017

@author: X-phile
"""
from __future__ import division, unicode_literals
import numpy as np
import functionsCollection as fc
import glob
import matplotlib.pyplot as plt
import sys
    #to exit script, code 'sys.exit()'


import os
computer = os.getenv('COMPUTERNAME')

def call_things(x, y, yerr, serial, status):
    time, magnitude, magError = fc.prelim(x, y, yerr)
    print 'finished prelim'
    rp_rstar, a_rstar, myt0, params, time, magnitude, magError, thismodel = fc.initialModel(time, magnitude, magError, serial, status)
    print 'finished initial model'
    samples = fc.emceeModel(rp_rstar, a_rstar, myt0, params, time, magnitude, magError, thismodel, serial, status)
    print 'finished emcee model'
#    fc.cornerPlot(samples, serial)
#    print 'finshed corner plot'
       

if computer == 'Juno':
	list_of_files = glob.glob('/Users/enewton/Dropbox (MIT)/k225_transits/data/epic*transits*.dat')
else:
    list_of_files = glob.glob('/Users/X-phile/Dropbox/k225_transits/data/epic*transits*.dat')
#    list_of_files = glob.glob('/Users/X-phile/Dropbox/k225_transits/data/K2-25*.ascii')
##LATER: UI that asks user to select either Las Cumbres or Spitzer 


#==============================================================================
# CHUNK 1: FOR INDIVIDUAL FILES    
#==============================================================================
#
#for file in list_of_files:
#    print file
#    
#    nlist = file.split('_')
#    #nlist = /Users/X-phile/Dropbox/k225, transits/data/epic2104, tel03, 2457366, transits, corrected.dat
#    #OR /Users/X-phile/Dropbox/k225, transits/data/K2-25, 189, t1, r2, CH1
#    if 'epic' in nlist[1]:
#        if 'un' in nlist[-1]:
#            status = 'uncorr'
#        else:
#            status = 'corr'
#        serial = nlist[2] + '_' + nlist[3] + '_' + status 
#        x, y, yerr = np.genfromtxt(file, delimiter=',', dtype=None, unpack=True, usecols=(0,1,2))
#    elif 'K2-25' in nlist[1]:
#        status = ''
#        serial = nlist[2] + '_' + nlist[3] + '_' + nlist[4] + '_' + nlist[5]
#        print serial
#        x, y, e1, e2 = np.genfromtxt(file, delimiter=' ', skip_header=1, dtype=None, unpack=True, usecols=(0,1,2,3))
#        i = 0
#        while i<len(x):
#            yerr[i] = e2[i]-e1[i]
#            i+=1
#            
#
###    call_things(x, y, yerr, serial, status)
##print 'finished!'
#    
#    
#==============================================================================
# CHUNK 2: FOR STACKED FILES 
#==============================================================================
#
#for file in list_of_files:
#    nlist = file.split('_')
#    if 'uncorrected' in nlist[5]:
#        list_of_files.remove(file)
#
#file_dict = {}
#for file in list_of_files:
#    nlist = file.split('_')
#    s = nlist[3]
##    print 'name got'
#    if s not in file_dict.keys():
#        file_dict[s] = [file]
##        print 'new sublist created'
#    else:
#        ##append file to list associated with 's' key
#        file_dict[s].append(file)
##        print 'file appended to sublist'   
##    print 'finished sorting'
#print 'finished making dictionary' 
#
#         
#
#for key in file_dict:
#    x, y, yerr = [], [], []
#    print ''
#    print key ##
#    
#    for file in file_dict[key]:
##        print file ##
#        t, m, mE = np.genfromtxt(file, delimiter=',', dtype=None, unpack=True, usecols=(0,1,2))
#        x.extend(t)
#        y.extend(m)
#        yerr.extend(mE)
#        
##        print len(x) ##
#    
#    x = np.array(x)
#    y = np.array(y)
#    yerr = np.array(yerr)
#    
#    serial = key + '_stacked'
#    status = ''
#    call_things(x, y, yerr, serial, status)
#    print 'another one down'
#print 'all finished!'  


#==============================================================================
# CHUNK 3: FOR MEGASTACK
#==============================================================================

#
#file_dict = {}  ## create dictionary of all files, with date as key 
#for file in list_of_files:
#    nlist = file.split('_')
#    s = nlist[3]
##    print 'name got'
#    if s not in file_dict.keys():
#        file_dict[s] = [file]
##        print 'new sublist created'
#    else:
#        ##append file to list associated with 's' key
#        file_dict[s].append(file)
##        print 'file appended to sublist'   
##    print 'finished sorting'
#print 'finished creating file dictionary' 
#
#
#x, y, yerr = [], [], []
## theta = [guess_rp, guess_a, guess_t0, guess_secosw, guess_sesinw, guess_u1, guess_u2, guess_slope, guess_offset] 
#
#
#for key in file_dict:   ## loading individual data from each file in dict into stacked lists by date 
#    print ''
#    print key ##
#    yind = []
#    for file in file_dict[key]:
##        print file ##
#        t, m, mE = np.genfromtxt(file, delimiter=',', dtype=None, unpack=True, usecols=(0,1,2))
#        x.extend(t)
#        yind.extend(m)
#        yerr.extend(mE)
#        
##        print len(x) ##
#    
#    print 'values gotten'
#
## remove slope and offset from individual fits
#    f=open('slope_file.txt',"r")
#    lines=f.readlines()
#    ser, slope_list, offset_list = [], [], []
#    for bit in lines:
#        ser.append(bit.split(';')[0])
#        slope_list.append(bit.split(';')[8])
#        offset_list.append(bit.split(';')[9]) 
#   
#    f.close()
#    print 'slope file parsed'
#
#    
#    for num in ser: ##iterating through list of serials
#        if key in num:  ##if the key is in the serial
#            foundit = key + '_stacked'
#            i = ser.index(foundit)
#            best_slope = float((slope_list[i][1:-2]).split(',')[0])
#            best_offset = float((offset_list[i][1:-2]).split(',')[0])
#            
#            print 'best values: ', best_slope, best_offset
#            
#    t0 = 2457062.57964
#    per = 3.484564 
#    numtransits = round((x[0]-t0)/per)
#    myt0 = t0 + per*numtransits   
#    tlin = np.linspace(myt0-0.2, myt0+0.2, 1000)
#    linmodel = best_slope*(tlin[0]-float(x[0])) + best_offset 
#    i = 0
#    while i<len(yind):
#        yind[i] = yind[i] - linmodel
#        i+=1
#    y.extend(yind)
#
#    print len(x)
#    print 'finished removing slope and offset from individual stack'
#    
#    
#
#x = np.array(x)
#y = np.array(y)
#yerr = np.array(yerr)
#
#key_list = file_dict.keys()
#    
#print 'all done!' 
#print ''   
#
#serial = 'Megastack'
#status = ''
#call_things(x, y, yerr, serial, status)
#print 'finished!!!'






#==============================================================================
# EXPLANATORY PLOTS
#==============================================================================
# secosw and sesinw vs. serial
# u1 and u2 vs. serial
# everything above but with error bars
# "with little elipses that show the correlations between the parameters"
# compare planet parameters as measured from diff telescopes from same transit


f=open('slope_file.txt',"r")
lines=f.readlines()
ser, rp_list, a_list, t0_list, ecw_list, esw_list, limb1_list, limb2_list, slope_list, offset_list, dur_list, depth_list = [], [], [], [], [], [], [], [], [], [], [], []
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
f.close()

i = 0
for num in ser:
    num = float(num.split('_')[0])
    ser[i] = num
    i+=1

rp_err, a_err, t0_err, ecw_err, esw_err, limb1_err, limb2_err, slope_err, offset_err, dur_err, dep_err = [], [], [], [], [], [], [], [], [], [], []

i = 0
while i<len(ser):    
    rp_list[i] = np.array(rp_list[i].strip('()').split(','), dtype=float)
    a_list[i] = np.array(a_list[i].strip('()').split(','), dtype='float')
    t0_list[i] = np.array(t0_list[i].strip('()').split(','), dtype='float')
    ecw_list[i] = np.array(ecw_list[i].strip('()').split(','), dtype='float')
    esw_list[i] = np.array(esw_list[i].strip('()').split(','), dtype='float')
    limb1_list[i] = np.array(limb1_list[i].strip('()').split(','), dtype='float')
    limb2_list[i] = np.array(limb2_list[i].strip('()').split(','), dtype='float')
    slope_list[i] = np.array(slope_list[i].strip('()').split(','), dtype='float')
    offset_list[i] = np.array(offset_list[i].strip('()').split(','), dtype='float')
    dur_list[i] = np.array(dur_list[i].strip('()').split(','), dtype='float')
    depth_list[i] = np.array(depth_list[i].strip('()').split(','), dtype='float')
    i+=1
    
i = 0
while i<len(ser):
    rp_err.append([rp_list[i][1], rp_list[i][2]])
    rp_list[i] = rp_list[i][0]
    a_err.append([a_list[i][1], a_list[i][2]])
    a_list[i] = a_list[i][0]
    t0_err.append([t0_list[i][1], t0_list[i][2]])
    t0_list[i] = t0_list[i][0]
    ecw_err.append([ecw_list[i][1], ecw_list[i][2]])
    ecw_list[i] = ecw_list[i][0]
    esw_err.append([esw_list[i][1], esw_list[i][2]])
    esw_list[i] = esw_list[i][0]
    limb1_err.append([limb1_list[i][1], limb1_list[i][2]])
    limb1_list[i] = limb1_list[i][0]
    limb2_err.append([limb2_list[i][1], limb2_list[i][2]])
    limb2_list[i] = limb2_list[i][0]
    slope_err.append([slope_list[i][1], slope_list[i][2]])
    slope_list[i] = slope_list[i][0]
    offset_err.append([offset_list[i][1], offset_list[i][2]])
    offset_list[i] = offset_list[i][0]
    dur_err.append([dur_list[i][0], dur_list[i][2]])
    dur_list[i] = dur_list[i][1]
    dep_err.append([depth_list[i][0], depth_list[i][2]])
    depth_list[i] = depth_list[i][1]
    i+=1

#### rp vs. date
#plt.figure()
#plt.scatter(ser, rp_list, c = 'b', alpha = .6, zorder = 1)
#plt.errorbar(ser, rp_list, np.array(rp_err).T,  fmt='.', c = 'k', alpha = .4, zorder = 0)
#plt.xlabel('Date of data collection')
#plt.ylabel('Radius of planet')
#plt.title('Planet radius v. date')
#plt.savefig('rp_plot.pdf')


### rp vs. transit duration
plt.figure()
plt.scatter(rp_list, dur_list, c = 'c', alpha =.6)
plt.errorbar(rp_list, dur_list, np.array(dur_err).T,  fmt='.', c = 'b', alpha = .4, zorder = 0)
plt.xlabel('Radius of planet')
plt.ylabel('Duration of transit')
plt.title('Duration of transit v. planet radius')
plt.savefig('rp-v-duration_plot.pdf')

#### t0 vs. serial
#plt.figure()
#plt.scatter(ser, t0_list, c = 'm', alpha = .6, zorder = 1)
#plt.errorbar(ser, t0_list, np.array(t0_err).T, fmt='.', c = 'k', alpha = .4, zorder = 0)
#plt.xlabel('Date of data collection')
#plt.ylabel('T0')
#plt.title('T0 v. date')
#plt.savefig('t0_plot.pdf')
#
#### ecc
#plt.figure()
#plt.scatter(ser, ecw_list, c = 'b', alpha = .6)
#plt.errorbar(ser, ecw_list, np.array(ecw_err).T,  fmt='.', c = 'c', alpha = .4, zorder = 0)
#plt.scatter(ser, esw_list, c = 'r', alpha = .6)
#plt.errorbar(ser, esw_list, np.array(esw_err).T,  fmt='.', c = 'm', alpha = .4, zorder = 0)
#plt.xlabel('Date of data collection')
#plt.ylabel('Eccentricity')
#plt.title('Eccentricity v. date')
#plt.savefig('ecc_plot.pdf')
#
#### limb
#plt.figure()
#u1 = plt.scatter(ser, limb1_list, c = 'r', alpha = .6, zorder = 1, label = 'u1')
#plt.errorbar(ser, limb1_list, np.array(limb1_err).T,fmt='.', c = 'm', alpha = .4,zorder = 0)
#u2 = plt.scatter(ser, limb2_list, c = 'b', alpha = .6, zorder = 1, label = 'u2')
#plt.errorbar(ser, limb2_list, np.array(limb2_err).T,fmt='.', c = 'c', alpha = .4,zorder = 0)
#plt.xlabel('Date of data collection')
#plt.ylabel('Limb darkening parameters')
#plt.legend(handles=[u1, u2], loc = 2)
#plt.title('Limb darkening v. date')
#plt.savefig('limb_plot.pdf')
#
#### slope, offset vs. date
#plt.figure()
#s = plt.scatter(ser, slope_list, c = 'm', alpha = .6, label = 'slope')
#plt.errorbar(ser, slope_list, np.array(slope_err).T, fmt='.', c = 'r', alpha = .4, zorder = 0)
#o = plt.scatter(ser, offset_list, c = 'c', alpha = .6, label = 'offset')
#plt.errorbar(ser, offset_list, np.array(offset_err).T, fmt='.', c = 'b', alpha = .4, zorder = 0)
#plt.xlabel('Date of data collection')
#plt.ylabel('Slope, offset')
#plt.legend(handles=[s, o], loc = 2)
#plt.title('Slope, offset v. date')
#plt.savefig('slope-offset_plot.pdf')


## rp v. depth


### depth v. duration
plt.figure()
de = plt.scatter(depth_list, dur_list, c = 'm', alpha = .6, label = 'slope')
plt.errorbar(depth_list, dur_list, np.array(dur_list).T, fmt='.', c = 'r', alpha = .4, zorder = 0)
axes = plt.gca()
axes.set_xlim([0.0, 0.002])
axes.set_ylim([-.005, .015])
plt.xlabel('Depth of transit')
plt.ylabel('Duration of transit')
plt.legend(handles=[de], loc = 2)
plt.title('Transit depth v. duration')
plt.savefig('depth-v-duration_plot.pdf')


### rp v. depth, duration
plt.figure()
s = plt.scatter(rp_list, depth_list, c = 'b', alpha = .6, label = 'depth')
plt.errorbar(rp_list, depth_list, np.array(dep_err).T, fmt='.', c = 'c', alpha = .4, zorder = 0)
o = plt.scatter(rp_list, dur_list, c = 'r', alpha = .6, label = 'duration')
plt.errorbar(rp_list, dur_list, np.array(dur_err).T, fmt='.', c = 'y', alpha = .4, zorder = 0)
plt.xlabel('Radius of planet')
plt.ylabel('Depth, duration of transit')
plt.legend(handles=[s, o], loc = 3)
plt.title('Depth, duration v. planet radius')
plt.legend()
plt.savefig('rp-v-depth-dur_plot.pdf')


    