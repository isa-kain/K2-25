#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 27 17:27:22 2018

@author: X-phile
"""

import numpy as np
from astropy.io import ascii


s_data = open('spitzer_slope_file.txt', 'r')
m_data = open('mearth_slope_file.txt', 'r')
l_data = open('LCO_slope_file.txt', 'r')
ktwo_data = open('K2_ktwo_slope_file.txt', 'r')
k2sc_data = open('K2_ktwo_slope_file.txt', 'r')

data = [s_data, m_data, ktwo_data, k2sc_data] # l_data
tel = ['spitz', 'mearth', 'ktwo', 'k2sc'] #, 'LCO'

#data = [l_data]
#tel = ['LCO']


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


## Change serial from '89_mearth' to 89.0 (float)
i = 0
for num in ser:
    num = float(num.split('_')[0])
    ser[i] = num
    i+=1

## Strip unreadables
i = 0
while i<len(ser):
    t0_list[i] = np.array(t0_list[i].strip('[]').split(','), dtype='float')
    i+=1

## Pull out errors
t0_err = []
i = 0
while i<len(ser):
    t0_err.append([t0_list[i][2], t0_list[i][1]])
    t0_list[i] = t0_list[i][0]
    i+=1

## Order parameters by transit number
index = np.argsort(ser)
ser = np.array(ser)[index]
t0_list = np.array(t0_list)[index]
t0_err = (np.array(t0_err)[index]).T
tel_list = np.array(tel_list)[index]


## Create 5-col file: tnum, t0, upper error, lower error, telescope/reduction
ascii.write([ser, t0_list, t0_err[0], t0_err[1], tel_list], 'transit_times.dat', names=['TransitNum', 'Time', 'UpperErr', 'LowerErr', 'Telescope'])

