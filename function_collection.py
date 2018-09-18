#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 15 09:49:10 2017

@author: X-phile
"""
import batman
import numpy as np
import math
import matplotlib.pyplot as plt
import emcee
import corner
import re
import copy
import time as timepackage
from zachopy import oned as bin
#from plotwalkers import walkers
import pickle

global dataset

def call_things(x, y, yerr, serial, conv_flux=True):
    time, magnitude, magError = prelim(x, y, yerr, conv_flux=conv_flux)
    print 'finished prelim'
    rp_rstar, a_rstar, midtransit, params, time, magnitude, magError, thismodel, transit_num = initialModel(time, magnitude, magError, serial)
    print 'finished initial model'
    samples, sampler, variables = emceeModel(rp_rstar, a_rstar, midtransit, params, time, magnitude, magError, thismodel, serial, transit_num)
    print 'finished emcee model'
    cornerPlot(samples, serial, variables)
    print 'finshed corner plot'
    #pw.walkers(sampler, serial, variables)
    #print 'finished plotting walkers'

# Kipping 2012
# https://academic.oup.com/mnras/article/435/3/2152/1024138/Efficient-uninformative-sampling-of-limb-darkening#18186712
# Eqs 15-18
def calc_u1u2(q1, q2):
	u1 = 2.*np.sqrt(q1)*q2
	u2 = np.sqrt(q1)*(1.-2.*q2)
	return u1, u2

def calc_q1q2(u1, u2):
	q1 = (u1 + u2)**2.
	q2 = u1/( 2*(u1+u2) )
	return q1, q2

# Eastman ??
def calc_eccw(sec, ses):
    ecc = sec**2. + ses**2.
    w = np.arctan2(ses, sec)
    return ecc, np.rad2deg(w)

def calc_escw(ecc, w):
    sec = np.sqrt(ecc)*np.cos(np.deg2rad(w))
    ses = np.sqrt(ecc)*np.sin(np.deg2rad(w))
    return sec, ses

# Dawson & Johnson
# requires very long term to be added to log likelihood
def calc_g(ecc, w):
    return 1.+ecc*np.sin(np.deg2rad(w))/np.sqrt(1.-ecc**2)

## find_key: combs through parameter file, returns true if key is found
def find_key(key):
    if 'Megastack' in key:
        with open('megaslope_file.txt', 'r') as slope_file:
            if re.search(key, slope_file.read()):
                return True
            else:
                return False
    elif 'mearth' in key:
        with open('mearth_slope_file.txt', 'r') as slope_file:
            key = key.split('_')[0] + '.0;'
            if re.search(key, slope_file.read()):
                return True
            else:
                return False
    elif 'spitzer' in key:
        with open('spitzer_slope_file.txt', 'r') as slope_file:
            key = key.split('_')[0] + '.0;'
            if re.search(key, slope_file.read()):
                return True
            else:
                return False
    elif 'LCO' in key:
        with open('LCO_slope_file.txt', 'r') as slope_file:
            key = key.split('_')[0] + '.0;'
            if re.search(key, slope_file.read()):
                return True
            else:
                return False


def lnprior(theta, minmax, a_mean=21.8, a_sigma=1.5572, rp_mean=0.11, rp_sigma=0.015):
    assert len(theta) == len(minmax)
    for i in range(len(theta)):
        if (theta[i] < minmax[i][0]) or (theta[i] > minmax[i][1]): # a parameter is out of range
            return -np.inf

    ## Gaussian prior on a/R*, rp
    lnprior_a = - (1.0/2.0)*((theta[1]-a_mean)/a_sigma)**2.0

    lnprior_r = - (1.0/2.0)*((theta[0]-rp_mean)/rp_sigma)**2.0

    ## Don't let the planet go inside the star
    ecc, w = calc_eccw(theta[3],theta[4])
    if (theta[1]*(1.-ecc) < 1):
        return -np.inf

    ## ecc must be 0 to 1
    if (w < -180.) or (w > 180.):
        return -np.inf
    if (ecc < 0) or (ecc > 1):
        return -np.inf

    return lnprior_r + lnprior_a

def lnlikelihood(theta, params, model, t, flux, err, myt0): #flux = adjusted magnitude, err = adjusted magError
    
    ## Initialize batman
    params.rp = copy.deepcopy(theta[0])
    params.a = copy.deepcopy(theta[1])
    params.t0 = copy.deepcopy(theta[2])
    ecc, w = calc_eccw(theta[3], theta[4])
    params.ecc = copy.deepcopy(ecc)
    params.w = copy.deepcopy(w)

    u1, u2 = calc_u1u2(theta[5], theta[6])
    params.u = copy.deepcopy([u1, u2])
    inclination = theta[9] # impact = math.acos(theta[9]/theta[1])*180./math.pi
    params.inc = copy.deepcopy(inclination)
    ltcurve = model.light_curve(params)

    ## Calc lnlikelihood
    residuals = flux/(theta[7]*(t-myt0) + theta[8]) - ltcurve
    ln_likelihood = -0.5*(np.sum((residuals/err)**2 + np.log(2*np.pi*(err)**2)))
    if np.random.random() < 1e-4:
     	print timepackage.ctime(timepackage.time()), "lnlikelihood:", ln_likelihood
        print theta

    ## Examine bad values
    if ln_likelihood<0:
        return -np.inf
    if not np.all(np.isfinite(residuals)):
        return -np.inf
    if not np.all(np.isfinite(err)):
        return -np.inf

    return ln_likelihood


def lnprobability(theta, params, model, t, flux, err, variables, myt0):
    order = ['rp', 'a', 't0', 'secosw', 'sesinw', 'q1', 'q2', 'slope', 'offset', 'inc']
    theta_all = []
    range_all = []

    ## Determine which variables are being varied
    i = 0
    for p in order:
        range_all.append( [variables[p]['min'],variables[p]['max']] )
        if variables[p]['vary']:
            theta_all.append(theta[i])
            i = i+1
        else:
            theta_all.append(variables[p]['value'])

    ## Find prior, likelihood
    assert (i) == len(theta)
    lp = lnprior(theta_all, range_all)
#    print 'lp: ', lp
    if not np.isfinite(lp):
        return -np.inf
    lli = lnlikelihood(theta_all, params, model, t, flux, err, myt0)
#    print 'lli: ', lli
    return lp + lli

##==================================================##
##==================================================##

def prelim(x, y, yerr, conv_flux=True):
    m = np.isfinite(x) & np.isfinite(y) & np.isfinite(yerr)
    time = x[m]
    magnitude = y[m]
    magError = yerr[m]

    ## If magnitude, convert to flux
    if conv_flux==True:
        magnitude = np.power(10, magnitude/-2.5)
        magError = magnitude*magError/-2.5*np.log(10)

#    print 'hopefully flux: ', magnitude

    ## Normalize flux for LCO data
    if 'L' in dataset:
        magError = magError/np.median(magnitude)
        magnitude = magnitude/np.median(magnitude)

    index = np.argsort(time)
    time = time[index]
    magnitude = magnitude[index]
    magError = magError[index]

    t0 = 2457062.57964
    per = 3.484552     #derived 3.484564
    numtransits = round((time[0]-t0)/per)
    myt0 = t0 + per*numtransits
    print 'myt0: ', myt0

    print len(time)

    ## Cut out-of-transit data   -- REDUCES TIME TO ZERO IF T0 EARLIER THAN TIME[0]
    ind = 0
    limit = np.ones(len(time), dtype = bool)
    for t in time:
        if t<=(myt0-.1) or t>=(myt0+.1):
            limit[ind] = False
        else:
            limit[ind] = True
        ind+=1

    time = time[limit]
    magnitude = magnitude[limit]
    magError = magError[limit]

    print len(time)

    ## Filter outliers
    for i in range(0,2):
        std = np.std(magnitude)
        avg = np.mean(magnitude)

        ind = 0
        check = np.ones(len(time), dtype = bool)
        for m in magnitude:
            if m>=(avg + (3*std)) or m<=(avg - (3*std)):
                check[ind] = False
            else:
                check[ind] = True
            ind+=1
        time = time[check]
        magnitude = magnitude[check]
        magError = magError[check]

    print len(time)

    return time, magnitude, magError


##==================================================##

def initialModel(time, magnitude, magError, serial):

    per = 3.484552     #derived 3.484564
    rp_rstar = 0.10
    #impact = .55
    a_rstar = 21.8
    #inclination = math.acos(impact/a_rstar)*180./math.pi
    inclination = 88.

    t0 = 2457062.57964
    numtransits = round((min(time)-t0)/per)
    myt0 = t0 + per*numtransits

    ## Fill in params for batman
    params = batman.TransitParams()
    params.t0 = myt0
    params.per = per
    params.rp = rp_rstar
    params.a = a_rstar
    params.inc = inclination
    params.ecc = 0.
    params.w = 90.
    params.ourlinearcoefficient = 0.
    params.ouroffset = 1.
    if 'S' in str(dataset):
        params.u = [0.0687, 0.1955]
#        params.u = [0.313, 0.154]  ##GJ 1132
    elif 'M' in str(dataset):
        params.u = [0.1, 0.3]
#        params.u = [0.215, 0.407]  ##GJ 1132
    elif 'L' in str(dataset):
        params.u = [0.1, 0.3] ##what should this be actually?
    params.limb_dark = "quadratic"

    ## Initialize batman
    thismodel = batman.TransitModel(params, time)
    flux = thismodel.light_curve(params)

    ## Bin data
    bt, by, bye = bin.binto(x=((time-myt0)*24), y=magnitude, yuncertainty=magError, binwidth=.01)
    print 'binning data, 1', len(time), len(bt)
    binned = True

    ## Plot initial fit
    plt.figure()
    if binned:
        plt.errorbar(bt, by, bye, fmt='o', c='k', alpha=.4)
    else:
        plt.errorbar((time-myt0)*24, magnitude, magError, fmt='o', c='k', alpha=.4)
    plt.plot((time-myt0)*24, flux, c='r', alpha = .7, linewidth=3, zorder = 10)
    plt.xlabel('Time from mid-transit in hours')
    plt.ylabel('Flux relative to baseline')
    if 'Megastack' in str(serial):
        plt.title('Initial fit of stacked transits')
    else:
        plt.title('Initial fit of transit ') #+ str(numtransits))
    title = str(serial) + '.pdf'
    print title
    plt.savefig('../plots/' + title)
#    plt.show()
    plt.close()

    return rp_rstar, a_rstar, myt0, params, time, magnitude, magError, thismodel, numtransits

##==================================================##

def emceeModel(rp_rstar, a_rstar, myt0, params, time, magnitude, magError, thismodel, serial, transit_num):

    ## Open parameter files
    mearth_slope_file = open('mearth_slope_file.txt', 'a+')
    megaslope_file = open('megaslope_file.txt', 'a+')
    spitzer_slope_file = open('spitzer_slope_file.txt', 'a+')
    LCO_slope_file = open('LCO_slope_file.txt', 'a+')

    ## Create dictionary of parameters
    oneparam = {'texname':'', 'vary': True, 'value':1.,
            'scale': 1., 'min':0., 'max':1.,
            'best': [np.nan, np.nan, np.nan]}
    variables = {'rp': copy.deepcopy(oneparam), 'a': copy.deepcopy(oneparam), 't0': copy.deepcopy(oneparam),
             'secosw': copy.deepcopy(oneparam), 'sesinw': copy.deepcopy(oneparam),
             'q1': copy.deepcopy(oneparam), 'q2': copy.deepcopy(oneparam),
             'slope': copy.deepcopy(oneparam), 'offset': copy.deepcopy(oneparam), 'inc': copy.deepcopy(oneparam)}

    guess_secosw, guess_sesinw  = calc_escw(params.ecc, params.w)
    guess_q1, guess_q2 = calc_q1q2(params.u[0], params.u[1])

    ## Set values, bounds for each parameter; determine if each will vary
    p = 'rp'
    variables[p]['texname'] = r'$rp$'
    variables[p]['value'] = params.rp
    variables[p]['vary'] = True
    variables[p]['scale'] = 1.e-3
    variables[p]['min'] = 0.
    variables[p]['max'] = 1.

    p = 'a'
    variables[p]['texname'] = r'$a$'
    variables[p]['value'] = params.a
    variables[p]['vary'] = True
    variables[p]['scale'] = 1.e-2
    variables[p]['min'] = 1.
    variables[p]['max'] = 30.

    p = 't0'
    variables[p]['texname'] = r'$t0$'
    variables[p]['value'] = myt0
    variables[p]['vary'] = True
    variables[p]['scale'] = 1.e-4
    variables[p]['min'] = 2456000.
    variables[p]['max'] = 2458000.

    p = 'secosw'
    variables[p]['texname'] = r'$secosw$'
    variables[p]['value'] = -0.42082681735021493
    variables[p]['vary'] = False
    variables[p]['scale'] = 0.1
    variables[p]['min'] = -1.
    variables[p]['max'] = 1.

    p = 'sesinw'
    variables[p]['texname'] = r'$sesinw$'
    variables[p]['value'] = 0.0023165374683225028
    variables[p]['vary'] = False
    variables[p]['scale'] = 0.1
    variables[p]['min'] = -1.
    variables[p]['max'] = 1.

    p = 'q1'
    variables[p]['texname'] = r'$q1$'
    variables[p]['value'] = guess_q1
    variables[p]['vary'] = False
    variables[p]['scale'] = 1.e-3
    variables[p]['min'] = 0.
    variables[p]['max'] = 1.

    p = 'q2'
    variables[p]['texname'] = r'$q2$'
    variables[p]['value'] = guess_q2
    variables[p]['vary'] = False
    variables[p]['scale'] = 1.e-3
    variables[p]['min'] = 0.
    variables[p]['max'] = 1.

    p = 'slope'
    variables[p]['texname'] = r'$slope$'
    variables[p]['value'] = 0.
    variables[p]['vary'] = True
    variables[p]['scale'] = 1.e-3
    variables[p]['min'] = -5.
    variables[p]['max'] = 5.

    p = 'offset'
    variables[p]['texname'] = r'$offset$'
    variables[p]['value'] = 1.
    variables[p]['vary'] = True
    variables[p]['scale'] = 1.e-3
    variables[p]['min'] = -5.
    variables[p]['max'] = 5.

    p = 'inc'
    variables[p]['texname'] = r'$i$'
    variables[p]['value'] = 88. #math.cos(params.inc*math.pi/180.)*params.a
    variables[p]['vary'] = True
    variables[p]['scale'] = 1.e-2
    variables[p]['min'] = 0.
    variables[p]['max'] = 90.

    ## Establish order & names of parameters:
    param_order = ['rp', 'a', 't0', 'secosw', 'sesinw', 'q1', 'q2', 'slope', 'offset', 'inc']

    ## Determine which parameters are being varied
    varyparams = [] # list of parameters that are being varied this run
    theta, scale, mins, maxs = [], [], [], [] # to be filled with parameter values

    for p in param_order: # record if this parameter is being varied
        if variables[p]['vary']:
            varyparams.append(p)
            theta.append(variables[p]['value'])
            scale.append(variables[p]['scale'])  ## how are we deciding scale?
            mins.append(variables[p]['min'])
            maxs.append(variables[p]['max'])
        else: # if parameter fixed, just record the starting value as the best value
            variables[p]['best'][0] = variables[p]['value']

#    print "Varying: ", varyparams
#    for p in variables.keys():
#        print p, variables[p]['value']

    ## Number of walkers, dimensions
    ndimensions = len(varyparams)
    nwalkers = 30

    ## Determine starting positions
    np.random.seed(42)
    pos = [theta + scale*np.random.randn(ndimensions) for i in range (nwalkers)]

    ## Replace each u with q
    for p in pos:
        q1, q2 = calc_q1q2(params.u[0]+1.e-3*np.random.randn(), params.u[1]+1.e-3*np.random.randn())
        p[varyparams=='q1']=q1
        p[varyparams=='q2']=q2

        ## Determine starting positions for secosw, sesinw based on ecc and w values
        if ('secosw' in varyparams) and ('sesinw' in varyparams):
            e = np.random.uniform(low=0., high=0.6)
            w = np.random.uniform(low=0., high=180.)
            sec, ses = calc_escw(e, w)
            p[varyparams=='secosw']=sec
            p[varyparams=='sesinw']=ses

    #------------------------------#
    # Run the MCMC chain
    #------------------------------#

    sampler = emcee.EnsembleSampler(nwalkers, ndimensions, lnprobability,
                                    args = (params, thismodel, time, magnitude, magError, variables, myt0))

    print "Starting MC", timepackage.ctime(timepackage.time())
    sampler.run_mcmc(pos, 10000) ##changed from pos_to_vary
    print "Finishing MC", timepackage.ctime(timepackage.time())

    #------------------------------#
    # Analyze results from MCMC
    #------------------------------#

    ## Remove the "burn-in" and calculate the best values
    samples = sampler.chain[:, 8000:, :].reshape((-1, ndimensions))
    best = map(lambda v: [v[1], v[2]-v[1], v[1]-v[0]], \
                       zip(*np.percentile(samples, [16, 50, 84], axis=0))) ## arranged: [50th, upper error, lower error]

    ## Pickle sampler
    filename = 'pickle_' + str(serial)
    sfile = open(filename, 'w+')
    pickle.dump(samples, sfile)
    sfile.close()
    print 'pickled'


    ## Update variables with best value
    i = 0
    for p in param_order:
        if variables[p]['vary']:
            variables[p]['best'] = best[i]
            i = i+1
            print p, variables[p]['best']
    assert (i) == len(theta)

    ## Calculate secondary best-fit parameters
    durations = []
    depths = []
    eccs = []
    omegas = []
    samp = {'rp':np.nan, 'a':np.nan, 't0':np.nan, 'secosw':np.nan, 'sesinw':np.nan, 'q1':np.nan, 'q2':np.nan, 'slope':np.nan, 'offset':np.nan, 'inc':np.nan}
    for i in samples:
        j = 0
        for p in param_order:
            if variables[p]['vary']:
                samp[p] = i[j]
                j+=1
            else:
                samp[p] = variables[p]['value']
        ecc, w = calc_eccw(samp['secosw'], samp['sesinw'])
#        inc = math.acos(samp['inc']/samp['a']) ##apparently these aren't single values?
        inc = np.deg2rad(samp['inc'])

        eccs.append(ecc)
        omegas.append(w)

        if ((1+samp['rp'])**2)>=(samp['a']*math.cos(inc)):
            if inc>=-math.pi/2 and inc<=math.pi/2 and inc!=0:
                # Calculating transit duration
                p1 = (1+samp['rp'])**2-(samp['a']*math.cos(inc))
                p2 = math.sin(inc)
                p3 = math.sqrt(p1)
                d1 = p3/p2

        d2 = math.asin((1/samp['a'])*(d1))
        dur = (params.per/math.pi)*d2*(math.sqrt(1-ecc**2)/(1+samp['sesinw']))
        durations.append(dur)

        ## Calculating transit depth
        dep = samp['rp']**2
        depths.append(dep)

    ## Listing best durations
    temp = np.percentile(durations, [16, 50, 84], axis=0)    ## arranged: [lower error, 50th, upper error]
    best_duration = [temp[1], temp[2]-temp[1], temp[1]-temp[0]]
    best_duration = str(best_duration)
    best_duration = '[' + best_duration[1:-1] + ']'
#    best_duration = '0'

    ## Listing best depths
    temp = np.percentile(depths, [16, 50, 84], axis=0)
    best_depth = [temp[1], temp[2]-temp[1], temp[1]-temp[0]]
    best_depth = str(best_depth)
    best_depth = '[' + best_depth[1:-1] + ']'

    ## Listing best eccentricities
    temp = np.percentile(eccs, [16, 50, 84], axis=0)
    best_ecc = [temp[1], temp[2]-temp[1], temp[1]-temp[0]]
    best_ecc = str(best_ecc)
    best_ecc = '[' + best_ecc[1:-1] + ']'

    ## Listing best omegas
    temp = np.percentile(omegas, [16, 50, 84], axis=0)
    best_w = [temp[1], temp[2]-temp[1], temp[1]-temp[0]]
    best_w = str(best_w)
    best_w = '[' + best_w[1:-1] + ']'

    ## Get the posterior, use max likelihood because median could be anywhere
    max = sampler.flatchain[np.argmax(sampler.flatlnprobability)]


    #-----------------------------------------#
    # Write out the best fit to slope file
    #-----------------------------------------#

    split = serial.split('_')
    # 89_mearth  OR  189_spitzer  OR  Megastack_spitzer  OR  Megastack_mearth

    if 'Megastack' in str(serial):
        tag = serial
    elif 'spitzer' in split[1] or 'mearth' in split[1] or 'LCO' in split[1] or 'GJ' in split[0]:
        tag = split[0]
    print tag

    thing = str(tag) + ';' + str(variables['rp']['best']) + ';' + str(variables['a']['best']) + \
                    ';' + str(variables['t0']['best']) + ';' + str(variables['secosw']['best']) + \
                    ';' + str(variables['sesinw']['best']) + ';' + str(variables['q1']['best']) + \
                    ';' + str(variables['q2']['best']) + ';' + str(variables['slope']['best']) + \
                    ';' + str(variables['offset']['best']) + ';' + best_duration + ';' + best_depth + \
                    ';' + str(variables['inc']['best']) + best_ecc + ';' + best_w + ';' + str(max) + '; \n'

    ## Write thing to file
#    print "serial:", serial
#    print 'thing:', thing
    if 'Megastack' in str(serial):
        print "found em!"
        if find_key(str(serial))==False:
            megaslope_file.write(thing)
            megaslope_file.close()
            print 'line inserted'
    elif 'mearth' in str(serial):
        print "found em!"
        if find_key(str(serial))==False:
            mearth_slope_file.write(thing)
            mearth_slope_file.close()
            print 'line inserted'
    elif 'spitzer' in str(serial):
        print "found em!"
        if find_key(str(serial))==False:
            spitzer_slope_file.write(thing)
            spitzer_slope_file.close()
            print 'line inserted'
    elif 'LCO' in str(serial):
        print "found em!"
        if find_key(str(serial))==False:
            LCO_slope_file.write(thing)
            LCO_slope_file.close()
            print 'line inserted'

    #--------------------#
    # Plot best fit
    #--------------------#

    ## Initialize and plot batman fits from all theta values
    theta_best =  np.percentile(samples, 50, axis=0)
    tlin = np.linspace(myt0-0.2, myt0+0.2, 1000)

    plt.figure()
    for theta in samples[np.random.randint(len(samples), size=100)]:
        lnprobability(theta, params, thismodel, time, magnitude, magError, variables, myt0)
        mlin = batman.TransitModel(params, tlin)  #initialize model

        flux = mlin.light_curve(params)*(variables['slope']['best'][0]*(tlin-time[0]) + variables['offset']['best'][0])
        plt.plot((tlin-myt0)*24, flux, c='b', linewidth=3, alpha=0.05, zorder=3)

    ## Initialize and plot batman fit from best values
    lnprobability(theta_best, params, thismodel, time, magnitude, magError, variables, myt0)
    mlin = batman.TransitModel(params, tlin)
    print 'method flux: ', np.median(flux)
    newflux = mlin.light_curve(params)*(variables['slope']['best'][0]*(tlin-time[0]) + variables['offset']['best'][0])
    print 'flux, where she be? ', np.median(newflux)

    ## Bin data
    binned = False
    if 'S' in str(dataset) or 'stacked' in serial or 'Megastack' in serial:
        bt, by, bye = bin.binto(x=time, y=magnitude, yuncertainty=magError, binwidth=.005) #binwidth=.0004
        print 'binning data, 2'
        binned = True
    if 'M' in str(dataset):
        bt, by, bye = bin.binto(x=time, y=magnitude, yuncertainty=magError, binwidth=.0004)
        print 'binning data, 2'
        binned = True

    if binned:
        plt.errorbar((bt-myt0)*24, by, bye, fmt='o', c='k', alpha=0.4, zorder = 0)
    else:
        plt.errorbar((time-myt0)*24, magnitude, magError, fmt='o', c='k', alpha=0.4, zorder = 0)
    plt.plot((tlin-myt0)*24, newflux, c='purple', linewidth=1, alpha=0.6, zorder = 3)
    plt.ylim(0.97,1.02)   #will this be constant across all datasets?

    plt.xlim((time[0]-myt0)*24, (time[-1]-myt0)*24)    #(time[0], time[-1])
    plt.xlabel('Time from mid-transit in hours')
    plt.ylabel('Flux relative to baseline')
    if 'Megastack' in str(serial):
        plt.title('Best fit of stacked transits')
    else:
        plt.title('Best fit of transit ' + str(serial))
    title = str(serial) + '_best.pdf'
    print title

    plt.savefig('../plots/'+title)
#    plt.show()
    plt.close()

    return samples, sampler, variables

##==================================================##


def cornerPlot(samples, serial, variables):
    ## Grab a, ecc, w sample values
    param_order = ['rp', 'a', 't0', 'secosw', 'sesinw', 'q1', 'q2', 'slope', 'offset', 'inc']
    i = 0
    labels = []
    do_corner = True
    for p in param_order:
        if variables[p]['vary']:
            labels.append(variables[p]['texname'])
            if p == 'a':
                a_ind = i
            if p == 'secosw':
                c_ind = i
            if p == 'sesinw':
                s_ind = i
            i = i+1
        else:
            if p in [ 'secosw', 'sesinw']:
                do_corner = False

    ## Corner plot for a, ecc, w
    if do_corner:
        ecc = np.zeros_like(samples[:,0])
        w = np.zeros_like(samples[:,0])
        for i, sam in enumerate(samples):
            e, omega = calc_eccw(sam[c_ind], sam[s_ind])
            ecc[i] = e
            w[i] = omega

        s = np.array([samples[:,1],ecc, w]).transpose()

        plt.figure()
        corner.corner(s, labels = ['a', 'ecc', 'w'])
        title = serial + '_corner_ecc.pdf'
        plt.savefig('../plots/' + title)
        #plt.show()
        plt.close()

    plt.figure()
    corner.corner(samples, labels = labels)
    title = serial + '_corner.pdf'
    #plt.show()
    plt.savefig('../plots/' + title)
    plt.close()



