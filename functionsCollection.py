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

def calc_eccw(sec, ses):  
    ecc = sec**2. + ses**2.
    w = np.arctan2(ses, sec)
    return ecc, np.rad2deg(w)

def calc_escw(ecc, w):
    sec = np.sqrt(ecc)*np.cos(np.deg2rad(w))
    ses = np.sqrt(ecc)*np.sin(np.deg2rad(w))
    return sec, ses    


def find_key(key):
    if 'stacked' in key:
        with open('slope_file.txt', 'r') as slope_file:
            if re.search(key, slope_file.read()):
                return True
            else:
                return False
    elif 'Megastack' in key:    
        with open('megaslope_file.txt', 'r') as slope_file:
            if re.search(key, slope_file.read()):
                return True
            else:
                return False

def lnprior(theta): 
    ecc, w = calc_eccw(theta[3], theta[4])
    if (ecc < 0.) or (ecc > 1.):
        return np.inf
    if (w < 0.) or (w > 180.):
        return np.inf
    if (theta[5] < 0) or (theta[6] < 0): #limb darkening
        return np.inf
    if (theta[6]>1.0):     
        return np.inf      
    if theta[6]<0 or theta[6]>1:     ## how to do soft limit?
        return np.inf
    return 0


#theta = [guess_rp, guess_a, guess_t0, guess_secosw, guess_sesinw, guess_u1, guess_u2, guess_slope]


def lnlikelihood(theta, params, model, t, flux, err): #flux = adjusted magnitude, err = adjusted magError
    params.rp = theta[0]
    params.a = theta[1]
    params.t0 = theta[2]
    
    ecc, w = calc_eccw(theta[3], theta[4])
    params.ecc = ecc
    params.w = w
    
    params.u = [theta[5], theta[6]]
    ltcurve = model.light_curve(params) 
    residuals = flux - ltcurve - theta[7]*(t-t[0]) - theta[8] 
    ln_likelihood = -0.5*(np.sum((residuals/err)**2 + np.log(2*np.pi*(err)**2)))  
    
    return ln_likelihood

def lnprobability(theta, params, model, t, flux, err):
    p = lnprior(theta)
    if not np.isfinite(p):   
        return -np.inf
    return p + lnlikelihood(theta, params, model, t, flux, err)  


##==================================================##


def prelim(x, y, yerr):
    m = np.isfinite(x) & np.isfinite(y) & np.isfinite(yerr)
    time = x[m]
    magnitude = y[m]
    magError = yerr[m]
    
    magnitude = np.power(10, magnitude/-2.5)
    magError = magnitude*magError/-2.5*np.log(10)
    
    index = np.argsort(time)
    time = time[index]
    magnitude = magnitude[index]
    magError = magError[index]
    
    #limit data to a 4hr window around transit
    t0 = 2457062.57964  ##how to make this auto?
    per = 3.484564
    numtransits = round((time[0]-t0)/per)
    myt0 = t0 + per*numtransits   
    
    ind = 0
    check = np.ones(len(time), dtype = bool)
    for t in time:
        if t<=(myt0-0.083) or t>=(myt0+0.1):
            check[ind] = False
        else:
            check[ind] = True
        ind+=1           
    
    time = time[check]
    magnitude = magnitude[check]
    magError = magError[check]   
    
    #filter outliers here:
    
    std = np.std(magnitude)
    avg = np.mean(magnitude)
    
    ind = 0
    check = np.ones(len(time), dtype = bool)
    for m in magnitude:
        if m>=(avg + (3*std)):
            check[ind] = False
        else: 
            check[ind] = True
        ind+=1
    time = time[check]
    magnitude = magnitude[check]
    magError = magError[check]
    
    std = np.std(magnitude)
    avg = np.mean(magnitude)
    
    #repeat to filter again:
    
    ind = 0
    check = np.ones(len(time), dtype = bool)
    for m in magnitude:
        if m>=(avg + (3*std)):
            check[ind] = False
        else: 
            check[ind] = True
        ind+=1
    time = time[check]
    magnitude = magnitude[check]
    magError = magError[check]
     
    return time, magnitude, magError


##==================================================##

def initialModel(time, magnitude, magError, serial, status):
    #how to get these values?
    per = 3.484564  
    rp_rstar = 0.1184 
    impact = .84  
    a_rstar = 21.9
    inclination = math.acos(impact/a_rstar)*180./math.pi
    #print "inclination", inclination
    
    t0 = 2457062.57964  ##how to make this auto?
    numtransits = round((time[0]-t0)/per)
    myt0 = t0 + per*numtransits   
    #print myt0, per, rp_rstar, a_rstar, inclination
    
    params = batman.TransitParams()
    params.t0 = myt0
    params.per = per
    params.rp = rp_rstar
    params.a = a_rstar
    params.inc = inclination
    ecc, w = calc_eccw(0., np.sqrt(0.4))
    params.ecc = ecc
    params.w = w
    #print "ecc:" , ecc , "  w:" , w
    params.u = [0.1, 0.3]  #example values
    params.limb_dark = "quadratic"
    
    thismodel = batman.TransitModel(params, time)
    flux = thismodel.light_curve(params)
    
    #PLOT
#    plt.figure()
#    plt.errorbar((time-myt0)*24, magnitude, magError, fmt='o', c='k', alpha=.4) 
#    plt.scatter((time-myt0)*24, flux, alpha = .4, zorder = 1)
#    plt.plot((time-myt0)*24, flux, c='r', alpha = .7, linewidth=3, zorder = 10) 
#    plt.xlabel('Time from mid-transit in hours')
#    plt.ylabel('Flux relative to baseline')
#    plt.title('Transit of K2-25b ' + serial)
#    title = serial + '.pdf'
#    plt.savefig(title)
#    plt.show()
#    plt.close()
    
    return rp_rstar, a_rstar, myt0, params, time, magnitude, magError, thismodel

##==================================================##


def emceeModel(rp_rstar, a_rstar, myt0, params, time, magnitude, magError, thismodel, serial, status):
    slope_file = open('slope_file.txt', 'a+')
    megaslope_file = open('megaslope_file.txt', 'a+')
    
    #how to automate these? is it ok to keep these guesses?
    guess_rp = 0.12
    guess_a = a_rstar
    guess_t0 = myt0
    guess_secosw, guess_sesinw  = calc_escw(0.1,90.) 
    guess_u1, guess_u2 = 0.1, 0.3 
    guess_slope = 0.0
    guess_offset = 0.0
    
    theta = [guess_rp, guess_a, guess_t0, guess_secosw, guess_sesinw, guess_u1, guess_u2, guess_slope, guess_offset] 
    
    ndimensions, nwalkers = len(theta), 20
    sampler = emcee.EnsembleSampler(nwalkers, ndimensions, lnprobability, args = (params, thismodel, time, magnitude, magError))
        
    scale = np.ones_like(theta)*1.e-8 
    pos = [theta + scale*np.random.randn(ndimensions) for i in range (nwalkers)]
    
    for p in pos: 
        e = np.random.uniform(low=0., high=0.6)
        w = np.random.uniform(low=50., high=130.)
        ses, sec = calc_escw(e, w)
        p[3]=sec
        p[4]=ses
        
   # print 'sec:', sec, 'ses:', ses    
    sampler.run_mcmc(pos,1000)
    print "done running emcee!"
    
    
    samples = sampler.chain[:, 800:, :].reshape((-1, ndimensions))  
    
    best_rp, best_a, best_t0, best_secosw, best_sesinw, best_u1, best_u2, best_slope, best_offset = \
                                 map(lambda v: 
                                     (v[1], v[2]-v[1], v[1]-v[0]),    ##numbers arranges: [50th, upper error, lower error]
                                     zip(*np.percentile(samples, [16, 50, 84], axis=0))) 
    durations = []     
    per = params.per
    inc = math.acos(.84/best_a[0])  ##sample_a or best_a?
    print 'creating durations list'                           
    for i in samples:
        sample_rp, sample_a, sample_t0, sample_secosw, sample_sesinw, sample_u1, sample_u2, sample_slope, sample_offset = i
        ecc, w = calc_eccw(sample_secosw, sample_sesinw) 
        if ecc>=0 and ecc<1:
            if ((1+sample_rp)**2)>=(sample_a*math.cos(inc)):  ##sincle a_rstar/rstar = a?
                if inc>=-math.pi/2 and inc<=math.pi/2 and inc!=0:
                    p1 = (1+sample_rp)**2-(sample_a*math.cos(inc))
                    p2 = math.sin(inc)
                    p3 = math.sqrt(p1)
                    d1 = p3/p2  ##equation split for debugging purposes; rejoin later?
                    
                    d2 = math.asin((1/sample_a)*(d1))
                    dur = (per/math.pi)*d2*(math.sqrt(1-ecc**2)/(1+sample_sesinw))
                    durations.append(dur)
        
    print np.array(durations).shape
#    print durations
    temp = np.percentile(durations, [16, 50, 84], axis=0)    ##numbers arranged: [lower error, 50th, upper error]
    best_duration = [temp[1], temp[2]-temp[1], temp[1]-temp[0]]
    best_duration = np.array(best_duration)
    print 'best durations', np.array(best_duration).shape, best_duration    
        
    
    depths = []
    for i in samples:
        sample_rp, sample_a, sample_t0, sample_secosw, sample_sesinw, sample_u1, sample_u2, sample_slope, sample_offset = i
        dep = sample_rp**2
        depths.append(dep)
    print np.asarray(depths).shape
    temp = np.percentile(depths, [16, 50, 84], axis=0)
    best_depth = [temp[1], temp[2]-temp[1], temp[1]-temp[0]]
    best_depth = np.array(best_depth)
    
    #print 'theta values: ', best_rp, best_a, best_t0, best_secosw, best_sesinw, best_u1, best_u2, best_slope, best_offset
    
    theta_best =  np.percentile(samples, 50, axis=0)

    tlin = np.linspace(myt0-0.2, myt0+0.2, 1000)
    
    plt.figure()
    for theta in samples[np.random.randint(len(samples), size=100)]:
        lnlikelihood(theta, params, thismodel, time, magnitude, magError)
        #initialize model
        mlin = batman.TransitModel(params, tlin)
        flux = mlin.light_curve(params) + theta[7]*(tlin-time[0]) + theta[8] ## OFFSET
        plt.plot( (tlin), flux, c='b', linewidth=3, alpha=0.05, zorder=3)
         
    lnlikelihood(theta, params, thismodel, time, magnitude, magError)
    mlin = batman.TransitModel(params, tlin)
    flux = mlin.light_curve(params) + theta_best[7]*(tlin-time[0]) + theta_best[8]  ##OFFSET
    
    #------------------------------#
    thing = str(serial) + ';' + str(best_rp) + ';' + str(best_a) + ';' + str(best_t0) + ';' + str(best_secosw) + ';' + str(best_sesinw) + ';' + str(best_u1) + ';' + str(best_u2) + ';' + str(best_slope) + ';' + str(best_offset) + ';' + str(best_duration) + ';' + str(best_depth) + '; \n'
    if 'stacked' in serial:
        if find_key(serial)==False:
            slope_file.write(thing)
    elif 'Megastack' in serial:
        if find_key(serial)==False:
            megaslope_file.write(thing)
    print 'line inserted'
    #------------------------------#
    
#    plt.plot( (tlin), flux, c='b', linewidth=3, alpha=0.05, zorder = 3)
#    plt.scatter((tlin), flux, alpha = .4, zorder = 0)
#    plt.errorbar(time, magnitude, magError, fmt='o', c='k', alpha=0.4, zorder = 0)
#    plt.ylim(0.97,1.02)   #will this be constant across all data sets?
#    plt.xlim(myt0-0.08, myt0+0.08)
#    plt.xlabel('Time in BJD')
#    plt.ylabel('Flux relative to baseline')
#    plt.title('Best fit transit of K2-25b for ' + serial)
#    title = serial + '_best.pdf'
#    plt.savefig(title)
#    plt.close()
#    plt.show()
#    
    return samples 

##==================================================##


def cornerPlot(samples, serial):
    figure = corner.corner(samples, labels = ["rp", "a", "t0", "secosw", "sesinw", 'u1', 'u2', "slope", 'offset'])
    title = serial + '_corner.pdf'
    plt.savefig(title)
    plt.close()
#    plt.show()