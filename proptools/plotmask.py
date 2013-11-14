# -*- coding: utf-8 -*-
"""
Created on Wed May  8 23:19:57 2013

@author: jpk
"""

import numpy as np
import pylab as pl
import pywcs
from matplotlib.patches import Rectangle,Circle
import sys,os

import PySpectrograph
from PySpectrograph.Models import RSSModel
from PySpectrograph.Spectra import Spectrum

import logging # only needed for debugging
import time
import guide_stars

def init_plot(cra, cdec):
    '''
    initialise the plotting window
    '''
    
    fig = pl.figure(figsize=(14,10))
    ax = fig.add_subplot(111)    
    pl.axis('equal')
    pl.xlabel('RSS x (unbinned pix)')
    pl.ylabel('RSS y (unbinned pix)')
    
    ax.set_xlim(cra + 0.121868, cra - 0.121868)
    ax.set_ylim(cdec - 0.07312, cdec + 0.07312)
    
    return ax
    

def test(ax, cra, cdec):
    '''
    setup the plotting window with the RSS FoV and chip locations
    '''
    
    pixscale = 0.12535
    
    # chip geometry in degrees
    chip_width = 2034 * pixscale / 3600.
    chip_height = 4102 * pixscale / 3600.
    chip_gap = 70 * pixscale / 3600.
    dcr = 4. / 60. # RSS FoV radius in degrees 
    
    # define the lower left corners of each chip
    chip1_corner = (cra + chip_gap + (chip_width / 2.), cdec - (chip_height / 2.))    
    chip2_corner = (cra - (chip_width / 2.), cdec - (chip_height / 2.))
    chip3_corner = (cra - chip_width - chip_gap - (chip_width / 2.), cdec - (chip_height / 2.)) 
     
    # draw the chips all relative to the central mask coordinates
    rect = Rectangle(chip1_corner, chip_width, chip_height, color='none', ec='y')
    ax.add_patch(rect)
    
    rect = Rectangle(chip2_corner, chip_width, chip_height, color='none',ec='y')
    ax.add_patch(rect)
    
    rect = Rectangle(chip3_corner, chip_width, chip_height, color='none',ec='y')
    ax.add_patch(rect)   
    
    circ = Circle((cra, cdec), dcr, color='none', ec='r')
    ax.add_patch(circ)
    
    
           
def plot_slitlets(ax, x):
    '''
    bla bla.... update when I'm done

    '''
    
    pixscale = 0.12535    
    
    in_FoV_ids = np.where((x['priority'] > 0) & (x['fov_flag'] == 1) & (x['collision_flag'] != 1))[0]
    out_FoV_ids = np.where((x['priority'] > 0) & (x['fov_flag'] == 0))[0]
    refstars_ids = np.where(x['refstar_flag'] == 1)[0]
    cols_ids = np.where(x['collision_flag'] ==1)[0]    
    

    slit_x0 = x['targ_ra'] - x['width'] / 2. / 3600.
    slit_y0 = x['targ_dec'] - x['len1'] / 3600.
    slit_width = x['width'] / 3600.
    slit_length = (x['len1'] + x['len2']) / 3600.
    
    spec_x0 = x['targ_ra'] - (3000 / 2.) *pixscale / 3600.
    spec_y0 =  x['targ_dec'] - x['len1'] / 3600.
    spec_width = np.array([3000 * pixscale / 3600.] * len(spec_x0))
    spec_length = slit_length
    
    # plot the reference stars first
    for i in refstars_ids:
        ax.plot(x[i]['targ_ra'], x[i]['targ_dec'], '*', color='k', ms=7)
        rect = Rectangle((spec_x0[i], spec_y0[i] - 2.5/3600), spec_width[i], 5./3600, color='g', ec='g', alpha=0.5)
        ax.add_patch(rect)
        
    # all the objects that lie outside the FoV
    for i in out_FoV_ids:
        ax.plot(x[i]['targ_ra'], x[i]['targ_dec'], 'o', color='r')
        
    # all the objects in the FoV that have to collisions
    for i in in_FoV_ids:
        ax.plot(x[i]['targ_ra'], x[i]['targ_dec'], '.', color='w', ms=3)
        rect = Rectangle((slit_x0[i], slit_y0[i]), slit_width[i], slit_length[i], color='k', ec='k', alpha=0.7)
        ax.add_patch(rect)     
        rect = Rectangle((spec_x0[i], spec_y0[i]), spec_width[i], spec_length[i], color='b', ec='b', alpha=0.3)
        ax.add_patch(rect)  
        
    # all the objsects in the FoV that have not collisions    
    for i in cols_ids:
        ax.plot(x[i]['targ_ra'], x[i]['targ_dec'], '.', color='w', ms=3)
        rect = Rectangle((slit_x0[i], slit_y0[i]), slit_width[i], slit_length[i], color='k', ec='k', alpha=0.7)
        ax.add_patch(rect)     
        rect = Rectangle((spec_x0[i], spec_y0[i]), spec_width[i], spec_length[i], color='r', ec='r', alpha=0.3)
        ax.add_patch(rect)  
        
#    ax.draw()
def get_guide_stars(cra, dec):
    # divide the cra by 15 so it's in decimal hours and not degrees
    ra = cra / 15.
    dec = cdec
    
    n_stars, data, header = guide_stars.findGuideStars(ra, dec, '', 'rss', targetRadius=3.)
    
    return n_stars, data, header
    
    
def plot_guide_stars(ax, cra, dec):
    
    n_stars, data, header = guide_stars.findGuideStars(cra, cdec, '', 'rss', targetRadius=3.)
    
    
    print n_stars, data, header

def clear_plot(ax):
    ax.cla()
    
if __name__ == "__main__":
    
    x = np.load('slitlets.npy')
    cra = 61.951306
    cdec = -12.193412
    pl.ion()
    
    ax = init_plot(cra, cdec)
    n_stars, gstar_data, gstar_header = get_guide_stars(61.951306, -12.193412)
#    drawCCDs(ax)
#    drawRSS_FoV(ax)
    test(ax, cra, cdec)
#    plot_refstars(ax, x)
    plot_slitlets(ax, x)
#    pl.draw()
#    time.sleep(5)
#    clear_plot(ax)
#    time.sleep(3)
#    init_plot(cra, cdec)
#    test(ax, cra, cdec)
    pl.draw()
    
