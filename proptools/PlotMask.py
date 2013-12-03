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



class MaskPlot():
    '''
    this is the main plotting class for the slit mask design
    '''
    
    def __init__(self, cra=0.0, cdec=0.0):
        '''
        initialise the plotting window and put the labels in place
        '''
                
        
        
        pl.ion()
        self.fig = pl.figure(figsize=(14,10))
        self.ax = self.fig.add_subplot(111)    
        pl.axis('equal')
        pl.xlabel('Right Ascension [deg]')
        pl.ylabel('Declination [deg]')
   
        pl.show(block=False)
        '''
        setup the plotting window with the RSS FoV and chip locations
        '''
    def draw_CCD_FoV(self, cra, cdec):
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
        rect = Rectangle(chip1_corner, chip_width, chip_height, color='none', ec='y', zorder=0, gid='chip1')
        self.ax.add_patch(rect)
        
        rect = Rectangle(chip2_corner, chip_width, chip_height, color='none', ec='y', zorder=0, gid='chip2')
        self.ax.add_patch(rect)
        
        rect = Rectangle(chip3_corner, chip_width, chip_height, color='none', ec='y', zorder=0, gid='chip3')
        self.ax.add_patch(rect)   
        
        circ = Circle((cra, cdec), dcr, color='none', ec='r', zorder=0, gid='FoV')
        self.ax.add_patch(circ)
        
        self.ax.set_xlim(cra + 0.121868, cra - 0.121868)
        self.ax.set_ylim(cdec - 0.07312, cdec + 0.07312)
        
        pl.draw()
        
    def plot_slitlets(self, x):
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
            self.ax.plot(x[i]['targ_ra'], x[i]['targ_dec'], '*', color='k', ms=7, gid='refstar')
            rect = Rectangle((spec_x0[i], spec_y0[i] - 2.5/3600), spec_width[i], 5./3600, color='g', ec='g', alpha=0.5, gid='refstar')
            self.ax.add_patch(rect)
            
        # all the objects that lie outside the FoV
        for i in out_FoV_ids:
            self.ax.plot(x[i]['targ_ra'], x[i]['targ_dec'], 'o', color='r', gid='out_point')
            
        # all the objects in the FoV that have no collisions
        for i in in_FoV_ids:
            self.ax.plot(x[i]['targ_ra'], x[i]['targ_dec'], '.', color='w', ms=3, gid='point')
            rect = Rectangle((slit_x0[i], slit_y0[i]), slit_width[i], slit_length[i], color='k', ec='k', alpha=0.7, gid='in_slit')
            self.ax.add_patch(rect)     
            rect = Rectangle((spec_x0[i], spec_y0[i]), spec_width[i], spec_length[i], color='b', ec='b', alpha=0.3, gid='in_spec')
            self.ax.add_patch(rect)  
            
        # all the objsects in the FoV that have collisions    
        for i in cols_ids:
            self.ax.plot(x[i]['targ_ra'], x[i]['targ_dec'], '.', color='w', ms=3, gid='point')
            rect = Rectangle((slit_x0[i], slit_y0[i]), slit_width[i], slit_length[i], color='k', ec='k', alpha=0.7, gid='col_slit')
            self.ax.add_patch(rect)     
            rect = Rectangle((spec_x0[i], spec_y0[i]), spec_width[i], spec_length[i], color='r', ec='r', alpha=0.3, gid='col_spec')
            self.ax.add_patch(rect)  
        
        pl.draw()
    
    def get_guide_stars(self, cra, cdec):
        '''
        Call the guide star script written by Keith Smith to obtain the 
        positions of candidate guide stars
        '''
        # divide the cra by 15 so it's in decimal hours and not degrees
        ra = cra / 15.
        dec = cdec
        
        n_stars, data, header = guide_stars.findGuideStars(ra, dec, '', 'rss', 20.)
        
        if n_stars > 0:
            print 'guide stars: ', n_stars
            gstar_RA = np.array([self.sex2RA(data[i,2]) for i in range(0,len(data))])
            gstar_Dec =  np.array([self.sex2Dec(data[i,3]) for i in range(0,len(data))]) 
            
            return gstar_RA, gstar_Dec
            
        else:
            return np.array([999.99]), np.array([999.99])
    
    def plot_guide_stars(self, cra, cdec):
        '''
        Plot the guide star postitions
        '''
        gstar_ra, gstar_dec = self.get_guide_stars(cra, cdec)
        self.ax.plot(gstar_ra, gstar_dec, 'c*', ms=15, gid='gstars')
        
        pl.draw()

    def sex2RA(self, x):
        x = x.split()
        if float(x[0])>=0:
            return (float(x[0])+float(x[1])/60.0+float(x[2])/3600.) * 15.
        else:
            return -(abs(float(x[0]))+float(x[1])/60.0+float(x[2])/3600.0) * 15.
        
    def sex2Dec(self, x):
        x = x.split()
        if float(x[0])>=0:
            return float(x[0])+float(x[1])/60.0+float(x[2])/3600.0
        else:
            return -(abs(float(x[0]))+float(x[1])/60.0+float(x[2])/3600.0)
        
         

    def clear_plot(self, cra, cdec):
        self.ax.cla()
        
    def clear_all(self):
        '''
        remove everything from the plot
        '''
        # remove the slits from the lines artist
        for i in plot.ax.get_lines():
            if (i.get_gid() == 'out_point') or \
                    (i.get_gid() == 'point') or \
                    (i.get_gid() == 'refstar') or \
                    (i.get_gid() == 'gstars'):
                i.remove()  

        # remove the slits from the children artist where the patches are        
        for i in plot.ax.get_children():
            if (i.get_gid() == 'in_slit') or \
                    (i.get_gid() == 'in_spec') or \
                    (i.get_gid() == 'col_slit') or \
                    (i.get_gid() == 'col_spec') or \
                    (i.get_gid() == 'refstar'):
                i.remove()
                
        pl.draw()                    
    
    
    
if __name__ == "__main__":

    x = np.load('slitlets.npy')
    cra = 61.951306
    cdec = -12.193412
#    gstar_RA, gstar_Dec = get_guide_stars(cra, cdec)    
#    ax = init_plot(cra, cdec)
    plot = MaskPlot(cra, cdec)
    plot.draw_CCD_FoV(cra, cdec)
    plot.plot_slitlets(x)
    plot.plot_guide_stars(cra, cdec)
    time.sleep(5)
#    plot.clear_plot(cra, cdec)
    plot.clear_all()
    time.sleep(3)
    plot.plot_slitlets(x)
#    draw_CCD_FoV(ax, cra, cdec)
#    plot_slitlets(ax, x)
#    glines = plot_guide_stars(ax, gstar_RA, gstar_Dec)
#    print glines
#    pl.show()

    

    
#    pl.draw()
#    time.sleep(5)
#    clear_plot(ax)
#    time.sleep(3)
#    init_plot(cra, cdec)
#    test(ax, cra, cdec)
#    pl.draw()
#    pl.show()
