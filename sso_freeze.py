#!/usr/bin/env python3

#sso_freeze python utility 
#by Brad Snios - 11/14/18
#
#Compiled to run on python 3+

import math                                     #required for trig functions
import numpy as np                              #define arrays
import operator                                 #subtract arrays
from astropy.io import fits                     #input .fits files
from astropy import wcs                         #convert to wcs coordinates
from astroquery.jplhorizons import Horizons     #automatically download ephemermis 
from astropy.time import Time                   #convert between different time coordinates
from astropy.time import TimeDelta              #add/subtract time intervals 
from scipy.interpolate import interp1d          #interpolate functions

if __name__ == "__main__":
    
    #Read inputs from separate file
    inputs = {}
    input_file = open("inputs.par")
    input_data = input_file.readlines()
    for line in input_data: 
        #Parse input, assign values to variables
        if not line.startswith('#'): 
            key, value = line.split(":")
            inputs[key.strip()] = value.strip()
    input_file.close()

    #Read in fits file, including relevant header information
    with fits.open(inputs['filename']) as fitsfile: 
        header = fitsfile[1].header
        eventlist = fitsfile[1].data
        w = wcs.WCS(fitsfile[1].header)
        
        #Read in observation start/stop times and convert to time format
        temp1 = fitsfile[1].header['TSTART']
        temp2 = fitsfile[1].header['TSTOP']
        tstart = Time(temp1, format='cxcsec')
        tstop = Time(temp2, format='cxcsec') 
        
        #Read in instrument used for observation
        instrument = fitsfile[1].header['INSTRUME']

        #Define telescope_id required for JPL-Horizons and wcs transformation keywords
        #required for wcs-to-pixel coordinate transformations, and vice versa 
        if instrument == 'ACIS': 
            telescope = '500@-151'
            w.wcs.ctype = [fitsfile[1].header['TCTYP11'], fitsfile[1].header['TCTYP12']]
            w.wcs.crval = [fitsfile[1].header['TCRVL11'], fitsfile[1].header['TCRVL12']]
            w.wcs.crpix = [fitsfile[1].header['TCRPX11'], fitsfile[1].header['TCRPX12']]
            w.wcs.cdelt = [fitsfile[1].header['TCDLT11'], fitsfile[1].header['TCDLT12']]
 
        if instrument == 'HRC': 
            telescope = '500@-151'
            w.wcs.ctype = [fitsfile[1].header['TCTYP8'], fitsfile[1].header['TCTYP9']]
            w.wcs.crval = [fitsfile[1].header['TCRVL8'], fitsfile[1].header['TCRVL9']]
            w.wcs.crpix = [fitsfile[1].header['TCRPX8'], fitsfile[1].header['TCRPX9']]
            w.wcs.cdelt = [fitsfile[1].header['TCDLT8'], fitsfile[1].header['TCDLT9']]

        if (instrument == 'EPN') or (instrument == 'MOS1') or (instrument == 'MOS2'):
            telescope = '500@-125989'
            w.wcs.ctype = [fitsfile[1].header['REFXCTYP'], fitsfile[1].header['REFYCTYP']]
            w.wcs.crval = [fitsfile[1].header['REFXCRVL'], fitsfile[1].header['REFYCRVL']]
            w.wcs.crpix = [fitsfile[1].header['REFXCRPX'], fitsfile[1].header['REFYCRPX']]
            w.wcs.cdelt = [fitsfile[1].header['REFXCDLT'], fitsfile[1].header['REFYCDLT']]

        #Define start and stop times for ephemeris data; since jpl does not accept seconds, 
        #all times are in YY:MM:DD hh:mm format;dt is added to stop time to ensure ephemeris 
        #data range extends beyond exposure time 
        eph_tstart = Time(tstart, out_subfmt='date_hm')
        dt = TimeDelta(0.125, format='jd') 
        eph_tstop = Time(tstop + dt, out_subfmt='date_hm')
        obj = Horizons(id=inputs['obj_id'],location=telescope,epochs={'start':eph_tstart.iso, 'stop':eph_tstop.iso, 'step':'5m'}, id_type=inputs['jpl_id_type'])
        eph = obj.ephemerides()
        
        #Create interpolation function for RA and DEC based on ephemeris data
        ra_int = interp1d(eph['datetime_jd'], eph['RA'], kind='cubic')
        dec_int = interp1d(eph['datetime_jd'], eph['DEC'], kind='cubic')
        #Define start time as time stamp of first recorded count during exposure
        t0 = Time(eventlist['time'][0], format='cxcsec')
        #Use starting time to define initial RA and DEC of the solar system object
        ra0 = ra_int(t0.jd)
        dec0 = dec_int(t0.jd)

        #Output starting time, RA, and DEC
        print('Start time: ', t0.iso) 
        print('Starting position for target: ', ra0, dec0)

        #Convert event list x,y coordinates to RA,DEC is pix2world transformation
        cnt_ra, cnt_dec = w.wcs_pix2world(eventlist['x'], eventlist['y'], 1)

        #Define event list time coordinate as array, then use to calculate
        #RA,DEC for each time stamp
        ttemp = Time(eventlist['time'], format='cxcsec')
        ra_temp = ra_int(ttemp.jd)
        dec_temp = dec_int(ttemp.jd) 

        #Calculate cnt_ra - (ra_temp - ra0) to determine coordinate shift for each count;
        #repeat process for DEC coordinates
        cnt_ra = list(map(operator.sub, cnt_ra, ra_temp))
        cnt_dec = list(map(operator.sub, cnt_dec, dec_temp))
        cnt_ra[:] = [x + ra0 for x in cnt_ra]
        cnt_dec[:] = [x + dec0 for x in cnt_dec]

        #Convert corrected RA,DEC to pixel coordinates and output to fits file
        eventlist['x'], eventlist['y'] = w.wcs_world2pix(cnt_ra,cnt_dec, 1)

        #Calculate angular radius, North pole angle, and North pole distance 
        #of solar system object based on JPL data
        ang_radius = 0.5*sum(eph['ang_width'])/len(eph['ang_width'])
        NPole_ang = sum(eph['NPole_ang'])/len(eph['NPole_ang'])
        NPole_dist = sum(eph['NPole_dist'])/len(eph['NPole_dist'])
        
        #Convert North Pole location to RA/DEC coordinates
        ra_pole = ra0 + (1/3600)*abs(NPole_dist)*math.sin(NPole_ang*math.pi/180)
        dec_pole = dec0 + (1/3600)*abs(NPole_dist)*math.cos(NPole_ang*math.pi/180)

        #Write region file that contains location of solar system object, average size
        #of object, and north pole location
        region_file = open(inputs['object_name']+'.reg','w') 
        region_file.write('# Region file format: DS9 version 4.1\nglobal color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\nfk5\n')
        region_file.write('circle('+str(ra0)+','+str(dec0)+','+str(ang_radius)+'")\n')
        region_file.write('line('+str(ra0)+','+str(dec0)+','+str(ra_pole)+','+str(dec_pole)+')')
        #If the North is on the opposite side of the object, the North Pole marker is highlighted magenta
        #rather than green 
        if NPole_dist < 0: 
            region_file.write('# line=0 0 color=magenta')
        region_file.close()
        
        #Read out fits file
        fitsfile.writeto(inputs['object_name']+'.fits') 
