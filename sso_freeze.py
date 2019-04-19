#!/usr/bin/env python3

#sso_freeze python utility 
#by Brad Snios - 11/14/18
#
#Compiled to run on python 3+

import math                                     #required for trig functions
import numpy as np                              #define arrays
import operator                                 #subtract arrays
import sys 			  							#import variables from command line
from astropy.io import fits                     #input .fits files
from astropy import wcs                         #convert to wcs coordinates
from astroquery.jplhorizons import Horizons     #automatically download ephemermis 
from astropy.time import Time                   #convert between different time coordinates
from astropy.time import TimeDelta              #add/subtract time intervals 
from scipy.interpolate import interp1d          #interpolate functions

#Define input parameters
input_file = ''
object_name = '' 
output_file = ''
obj_id = '' 
jpl_id_type = ''

#Compiled target list of JPL targets, including required object IDs and target types
#NOTE: This list is a work in progress and will continue to expand
target_list = ['mercury','venus','earth','mars','io','europa','jupiter','saturn','neptune','pluto']
obj_id_list = [199,299,399,499,501,502,599,699,799,899,999]
jpl_id_type_list = ['majorbody','majorbody','majorbody','majorbody','majorbody','majorbody','majorbody','majorbody','majorbody']

#Verify all required input parameters are included
if len(sys.argv) != 4: 
	print('\nPlease check input format: sso_freeze.py <input_file> <object_name> <output_file>\n')
	sys.exit()
else: 
	input_file = sys.argv[1]
	object_name = sys.argv[2]
	output_file = sys.argv[3]

#Find the user-selected target from the list and update JPL parameters
if object_name.lower() in (name for name in target_list):
	obj_id = obj_id_list[target_list.index(object_name.lower())]
	jpl_id_type = jpl_id_type_list[target_list.index(object_name.lower())]
	print('\nObject recognized as "'+target_list[target_list.index(object_name.lower())]+'"')
else: 
	print('\nTarget not found. Please check spelling.\n')
	sys.exit()

if __name__ == "__main__":
    
    #Read in fits file, including relevant header information
    with fits.open(input_file) as fitsfile: 
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
        print('Telescope configuration recognized as "'+instrument+'"\n')
        if instrument == 'ACIS': 
            telescope = '500@-151'
            ctypex, ctypey = fitsfile[1].header['TCTYP11'], fitsfile[1].header['TCTYP12']
            crvalx, crvaly = fitsfile[1].header['TCRVL11'], fitsfile[1].header['TCRVL12']
            crpixx, crpixy = fitsfile[1].header['TCRPX11'], fitsfile[1].header['TCRPX12']
            cdeltx, cdelty = fitsfile[1].header['TCDLT11'], fitsfile[1].header['TCDLT12']
 
        elif instrument == 'HRC': 
            telescope = '500@-151'
            ctypex, ctypey = fitsfile[1].header['TCTYP8'], fitsfile[1].header['TCTYP9']
            crvalx, crvaly = fitsfile[1].header['TCRVL8'], fitsfile[1].header['TCRVL9']
            crpixx, crpixy = fitsfile[1].header['TCRPX8'], fitsfile[1].header['TCRPX9']
            cdeltx, cdelty = fitsfile[1].header['TCDLT8'], fitsfile[1].header['TCDLT9']

        elif (instrument == 'EPN') or (instrument == 'EMOS1') or (instrument == 'EMOS2'):
            telescope = '500@-125989'
            ctypex, ctypey = fitsfile[1].header['REFXCTYP'], fitsfile[1].header['REFYCTYP']
            crvalx, crvaly = fitsfile[1].header['REFXCRVL'], fitsfile[1].header['REFYCRVL']
            crpixx, crpixy = fitsfile[1].header['REFXCRPX'], fitsfile[1].header['REFYCRPX']
            cdeltx, cdelty = fitsfile[1].header['REFXCDLT'], fitsfile[1].header['REFYCDLT']

        else: 
            print('Telescope configuration not recognized. Please verify information is included in .fits file header\n')
            sys.exit()

        #Input header keywords into wcs transformation matrix
        w.wcs.ctype = [ctypex, ctypey]
        w.wcs.crval = [crvalx, crvaly]
        w.wcs.crpix = [crpixx, crpixy]
        w.wcs.cdelt = [cdeltx, cdelty]

        #Define start and stop times for ephemeris data; since jpl does not accept seconds, 
        #all times are in YY:MM:DD hh:mm format;dt is added to stop time to ensure ephemeris 
        #data range extends beyond exposure time 
        eph_tstart = Time(tstart, out_subfmt='date_hm')
        dt = TimeDelta(0.125, format='jd') 
        eph_tstop = Time(tstop + dt, out_subfmt='date_hm')
        obj = Horizons(id=obj_id,location=telescope,epochs={'start':eph_tstart.iso, 'stop':eph_tstop.iso, 'step':'5m'}, id_type=jpl_id_type)
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
        print('Starting R.A., Decl. for target: ', ra0, dec0,'\n')

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
        #eventlist['x'], eventlist['y'] = w.wcs_world2pix(cnt_ra,cnt_dec, 1)
        ocx, ocy = w.wcs_world2pix(cnt_ra,cnt_dec, 1)
        ocx_col = fits.Column(name='ocx', format='1E', unit='pixel', coord_type=ctypex, coord_unit='deg', coord_ref_point=crpixx, coord_ref_value=crvalx, coord_inc=cdeltx, array=ocx)
        ocy_col = fits.Column(name='ocy', format='1E', unit='pixel', coord_type=ctypey, coord_unit='deg', coord_ref_point=crpixy, coord_ref_value=crvaly, coord_inc=cdelty, array=ocy)

        #update fits table information with new object-centered coordinate columns
        table = fitsfile["EVENTS"]
        newtable = fits.BinTableHDU.from_columns(table.columns + fits.ColDefs([ocx_col]) + fits.ColDefs([ocy_col]))
        fitsfile["EVENTS"].data = newtable.data
        fitsfile["EVENTS"].header.update(newtable.header)

        #output fits file
        fitsfile.writeto(output_file+'.fits')

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
        region_file = open(output_file+'.reg','w') 
        region_file.write('# Region file format: DS9 version 4.1\nglobal color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\nfk5\n')
        region_file.write('circle('+str(ra0)+','+str(dec0)+','+str(ang_radius)+'")\n')
        region_file.write('line('+str(ra0)+','+str(dec0)+','+str(ra_pole)+','+str(dec_pole)+')')
        #If the North is on the opposite side of the object, the North Pole marker is highlighted magenta
        #rather than green 
        if NPole_dist < 0: 
            region_file.write('# line=0 0 color=magenta')
        region_file.close()
        

