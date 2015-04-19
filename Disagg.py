# -*- coding: utf-8 -*-
"""
Created on Sun Apr 12 12:11:46 2015

@author: admin
"""

import numpy as np
#import matplotlib.pyplot as plt
#from scipy.optimize import minimize
from InterpolateRadii import get_radius
from DistanceBearing import get_r_theta
from Holland2010 import fitHolland2010, Holland2010
import csv
import os
import datetime

print(datetime.datetime.now())

#  Read each line of HurPlatLatLon.csv
plat_file  = open('HurPlatLatLon.csv', "rb")
plat_reader = csv.reader(plat_file)


out_file = open('Kara.csv', "wb")
out_file.write('"ID",       "PLATFORM", "LAT",    "LON",   "Date", "Time", "Platform Distance", "WindSpeed"\n')

max_deg = 765
# Vector of angles for which wind radii are defined
# Extended to nine to get smooth periodic contour
x=[45,135,225,315,405,495,585,675,765] 

for plat_row in plat_reader:
    hurr_id = plat_row[0]
    rig_id = plat_row[1]
    dest_lat = float(plat_row[2])
    dest_lon = float(plat_row[3])

    #print "---------------- Hurricane id:", hurr_id
    # Search for each hurricane id in HurPlatTest.csv and save the records in a temp file
    r_c = os.system("grep {id} hurdat4.csv > infile.csv".format(id=hurr_id))
    #print "grep return_code:{}",r_c
    in_file  = open('infile.csv', "rb")
    in_reader = csv.reader(in_file)

    # Read each line of hurdat corresponding to our hurricane id
    for row in in_reader:
        date = float(row[4])
        time = float(row[5])
        origin_lat = float(row[8])
        origin_lon = float(row[9])
        v_max = float(row[10])
        min_pressure = float(row[11])
        wind_rad34 = [float(i) for i in row[12:16]]
        wind_rad50 = [float(i) for i in row[16:20]]
        wind_rad64 = [float(i) for i in row[20:24]]
        #print("origin_lat={} \norigin_lon={}".format(origin_lat,origin_lon))
        #print("dest_lat={} \ndest_lon={}".format(dest_lat,dest_lon))
        #print "wr34:" ,(wind_rad34)
        plat_dist,theta = get_r_theta((origin_lat, origin_lon), (dest_lat, dest_lon))
    
        r = wind_rad34
        extended_r = [r[0], r[1],r[2],r[3],r[0],r[1],r[2],r[3],r[0]]
    
        # Interlopate radii quadrant data points for given theta
        rad34 = get_radius(theta, extended_r)
    
        r = wind_rad50
        extended_r = [r[0], r[1],r[2],r[3],r[0],r[1],r[2],r[3],r[0]]
    
        # Interlopate radii quadrant data points for given theta
        rad50 = get_radius(theta, extended_r)
    
        r = wind_rad64
        extended_r = [r[0], r[1],r[2],r[3],r[0],r[1],r[2],r[3],r[0]]
    
        # Interlopate radii quadrant data points for given theta
        rad64 = get_radius(theta, extended_r)
       #print "     rig_dist:%d, theta:%d \n"% (rig_dist,theta)
        #print "wind_rad34:" , wind_rad34
        #print "rad34:" , rad34
    
        # Max speed is too low to fit to model
        if(v_max < 40):
            plat_wind_speed = -1
        
        # No speed data
        elif(np.min(wind_rad34)<0 or np.sum(np.abs(wind_rad34))==0):
            plat_wind_speed = -1
        
        # All 50kt radii are zero: two data points
        elif(np.sum(np.abs(wind_rad50))==0):
            case = 1
    
        # All 64kt radii are zero: three data points    
        elif(np.sum(np.abs(wind_rad64))==0):
            case = 2
    
        #  fit Holland equation for a, b, r_max: Four data points
        else:
            case = 3
    
        r_max0 = 8.
        a0 = 0.01
        b0 = 1.
        case = 3    
        # solve Holland equation for b, r_max
        a, b, r_max = fitHolland2010(a0, b0, r_max0, rad34, rad50, rad64, v_max, case)
        plat_wind_speed = -1
        #if(Success):
        plat_wind_speed = Holland2010(plat_dist, a, b, r_max, v_max)
        if(plat_wind_speed<0 or plat_wind_speed > v_max):
            plat_wind_speed = 0.
        # Create a record to write to the output file
        out_row = '"%s", "%s", %.4f, %.4f, %d, %4d, %.1f, %.2f\n' % \
          (hurr_id, rig_id, dest_lat, dest_lon, date, time, plat_dist, plat_wind_speed)
        out_file.write(out_row)

    in_file.close()

plat_file.close()
out_file.close()
print(datetime.datetime.now())

