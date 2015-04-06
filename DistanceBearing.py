#!/usr/bin/env python
#
# Surface distance and bearing using Haversine formula
# Kamran Ossia 3/14/2015
# Based on http://www.movable-type.co.uk/scripts/latlong.html
#
# Inputs: 
#    Origin (lat1, lon1), destination (lat2,lon2) in degrees, 
#    Boolean: if bearing is needed in radians or degrees
#
# Outputs:
#    Great Circle Distance: nautical miles
#    Bearing: radians (default) or degrees, CCW from due east
# Usage: 
#    import DistanceBearing
#    r,theta = DistanceBearing.get((lat1, lon1), (lat2, lon2))
#
#  to have bearing returned in degrees:   
#
#    r,theta = dist_bear.d_b((lat1, lon1), (lat2, lon2), True)
#

import math
 
def get(origin, destination, giveDegrees=False):
    lat1, lon1 = origin
    lat2, lon2 = destination
    earth_radius = 3440.0648  # nm or 6371.0 km 

    r_lat1 = math.radians(lat1)
    r_lat2 = math.radians(lat2)
    r_lon1 = math.radians(lon1)
    r_lon2 = math.radians(lon2)

    d_lon = math.radians(lon2-lon1)
    d_lat = math.radians(lat2-lat1)

    a = math.sin(d_lat/2.)**2 + math.cos(r_lat1)*math.cos(r_lat2)*math.sin(d_lon/2.)**2

    distance = 2.0 * earth_radius * math.atan2(math.sqrt(a), math.sqrt(1-a))

    b = math.sin(d_lon)*math.cos(r_lat2)
    c = math.cos(r_lat1)*math.sin(r_lat2)-math.sin(r_lat1)*math.cos(r_lat2)*math.cos(d_lon)

    bearing = math.atan2(c,b)

    if(giveDegrees): bearing *= 180.0 / math.pi
    return(distance, bearing)
