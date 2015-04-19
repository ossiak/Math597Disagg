# -*- coding: utf-8 -*-

import numpy as np
from scipy.optimize import root

"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
 Holland 2012 formula for wind speed
 We use the log of speed for simpler calculation
"""
def Holland2010(radius, a, b, r_max, v_max):
    # The x parameter is constant near the eye, 
    #   then linear function of r, not a function of direction
    if(radius==0):
        return 0

    radius = float(radius)    
    if(radius<r_max):
        x = 0.5
    else:
        x = (0.5+a*(radius-r_max))

    t = np.power(r_max/radius,b)
    wind_speed = np.log(v_max) + x*(np.log(t) + 1-t)

    #print t, radius, r_max, b, a, v_max,np.exp(wind_speed)
    return np.exp(wind_speed)

"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
 Holland 2012 root finding function
"""
def HollandRoot1(x, a, r_max, wr34, v_max):
    b  = x[0]
    
    f = Holland2010(wr34, a, b, r_max, v_max)- 34.0
   
    return f

"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
 Holland 2012 root finding function
"""
def HollandRoot2(x, a, wr34, wr50, v_max):
    b  = x[0]
    r_max = x[1]
    
    f = np.array([Holland2010(wr34, a, b, r_max, v_max)- 34.0, 
         Holland2010(wr50, a, b, r_max, v_max)- 50.0])
   
    return f
    
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
 Holland 2012 root finding function
"""
def HollandRoot3(x,wr34, wr50, wr64, v_max):
    a  = x[0]
    b  = x[1]
    r_max = x[2]
    
    f = np.array([Holland2010(wr34, a, b, r_max, v_max)- 34.0, 
         Holland2010(wr50, a, b, r_max, v_max)- 50.0,
         Holland2010(wr64, a, b, r_max, v_max)- 64.0])
    
    return f

"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
 Fit Holland 2012 formula 
"""
def fitHolland2010(a0, b0, r_max0, wr34, wr50, wr64, v_max, case):
    if case==1:
        a = a0
        x0 = b0
        r_max = r_max0
        sol = root(HollandRoot1, x0, args=(a0, r_max0, wr34, v_max))
        b = sol.x[0]
    elif case==2:
        a = a0
        x0 = [b0, r_max0]
        sol = root(HollandRoot2, x0, args=(a0, wr34, wr50, v_max))
        b = sol.x[0]
        r_max = sol.x[1]
    elif case==3:
        x0 = [a0, b0, r_max0]
        sol = root(HollandRoot3, x0, args=(wr34, wr50, wr64, v_max))
        a = sol.x[0]
        b = sol.x[1]
        r_max = sol.x[2]
    #print sol.success
    return a, b, r_max
