# -*- coding: utf-8 -*-
"""
Kamran Ossia

splrep(x, y, weights=None, xbegin=None, xend=None, order_k=3, task=0,
       smooth=None, t=None, full_output=0, periodic=0, quiet=1)

"""

from scipy.interpolate import splrep,splev

x  = [45,135,225,315,405,495,585,675,765] 
x_begin = 45
x_end   = 765

def get_radius(theta, quad_radii):
    #theta %= 360
    t = splrep(x, quad_radii, xb=x_begin, xe=x_end,k=5,s=10,per=1)
    r = splev(theta,t)
    return max(r,0)