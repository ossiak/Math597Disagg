"""
Code for Disaggregation of Wind Speed Usinf NHC Wind Radii
Kamran Ossia
MATH 597
March 2015

"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize

# Hard-coded numbers for now
# Hurricane Gordon AL082012 20120819 0600
v_max = 95
wr34 = np.array([80, 110, 110, 80])
wr50 = np.array([40,  60,  60, 40])
wr64 = np.array([20,  30,  30, 20])
minrm = np.min(wr64)/4.

# Initial guess
# A rough guess for r_max
x0 = [0, 64*wr64[1]/v_max]   # a, r_max

#theta = np.array(range(16))*np.pi/8.0

# b is the shape parameter, varies with direction
x0.extend(1.75*np.ones(16)) # b(theta)

# v_max also varies with direction
x0.extend(np.ones(16)*v_max) # v_max(theta)

# Initialize the radius data, one for each speed
# The NHC provides speed thresholds for 34,50,64 kt
# We divide the circle into 16 segments
r34 = np.array(np.zeros(16))
r50 = np.array(np.zeros(16))
r64 = np.array(np.zeros(16))

#====================================================================
# Holland 2012 formula for wind speed
# W use the log of speed for simpler calculation
def Holland2010(radius, r_max, b, a, v_max):
    # The x parameter is constant near the eye, 
    #   then linear function of r, not a function of direction
    if(radius<r_max):
        x = 0.5
    else:
        x = (0.5+a*(radius-r_max))

    t = pow(r_max/radius,b)
    wind_speed = np.log(v_max) + x*(np.log(t) + 1-t)
    return np.exp(wind_speed)
    
#=====================================================================
def Holland2010Cost(x):
    a  = x[0]
    rm = x[1]
    b  = np.array(np.zeros(16))
    vm = np.array(np.zeros(16))
    for i in range(16): b[i]  = x[i+2]
    for i in range(16): vm[i] = x[i+18] 

    # 16 points on the unit circle
    #theta = np.array(range(16))*np.pi/8.0

    # Fill in the radii and interpolate
    r34[2]  = wr34[0]
    r34[6]  = wr34[1]
    r34[10] = wr34[2]
    r34[14] = wr34[3]
    r34[0]  = (r34[14]+r34[2])/2.
    r34[4]  = (r34[2]+r34[6])/2.
    r34[8]  = (r34[6]+r34[10])/2.
    r34[12] = (r34[10]+r34[14])/2.
    r34[1]  = (r34[0]+r34[2])/2.
    r34[3]  = (r34[2]+r34[4])/2.
    r34[5]  = (r34[4]+r34[6])/2.
    r34[7]  = (r34[6]+r34[8])/2.
    r34[9]  = (r34[8]+r34[10])/2.
    r34[11] = (r34[10]+r34[12])/2.
    r34[13] = (r34[12]+r34[14])/2.
    r34[15] = (r34[14]+r34[0])/2.

    #ax1 = plt.subplot(111, polar=True)
    #ax1.plot(theta, r34, color='r', linewidth=2)
    #ax1.grid(True)
    
    r50[2]  = wr50[0]
    r50[6]  = wr50[1]
    r50[10] = wr50[2]
    r50[14] = wr50[3]
    r50[0]  = (r50[14]+r50[2])/2.
    r50[4]  = (r50[2]+r50[6])/2.
    r50[8]  = (r50[6]+r50[10])/2.
    r50[12] = (r50[10]+r50[14])/2.
    r50[1]  = (r50[0]+r50[2])/2.
    r50[3]  = (r50[2]+r50[4])/2.
    r50[5]  = (r50[4]+r50[6])/2.
    r50[7]  = (r50[6]+r50[8])/2.
    r50[9]  = (r50[8]+r50[10])/2.
    r50[11] = (r50[10]+r50[12])/2.
    r50[13] = (r50[12]+r50[14])/2.
    r50[15] = (r50[14]+r50[0])/2.

    #ax2 = plt.subplot(111, polar=True)
    #ax2.plot(theta, r50, color='g', linewidth=2)
    #ax2.grid(True)

    r64[2]  = wr64[0]
    r64[6]  = wr64[1]
    r64[10] = wr64[2]
    r64[14] = wr64[3]
    r64[0]  = (r64[14]+r64[2])/2.
    r64[4]  = (r64[2]+r64[6])/2.
    r64[8]  = (r64[6]+r64[10])/2.
    r64[12] = (r64[10]+r64[14])/2.
    r64[1]  = (r64[0]+r64[2])/2.
    r64[3]  = (r64[2]+r64[4])/2.
    r64[5]  = (r64[4]+r64[6])/2.
    r64[7]  = (r64[6]+r64[8])/2.
    r64[9]  = (r64[8]+r64[10])/2.
    r64[11] = (r64[10]+r64[12])/2.
    r64[13] = (r64[12]+r64[14])/2.
    r64[15] = (r64[14]+r64[0])/2.

    #ax3 = plt.subplot(111, polar=True)
    #ax3.plot(theta, r64, color='b', linewidth=2)
    #ax3.grid(True)

    # Calculate cost function, using the formula directly
    # The sum of squares is used here
    cost = 0.0
    for i in range(15):
        t = pow(rm/r34[i],b[i])
        res = np.log(34) - np.log(vm[i]) - (0.5+a*(r34[i]-rm))*(np.log(t) + 1-t)
        cost += res**2

    for i in range(15):
        t = pow(rm/r50[i],b[i])
        res = np.log(50) - np.log(vm[i]) - (0.5+a*(r50[i]-rm))*(np.log(t) + 1-t)
        cost += res**2

    for i in range(15):
        t = pow(rm/r64[i],b[i])
        res = np.log(64) - np.log(vm[i]) - (0.5+a*(r64[i]-rm))*(np.log(t) + 1-t)
        cost += res**2

    return (cost)

# Bounds for the parameters: without these many local minima are found
bou = [(-1,1),(minrm, 2*minrm), \
      (0.5,3), (0.5,3), (0.5,3), (0.5,3), (0.5,3), (0.5,3), (0.5,3), (0.5,3), \
      (0.5,3), (0.5,3), (0.5,3), (0.5,3), (0.5,3), (0.5,3), (0.5,3), (0.5,3), \
      (80,110), (80,110), (80,110), (80,110), (80,110), (80,110), (80,110), (80,110), 
      (80,110), (80,110), (80,110), (80,110), (80,110), (80,110), (80,110), (80,110)]

res = minimize(Holland2010Cost, x0, method='SLSQP', bounds=bou, \
      options={'maxiter': 1000, 'disp': True})

print("Final x:")
print(res.x)
#print(res.success)
print(res.message)

# Plot the 16 wind profiles up to 200 kt distance
a = res.x[0]
r_max = res.x[1]
x_axis = range(1,201)
ws = np.array(np.zeros(200))
for j in range(15):
    b = res.x[j+2]
    v_max = res.x[j+18]
    for i in x_axis:
        ws[i] = Holland2010(i, r_max, b, a, v_max)

    ax0 = plt.subplot(111, polar=False)
    ax0.plot(x_axis, ws, color='r', linewidth=1)
    ax0.grid(True)
 