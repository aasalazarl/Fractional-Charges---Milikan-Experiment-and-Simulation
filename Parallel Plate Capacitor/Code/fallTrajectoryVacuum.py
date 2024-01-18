# -*- coding: utf-8 -*-
"""
Created on Fri Jun 21 16:08:38 2019

@author: alejandrosalazar

Fall trajectory in vacuum.
"""

import matplotlib.pyplot as plt
import numpy as np
import math
from numpy import zeros

def Int(i, variable):
    return (str(int(variable[i])))  

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
# Prepare variables.
#------------------------------------------------------------------------------
E = 1000000                                         # volts / meter; electric field.
g = 9.806                                           # m / s^2; acceleration due to gravity.
e = 1.602176634 * 10**(-19)                         # coulombs; electron charge.
#------------------------------------------------------------------------------
# Air:
Rair = 287.05                                       # joules / (kg mol); individual gas constant of air.
#------------------------------------------------------------------------------
# Glass sphere:
a = (25 * 10**(-6))/2                               # meters; radius of glass sphere.
sigma = 8000                                        # k / m^3; density of glass.
m = (4/3) * np.pi * a**3 * sigma                       # kg; mass of glass sphere.
V = (4 / 3) * np.pi * a**3                          # m^3; radius of sphere.
q = 1 * e
#------------------------------------------------------------------------------
# Time:
t0 = 0                                                                                 # Initial time.
tFall = np.sqrt(20 / g)
steps = tFall                                              # Number of time intervals.
time = np.linspace(t0, math.ceil(tFall), 10000)        # Times.      
#------------------------------------------------------------------------------
# Storage lists:
x_coordinate = np.linspace(t0, int(np.max(time)), len(time) + 1)     # Variable to store x-coordinates.
y_coordinate = np.linspace(t0, int(np.max(time)), len(time) + 1)     # Variable to store y-coordinates.

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
# Determine landing positions and times of fall.
#------------------------------------------------------------------------------
for l in range(0, len(time)):
    #------------------------------------------------------------------
    # x-coordinates:
    x_coordinate[l] = ((q * E) / (2 * m)) * time[l]**2
    #------------------------------------------------------------------
    # y-coordinates:
    y_coordinate[l] = - (g / 2) * time[l]**2
        
    if y_coordinate[l] < -10:
        size = l
        break
#------------------------------------------------------------------------------
# Store definite values of x- and y-coordinates.
r = zeros([2, size + 1])                    # Create a variable to store x- and y-coordinates.
    
for i in range(0, size + 1):                # Store them.
    r[0, i] = x_coordinate[i] * 1000
    r[1, i] = y_coordinate[i]
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
# Plot the trajectory.
plt.plot(r[0, :], r[1, :], linewidth = 0.9, label = 'Analytical solution') # Plots y over x analytically.
plt.grid()
plt.xlabel('Horizontal axis [millimeters]')
plt.ylabel('Vertical axis [meters]')
plt.text(1.5, 1, 'E-field = 1000 V/mm; vacuum:', va = 'center', ha = 'center', fontsize = 11)
plt.legend()
plt.savefig('fallTrajectoryVacuum.png', format = 'png', dpi = 300, bbox_inches = 'tight')

# Print some information.
print ("Mass of sphere:", m, "kg")