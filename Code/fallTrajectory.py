# -*- coding: utf-8 -*-
"""
Created on Fri Jun 21 15:31:19 2019

@author: alejandrosalazar

Fall trajectory.
"""
import matplotlib.pyplot as plt
import numpy as np
import math
import random
from numpy import zeros
from scipy.optimize import fsolve

def Int(i, variable):
    return (str(int(variable[i])))  

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
# Prepare variables.
#------------------------------------------------------------------------------
E = 1000000                                         # volts / meter; electric field.
g = 9.806                                           # m / s^2; acceleration due to gravity.
e = 1.602176634 * 10**(-19)                         # coulombs; electron charge.
h = 10                                              # meters; height of fall.
#------------------------------------------------------------------------------
# Air:
pressure = 10 * 100000                              # pascals; pressure.
Rair = 287.05                                       # joules / (kg mol); individual gas constant of air.
T = 293.15                                          # kelvin = 20 C; temperature.
rho = pressure / (Rair * T)                               # kg / m^3; air density.
eta = 1.8 * 10**(-5)                                # kg / (m s); (dynamic) viscosity @ 20 C and 10 atm.
#------------------------------------------------------------------------------
# Glass sphere:
a = (25 * 10**(-6))/2                               # meters; radius of glass sphere.
sigma = 8000                                        # k / m^3; density of glass.
m = (4/3) * np.pi * a**3 * sigma                       # kg; mass of glass sphere.
V = (4 / 3) * np.pi * a**3                          # m^3; radius of sphere.
m_eff = V * (sigma - rho)                           # kg; effective mass.
#q = [1 * e, 0.1 * e, 0.01 * e, 0.001 * e]           # coulombs; charge of the sphere.
q = 1 * e
#------------------------------------------------------------------------------
# Other constants:
B = 0.100                                           # Experimental constant. Taken from McKeehan, p. 37: B = 0.00753 mm of Hg * cm.
k = (1 + (B / (a * pressure)))**(-1)
K = ((9 * eta) / (2 * a**2 *(sigma - rho))) * k   
#------------------------------------------------------------------------------
# Time:
t0 = 0                                                                                 # Initial time.

# Define function to be solved for fall time.
def timeFall(tFall, pressure):  
    rho = pressure / (Rair * T)                     # Density of the air. 
    K = ((9 * eta) / (2 * a**2 *(sigma - rho))) * (1 + (B / (a * pressure)))**(-1)  
    return ((10 * K) / g) * (1 - (rho / (sigma - rho)))**(-1) - (1 / K) * (np.exp(- K * tFall) - 1) - tFall

tFall = fsolve(timeFall, 60, pressure)
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
    x_coordinate[l] = ((q * E) / ((4 / 3) * np.pi * a**3 * (sigma - rho) * K)) * (1 - (rho / (sigma - rho))) * (time[l] + (1 / K) * (np.exp(-K * time[l]) - 1))
    #------------------------------------------------------------------
    # y-coordinates:
    y_coordinate[l] = - (g / K) * (time[l] + (1 / K) * (np.exp(- K * time[l]) - 1)) * (1 - (rho / (sigma - rho)))
        
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
plt.text(1.5, 1.3, 'E-field = 1000 V/mm; pressure = 10 bar:', va = 'center', ha = 'center', fontsize = 11)
plt.legend()
plt.savefig('fallTrajectory.png', format = 'png', dpi = 300, bbox_inches = 'tight')