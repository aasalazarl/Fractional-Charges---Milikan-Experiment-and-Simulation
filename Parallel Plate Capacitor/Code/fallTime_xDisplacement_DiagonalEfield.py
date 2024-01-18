# -*- coding: utf-8 -*-
"""
Created on Fri Jun 21 10:58:31 2019

@author: alejandrosalazar

Determines the fall time of the glass sphere of charge e as a function of time 
and maximum horizontal displacement when E-field is diagonally oriented.
"""


import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
from scipy import interpolate


#------------------------------------------------------------------------------
# Assign values to parameters.
Ey = 1000000 / np.sqrt(2)
q = 1.602176634 * 10**(-19)     # coulombs; charge of glass sphere.
eta = 1.8 * 10**(-5)            # kg/(m s); dynamic viscosity of air
sigma = 8000                    # kg/m^2; density of glass.
g = 9.80665                     # m/s^2; acceleration due to gravity.
Rair = 287.05                   # joules / (kg mol); individual gas constant of air.
T = 293.15                      # kelvin; temperature.
B = 0.100                       # Experimental constant (see McKeehan).
a = (25 * 10**(-6)) / 2         # meters; radius of sphere.
V = (4 / 3) * np.pi * a**3      #m^3; volume of sphere.

#------------------------------------------------------------------------------
# Find fall time.

# Define function to be solved for fall time.
def timeFall(tFall, pressure):  
    rho = pressure / (Rair * T)                     # Density of the air. 
    m_eff = V * (sigma - rho)
    g_eff = g - ((q * Ey) / m_eff)             # m/s^2; effective acceleration of gravity due to vertical component of E-field.
    K = ((9 * eta) / (2 * a**2 *(sigma - rho))) * (1 + (B / (a * pressure)))**(-1)  
    return ((10 * K) / g_eff) * (1 - (rho / (sigma - rho)))**(-1) - (1 / K) * (np.exp(- K * tFall) - 1) - tFall

# Pressure.
pressure = np.multiply([1, 2, 3, 4, 5, 10, 15, 20, 25, 30, 40, 50, 100], 100000) # pascals; pressure.

# Compute fall times.
tFall = []

for i in range(0, len(pressure)):
    tFall.append(fsolve(timeFall, 60, pressure[i]))
    
# Perform cubic spline interpolation.
tck = interpolate.splrep(pressure, tFall, s = 0)
pressureInterpol = pressure
timeInterpol = interpolate.splev(pressureInterpol, tck, der = 0)

#------------------------------------------------------------------------------
# Maximize horizontal displacement.
Ex = np.divide([1000000], np.sqrt(2))                # volts/meter; electric field (horizontal).

def xDisplacement(tFall, pressure, Exfield):
    rho = pressure / (Rair * T)                     # Density of the air. 
    K = ((9 * eta) / (2 * a**2 *(sigma - rho))) * (1 + (B / (a * pressure)))**(-1)  
    return ((q * Exfield) / ((4 / 3) * np.pi * a**3 * (sigma - rho) * K)) * (1 - (rho / (sigma - rho))) * (tFall + (1 / K) * (np.exp(-K * tFall) - 1))

#------------------------------------------------------------------------------    
# Determine horizontal displacements.

def Num(i):
    return (str(int(Ex[i]/1000)))    

#Create dictionary.
x = {}
xDummy = {}

for i in range(0, len(Ex)):
    xDummy['x{}'.format(Num(i))] = []
    x.update(xDummy)
    
# Compute horizontal displacements.
for i in range(0, len(pressure)):
    for j in range(0, len(Ex)):
        x['x' + Num(j)].append(xDisplacement(tFall[i], pressure[i], Ex[j]))
    
# Perform cubic spline interpolation.
tck = {}
xInterpol = {}

for i in range(0, len(Ex)):
    tck['tck{}'.format(Num(i))] = interpolate.splrep(pressure, x['x' + Num(i)], s = 0)
    pressureInterpol = pressure
    xInterpol['x{}Interpol'.format(Num(i))] = interpolate.splev(pressureInterpol, tck['tck' + Num(i)], der = 0) 

#------------------------------------------------------------------------------
# Plot.
f, axarr = plt.subplots(2, sharex = True, figsize = (7, 7))

# Time plot.
axarr[0].plot(pressure/100000, tFall, 'bo', markersize = 5)
axarr[0].plot(pressureInterpol/100000, timeInterpol, 'k-', markersize = 1)
#plt.text(100, 64, 'Max. fall time:\n{:f} sec at {:f} bar'.format(np.max(tFall), np.max(pressure) / 100000), fontsize = 9)
axarr[0].set_ylabel('Fall time \n (seconds)', fontsize = 13)

# Horizontal displacement plot.

for i in range(0, len(Ex)):
    axarr[i + 1].plot(pressure/100000, np.multiply(x['x' + Num(i)], 1000), 'bo', markersize = 5)
    axarr[i + 1].plot(pressureInterpol/100000, xInterpol['x' + Num(i) + 'Interpol'] * 1000, 'k-', markersize = 1)
    axarr[i + 1].set_ylabel('Horizontal displacement \n (mm)', fontsize = 13)

# Annotate.
f.text(0.5, 0.9, 'E-field =' + str(int((Ey * np.sqrt(2))/1000)) + ' volts/mm:', va = 'center', ha = 'center', fontsize = 11)

# Axes labels.
f.text(0.5, 0.07, 'Pressure (bar)', va='center', ha='center', fontsize = 15)

plt.hold(True)

plt.savefig('fallTime_xDisplacement_diagonalEfield.png', format = 'png', dpi = 300, bbox_inches = 'tight')

# Print some information.
rho = pressure / (Rair * T)                     
m_eff = V * (sigma - rho)
g_eff = g - ((q * Ey) / m_eff)
print("Effective mass:", m_eff, "kg")
print("Effective acceleration of gravity:", g_eff, "m/s^2")
