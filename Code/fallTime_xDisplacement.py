# -*- coding: utf-8 -*-
"""
Created on Wed Jun 19 15:34:26 2019

@author: alejandrosalazar

Determines the fall time of the glass sphere of charge e as a function of time and maximum horizontal displacement.
"""

import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
from scipy import interpolate


#------------------------------------------------------------------------------
# Assign values to parameters.
eta = 1.8 * 10**(-5)            # kg / (m s); dynamic viscosity of air
sigma = 8000                    # kg / m^2; density of glass.
g = 9.80665                     # m / s^2; acceleration due to gravity.
Rair = 287.05                   # joules / (kg mol); individual gas constant of air.
T = 293.15                      # kelvin; temperature.
B = 0.100                       # Experimental constant (see McKeehan).
a = 12.5 * 10**(-6)             # meters; radius of sphere.

#------------------------------------------------------------------------------
# Find fall time.

# Define function to be solved for fall time.
def timeFall(tFall, pressure):  
    rho = pressure / (Rair * T)                     # Density of the air. 
    K = ((9 * eta) / (2 * a**2 *(sigma - rho))) * (1 + (B / (a * pressure)))**(-1)  
    return ((10 * K) / g) * (1 - (rho / (sigma - rho)))**(-1) - (1 / K) * (np.exp(- K * tFall) - 1) - tFall

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
q = 1.602176634 * 10**(-19)           # coulombs; charge of glass sphere.
E = [500000, 1000000]                      # volts/meter; electric field (horizontal).

def xDisplacement(tFall, pressure, Efield):
    rho = pressure / (Rair * T)                     # Density of the air. 
    K = ((9 * eta) / (2 * a**2 *(sigma - rho))) * (1 + (B / (a * pressure)))**(-1)  
    return ((q * Efield) / ((4 / 3) * np.pi * a**3 * (sigma - rho) * K)) * (1 - (rho / (sigma - rho))) * (tFall + (1 / K) * (np.exp(-K * tFall) - 1))

#------------------------------------------------------------------------------    
# Determine horizontal displacements.

def Num(i):
    return (str(int(E[i]/1000)))    

#Create dictionary.
x = {}
xDummy = {}

for i in range(0, len(E)):
    xDummy['x{}'.format(Num(i))] = []
    x.update(xDummy)
    
# Compute horizontal displacements.
for i in range(0, len(pressure)):
    for j in range(0, len(E)):
        x['x' + Num(j)].append(xDisplacement(tFall[i], pressure[i], E[j]))
    
# Perform cubic spline interpolation.
tck = {}
xInterpol = {}

for i in range(0, len(E)):
    tck['tck{}'.format(Num(i))] = interpolate.splrep(pressure, x['x' + Num(i)], s = 0)
    pressureInterpol = pressure
    xInterpol['x{}Interpol'.format(Num(i))] = interpolate.splev(pressureInterpol, tck['tck' + Num(i)], der = 0) 

#------------------------------------------------------------------------------
# Plot.
f, axarr = plt.subplots(3, sharex = True, figsize = (7, 7))

# Time plot.
axarr[0].plot(pressure/100000, tFall, 'bo', markersize = 5)
axarr[0].plot(pressureInterpol/100000, timeInterpol, 'k-', markersize = 1)
#plt.text(100, 64, 'Max. fall time:\n{:f} sec at {:f} bar'.format(np.max(tFall), np.max(pressure) / 100000), fontsize = 9)
axarr[0].set_ylabel('Fall time \n (seconds)', fontsize = 13)

# Horizontal displacement plot.

for i in range(0, len(E)):
    axarr[i + 1].plot(pressure/100000, np.multiply(x['x' + Num(i)], 1000), 'bo', markersize = 5)
    axarr[i + 1].plot(pressureInterpol/100000, xInterpol['x' + Num(i) + 'Interpol'] * 1000, 'k-', markersize = 1)

# Annotate.
f.text(0.5, 0.625, 'E-field = 500 volts/mm:', va = 'center', ha = 'center', fontsize = 11)
f.text(0.5, 0.36, 'E-field = 1000 volts/mm:', va = 'center', ha = 'center', fontsize = 11)

# Axes labels.
f.text(0.5, 0.07, 'Pressure (bar)', va='center', ha='center', fontsize = 15)
f.text(0.04, 0.37, 'Horizontal displacement \n (mm)', va='center', ha='center', rotation='vertical', fontsize= 13)

plt.hold(True)

plt.savefig('fallTime_xDisplacement.png', format = 'png', dpi = 300, bbox_inches = 'tight')

# Print some information:
rho = pressure / (Rair * T)
V = (4 / 3) * np.pi * a**3
m = V * sigma
m_eff = V * (sigma - rho)
print("Volume of sphere:", V, "m^3")
print("Mass of sphere:", m, "kg")
print("Effective mass of sphere:", m_eff, "kg")