# -*- coding: utf-8 -*-
"""
Created on Thu Jun 20 13:31:04 2019

@author: alejandrosalazar
"""

"""
Created on Wed Jun 19 15:34:26 2019

@author: alejandrosalazar

Determines the fall time of the glass sphere of charge e as a function of time and maximum horizontal displacement.
"""

import numpy as np
from scipy.optimize import fsolve
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
e = 1.602176634 * 10**(-19)                               # coulombs; electron charge.
q = [1 * e, 0.1 * e, 0.01 * e, 0.001 * e]           # coulombs; charge of glass sphere.
E = [500000, 1000000]                               # volts/meter; electric field (horizontal).

def xDisplacement(tFall, pressure, Efield, charge):
    rho = pressure / (Rair * T)                     # Density of the air. 
    K = ((9 * eta) / (2 * a**2 *(sigma - rho))) * (1 + (B / (a * pressure)))**(-1)  
    return ((charge * Efield) / ((4 / 3) * np.pi * a**3 * (sigma - rho) * K)) * (1 - (rho / (sigma - rho))) * (tFall + (1 / K) * (np.exp(-K * tFall) - 1))
#------------------------------------------------------------------------------    
# Determine horizontal displacements.

def Int(i, variable):
    return (str(int(variable[i])))  

def Float(i, variable):
    return (str(float(variable[i])))
    
#Create dictionary.
x = {}
xDummy = {}

for i in range(0, len(q)):  
    
    for j in range(0, len(E)):
        xDummy['x{}_charge{}'.format(Int(j, np.divide(E, 1000)), Float(i, np.divide(q, e)) + 'e')] = []
        x.update(xDummy)
        
    # Compute horizontal displacements.
    for j in range(0, len(pressure)):
        
        for k in range(0, len(E)):
            x['x' + Int(k, np.divide(E, 1000)) + '_charge' + Float(i, np.divide(q, e)) + 'e'].append(xDisplacement(tFall[j], pressure[j], E[k], q[i]))

#------------------------------------------------------------------------------
# Create file with results.
f = open("Maximal horizontal displacements.txt",'w')
f.write("Maximal horizontal displacements \n")

for i in range(0, len(E)):
    f.write('\n-----------------------------------------------------\n')
    f.write('E-field = {} volts/mm \n'.format(Int(i, np.divide(E, 1000))))
    
    for j in range(0, len(q)):
        f.write('\nCharge {}e \n'.format(Float(j, np.divide(q, e))))
        f.write('Pressure (bar) | Horizontal displacement (mm) | Horizontal displacements [sphere diameters] \n')
        
        for k in range(0, len(pressure)):
            f.write('{} | {} | {} \n'.format(str(pressure[k]/100000), str(float('%.8g' % x['x' + Int(i, np.divide(E, 1000)) + '_charge' + Float(j, np.divide(q, e)) + 'e'][k]) * 1000), str(float('%.4g' % ((x['x' + Int(i, np.divide(E, 1000)) + '_charge' + Float(j, np.divide(q, e)) + 'e'][k]) / (a *2))))))

f.close()