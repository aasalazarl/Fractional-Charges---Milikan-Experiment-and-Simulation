# -*- coding: utf-8 -*-
"""
Created on Thu Jun 27 15:20:31 2019

@author: alejandrosalazar
"""

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
#from scipy.optimize import fsolve

#------------------------------------------------------------------------------
# Variables.

# Constants.
e = 1.602176634 * 10**(-19)                 # coulombs; elementary charge of the electron.
g = 9.806                                   # m/s^2; acceleration due to gravity.

# Air.
pressure = 1 * 100000                      # bar; air pressure.
Rair = 287.05                               # joules / (kg mol); individual gas constant of air.
T = 293.15                                  # kelvin; temperature.
rho = pressure / (Rair * T)                 # kg/m^3; air density.
eta = 1.8 * 10**(-5)                        # kg/(m s); (dynamic) viscosity of air.

# Glass sphere.
a = (25 * 10**(-6)) / 2                     # m; radius.
sigma = 8000                                # kg/m^3; density of glass.
V = (4 / 3) * np.pi * a**3                  # m^3; volume.
m = V * sigma                               # kg; mass.
m_eff = V * (sigma - rho)                   # kg; effective mass.
q = [1* e, 0.1 * e, 0.01 * e, 0.001 * e]               # coulombs; charge.

# Other constants.
B = 0.100                                           # Experimental constant. Taken from McKeehan, p. 37: B = 0.00753 mm of Hg * cm.
k = (1 + (B / (a * pressure)))**(-1)
K = ((9 * eta) / (2 * a**2 *(sigma - rho))) * k   

# Electrostatics.
sigma1 = 0.000025                           # coulombs/m^3; surface charge density of central pillar.
r1 = 0.005                                   # m; raidus of the central pillar.
r2 = 0.10
epsilon = 8.8541878128 * 10**(-12)           # Vacuum electric permittivity.
R = 0.05                                      # m; initial radius.
Epotential = - ((r1 * sigma1) / epsilon) * np.log(r1 / r2)    # volts; electric potential between central pillar and outer cylinder.

#------------------------------------------------------------------------------
# To compute height of fall.
def height(tfall):
    return (g / K) * (tfall + (1 / K) * (np.exp(-K * tfall) -1)) * (1 - (rho / (sigma - rho)))

#------------------------------------------------------------------------------
# Determine changes in radius.

f = open('Radius change under pressure in free fall.txt', 'w')
f.write('Change in radius for different charges in air \n \n')
f.write('Pressure: {} bar \n'.format(str(pressure/100000)))
f.write('Surface charge density on central pillar: {} coulombs/m^2 \n'.format(str(sigma1)))
f.write('Electric potential: {} volts \n'.format(str('%.4g' % float(Epotential))))

print('Pressure:', pressure/100000, 'bar')
print('Surface charge density on central pillar: {} coulombs/m^2'.format(str(sigma1)))
print('Electric potential:', Epotential, 'volts')

for i in range(0, len(q)):
    
    f.write('\nCharge of {} e: \n'.format(str('%.3g' % float(q[i]/e))))
    
    print()
    print('Charge of', q[i]/e, 'e:') 
    
    # Define parameters.
    C = (q[i] * r1 * sigma1) / (epsilon * m_eff)

    # System of ordinary differential equations.

    def sys(X, t):
        # I want to solve the system:
        # d2rdt2 = vtheta^2/r - qE/m_eff
        # dvxdt = -Kvtheta - vtheta (drdt/r) 
        
        # Define parameters.
        C = (q[i] * r1 * sigma1) / (epsilon * m_eff)
        
        # Assign values to variables.
        r = X[0]
        z = X[1]
        
        # Differential equztions.
        dzdt = (C / r) * (((rho * V) / (m_eff)) - 1) - (K * z)
        drdt = z
        
        return [drdt, dzdt]
    
    # Initial values for r, vtheta, z.
    X0 = [R, 0]  # inital radius = 1 m, inital tangential speed = 0 m/s, inital radial speed = 0.
    
    # Time.
    t = np.linspace(0, 60, 10000)
    
    # Solve the system of ODEs.
    sol = odeint(sys, X0, t)
    
    r = sol[:, 0]
    z = sol[:, 1]
    
    rfinal = np.min(r)          # Final radius.
    tfall = np.max(t)           # Fall time.
    yfall = height(tfall)       # Fall height.
    
    f.write('Change in radius: {} mm \n'.format(str('%.4g' % float((R - rfinal) * 1000))))
    f.write('Radial displacement in diameters: {} diameters \n'.format(float((R - rfinal) / (2 * a))))
    f.write('Fall time: {} seconds \n'.format(str('%.4g' % float(tfall))))
    f.write('Fall height: {} meters \n'.format(str('%.4g' % float(yfall))))
    
    print('Change in radius:', (R - rfinal) * 1000, 'mm')
    print('Fall time:', tfall, 'seconds')
    print('Fall height:', yfall, 'meters')
        
'''    for j in range(0, len(r)):
        if r[j] < 0.03:
            rfinal = r[j]                   # Final radius.
            tfall = t[j]                    # Fall time.
            yfall = height(tfall)           # Fall height.
            
            f.write('Change in radius: {} meters \n'.format(str('%.4g' % float(R - rfinal))))
            f.write('Fall time: {} seconds \n'.format(str('%.4g' % float(tfall))))
            f.write('Fall height: {} meters \n'.format(str('%.4g' % float(yfall))))
            
            print('Change in radius:', R - rfinal, 'meters')
            print('Fall time:', tfall, 'seconds')
            print('Fall height:', yfall, 'meters')
            print() 
            
            break'''
    
'''    time = np.linspace(0, tfall, 100)
    
    # Solve the ODE system to plot with new time.
    sol = odeint(sys, X0, time)
    radius = sol[:, 0]
    drdt = sol[:, 1]'''
    
    #plt.plot(t, r, label = 'Sphere of charge ' + str('%.3g' % float(q[i]/e)) + 'e') 
    #plt.legend(loc = 'best')
    #plt.xlabel('Time (seconds)')
    #plt.ylabel('Radius (meters)')

    #if i == len(q) - 1:
     #   plt.text(0.0, 0.320, 'Surface charge density $\sigma_1$ = ' + str(sigma1) + ' coulombs/m$^2$;\nelectric potential V = '\
      #           + str('%.4g' % float(Epotential)) + ' volts;\npressure = ' + str(pressure/100000) + ' bar; free fall:')

    #plt.hold(False)
    
#plt.savefig('freeFallAir_differentCharges_sameSigma1.png', format = 'png', dpi = 300, bbox_inches = 'tight')


f.close()