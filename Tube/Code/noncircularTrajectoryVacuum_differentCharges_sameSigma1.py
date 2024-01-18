# -*- coding: utf-8 -*-
"""
Created on Thu Jun 27 14:37:10 2019

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
E = 0                                       # volts/m; vertical electric field from parallel plate capacitors

# Glass sphere.
a = (25 * 10**(-6)) / 2                     # m; radius.
sigma = 8000                                # kg/m^3; density of glass.
V = (4 / 3) * np.pi * a**3                  # m^3; volume.
m = V * sigma                               # kg; mass.
m_eff = m                                   # kg; effective mass.
q = [1* e, 0.1 * e, 0.01 * e, 0.001 * e]    # coulombs; charge.

# Electrostatics.
sigma1 =  0.00003                                 # coulombs/m^3; surface charge density of central pillar.
r1 = 0.005                                         # m; radius of the central pillar.
r2 = 0.10                                          # m; radius of the outer cylinder.
epsilon = 8.8541878128 * 10**(-12)                 # Vacuum electric permittivity.
R = 0.05                                           # m; initial radius.
Epotential = - ((r1 * sigma1) / epsilon) * np.log(r1 / r2)    # volts; electric potential between central pillar and outer cylinder.

# Effective acceleration of gravity.
g_eff = np.divide(np.subtract((m_eff * g), (np.multiply(q, E))), m_eff)


#------------------------------------------------------------------------------
# To compute height of fall.
def height(tfall, i):
    return (0.5) * g_eff[i] * tfall**2

#------------------------------------------------------------------------------
# Determine changes in radius.

f = open('Vacuum radius change in free fall.txt', 'w')
f.write('Changes in radius for different charges in vacuum \n \n')
f.write('Vacuum \n')
f.write('Surface charge density on central pillar: {} coulombs/m^2 \n'.format(str(sigma1)))
f.write('Electric potential: {} kilo-volts \n \n'.format(str('%.4g' % float(Epotential/1000))))

print('Electric potential:', Epotential, 'volts')

for i in range(0, len(q)):
    
    f.write('\nCharge of {} e: \n'.format(str('%.3g' % float(q[i]/e))))
    
    print()
    print('Charge of', q[i]/e, 'e:') 
    print('Effective acceleration of gravity', g_eff[i], 'm/s^2')
    
    # Define parameters.
    C = (q[i] * r1 * sigma1) / (epsilon * m_eff)

    # System of ordinary differential equations.

    def sys(X, t):
        # I want to solve the system:
        # d2rdt2 = vtheta^2/r - qE/m_eff
        # dvxdt = - vtheta (drdt/r) 
        
        # Assign values to variables.
        r = X[0]
        z = X[1]
        
        # Differential equztions.
        dzdt = - C / r
        drdt = z
        
        return [drdt, dzdt]
    
    # Initial values for r, vtheta, z.
    X0 = [R, 0]  # inital radius = 1 m, inital radial speed = 0.
    
    # Time.
    t = np.linspace(0, 10, 100000)
    
    # Solve the system of ODEs.
    sol = odeint(sys, X0, t)
    
    r = sol[:, 0]
    z = sol[:, 1]
    
    plt.plot(t, r, label = 'Sphere of charge ' + str('%.1g' % float(q[i]/e)) + 'e') 
    plt.legend(loc = 'best')
    plt.xlabel('Time (seconds)')
    plt.ylabel('Radius (meters)')

    if i == len(q) - 1:
        plt.text(0.0, 0.055, 'Surface charge density $\sigma_1$ = ' + str(sigma1) + ' coulombs/m$^2$;\nelectric potential V = '\
                 + str('%.4g' % float(Epotential/1000)) + ' kilo-volts;\nvacuum; free fall:')
    plt.hold(True)

    plt.savefig('freeFallVacuum_differentCharges_sameSigma1.png', format = 'png', dpi = 300, bbox_inches = 'tight')   

    rfinal = np.min(r)
    tfall = np.max(t)
    yfall = height(tfall,i)
    
    f.write('Change in radius: {} mm \n'.format(str('%.4g' % float((R - rfinal)*1000))))
    f.write('Radial displacement in diameters: {} diameters \n'.format(float((R - rfinal) / (2 * a))))
    f.write('Fall time: {} seconds \n'.format(str('%.4g' % float(tfall))))
    f.write('Fall height: {} meters \n'.format(str('%.4g' % float(yfall))))
    
    print('Change in radius:', (R - rfinal)*1000, 'mm')
    print('Fall time:', tfall, 'seconds')
    print('Fall height:', yfall, 'meters')
    
'''    for j in range(0, len(r)):
        if r[j] < 0.02:
            rfinal = r[j]                   # Final radius.
            tfall = t[j]                    # Fall time.
            yfall = height(tfall, i)        # Fall height.
            
            f.write('Change in radius: {} meters \n'.format(str('%.4g' % float(R - rfinal))))
            f.write('Fall time: {} seconds \n'.format(str('%.4g' % float(tfall))))
            f.write('Fall height: {} meters \n'.format(str('%.4g' % float(yfall))))
            
            print('Change in radius:', R - rfinal, 'meters')
            print('Fall time:', tfall, 'seconds')
            print('Fall height:', yfall, 'meters')
            
            break'''
    
'''    time = np.linspace(0, tfall, 100)
    
    # Solve the ODE system to plot with new time.
    sol = odeint(sys, X0, time)
    radius = sol[:, 0]
    z = sol[:, 1]'''
    
f.close()