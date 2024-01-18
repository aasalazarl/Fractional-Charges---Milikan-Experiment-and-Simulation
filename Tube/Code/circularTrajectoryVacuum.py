# -*- coding: utf-8 -*-
"""
Created on Wed Jun 26 11:38:42 2019

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

# Glass sphere.
a = (25 * 10**(-6)) / 2                     # m; radius.
sigma = 8000                                # kg/m^3; density of glass.
V = (4 / 3) * np.pi * a**3                  # m^3; volume.
m = V * sigma                               # kg; mass.
m_eff = m                                   # kg; effective mass.
q = 0.1 * e                                   # coulombs; charge.
v0 = 0.1
#v0 = 0.001                                      # m/s; initial velocity (ramp at theta degrees from horizontal).
theta = 0                                  # degrees; ramp angle of inclination from the horizontal.
vx0 = v0 * np.cos(theta)                    # m/s; initial velocity horizontal component.
vy0 = v0 * np.sin(theta)                    # m/s; initial velocity vertical component.
vtheta0 = vx0                               # m/s; initial tangential velocity to circular motion.

# Electrostatics.
sigma1 = 0.2                                   # coulombs/m^3; surface charge density of central pillar.
#sigma1 = 0.00003
#sigma1 = 0.02
r1 = 0.005                                     # m; raidus of the central pillar.
r2 = 0.10                                      # m; radius of the outer tube.
epsilon = 8.8541878128 * 10**(-12)             # Vacuum electric permittivity.
Epotential = - ((r1 * sigma1) / epsilon) * np.log(r1 / r2)
R = 0.05                                        # m; initial radius.

# Define parameters.
C = (q * r1 * sigma1) / (epsilon * m_eff)

#------------------------------------------------------------------------------
# System of ordinary differential equations.

def sys(X, t):
    # I want to solve the system:
    # d2rdt2 = vtheta^2/r - qE/m_eff
    # dvxdt = - vtheta (drdt/r) 
    
    # Assign values to variables.
    r = X[0]
    vtheta = X[1]
    z = X[2]
    
    # Differential equztions.
    dzdt = (1 / r) * (vtheta**2 - C)
    drdt = z
    dvthetadt = -vtheta * (z / r)
    
    return [drdt, dvthetadt, dzdt]

# Initial values for r, vtheta, z.
X0 = [R, vtheta0, 0]  # inital radius = 1 m, inital tangential speed = 1 m/s, inital radial speed = 0.

# Time.
time = np.linspace(0, 1, 10000)
#time = np.linspace(0, 80, 10000)

# Solve the system of ODEs.
sol = odeint(sys, X0, time)

r = sol[:, 0]
vtheta = sol[:, 1]
z = sol[:, 2]

plt.plot(time, r, label = 'Radius trajectory (meters)') 
plt.plot(time, vtheta, label = 'Tangential velocity to circular motion (meters/second)')
plt.legend(loc = 'best')
plt.xlabel('Time (seconds)')

plt.text(0, 0.27, 'Surface charge density $\sigma_1$ = ' + str(float('%.3g' % sigma1)) + ' C/m$^2$;'\
         '\nelectric potential = ' + str('%.4g' % float(Epotential/1000000)) + \
         ' mega-volts; \ncharge of sphere = '+ str('%.1g' % float(q/e)) + 'e; vacuum; \nv$_0$ = ' + str(v0) + ' m/s:')

#plt.text(0, 0.055, 'Surface charge density $\sigma_1$ = ' + str(float('%.3g' % sigma1)) + ' C/m$^2$;'\
 #        '\nelectric potential = ' + str('%.4g' % float(Epotential/1000)) + \
  #       ' kilo-volts; \ncharge of sphere = '+ str('%.1g' % float(q/e)) + 'e; vacuum; \nv$_0$ = ' + str(v0) + ' m/s:')


plt.hold(True)

plt.savefig('circularTrajectoryVacuum_ellipse1.png', format = 'png', dpi = 300, bbox_inches = 'tight')
#plt.savefig('circularTrajectoryVacuum_ellipse2.png', format = 'png', dpi = 300, bbox_inches = 'tight')

# To cancel gravity:
E = (m_eff * g) / e
print('Evertical:', E/1000, 'volts/mm')