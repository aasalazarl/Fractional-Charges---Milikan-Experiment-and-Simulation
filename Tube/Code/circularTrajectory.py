# -*- coding: utf-8 -*-
"""
Created on Wed Jun 26 11:07:33 2019

@author: alejandrosalazar
"""

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

#------------------------------------------------------------------------------
# Variables.

# Constants.
e = 1.602176634 * 10**(-19)                 # coulombs; elementary charge of the electron.
g = 9.806                                   # m/s^2; acceleration due to gravity.

# Air.
pressure = 1 * 100000                      # pascals; air pressure.
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
q = 0.1 * e                                   # coulombs; charge.
v0 = 0.1                                     # m/s; initial velocity (ramp at theta degrees from horizontal).
theta = 0                                  # degrees; ramp angle of inclination from the horizontal.
vx0 = v0 * np.cos(theta)                    # m/s; initial velocity horizontal component.
vy0 = v0 * np.sin(theta)                    # m/s; initial velocity vertical component.
vtheta0 = vx0                               # m/s; initial tangential velocity to circular motion.

# Other constants.
B = 0.100                                           # Experimental constant. Taken from McKeehan, p. 37: B = 0.00753 mm of Hg * cm.
k = (1 + (B / (a * pressure)))**(-1)
K = ((9 * eta) / (2 * a**2 *(sigma - rho))) * k   
R = 0.05                                               # m; initial radius.

# Electrostatics.
sigma1 = 0.000025                                 # coulombs/m^3; surface charge density of central pillar.
r1 = 0.005                                     # m; radius of the central pillar.
r2 = 0.10                                      # m; radius of the outer tube.
epsilon = 8.8541878128 * 10**(-12)             # Vacuum electric permittivity.
Epotential = - ((r1 * sigma1) / epsilon) * np.log(r1 / r2)

# Define parameters.
C = (q * r1 * sigma1) / (epsilon * m_eff)

#------------------------------------------------------------------------------
# System of ordinary differential equations.

def sys1(X, t):
    # I want to solve the system:
    # d2rdt2 = vtheta^2/r - qE/m_eff
    # dvxdt = -Kvtheta - vtheta (drdt/r) 
    
    # Assign values to variables.
    r = X[0]
    vtheta = X[1]
    z = X[2]
    
    # Differential equztions.
    dzdt = (1 / r) * (vtheta**2 + C * (((rho * V) / (m_eff)) - 1)) - (K * z)
    drdt = z
    dvthetadt = -(K + (z / r)) * vtheta
   
    return [drdt, dvthetadt, dzdt]

# Initial values for r, vtheta, z.
X0 = [R, vtheta0, 0]  # inital radius = 1 m, inital tangential speed = 1 m/s, inital radial speed = 0.

# Time.
time = np.linspace(0, 0.5, 100)

# Solve the system of ODEs.
sol = odeint(sys1, X0, time)

r = sol[:, 0]
vtheta = sol[:, 1]
z = sol[:, 2]

plt.plot(time, r, label = 'Radius trajectory (meters)') 
plt.plot(time, vtheta, label = 'Tangential velocity to circular motion\n(meters/second)')
#plt.plot(time, z, label = 'Radial velocity')
plt.legend(loc = 'best')
plt.xlabel('Time (seconds)')
#plt.text(0.01, 0.11, 'Surface charge density $\sigma_1$ = ' + str(sigma1) + ' C/m$^2$;\nelectric potential = ' + str('%.4g' % float(Epotential/1000))\
 #        + ' kilo-volts;\ncharge of sphere =' + str('%.1g' % float(q/e)) + 'e; pressure = ' + str(pressure/100000) + ' bar \nv$_0$ = 0.1 m/s:', ha = 'left')

plt.text(0.01, 0.11, 'Surface charge density $\sigma_1$ = ' + str(sigma1) + ' C/m$^2$;\nelectric potential = ' + str('%.4g' % float(Epotential/1000))\
         + ' kilo-volts;\ncharge of sphere = ' + str('%.1g' % float(q/e)) + 'e; pressure = ' + str(pressure/100000) + ' bar; \nv$_0$ = ' + str(v0) + ' m/s:', ha = 'left')

#plt.savefig('circularTrajectory_1e.png', format = 'png', dpi = 300, bbox_inches = 'tight')
plt.savefig('circularTrajectory_01e.png', format = 'png', dpi = 300, bbox_inches = 'tight')
