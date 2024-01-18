# -*- coding: utf-8 -*-
"""
Created on Thu Jun 20 15:25:23 2019

@author: alejandrosalazar

Determines the deviation from expected landing position of a glass sphere due to Brownian motion.
"""
''' Perform all the calculations in the frame of reference of the accelerated sphere unaffected
by Brownian motion'''

import numpy as np
import math
import random
from scipy.optimize import fsolve

def Int(i, variable):
    return (str(int(variable[i])))  

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
# Prepare variables.
#------------------------------------------------------------------------------
E = [500000, 1000000]                               # volts / meter; electric field.
g = 9.806                                           # m / s^2; acceleration due to gravity.
e = 1.602176634 * 10**(-19)                         # coulombs; electron charge.
h = 10                                              # meters; height of fall.
pi = math.pi
exp = math.exp
#------------------------------------------------------------------------------
# Air:
pressure = np.multiply([1, 2, 5, 10, 20, 50, 100], 100000) # pascals; pressure.
Rair = 287.05                                      # joules / (kg mol); individual gas constant of air.
T = 293.15                                          # kelvin = 20 C; temperature.
rho = pressure / (Rair * T)                               # kg / m^3; air density.
eta = 1.8 * 10**(-5)                                # kg / (m s); (dynamic) viscosity @ 20 C and 10 atm.
#------------------------------------------------------------------------------
# Glass sphere:
a = (25 * 10**(-6))/2                               # meters; radius of glass sphere.
sigma = 8000                                        # k / m^3; density of glass.
m = (4/3) * pi * a**3 * sigma                       # kg; mass of glass sphere.
V = (4 / 3) * np.pi * a**3                          # m^3; radius of sphere.
m_eff = V * (sigma - rho)                           # kg; effective mass.
#q = [1 * e, 0.1 * e, 0.01 * e, 0.001 * e]           # coulombs; charge of the sphere.
q = [1 * e, 0.1*e, 0.01*e, 0.001*e]
#------------------------------------------------------------------------------
# Other constants:
B = 0.100                                           # Experimental constant. Taken from McKeehan, p. 37: B = 0.00753 mm of Hg * cm.
k = (1 + (B / (a * pressure)))**(-1)
K = ((9 * eta) / (2 * a**2 *(sigma - rho))) * k   
#------------------------------------------------------------------------------
# Variables for Brownian motion:
R = 8.314462618                                     # joules / (mol kelvin); universal gas constant.
kB =  1.380649 *10**(-23)                          # m^2 kg s^-2 K^-1; Boltzmann's constant.
N = R/kB                                            # number of air molecules in one gram.
#------------------------------------------------------------------------------
# Time:
t0 = 0                                                                                 # Initial time.

# Define function to be solved for fall time.
def timeFall(tFall, pressure):  
    rho = pressure / (Rair * T)                     # Density of the air. 
    K = ((9 * eta) / (2 * a**2 *(sigma - rho))) * (1 + (B / (a * pressure)))**(-1)  
    return ((10 * K) / g) * (1 - (rho / (sigma - rho)))**(-1) - (1 / K) * (np.exp(- K * tFall) - 1) - tFall

dt = 1/400

# Fall time is the same for all charges.
tL = {}
for i in range(0, len(pressure)):
    tFall = fsolve(timeFall, 60, pressure[i])
    steps = tFall/dt                                        # Number of time intervals.
    tL['tL{}'.format(Int(i, pressure/100000))] = np.linspace(t0, math.ceil(tFall), steps + 1)        # Times.  
    
#------------------------------------------------------------------------------

# For 300 spheres:
size = np.linspace(0, 101, 100)                                 # Variable to store size of lists when sphere touches the ground.
landing_positions = np.linspace(0, 101, 100)                    # Variable to store landing positions (x-positions).
times_fall = np.linspace(0, 101, 100)                           # Variable to store times of fall.

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
# Determine landing positions and times of fall.
#------------------------------------------------------------------------------
# Landing positions statistics:
mean_landing_positions = [0, 0, 0, 0]           # Variable to store mean landing positions for different charges.
variances = [0, 0, 0, 0]
std_deviations = [0, 0, 0, 0]
mean_fall_times = [0, 0, 0, 0]

# Open txt file. (A)
f = open("Brownian motion contribution (2).txt",'w')
f.write("Brownian motion contribution \n")
f.write("100 spheres and {} collisions per second".format(str(dt**(-1))))

for iii in range(0, len(E)):
    
    # Record results in txt file.
    f.write("\n-------------------------------------------\n")
    f.write("-------------------------------------------\n")
    f.write("Electric field: {} volts/mm \n".format(str(E[iii]/1000)))
    f.write("---------Results---------\n")
    
    # Print results.
    print()
    print('-----------------------------------------------')
    print('-----------------------------------------------')
    print("Electric field", E[iii]/1000, 'volts/mm')
    print("---------Results---------")
    
    for i in range(0, len(q)):    
        #--------------------------------------------------------------------------
        # Buoyancy force.
        Fbx = rho * V * ((q[i] * E[iii]) / m_eff)
        Fby = rho * V * g
    
        # Record results in txt file.
        f.write("-------------------------------------------\n")
        f.write("For a sphere of charge {} e: \n \n".format(str(q[i]/e)))
        
        # Print results
        print('------------------------------------------')
        print("For a sphere of charge", q[i] / e, "e:")
        print()
    
        for ii in range(0, len(pressure)):
            time = tL['tL' + Int(ii, pressure/100000)]
            #print('charge' + str(q[i]/e) + 'pressure' + Int(ii, pressure/100000), time)
            #print(time[2])
            # Storage lists:
            x_coordinates = np.linspace(t0, int(np.max(time)), len(time) + 1)     # Variable to store x-coordinates.
            y_coordinates = np.linspace(t0, int(np.max(time)), len(time) + 1)     # Variable to store y-coordinates.
            
            for j in range(0, len(size)): # 500 spheres are let fall.
                
                for l in range(0, len(time)):
                    #------------------------------------------------------------------
                    #------------------------------------------------------------------
                    # Average displacements from Brownian motion:
                    #------------------------------------------------------------------
                    # x-displacement:
                    pos_neg = [-1, 1]           # +/- 1 for right/up and left/down displacement due to Brownian motion.
                    random.shuffle(pos_neg)     # Randomly, generate a right- or left-oriented collision. 
                    fx = pos_neg[0]
                    delta_x = np.sqrt((4 / (np.pi * K[ii])) * (((R * T) / N) - (Fbx[ii] - q[i] * E[iii]) *(dt - (m_eff[ii] / K[ii])) * np.sqrt((R * T) / (N * m_eff[ii]))) * dt)
                    #------------------------------------------------------------------
                    # y-displacement:
                    random.shuffle(pos_neg)     # Randomly, generate upward- or downward-oriented collision.
                    fy = pos_neg[0]
                    delta_y = np.sqrt((4 / (np.pi * K[ii])) * (((R * T) / N) - (Fby[ii] - (m_eff[ii] * g)) * (dt - (m_eff[ii] / K[ii])) * np.sqrt((R *T) / (N * m_eff[ii]))) * dt)
                    #------------------------------------------------------------------
                    # x-coordinates:
                    #x_coordinates[l] = ((q[i] * E) / (6 * pi * eta * a * k)) * (tL[l] + (m / (6 * pi * eta * a * k)) * (exp(- ( (6 * pi * eta * a * k) / m) * tL[l]) - 1)) + delta_x    # Replace values for x-coordinates.
                    # Turn E-field off:
                    x_coordinates[l] = 0 + delta_x
                    
                    #------------------------------------------------------------------
                    # y-coordinates:
                    y_coordinates[l] = ((g / K[ii]) * (time[l] + (1 / K[ii]) * (np.exp(- K[ii] * time[l]) - 1)) * (1 - (rho[ii] / (sigma - rho[ii]))))      # Replace non-terminal y-displacements.
                    
                    if y_coordinates[l] > 10:
                        landing_positions[j] = x_coordinates[l] * 1000          # Update landing positions.
                        times_fall[j] = time[l]                                 # Update fall times.
                        break
                    
            #--------------------------------------------------------------------------
            #--------------------------------------------------------------------------
            # Determine statistics.
            #--------------------------------------------------------------------------
            # Mean:
            mean_landing_positions[i] = np.mean(landing_positions)
            mean_fall_times[i] = np.mean(times_fall)
            #--------------------------------------------------------------------------
            # Variance:
            di = np.linspace(0, len(landing_positions) + 1, len(landing_positions))
            
            for j in range(0, len(landing_positions)):
                di[j] = landing_positions[j] - mean_landing_positions[i]
                
            variances[i] = (1 / (len(landing_positions) - 1)) * np.sum(di**2)
            #--------------------------------------------------------------------------
            # Std. deviation:
            std_deviations[i] = variances[i]**0.5
    
            # Record results in txt file.
            f.write("Pressure: {} bar \n".format(str(pressure[ii]/100000)))
            f.write("Mean landing position: {} mm \n".format(str(mean_landing_positions[i])))
            f.write('Variance: {} mm^2. \n'.format(str(variances[i])))
            f.write('Std. deviation: {} mm. \n'.format(str(std_deviations[i])))
            f.write('Mean fall time: {} seconds.\n \n'.format(str(mean_fall_times[i])))
    
            # Print results.
            print("Pressure", pressure[ii]/100000, "bar:")
            print("Mean landing position:", mean_landing_positions[i], "mm.")
            print("Variance:", variances[i], "mm^2.")
            print("Std. deviation:", std_deviations[i], "mm.")
            print("Mean fall time:", mean_fall_times[i], "seoconds.")
            print()

# (A) Close txt file.
f.close()