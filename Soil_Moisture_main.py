# -*- coding: utf-8 -*-
"""
Created on Thu Feb 23 15:47:54 2017

@author: Qina Yan and Phong Le

This program simulates vertical 1D soil moisture. The model resolves the mixed form variably saturated Richard equation. The method is based on finite difference following Celia et al (1990).

The current version has:
    1. uniform grid size
    2. uniform soil properties, e.g. builk density rho_s, Ksat, 
    3. BC: Constant head at two ends

Key state variables:
    Psi   : pressure head [L]
    theta : soil moisture [-]

Other variables:
    C    : Specific moisture capacity [1/L]
    K    : Hydraulic conductivity [L/T]

Notation:
    _ini    : Initial condition (time 0)
    _n      : At time step n
    _np1    : At time n+1 (known)
    _np1m   : At time n+1, iteration m (known)
    _np1mp1 : At time n+1, iteration m+1 (unknown)

"""

from parameters import *
import numpy as np
from F_SoilMoisture import *
import matplotlib.pyplot as plt


#------------------------------------------------------------------------
# ADD ATMOSPHERIC FORCING
#------------------------------------------------------------------------
#...................Rainfall
RF  = 1.5*1e-4*60.0             # Rainfall rate in unit [m/hr]
PPT = RF * dt                   # Rainfall depth per time step. [m]

#...................Evapotranspiration 
Evap = 8.4*1e-5*60.0            # Rainfall rate in unit [m/hr]
ET   = 0*Evap *dt                 # Evaporation per time step dt [m] 

#------------------------------------------------------------------------
# Initial Condition and grid size
#------------------------------------------------------------------------
nz       = int(total_soil_depth/dz)     # number of vertical grid size
Psi_ini  = np.linspace(-3,2,nz)         # Initial pressure head
Psi      = Psi_ini

C_ini,K_ini,theta_ini = vanGenuchten(Psi_ini,alpha, theta_S, theta_R, n, m, Ksat)
 
C        = C_ini
K        = K_ini
theta    = theta_ini

plt.plot(theta_ini, -dz*np.arange(nz))

#......... Model Stopping tolerance for convegence 
stop_tol = 0.01
max_iter = 100
Psimin   = 0.0001

#------------------------------------------------------------------------
# Boundary Conditions (constant head at both end)
#------------------------------------------------------------------------
Psi_top = -2.0 #PPT-ET
Psi_bot = 1.0

#------------------------------------------------------------------------
# The Main loop
#------------------------------------------------------------------------
time_steps = 5
Psi_store = np.zeros([nz,time_steps])

for i in xrange(time_steps): 
    # Initialize the Picard iteration solver
    
    Psi_n = Psi
    C_n,K_n,theta_n = vanGenuchten(Psi_n,alpha, theta_S, theta_R, n, m, Ksat)
    Psi_np1m = Psi_n
    loop_num = 0
    deltam   = 10000.0
    while deltam > stop_tol and loop_num < max_iter:
        C_np1m, K_np1m, theta_np1m = vanGenuchten(Psi_np1m,alpha, theta_S, theta_R, n, m, Ksat)
        K_np1m_ph, K_np1m_mh = Khalf(nz,K_np1m)
            
        Psi_np1mp1 = ImplicitSolver_Psi_npimp1(Psi_n, Psi_np1m, theta_n, theta_np1m, C_np1m,
                                               K_np1m_ph, K_np1m_mh, nz, Psi_top, Psi_bot)
        deltam = np.max(np.abs(Psi_np1mp1 - Psi_np1m))
        Psi_np1m = Psi_np1mp1
        loop_num = loop_num+1        
    
    print loop_num
    Psi = Psi_np1mp1
    Psi_store[:,i] = Psi
    
    if (i)%10 == 0:    
        plt.figure()
        plt.plot(Psi, -dz*np.arange(nz))

plt.figure()
plt.imshow(Psi_store)
plt.show()

#    Psi_np1 = 
#    Psi = Psi_np1
#   C,K,theta = vanGenuchten(Psi,alpha, theta_S, theta_R, n, m, Ksat)





