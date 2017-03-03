# -*- coding: utf-8 -*-
"""
Created on Thu Feb 23 15:47:54 2017

@author: Qina and Phong

This code is to simulate soil  moisture in 1D soil column

The current version has:
1. uniform grid size
2. uniform soil properties, e.g. builk density \rho_s, Ksat, 
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
RF  = 1.5*1e-4*60.0            #Rainfall rate in unit [m/hr]
PPT = RF * dt                 # Rainfall depth per time step. [m]

#...................Evapotranspiration 
Evap = 8.4*1e-5*60.0         #Rainfall rate in unit [m/hr]
ET   = Evap *dt                 #Evaporation per time step dt [m] 

#------------------------------------------------------------------------
# Initial Condition and grid size
#------------------------------------------------------------------------
nz       = int(total_soil_depth/dz) # number of vertical grid size
Psi_ini  = np.linspace(-3,2,nz)
Psi      = Psi_ini

C_ini,K_ini,theta_ini = vanGenuchten(Psi_ini,alpha, theta_S, theta_R, n, m, Ksat)
 
C        = C_ini
K        = K_ini
theta    = theta_ini

plt.plot(Psi_ini, -dz*np.arange(nz))
#......... Model Stopping tolerance for convegence 
stop_tol = 0.01
max_iter = 100
Psimin   = 0.0001

#------------------------------------------------------------------------
# Boundary Conditions (constant head at both end)
#------------------------------------------------------------------------
Psi_top = 0.0 #PPT-ET
Psi_bot = 1.0

#------------------------------------------------------------------------
# The Main loop
#------------------------------------------------------------------------
time_steps = 100
diff_mp1=np.zeros((nz,max_iter))

psi_store =np.zeros((nz,time_steps))

for i in xrange(time_steps): 
    # Initialize the Picard iteration solver
    
    Psi_n = Psi
    C_n,K_n,theta_n = vanGenuchten(Psi_n,alpha, theta_S, theta_R, n, m, Ksat)
    Psi_np1m = Psi_n
    
    loop_num = 0
    deltam   = 10000.0
    while deltam > stop_tol and loop_num<max_iter:
        loop_num = loop_num+1
        print loop_num
        C_np1m,K_np1m,theta_np1m = vanGenuchten(Psi_np1m,alpha, theta_S, theta_R, n, m, Ksat)
        K_np1m_ph,K_np1m_mh =  Khalf(nz,K_np1m)
            
        A1d, Psi_np1mp1 = ImplicitSolver_Psi_np1mp1(Psi_n,Psi_np1m,theta_n, theta_np1m,C_np1m,K_np1m_ph,K_np1m_mh,nz,Psi_top,Psi_bot)
        
        diff_mp1[:,loop_num-1] = (Psi_np1mp1 - Psi_np1m)
        deltam = np.max(np.abs(Psi_np1mp1 - Psi_np1m))

        Psi_np1m = Psi_np1mp1
#        plt.plot(Psi_np1m, -dz*np.arange(nz))
    psi_store[:,i]= Psi_np1mp1
    Psi = Psi_np1mp1
    
    plt.plot(Psi, -dz*np.arange(nz))
    
plt.figure()
plt.imshow(psi_store, interpolation='bilinear')
plt.colorbar()
plt.show()
    
#    Psi_np1 = 
#    Psi = Psi_np1
#   C,K,theta = vanGenuchten(Psi,alpha, theta_S, theta_R, n, m, Ksat)





