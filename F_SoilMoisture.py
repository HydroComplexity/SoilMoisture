# -*- coding: utf-8 -*-
"""
Created on Thu Feb 23 15:47:54 2017

@author: Qina and Phong

This code is to simulate soil  moisture in 1D soil column

The current version has:
1. uniform grid size
2. uniform soil properties, e.g. builk density \rho_s, Ksat, 
3. BC: Constant head at two ends
"""

from parameters import *
from scipy.sparse import diags
from scipy.sparse import dia_matrix
from scipy.sparse.linalg import spsolve

def vanGenuchten(Psi,alpha, theta_S, theta_R, n, m, Ksat):

    #Convert unit of h from [m] to [cm] to match with alpha
    Psi = Psi*100.0
    
    #Compute the volumetric moisture content [eqn 21]
    theta = theta_R + (theta_S - theta_R)/(1.0 + (alpha*np.abs(Psi))**n)**m
    
    theta[Psi>=0.0] = theta_S
    
    #Compute the effective saturation [eqn 2]
    Se = ((theta - theta_R)/(theta_S - theta_R))
     
    #Compute the hydraulic conductivity [eqn 8]
    K = Ksat*Se**(0.5)*(1.0 - (1.0 - Se**(1.0/m))**m)**2
     
    #Compute the specific moisture storage: 
    #C = d(theta)/dh
    #C = alpha*n*(1.0/n - 1.0)*(alpha*np.abs(Psi))**(n - 1.0)*(theta_R - theta_S)*((alpha*np.abs(Psi))**n + 1.0)**(1.0/n - 2.0)*100.0
    C = -100.0/Psi*m*(theta_S - theta_R)*Se*(1-Se**(1.0/m))/(1.0-m)
    C[Psi>=0.0] = 0.0
    
    return C,K,theta
    
def Khalf(nz,K):
    
    # K_i+1/2, Kp1
    diag_p = np.ones(nz,dtype = 'int')
    diag_p[nz-1] = 2
    A_p = diags( [diag_p,1], [0,1],shape=(nz, nz))
    K_ph = A_p.dot(K)/2.0 
    
    # K_i-1/2, Km1
    diag_m = np.ones(nz,dtype = 'int')
    diag_m[0] = 2
    A_m = diags([diag_m,1], [0,-1],shape=(nz,nz))
    K_mh = A_m.dot(K)/2.0

    return K_ph,K_mh
    

def ImplicitSolver_Psi_npimp1(Psi_n,Psi_np1m,theta_n, theta_np1m,C_np1m,K_np1m_ph,K_np1m_mh,nz,Psi_top,Psi_bot ):
    
    # Build the diagal and up,low diagnal for matrix A, where Ax = b
    Mdia = Ss*theta_np1m/poros/dt + C_np1m/dt+(K_np1m_ph+K_np1m_mh)/dz**2
    Mdia[0] = 1.0
    Mdia[nz-1] = 1.0
    Udia = -K_np1m_ph/dz**2
    Udia[0] = 0.0
    Ldia = -K_np1m_mh/dz**2   
    Ldia[nz-1] = 0.0
    
    # known values on the RHS
    RHS = Ss*theta_np1m/poros*Psi_n/dt+(theta_n-theta_np1m)/dt+C_np1m*Psi_np1m/dt + (K_np1m_ph-K_np1m_mh)/dz
    RHS[0] = Psi_top 
    RHS[nz-1] = Psi_bot
    
    # Do not change after this line
    Ldia_r  = np.roll(Ldia,-1)
    Udia_r  = np.roll(Udia,1)
    
    ### build the matrix A, Ax = b
    A1d = dia_matrix( ([Ldia_r,Mdia,Udia_r],[-1,0,1]),shape=(nz, nz))	#
    A1d = A1d.tocsr()
     
    ### Solve x, Ax = b
    Psi_np1mp1 = spsolve(A1d,RHS)
    
    return Psi_np1mp1
    



    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    


