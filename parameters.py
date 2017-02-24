#Parameters
import numpy as np

#####################################
## soil moisture 
######################################
total_soil_depth = 2.0 	# [m]
dt = 1.0/60.0 			#
dz = 0.2 			   	# [m]
poros = 0.4       		# soil porosity
Ksat=6.94*1e-2*60.0 	# Sat. hydr.conductivity [m/hr]
alpha   = 0.01          # parameter related to the inverse of the air entry suction [1/cm]
theta_S = poros         # Saturated water content [-]
theta_R = 0.08          # Residual water content [-]
n       = 2.0           # Pore-size distributions [-]
m       = 1.0-1.0/n;    
Ss = 5.0*1e-4 			# Specific storage
