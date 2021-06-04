'''
Created on Sep 19, 2011 

@Author: Huynh Anh Nguyen LE and Soojong PAK
'''

# import math

from math import *
from numpy import *
import numpy as np

# from scipy import *

ETC_version = 'IGRINS-ETC v2.1.1'

# Constant 
C = 299792500           # [m s-1]
h = 6.6260755E-34         # [J s]

# Parameters: Basic
D_telescope   = 2.7     # 2.7m for the full aperture. 
#D_telescope_e = 21.5    # 21.5m for the effective aperture
T_ambient     = 275     # [K]
PI = math.pi
W_slit = 1.0
n_slit = 3.66

R = [40000, 40000]   # [H-band, K-band]

# Read noise
n_read = 5               # [electron, Read Noise, IR Array parameters]
# dark current
n_dark = 0.02            # [electron sec-1]


#================================================
# Reading Telluric absorption and emission data table

print('>>>>>>>>>>>>>>> ' + ETC_version + ' <<<<<<<<<<<<<<<')

# Tau atmosphere
Tau_atmosphere_band = [0.95, 0.95]  

print('...... Reading Le_IGRINS_AM1_H_40000.dat ......')
ATMO_LE_H = genfromtxt("Le_IGRINS_AM1_H_40000.dat")
print('...... Reading Le_IGRINS_AM1_K_40000.dat ......')
ATMO_LE_K = genfromtxt("Le_IGRINS_AM1_K_40000.dat")

print('...... Reading Kim_IGRINS_OH_H_40000.dat ......')
OH_H = genfromtxt("Kim_IGRINS_OH_H_40000.dat")
print('...... Reading Kim_IGRINS_OH_K_40000.dat ......')
OH_K = genfromtxt("Kim_IGRINS_OH_K_40000.dat") 
#================================================

# Parameters: Background emission

# [um]
wave_band  	= [1.63,	2.22]
delta_wave	= [0.35,	0.45]
band_min 	= [1.4,	1.9]
band_max 	= [1.9,	2.5]

# Parameters: Convolution calibrated emission line by Gauss function with R=40000
sigma = [4.13E-05 / 2.35482, 5.50E-05 / 2.35482]
delta = [1.00E-05, 1.50E-05]

# [W m-2 Hz-2], Allen's AQ (p150)
S_ZM        = [1.05E-23,   6.55E-24]
# [ph s-1 m-2 arcsec-2], Zodiacal Light
phi_zod     = [77.,        0.]

# phi OH
# [ph s-1 m-2 arcsec-2], Use H-band data for J and K-bands.
phi_OH_band = [7400.,      7400.] 

			  
#================================================
# Parameters: in Cal_emissivity_total.py
Tau_M_1                  = 0.95
Tau_M_1_e                = [0.55, 0.55]                 # Tau_M_1 * np.square(D_telescope_e / D_telescope)
Tau_M_2                  = [1.00, 1.00]
Tau_M_3                  = [1.00, 1.00]
Tau_Window               = [0.95, 0.95]                 # (and AO)
Tau_slit_loss            = [0.64, 0.64]	        # only for point source
Pupil_Stop = [0.90, 0.90]
Doublet    = [0.96, 0.96]
Collimator = [0.58, 0.58]
Immersion  = [0.88, 0.88]
X_disp     = [0.80, 0.75]
camera     = [0.89, 0.89]
detector   = [0.80, 0.80]

Tau_Optics = np.zeros(2)
for i in range(0,2):
    Tau_Optics[i] = Tau_M_1_e[i]*Tau_Window[i]*Tau_slit_loss[i]*Pupil_Stop[i] \
	*Doublet[i]*Collimator[i]*Immersion[i]*X_disp[i]*camera[i]*detector[i]

#Tau_Optics               = [0.08399639701094, 0.07874662219776]









