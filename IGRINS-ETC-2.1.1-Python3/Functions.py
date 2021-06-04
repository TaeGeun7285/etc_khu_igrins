"""
Updated on 2013.05.21 by Huynh Anh and Soojong

@Author: Huynh Anh Nguyen LE, Jongmin LEE, Soojong PAK
"""

from math import *

from numpy import *
import numpy as np
from numpy import convolve

from pylab import *
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter, OldScalarFormatter

from Constants  import *


# Function calculations

#===============================================================================

# This function calculates single Signal-to-Noise values from input magnitude.
def Mode_SN_band(t_exp, n_exp, K_template, T_template, lineflux, linewidth, restwave, vshift):
    mag = Cal_Template_Mag_band(K_template, T_template)
    SN = np.zeros(2)
    for BAND in np.arange(0,2):
        y = Cal_SN(Tau_atmosphere_band[BAND], wave_band[BAND], phi_OH_band[BAND], mag[BAND], t_exp, n_exp, \
                   Cal_LSignal(lineflux, wave_band[BAND], linewidth, restwave, vshift), BAND)
        SN[BAND] = round(y, 2)
    return SN

# Test with calibrated sginal
def Mode_MS_band(t_exp, n_exp, K_template, T_template):
    mag = Cal_Template_Mag_band(K_template, T_template)
    SN = np.zeros(2)
    for BAND in np.arange(0,2):
        #y = Cal_SN(Tau_atmosphere_band[BAND], wave_band[BAND], phi_OH_band[BAND], mag[BAND], t_exp, n_exp, BAND)
        y = Cal_MSignal(Tau_atmosphere_band[BAND], wave_band[BAND], mag[BAND], t_exp , n_exp, BAND)
        SN[BAND] = y #round(y, 2)
    return SN

   
# This function calculate Signal-to-Noise vs. Magnitude values.
def Mode_SN_Mag_Plot(t_exp, n_exp, mag_max, mag_min, lineflux, linewidth, restwave, vshift):
    mag = np.arange(mag_min, mag_max + 0.5, 0.5)
    S_N = [np.arange(mag_min, mag_max + 0.5, 0.5), np.arange(mag_min, mag_max + 0.5, 0.5)]

    for BAND in range(0, 2):
        for i in range(0, len(mag)):
            S_N[BAND][i] = Cal_SN(Tau_atmosphere_band[BAND], wave_band[BAND], phi_OH_band[BAND], mag[i], t_exp , n_exp, \
                           Cal_LSignal(lineflux, wave_band[BAND], linewidth, restwave, vshift), BAND)
    #Plot
    Plot_SN_Mag(mag, mag_min, mag_max, S_N[0], S_N[1], n_exp, t_exp)
	
# This function calculate Signal-to-Noise vs wavelength from "Get_Tau_atmo()" function.	
def Mode_SN_wave_Plot(t_exp, n_exp, pwv_type, min_x, max_x, K_template, T_template, lineflux, linewidth, restwave, vshift):

    if min_x >= band_min[0] and max_x <= band_max[0]:
        Tau_atmosphere = Get_Tau_atmo(ATMO_LE_H[:,0], ATMO_LE_H[:,1], ATMO_LE_H[:,2], ATMO_LE_H[:,3], pwv_type)
        S_N = np.zeros(len(ATMO_LE_H[:,0]))
        for i in range(0,len(ATMO_LE_H[:,0])):
            BAND = 0
            S_N[i] = Cal_SN(Tau_atmosphere[i], ATMO_LE_H[:,0][i], OH_H[:,1][i], \
            Cal_Template_Mag_wave(ATMO_LE_H[:,0][i], K_template, T_template, BAND), t_exp , n_exp, \
                            Cal_LSignal(lineflux, ATMO_LE_H[:,0][i], linewidth, restwave, vshift), BAND)

    if min_x >= band_min[1] and max_x <= band_max[1]:
        Tau_atmosphere = Get_Tau_atmo(ATMO_LE_K[:,0], ATMO_LE_K[:,1], ATMO_LE_K[:,2], ATMO_LE_K[:,3], pwv_type)
        S_N = np.zeros(len(ATMO_LE_K[:,0]))
        for i in range(0,len(ATMO_LE_K[:,0])):
            BAND = 1
            S_N[i] = Cal_SN(Tau_atmosphere[i], ATMO_LE_K[:,0][i], OH_K[:,1][i], \
            Cal_Template_Mag_wave(ATMO_LE_K[:,0][i], K_template, T_template, BAND), t_exp , n_exp, \
                            Cal_LSignal(lineflux, ATMO_LE_K[:,0][i], linewidth, restwave, vshift), BAND)
        
    #Plot
    if min_x >= band_min[0] and max_x <= band_max[0]:
        Title_para = ' (t=' + '%d' %t_exp + 'sec, N=' + '%d'%n_exp + ', PWV=' + '%d' % pwv_type + ')'
        Plot_SN_wave(ATMO_LE_H[:,0], S_N, min_x, max_x, Title_para) 	
    if min_x >= band_min[1] and max_x <= band_max[1]:
        Title_para = ' (t=' + '%d' %t_exp + 'sec, N=' + '%d'%n_exp + ', PWV=' + '%d' % pwv_type + ')'
        Plot_SN_wave(ATMO_LE_K[:,0], S_N, min_x, max_x, Title_para)	

# This function calculate measured-signal vs wavelength from "Get_Tau_atmo()" function.	
def Mode_MS_wave_Plot(t_exp, n_exp, pwv_type, min_x, max_x, K_template, T_template):

    if min_x >= band_min[0] and max_x <= band_max[0]:
        Tau_atmosphere = Get_Tau_atmo(ATMO_LE_H[:,0], ATMO_LE_H[:,1], ATMO_LE_H[:,2], ATMO_LE_H[:,3], pwv_type)
        M_S = np.zeros(len(ATMO_LE_H[:,0]))
        for i in range(0,len(ATMO_LE_H[:,0])):
            BAND = 0
            M_S[i] = Cal_MSignal(Tau_atmosphere[i], ATMO_LE_H[:,0][i], \
            Cal_Template_Mag_wave(ATMO_LE_H[:,0][i], K_template, T_template, BAND), t_exp , n_exp, BAND)
            
    if min_x >= band_min[1] and max_x <= band_max[1]:
        Tau_atmosphere = Get_Tau_atmo(ATMO_LE_K[:,0], ATMO_LE_K[:,1], ATMO_LE_K[:,2], ATMO_LE_K[:,3], pwv_type)
        M_S = np.zeros(len(ATMO_LE_K[:,0]))
        for i in range(0,len(ATMO_LE_K[:,0])):
            BAND = 1
            M_S[i] = Cal_MSignal(Tau_atmosphere[i], ATMO_LE_K[:,0][i], \
            Cal_Template_Mag_wave(ATMO_LE_K[:,0][i], K_template, T_template, BAND), t_exp , n_exp, BAND)
        
    #Plot
    if min_x >= band_min[0] and max_x <= band_max[0]:
        Title_para = ' (t=' + '%d' %t_exp + 'sec, N=' + '%d'%n_exp + ', PWV=' + '%d' % pwv_type + ')'
        Plot_MS_wave(ATMO_LE_H[:,0], M_S, min_x, max_x, Title_para) 	
    if min_x >= band_min[1] and max_x <= band_max[1]:
        Title_para = ' (t=' + '%d' %t_exp + 'sec, N=' + '%d'%n_exp + ', PWV=' + '%d' % pwv_type + ')'
        Plot_MS_wave(ATMO_LE_K[:,0], M_S, min_x, max_x, Title_para)	

# This function calculate measured-signal-lines vs wavelength from "Get_Tau_atmo()" function.	
def Mode_MSL_wave_Plot(t_exp, n_exp, pwv_type, min_x, max_x, K_template, T_template, lineflux, linewidth, restwave, vshift):

    if min_x >= band_min[0] and max_x <= band_max[0] and restwave >= band_min[0] and restwave <= band_max[0] :
        Tau_atmosphere = Get_Tau_atmo(ATMO_LE_H[:,0], ATMO_LE_H[:,1], ATMO_LE_H[:,2], ATMO_LE_H[:,3], pwv_type)
        M_S = np.zeros(len(ATMO_LE_H[:,0]))
        M_N = np.zeros(len(ATMO_LE_H[:,0]))
        R_N = np.zeros(len(ATMO_LE_H[:,0]))

        for i in range(0,len(ATMO_LE_H[:,0])):
            BAND = 0
            M_S[i] = Cal_MSignal(Tau_atmosphere[i], ATMO_LE_H[:,0][i], \
            Cal_Template_Mag_wave(ATMO_LE_H[:,0][i], K_template, T_template, BAND), t_exp , n_exp, BAND)

            M_N[i] = Cal_MNoise(Tau_atmosphere[i], ATMO_LE_H[:,0][i], OH_H[:,1][i], \
            Cal_Template_Mag_wave(ATMO_LE_H[:,0][i], K_template, T_template, BAND), t_exp , n_exp, \
                                Cal_LSignal(lineflux, ATMO_LE_H[:,0][i], linewidth, restwave, vshift), BAND)
            R_N[i] = Cal_RNoise(0, M_N[i])

        E_L =  Cal_LSignal(lineflux, ATMO_LE_H[:,0], linewidth, restwave, vshift)
        BAND = 0 
        E_L_C = convolve(E_L, f_gauss(sigma[BAND], delta[BAND]), mode = 1) * delta[BAND]
        S_N_L = Cal_SNL(ATMO_LE_H[:,0], M_S + E_L_C + R_N, linewidth, restwave, delta[BAND])

        print('S/N of The Emission Line = %d' %S_N_L)
        '''
        # Test S/N calculation emission line
        
        f = open('signal_noise_h_range.dat','w')
        lamda = ATMO_LE_H[:,0][(ATMO_LE_H[:,0]>=min_x) & (ATMO_LE_H[:,0]<=max_x)]
        signoise = M_S + E_L_C + R_N
        lamsig = signoise[(ATMO_LE_H[:,0]>=min_x) & (ATMO_LE_H[:,0]<=max_x)]
                
        A = np.arange(len(lamda))
        for i in A:
            temp = "%.6f %.6e \n" %(lamda[i], lamsig[i])
            f.write(temp)
        f.close()
        '''
    if min_x >= band_min[1] and max_x <= band_max[1] and restwave >= band_min[1] and restwave <= band_max[1]:
        Tau_atmosphere = Get_Tau_atmo(ATMO_LE_K[:,0], ATMO_LE_K[:,1], ATMO_LE_K[:,2], ATMO_LE_K[:,3], pwv_type)
        M_S = np.zeros(len(ATMO_LE_K[:,0]))
        M_N = np.zeros(len(ATMO_LE_K[:,0]))
        R_N = np.zeros(len(ATMO_LE_K[:,0]))

        for i in range(0,len(ATMO_LE_K[:,0])):
            BAND = 1
            M_S[i] = Cal_MSignal(Tau_atmosphere[i], ATMO_LE_K[:,0][i], \
            Cal_Template_Mag_wave(ATMO_LE_K[:,0][i], K_template, T_template, BAND), t_exp , n_exp, BAND)

            M_N[i] = Cal_MNoise(Tau_atmosphere[i], ATMO_LE_K[:,0][i], OH_K[:,1][i], \
            Cal_Template_Mag_wave(ATMO_LE_K[:,0][i], K_template, T_template, BAND), t_exp , n_exp, \
                                Cal_LSignal(lineflux, ATMO_LE_K[:,0][i], linewidth, restwave, vshift), BAND)
            R_N[i] = Cal_RNoise(0, M_N[i])

        E_L = Cal_LSignal(lineflux, ATMO_LE_K[:,0], linewidth, restwave, vshift)
        BAND = 1 
        E_L_C = convolve(E_L, f_gauss(sigma[BAND], delta[BAND]), mode = 1) * delta[BAND]
        S_N_L = Cal_SNL(ATMO_LE_K[:,0], M_S + E_L_C + R_N, linewidth, restwave, delta[BAND])

        print('S/N of The Emission Line = %d' %S_N_L)
        '''
        # Test S/N calculation emission line
        f = open('signal_noise_k_range.dat','w')
        lamda = ATMO_LE_K[:,0][(ATMO_LE_K[:,0]>=min_x) & (ATMO_LE_K[:,0]<=max_x)]
        signoise = M_S + E_L_C + R_N
        lamsig = signoise[(ATMO_LE_K[:,0]>=min_x) & (ATMO_LE_K[:,0]<=max_x)]
               
        A = np.arange(len(lamda))
        for i in A:
            temp = "%.6f %.6e \n" %(lamda[i], lamsig[i])
            f.write(temp)
        f.close()
        #f = open('signal_noise_k.dat','w')
        #A = np.arange(len(ATMO_LE_K[:,0]))
        #for i in A:
        #    temp = "%.6e %.6e %.6e %.6e \n" %(ATMO_LE_K[:,0][i], M_S[i], E_L_C[i], M_N[i])
        #    f.write(temp)
        #f.close()
        '''
     
    #Plot
    if min_x >= band_min[0] and max_x <= band_max[0]:
        Title_para = ' (t=' + '%d' %t_exp + 'sec, N=' + '%d'%n_exp + ', PWV=' + '%d' % pwv_type + ')'
        Plot_MS_wave(ATMO_LE_H[:,0], M_S + E_L_C + R_N, min_x, max_x, Title_para) 	
    if min_x >= band_min[1] and max_x <= band_max[1]:
        Title_para = ' (t=' + '%d' %t_exp + 'sec, N=' + '%d'%n_exp + ', PWV=' + '%d' % pwv_type + ')'
        Plot_MS_wave(ATMO_LE_K[:,0], M_S + E_L_C + R_N, min_x, max_x, Title_para)	

# This function calculate measured-noise vs wavelength from "Get_Tau_atmo()" function.	
def Mode_MN_wave_Plot(t_exp, n_exp, pwv_type, min_x, max_x, K_template, T_template, lineflux, linewidth, restwave, vshift):

    if min_x >= band_min[0] and max_x <= band_max[0]:
        Tau_atmosphere = Get_Tau_atmo(ATMO_LE_H[:,0], ATMO_LE_H[:,1], ATMO_LE_H[:,2], ATMO_LE_H[:,3], pwv_type)
        M_N = np.zeros(len(ATMO_LE_H[:,0]))
        for i in range(0,len(ATMO_LE_H[:,0])):
            BAND = 0
            M_N[i] = Cal_MNoise(Tau_atmosphere[i], ATMO_LE_H[:,0][i], OH_H[:,1][i], \
            Cal_Template_Mag_wave(ATMO_LE_H[:,0][i], K_template, T_template, BAND), t_exp , n_exp, \
                                Cal_LSignal(lineflux, ATMO_LE_H[:,0][i], linewidth, restwave, vshift), BAND)
        #print 'stdev noise', np.std(M_N)            

    if min_x >= band_min[1] and max_x <= band_max[1]:
        Tau_atmosphere = Get_Tau_atmo(ATMO_LE_K[:,0], ATMO_LE_K[:,1], ATMO_LE_K[:,2], ATMO_LE_K[:,3], pwv_type)
        M_N = np.zeros(len(ATMO_LE_K[:,0]))
        for i in range(0,len(ATMO_LE_K[:,0])):
            BAND = 1
            M_N[i] = Cal_MNoise(Tau_atmosphere[i], ATMO_LE_K[:,0][i], OH_K[:,1][i], \
            Cal_Template_Mag_wave(ATMO_LE_K[:,0][i], K_template, T_template, BAND), t_exp , n_exp, \
                                Cal_LSignal(lineflux, ATMO_LE_K[:,0][i], linewidth, restwave, vshift), BAND)

        #print 'stdev noise', np.std(M_N)
    #Plot
    if min_x >= band_min[0] and max_x <= band_max[0]:
        Title_para = ' (t=' + '%d' %t_exp + 'sec, N=' + '%d'%n_exp + ', PWV=' + '%d' % pwv_type + ')'
        Plot_MN_wave(ATMO_LE_H[:,0], M_N, min_x, max_x, Title_para) 	
    if min_x >= band_min[1] and max_x <= band_max[1]:
        Title_para = ' (t=' + '%d' %t_exp + 'sec, N=' + '%d'%n_exp + ', PWV=' + '%d' % pwv_type + ')'
        Plot_MN_wave(ATMO_LE_K[:,0], M_N, min_x, max_x, Title_para)

    
# This function calculate limiting mangitude from input signal-to-noise value
def Mode_Mag_Limiting_band(t_exp, n_exp, S_N):
    magnitude = np.zeros(2)
    for BAND in np.arange(0,2):
        y = Cal_Mag_Limit(Tau_atmosphere_band[BAND], wave_band[BAND], phi_OH_band[BAND], S_N, t_exp , n_exp, BAND)
        magnitude[BAND] = round(y,2)
    return magnitude

def Cal_Mag_Limit(Tau_atmosphere, wave_band, phi_OH_band, S_N_value, t_exp , n_exp, BAND):
    Tau_band = Cal_Tau_total(Tau_atmosphere, BAND)
    Emissivity_band = Cal_emissivity_total(Tau_atmosphere, Tau_band, BAND)
    n_zod_OH_band = Cal_n_zod_OH(wave_band, Tau_band, phi_OH_band, BAND)
    B_band = Cal_Background(wave_band, T_ambient)
    n_thermal_band = Cal_n_thermal(Emissivity_band, B_band, BAND)
    sig_band = (pow(S_N_value,2) + np.sqrt(pow(S_N_value,4) \
            + 4 * pow(S_N_value,2) * (2 * pow(n_slit,2) * n_exp \
            * (t_exp * (n_zod_OH_band + n_thermal_band + n_dark)+ pow(n_read,2))))) / 2
    mag_band =  (np.log10(((t_exp * n_exp) *(PI * pow((D_telescope/2),2)) \
            * Tau_band * S_ZM[BAND]) / ((sig_band * (h * R[BAND]))))) / 0.4
    return mag_band

# This is Cal_Signal-to-Noise fucntion.
# This function is called "Mode_SN_band()", "Plot()", Plot_wavelength()" fucntions.
def Cal_SN(Tau_atmosphere, wave_um, phi_OH_band, mag, t_exp , n_exp, linesignal, BAND):
    Tau_band = Cal_Tau_total(Tau_atmosphere, BAND)
    Emissivity_band = Cal_emissivity_total(Tau_atmosphere, Tau_band, BAND)
    n_zod_OH_band = Cal_n_zod_OH(wave_um, Tau_band, phi_OH_band, BAND)
    B_band = Cal_Background(wave_um, T_ambient)
    n_thermal_band = Cal_n_thermal(Emissivity_band, B_band, BAND)
    y_signal = Cal_Signal(mag, Tau_band, t_exp, n_exp, BAND)
    #print 'signal', y_signal
    y_line = (linesignal *  (n_exp * t_exp * (PI * pow((D_telescope/2),2)) * Tau_band * pow(wave_um,2) \
                                     *  pow(10, -6))) / (h * C * R[i])
    y_signal_c = y_signal + y_line
    y_noise  = Cal_Noise(n_zod_OH_band, n_thermal_band, y_signal_c, t_exp, n_exp)
    y = y_signal_c / y_noise
    return y

# This function calcalates atmostphere transmission vs. wavelength.
# This function is called by "Plot_wavelength" function.
def Get_Tau_atmo(wavelength, transmission2, transmission4, transmission8, pwv):
    N_data = len(wavelength)
    y = np.zeros(N_data)
    #if pwv >= 1 and pwv <= 2:
    #    for i in np.arange(0, N_data):
#	    y[i] = transmission1[i] + (pwv - 1) * (transmission2[i] - transmission1[i])
    if pwv >= 2 and pwv <= 4:
        for i in np.arange(0, N_data):
            y[i] = transmission2[i] + (pwv - 2) * (transmission4[i] - transmission2[i]) / (4 - 2)
    elif pwv > 4 and pwv <= 8:
        for i in np.arange(0, N_data):
            y[i] = transmission4[i] + (pwv - 4) * (transmission8[i] - transmission4[i]) / (8 - 4)
    return y

# This fucntion plots Signal-to-Noise vs. Magnitude plots.
# This fucntion is called by "Plot()" function.
def Plot_SN_Mag(mag, mag_min, mag_max, S_N_H, S_N_K, n_exp, t_exp):
    plt.figure(num=None,figsize=(12, 8), dpi=80, facecolor='w', edgecolor='k')
    ax = plt.subplot(111)
    ax.plot(mag, S_N_H, 'm', mag, S_N_K, 'g', linewidth = 3)
    plt.title(ETC_version + ' (t=' + '%d' % t_exp + 's, N=' + '%d' % n_exp + ')',fontsize=16)
    plt.xlabel('Target Vega Magnitude', fontsize=15)
    plt.ylabel('Signal-to-Noise',fontsize=15)
    plt.legend(['H', 'K'])

    locs, labels = xticks()
    plt.setp(labels, 'fontsize', 'large')
    locs, labels = yticks()
    plt.setp(labels,'fontsize','large')
    '''
    x_majorLocator = MultipleLocator(1)
    x_majorFormatter = FormatStrFormatter('%d')
    x_minorLocator = MultipleLocator(0.5)
    y_majorLocator = MultipleLocator(10)
    y_majorFormatter = FormatStrFormatter('%0.1f')
    y_minorLocator = MultipleLocator(1)

    ax.xaxis.set_major_locator(x_majorLocator)
    ax.xaxis.set_major_formatter(x_majorFormatter)
    ax.xaxis.set_minor_locator(x_minorLocator)
    ax.yaxis.set_major_locator(y_majorLocator)
    ax.yaxis.set_major_formatter(y_majorFormatter)
    ax.yaxis.set_minor_locator(y_minorLocator)
    '''
    ax.set_yscale('log')
    ax.axis([mag_min, mag_max, 1, 1000])
    ax.grid(color='k', linestyle='-', which='minor', linewidth=0.5)
    ax.grid(color='k', linestyle='-', which='major', linewidth=1)
    plt.show()

# This function plots Signal-to-Noise vs. wavelength plots.
# This function is called by "Plot_wavelength()" fucntion.
def Plot_SN_wave(wavelength, S_N, min_x, max_x, Title_para):
    plt.figure(num=None,figsize=(12, 8), dpi=80, facecolor='w', edgecolor='k')
    ax = plt.subplot(111)
    ax.plot(wavelength, S_N, 'k', linewidth = 1)
    plt.xlabel('Wavelength (microns)', fontsize=15)
    plt.ylabel('S/N per Resolution Element', fontsize=15)
    plt.title(ETC_version + ' ' + Title_para, fontsize=15)

    locs, labels = xticks()
    plt.setp(labels, 'fontsize', 'large')
    locs, labels = yticks()
    plt.setp(labels,'fontsize','large')

    #x_majorLocator = MultipleLocator((max_x - min_x)/5)

    if (max_x - min_x) >= 0.1:
        x_majorFormatter = FormatStrFormatter('%0.2f')
    if (max_x - min_x) >= 0.01 and (max_x - min_x) < 0.1 :
        x_majorFormatter = FormatStrFormatter('%0.3f')
    if (max_x - min_x) >= 0.001 and (max_x - min_x) < 0.01 :
        x_majorFormatter = FormatStrFormatter('%0.4f')
    if (max_x - min_x) >= 0.0001 and (max_x - min_x) < 0.001 :
        x_majorFormatter = FormatStrFormatter('%0.5f')
    if (max_x - min_x) >= 0.00001 and (max_x - min_x) <= 0.0001 :
        x_majorFormatter = FormatStrFormatter('%0.6f')
    if (max_x - min_x) >= 0.000001 and (max_x - min_x) <= 0.00001:
        x_majorFormatter = FormatStrFormatter('%0.7f')
    if (max_x - min_x) >= 0.0000001 and (max_x - min_x) <= 0.000001:
        x_majorFormatter = FormatStrFormatter('%0.8f')
    if (max_x - min_x)  < 0.0000001:
        x_majorFormatter = FormatStrFormatter('%0.9f')
    '''   
    x_minorLocator = MultipleLocator((max_x - min_x)/25)
    y_majorLocator = MultipleLocator((max(log10(S_N)) - min(log10(S_N)))/10)
    y_majorFormatter = FormatStrFormatter('%0.1f')
    y_minorLocator = MultipleLocator((max(log10(S_N)) - min(log10(S_N)))/50)

    ax.xaxis.set_major_locator(x_majorLocator)
    ax.xaxis.set_major_formatter(x_majorFormatter)
    ax.xaxis.set_minor_locator(x_minorLocator)
    ax.yaxis.set_major_locator(y_majorLocator)
    ax.yaxis.set_major_formatter(y_majorFormatter)
    ax.yaxis.set_minor_locator(y_minorLocator)
    '''
    #ax.set_yscale('log')
    #if max(S_N) >= 20:
    ax.axis([min_x, max_x, 0, max(S_N) + 0.1 * max(S_N)])
    #else:
    #ax.axis([min_x, max_x, 0, max(log10(S_N)) + 0.05 * max(log10(S_N))])
    #ax.grid(color='k', linestyle='-', which='minor', linewidth=0.5)
    #ax.grid(color='k', linestyle='-', which='major', linewidth=1)
    plt.show()

def Plot_MS_wave(wavelength, S_N, min_x, max_x, Title_para):
    plt.figure(num=None,figsize=(12, 8), dpi=80, facecolor='w', edgecolor='k')
    ax = plt.subplot(111)
    ax.plot(wavelength, S_N, 'k', linewidth = 1)
    plt.xlabel('Wavelength (microns)', fontsize=15)
    plt.ylabel('Calibrated-Signal (W m-2 um-1)', fontsize=15)
    plt.title(ETC_version + ' ' + Title_para, fontsize=15)

    locs, labels = xticks()
    plt.setp(labels, 'fontsize', 'large')
    locs, labels = yticks()
    plt.setp(labels,'fontsize','large')

    #x_majorLocator = MultipleLocator((max_x - min_x)/5)
    
    if (max_x - min_x) >= 0.1:
        x_majorFormatter = FormatStrFormatter('%0.2f')
    if (max_x - min_x) >= 0.01 and (max_x - min_x) < 0.1 :
        x_majorFormatter = FormatStrFormatter('%0.3f')
    if (max_x - min_x) >= 0.001 and (max_x - min_x) < 0.01 :
        x_majorFormatter = FormatStrFormatter('%0.4f')
    if (max_x - min_x) >= 0.0001 and (max_x - min_x) < 0.001 :
        x_majorFormatter = FormatStrFormatter('%0.5f')
    if (max_x - min_x) >= 0.00001 and (max_x - min_x) <= 0.0001 :
        x_majorFormatter = FormatStrFormatter('%0.6f')
    if (max_x - min_x) >= 0.000001 and (max_x - min_x) <= 0.00001:
        x_majorFormatter = FormatStrFormatter('%0.7f')
    if (max_x - min_x) >= 0.0000001 and (max_x - min_x) <= 0.000001:
        x_majorFormatter = FormatStrFormatter('%0.8f')
    if (max_x - min_x)  < 0.0000001:
        x_majorFormatter = FormatStrFormatter('%0.9f')
    '''   
    x_minorLocator = MultipleLocator((max_x - min_x)/25)
    y_majorLocator = MultipleLocator((max(S_N) - min(S_N))/5)
    y_majorFormatter = FormatStrFormatter('%0.1f')
    y_minorLocator = MultipleLocator((max(S_N) - min(S_N))/25)

    ax.xaxis.set_major_locator(x_majorLocator)
    ax.xaxis.set_major_formatter(x_majorFormatter)
    ax.xaxis.set_minor_locator(x_minorLocator)
    ax.yaxis.set_major_locator(y_majorLocator)
    ax.yaxis.set_major_formatter(y_majorFormatter)
    ax.yaxis.set_minor_locator(y_minorLocator)
    '''
    #ax.xaxis.set_major_formatter(OldScalarFormatter())
    #ax.yaxis.set_major_formatter(OldScalarFormatter())
    #ax.set_yscale('log')
    
    #ax.axis([min_x, max_x, min(S_N), max(S_N) + 0.1 * max(S_N)])
    ax.axis([min_x, max_x, min(S_N[(wavelength >= min_x) & (wavelength <= max_x)]) - 0.1 * \
             min(S_N[(wavelength >= min_x) & (wavelength <= max_x)]), 
             max(S_N[(wavelength >= min_x) & (wavelength <= max_x)]) + 0.1 * \
             max(S_N[(wavelength >= min_x) & (wavelength <= max_x)])])
    
    #ax.grid(color='k', linestyle='-', which='minor', linewidth=0.5)
    #ax.grid(color='k', linestyle='-', which='major', linewidth=1)
    plt.show()

def Plot_MN_wave(wavelength, S_N, min_x, max_x, Title_para):
    plt.figure(num=None,figsize=(12, 8), dpi=80, facecolor='w', edgecolor='k')
    ax = plt.subplot(111)
    ax.plot(wavelength, S_N, 'k', linewidth = 1)
    #ax.plot(wavelength, S_N_C, 'r', linewidth = 1)
    plt.xlabel('Wavelength (microns)', fontsize=15)
    plt.ylabel('Calibrated-Noise (W m-2 um-1)', fontsize=15)
    plt.title(ETC_version + ' ' + Title_para, fontsize=15)

    locs, labels = xticks()
    plt.setp(labels, 'fontsize', 'large')
    locs, labels = yticks()
    plt.setp(labels,'fontsize','large')
    
    #x_majorLocator = MultipleLocator((max_x - min_x)/5)
    
    if (max_x - min_x) >= 0.1:
        x_majorFormatter = FormatStrFormatter('%0.2f')
    if (max_x - min_x) >= 0.01 and (max_x - min_x) < 0.1 :
        x_majorFormatter = FormatStrFormatter('%0.3f')
    if (max_x - min_x) >= 0.001 and (max_x - min_x) < 0.01 :
        x_majorFormatter = FormatStrFormatter('%0.4f')
    if (max_x - min_x) >= 0.0001 and (max_x - min_x) < 0.001 :
        x_majorFormatter = FormatStrFormatter('%0.5f')
    if (max_x - min_x) >= 0.00001 and (max_x - min_x) <= 0.0001 :
        x_majorFormatter = FormatStrFormatter('%0.6f')
    if (max_x - min_x) >= 0.000001 and (max_x - min_x) <= 0.00001:
        x_majorFormatter = FormatStrFormatter('%0.7f')
    if (max_x - min_x) >= 0.0000001 and (max_x - min_x) <= 0.000001:
        x_majorFormatter = FormatStrFormatter('%0.8f')
    if (max_x - min_x)  < 0.0000001:
        x_majorFormatter = FormatStrFormatter('%0.9f')
    '''  
    x_minorLocator = MultipleLocator((max_x - min_x)/25)
    y_majorLocator = MultipleLocator((max(S_N) - min(S_N))/5)
    y_majorFormatter = FormatStrFormatter('%0.1f')
    y_minorLocator = MultipleLocator((max(S_N) - min(S_N))/25)

    ax.xaxis.set_major_locator(x_majorLocator)
    ax.xaxis.set_major_formatter(x_majorFormatter)
    ax.xaxis.set_minor_locator(x_minorLocator)
    ax.yaxis.set_major_locator(y_majorLocator)
    ax.yaxis.set_major_formatter(y_majorFormatter)
    ax.yaxis.set_minor_locator(y_minorLocator)
    '''
    #ax.xaxis.set_major_formatter(OldScalarFormatter())
    #ax.yaxis.set_major_formatter(OldScalarFormatter())
    #ax.set_yscale('log')
    
    #ax.axis([min_x, max_x, min(S_N) - 0.25 * min(S_N), max(S_N) + 0.25 * max(S_N)])
    ax.axis([min_x, max_x, min(S_N) - 0.25 * min(S_N), max(S_N[(wavelength >= min_x) & (wavelength <= max_x)]) + \
             0.1 * max(S_N[(wavelength >= min_x) & (wavelength <= max_x)])])
    #ax.axis([min_x, max_x, min(S_N) - 0.25 * min(S_N), 1e-12])
    #ax.grid(color='k', linestyle='-', which='minor', linewidth=0.5)
    #ax.grid(color='k', linestyle='-', which='major', linewidth=1)
    plt.show()

# This fucntion calculates Template magnitude values.
# This function is called by "Mode_SN_band()", "Plot_wavelegnth()" function.
def Cal_Template_Mag_band(K_template, T_template):
    mag = np.zeros(2)
    Omega_template = S_ZM[1] / Cal_Background(wave_band[1], T_template)
    mag = np.zeros(2)
    for BAND in np.arange(0,2):
        y_flux = Omega_template * Cal_Background(wave_band[BAND], T_template)
        y_mag  = K_template - 2.5*log10(y_flux/S_ZM[BAND])
        mag[BAND] = round(y_mag, 16)
    return mag

# This fucntion calculate Magnitude values vs wavelength.
# This function is called by "Mode_SN_band()", "Plot_wavelegnth()" function.
def Cal_Template_Mag_wave(x, K_template, T_template, i):
    Omega_template = S_ZM[1] / Cal_Background(wave_band[1], T_template)
    flux = Omega_template * Cal_Background(x, T_template)
    magnitude = K_template - 2.5*log10(flux/S_ZM[i])
    return magnitude

# This is Signal fucntion.
def Cal_Signal(mag, Tau_total, t_exp, n_exp, i):
    x = (t_exp * n_exp) *(PI * pow((D_telescope/2),2))* Tau_total * S_ZM[i] * pow(10,(-0.4 * mag))/ (h * R[i])   
    return x

# This is Noise function.
def Cal_Noise(t_n_zod, t_thermal, t_signal, t_exp, n_exp):
    x = np.sqrt(2 * pow(n_slit,2) * n_exp * (t_exp * (t_n_zod + t_thermal + n_dark)+ pow(n_read,2))+ t_signal)
    return x

# This function calculated for measure signal based al Cal_Signal() function

def Cal_MSignal(Tau_atmosphere, wave_um, mag, t_exp , n_exp, BAND):
    Tau_band = Cal_Tau_total(Tau_atmosphere, BAND)
    #print 'tau', Tau_band
    y_signal = Cal_Signal(mag, Tau_band, t_exp, n_exp, BAND) * pow(10, 6)
    y_m      = (h * C * R[BAND]) / ((n_exp * t_exp * (PI * pow((D_telescope/2),2)) * Tau_band * pow(wave_um,2)))
    y = y_signal * y_m 
    return y

# This function calculated for measure noise based al Cal_Noise() function

def Cal_MNoise(Tau_atmosphere, wave_um, phi_OH_band, mag, t_exp , n_exp, linesignal, BAND):
    Tau_band = Cal_Tau_total(Tau_atmosphere, BAND)
    Emissivity_band = Cal_emissivity_total(Tau_atmosphere, Tau_band, BAND)
    n_zod_OH_band = Cal_n_zod_OH(wave_um, Tau_band, phi_OH_band, BAND)
    B_band = Cal_Background(wave_um, T_ambient)
    n_thermal_band = Cal_n_thermal(Emissivity_band, B_band, BAND)
    y_signal = Cal_Signal(mag, Tau_band, t_exp, n_exp, BAND)
    #print 'signal', wave_um, y_signal
    y_linesignal   = (linesignal *  (n_exp * t_exp * (PI * pow((D_telescope/2),2)) * Tau_band * pow(wave_um,2) \
                                     *  pow(10, -6))) / (h * C * R[BAND])
    #print 'line signal', wave_um, y_linesignal
    y_signal_c = y_signal + y_linesignal
    y_noise  = Cal_Noise(n_zod_OH_band, n_thermal_band, y_signal_c, t_exp, n_exp) *  pow(10, 6)
    y_m      = (h * C * R[BAND]) / (n_exp * t_exp * (PI * pow((D_telescope/2),2)) * Tau_band * pow(wave_um,2))
    y = y_noise * y_m
    return y

# This function for test the measure noise with the intensity of OH lines based al Cal_Noise() function
def Cal_MNoise_test(Tau_atmosphere, wave_um, phi_OH_band, mag, t_exp , n_exp, BAND):
    Tau_band = Cal_Tau_total(Tau_atmosphere, BAND)
    Emissivity_band = Cal_emissivity_total(Tau_atmosphere, Tau_band, BAND)
    n_zod_OH_band = Cal_n_zod_OH(wave_um, Tau_band, phi_OH_band, BAND)
    B_band = Cal_Background(wave_um, T_ambient)
    n_thermal_band = Cal_n_thermal(Emissivity_band, B_band, BAND)
    y_signal = Cal_Signal(mag, Tau_band, t_exp, n_exp, BAND)
    y_noise  = Cal_Noise(n_zod_OH_band, n_thermal_band, y_signal, t_exp, n_exp) 
    y_m      = (h * C * R[BAND]) / (n_exp * t_exp * (PI * pow((D_telescope/2),2)) * Tau_band * pow(wave_um,2))
    y = y_noise * y_noise * y_m *  pow(10, 6)
    return y

def Cal_RNoise(signal, noise):
    s = np.random.normal(signal, noise, 1)
    return s

# This Gaussian function calculate for the emission lines
def GaussFunction(wave_um, linewidth, restwave):
    sigma = linewidth / 2.35482
    A = (wave_um - restwave) / sigma
    g = (1/(sigma*np.sqrt(2.0*np.pi)))* np.exp(-0.5 * A * A)
    return g

# This Gaussian function is used to do convolution calibrated emission line with R=40000
def f_gauss(sigma, delta):
    x = np.linspace(-0.001, 0.001, ((0.001 + 0.001) / delta) + 1)
    A = 1.0/(sigma*sqrt(2.0*pi))
    B = exp((-1.0/2.0)*((x/sigma)**2))
    f = A*B
    return f

# This function calculate for the calibrated signal emission lines based on GaussFunction()
def Cal_LSignal(lineflux, wave_um, linewidth, restwave, vshift):
    line = lineflux * GaussFunction(wave_um, linewidth, Doppler(restwave, vshift))
    #print 'emission', line
    return line

# This function calculate for doppler shift
def Doppler(restwave, vshift):
    x = ((restwave * vshift ) / (C * pow(10, -3))) + restwave
    return x

# This function calculate for signal-to-noise of the emission line
def Cal_SNL(wave_um, signal, linewidth, restwave, delta_wave):
    signal_line = signal[(wave_um >= (restwave - linewidth)) & (wave_um <= (restwave + linewidth))]
    continuum_left = signal[(wave_um >= (restwave - (3 * linewidth))) & (wave_um <= (restwave - linewidth))]
    continuum_right = signal[(wave_um >= (restwave + linewidth)) & (wave_um <= (restwave + (3 * linewidth)))]
    continuum = np.concatenate([continuum_left, continuum_right])
    aver_cont = np.average(continuum) 
    #print 'average continuum', aver_cont
    signal_flux = delta_wave * np.sum(signal_line - aver_cont)
    #print 'signal flux', signal_flux
    noise_stdev = np.std(continuum)
    #print 'stdev continuum', noise_stdev
    noise_flux = delta_wave * noise_stdev * np.sqrt(len(signal_line) * (1 + (1 / len(continuum))))
    #print 'noise flux', noise_flux
    SN = signal_flux / noise_flux
    print('Emission Line Range [um]: ' + '%0.6f' % (restwave - linewidth) + ' - ' + '%0.6f' % (restwave + linewidth))
    print('Continuum Range [um]: ' + '%0.6f' % (restwave - (3 * linewidth)) + ' - ' + '%0.6f' % (restwave - linewidth) \
          + ', ' + '%0.6f' % (restwave + linewidth) + ' - ' + '%0.6f' % (restwave + (3 * linewidth)))

    return SN
    

# This function calculate the photo-electrons from background emissions
# This function is called by "Signal_to_Noise()", "Magnitude()" functions.
def Cal_n_zod_OH(wave_band, Tau_total, phi_OH_band, BAND):
    # A Omega calculation
    # [arcsec], We are planning on around 85 mas for all wavelengths.
    # [pixel], 4 pixels in the spectral and spatial resolution elements
    # We assume that the image in the cross-dispersion direction has a FWHM equal to the slit width.
    A_Omega = PI * np.square(D_telescope/2) * np.square((W_slit/n_slit)/206265)  # [m^2 sr]
    if phi_OH_band > 5000:
        y = (wave_band / delta_wave[BAND] / R[BAND]) * Tau_total * (A_Omega * pow(206265,2)) * (phi_zod[BAND] + 0.2 * phi_OH_band)
    else:
        y_zod = (wave_band / delta_wave[BAND] / R[BAND]) * Tau_total * (A_Omega * pow(206265,2)) * (phi_zod[BAND])
        y_OH  = Tau_total * (A_Omega * pow(206265,2)) * (phi_OH_band)
        y = y_zod + y_OH
    return y

# This function calculated background emission.
# This function is called by "n_thermal()", "Mode_SN_band()", "Plot()", "Plot_wavelegnth()", "Signal_to_Noise()", "Magnitude()" fucntions.
def Cal_Background(wave_um, Temperature):
    # Background Emission Calculation
    # Units in bands, [um], from <gmtnirs_AO_strawman_070706b.doc>
    # Background for Zodiacal and OH emissions
    y = 1.4745E-50 * pow((C / (wave_um * 0.000001)),3)/ (exp(0.0000000000479922 * C / (wave_um * 0.000001)/ Temperature)-1)
    return y

# This function calculate the thermal-electrons.
# This function is called by "Signal_to_Noise()", "Mode_SN_band()", "Plot()", "Plot_wavelegnth()", "Magnitude()" fucntions.
def Cal_n_thermal(emissivity, background, BAND):
    A_Omega = PI * np.square(D_telescope/2) * np.square((W_slit/n_slit)/206265)  # [m^2 sr]
    y = (A_Omega * emissivity * background) / (h * R[BAND])
    return y

# This function calculate tau from optical parameters. 
# This function is called by "Signal_to_Noise()", "Magnitude()" functions.
# Tau_Optics[] values are calculated seperately (see Cal_Tau_Optics())
def Cal_Tau_total(Tau_atmosphere, BAND):
    #Total transmission for point source
	y = Tau_atmosphere * Tau_Optics[BAND]
	return y

# This function calculated emissivity from optical parameters
# This function is called by "Cal_n_thermal()", "Signal_to_Noise()", "Magnitude()" functions. 
# The emissivities at different wavelengths are different enough that this should be folded into the real calculator
def Cal_emissivity_total(Tau_atmosphere, Tau_total, BAND):
    #Total emissivity        
    Eta_M1 = 0.25
    Eta_M2 = 0.00
    Eta_M3 = 0.00 
    Eta_Window     = 0.05
    Eta_atmosphere = 0.05  
    y = (Eta_atmosphere + (Eta_M1 + (Eta_M2 + (Eta_M3 + Eta_Window / Tau_Window[BAND]) / Tau_M_3[BAND]) / Tau_M_2[BAND]) / Tau_M_1_e[BAND]) \
         / Tau_atmosphere * (Tau_total / Tau_slit_loss[BAND])
    return y 
