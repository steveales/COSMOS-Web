

# This program calculates the ratio of dust-to-stellar mass versus stellar
# mass for the stacking data from COSMOS Webb using a new method that 
# averages over each bin of redshift and stellar mass. This program
# uses the results from SIMSTACK, including the estimates
# of the effect of the redshift and mass errors.

# This version includes a correction for the effect of the CMB.
# This version also includes the T-z relation from Liang et al. (2019)
# and a two-temperature model from Bakx et al. (2018) 

# Last edited: 25th November 2025

import numpy as np
import matplotlib.pyplot as plt
from astropy.cosmology import Planck15 as cosmo
from lmfit import minimize,Parameters

# Constants
# ****************************************************************************

# CMB temperature

Tcmb = 2.7

# Dust mass opacity coefficient

kappa = 0.077

# constant for mass equation

con1 = 0.326

# Constant in modified blackbody

con_bb = 48.01

# Temperature and dust emissivity index in dust calculation

Tdust = 22.0
beta = 2.0

# Reference frequency in units of 10^12 Hz

freref = 0.353

# Set up paramaters for fit
# ****************************************************************************

params = Parameters()
params.add('a',value=-0.001,min=-0.5,max=0.0)
params.add('b',value = 0.018,min=0.0,max=0.5)

# Set up arrays and parameters 
# ****************************************************************************

Mlow = np.array([])
Mup = np.array([])
zlow = np.array([])
zup = np.array([])
zmed = np.array([])
S850 = np.array([])
Error = np.array([])
art_noise = np.array([])

# Function to calculate residuals
# ***************************************************************************

def residuals(params,x,y,sigma):
    a = params['a'].value
    b = params['b'].value
    
    model = a * x + b
    return (y-model)/sigma

# Function for calculating the dust-to-mass ratio a redshift, flux density
# and stellar mass
# ****************************************************************************

def calc_ratio(z,s850,Mstar):
    dislum = cosmo.luminosity_distance(z).value

# Calculate new dust temperature

    tp1 = Tdust**(4.0+beta)
    tp2 = Tcmb**(4.0+beta)
    tp3 = ( (1+z)**(4.0+beta) - 1.0)
    tp4 = tp1 + tp2 * tp3
    Tdust_corr = tp4**(1.0/(4.0+beta))

# Now calculate the fractional change in the flux

    Tcmb_new = Tcmb * (1.0+z)
    freq_em = freref * (1.0+z)
    tp1 = (con_bb * freq_em) / Tdust_corr
    tp2 = (con_bb * freq_em) / Tcmb_new

    tp3 = 1.0 / (np.exp(tp1) - 1.0)
    tp4 = 1.0 / (np.exp(tp2) - 1.0)

    corr_cmb = 1.0 - tp4/tp3
    corr_cmb = 1.0 / corr_cmb

# Calculate rest-frame 850-micron flux density

    fre_em = freref * (1.0 + z)
    tp1 = fre_em**(3.0+beta)
    tp2 = (con_bb * fre_em) / Tdust_corr
    tp1 = tp1 / (np.exp(tp2) - 1.0)
    
    tp3 = freref**(3.0+beta)
    tp4 = (con_bb * freref) / Tdust_corr
    tp3 = tp3 / (np.exp(tp4) - 1.0)
    
    s850_new = s850 * corr_cmb
    s850_rf = (tp3 / tp1) * s850_new

# calculate dust mass

    bb = freref**3.0 / (np.exp( (con_bb*freref)/Tdust_corr  ) - 1.0)

    tp1 = s850_rf * dislum * dislum
    tp2 = (1.0+z) * kappa * bb
    
    dust_mass_total = con1 * (tp1/tp2)
    dust_to_gas = dust_mass_total/Mstar

    return dust_to_gas

# Function for calculating the dust-to-mass ratio a redshift, flux density
# and stellar mass, using the T-z relationship from Liang et al. (2019)
# ****************************************************************************

def calc_ratio_temp_var(z,s850,Mstar):
    dislum = cosmo.luminosity_distance(z).value

# Calculate new dust temperature, using the relationship from Liang et al. 
# (2019), with the redshift 0 temperature changed to 22 K

    tp1 = -0.02 + 0.25 * np.log10(1.0+z)

# Change Liang formula to make everything consistent at z=0

    Tdust = 23.04 * np.power(10.0,tp1)
       
    tp1 = Tdust**(4.0+beta)
    tp2 = Tcmb**(4.0+beta)
    tp3 = ( (1+z)**(4.0+beta) - 1.0)
    tp4 = tp1 + tp2 * tp3
    Tdust_corr = tp4**(1.0/(4.0+beta))

# Now calculate the fractional change in the flux

    Tcmb_new = Tcmb * (1.0+z)
    freq_em = freref * (1.0+z)
    tp1 = (con_bb * freq_em) / Tdust_corr
    tp2 = (con_bb * freq_em) / Tcmb_new

    tp3 = 1.0 / (np.exp(tp1) - 1.0)
    tp4 = 1.0 / (np.exp(tp2) - 1.0)

    corr_cmb = 1.0 - tp4/tp3
    corr_cmb = 1.0 / corr_cmb

# Calculate rest-frame 850-micron flux density

    fre_em = freref * (1.0 + z)
    tp1 = fre_em**(3.0+beta)
    tp2 = (con_bb * fre_em) / Tdust_corr
    tp1 = tp1 / (np.exp(tp2) - 1.0)
    
    tp3 = freref**(3.0+beta)
    tp4 = (con_bb * freref) / Tdust_corr
    tp3 = tp3 / (np.exp(tp4) - 1.0)
    
    s850_new = s850 * corr_cmb
    s850_rf = (tp3 / tp1) * s850_new

# calculate dust mass

    bb = freref**3.0 / (np.exp( (con_bb*freref)/Tdust_corr  ) - 1.0)

    tp1 = s850_rf * dislum * dislum
    tp2 = (1.0+z) * kappa * bb
    
    dust_mass_total = con1 * (tp1/tp2)
    dust_to_gas = dust_mass_total/Mstar

    return dust_to_gas

# Function for calculating the dust-to-mass ratio a redshift, flux density
# and stellar mass, using the two-temperature model from Bakx et al. (2018)
# ****************************************************************************

def calc_ratio_bakx(z,s850,Mstar):
    dislum = cosmo.luminosity_distance(z).value

# Basic parameters of model from Bakx et al. (2018)

    Tcold = 21.29
    Thot = 45.80
    mass_ratio = 26.62

# correct dust temperatures for the effect of the CMB
       
    tp1 = Tcold**(4.0+beta)
    tp2 = Tcmb**(4.0+beta)
    tp3 = ( (1+z)**(4.0+beta) - 1.0)
    tp4 = tp1 + tp2 * tp3
    Tcold_corr = tp4**(1.0/(4.0+beta))
    
    tp1 = Thot**(4.0+beta)
    tp2 = Tcmb**(4.0+beta)
    tp3 = ( (1+z)**(4.0+beta) - 1.0)
    tp4 = tp1 + tp2 * tp3
    Thot_corr = tp4**(1.0/(4.0+beta))

# Now calculate the fractional change in the flux, separately for each 
# temperature

# The cold dust

    Tcmb_new = Tcmb * (1.0+z)
    freq_em = freref * (1.0+z)
    
    tp1 = (con_bb * freq_em) / Tcold_corr
    tp2 = (con_bb * freq_em) / Tcmb_new

    tp3 = 1.0 / (np.exp(tp1) - 1.0)
    tp4 = 1.0 / (np.exp(tp2) - 1.0)

    corr_cold_cmb = 1.0 - tp4/tp3
    corr_cold_cmb = 1.0 / corr_cold_cmb
    
# The hot dust
    
    tp1 = (con_bb * freq_em) / Thot_corr
    tp2 = (con_bb * freq_em) / Tcmb_new

    tp3 = 1.0 / (np.exp(tp1) - 1.0)
    tp4 = 1.0 / (np.exp(tp2) - 1.0)

    corr_hot_cmb = 1.0 - tp4/tp3
    corr_hot_cmb = 1.0 / corr_hot_cmb
    
# Calculate 850-micron flux density at emitted frequency

    fre_em = freref * (1.0 + z)

    tp1 = fre_em**(3.0+beta)
    
    tp2 = (con_bb * fre_em) / Tcold_corr
    tp3 = tp1 / (np.exp(tp2) - 1.0)
    
    tp4 = (con_bb * fre_em) / Thot_corr
    tp5 = tp1 / (np.exp(tp4) - 1.0)
    
# calculate the correction factor for the 850-micron flux
# density

# model without CMB correction

    model_flux_no_cmb = tp5 + mass_ratio * tp3
    
# model flux with CMB correction

    model_flux_cmb = tp5 * corr_hot_cmb + mass_ratio * tp3 * corr_cold_cmb
    
    cmb_corr = model_flux_cmb / model_flux_no_cmb
    
# calculate 850-micon emission at reference frequency
   
    tp1 = freref**(3.0+beta)
    
    tp2 = (con_bb * freref) / Tcold_corr
    tp3 = tp1 / (np.exp(tp2) - 1.0)
    
    tp4 = (con_bb * freref) / Thot_corr
    tp5 = tp1 / (np.exp(tp4) - 1.0) 
    
    model_rf_flux_cold = mass_ratio * tp3
    model_rf_flux_hot = tp5
    
    
    s850_new = s850 * cmb_corr
    
# calculate flux in the rest frame from the cold and hot components

    s850_rf_cold = (model_rf_flux_cold / model_flux_no_cmb) * s850_new
    s850_rf_hot = (model_rf_flux_hot / model_flux_no_cmb) * s850_new

# calculate dust mass of two components

    bb1 = freref**3.0 / (np.exp( (con_bb*freref)/Tcold_corr  ) - 1.0)
    bb2 = freref**3.0 / (np.exp( (con_bb*freref)/Thot_corr  ) - 1.0)

    tp1 = s850_rf_cold * dislum * dislum
    tp2 = (1.0+z) * kappa * bb1
    dust_mass_cold = con1 * (tp1/tp2)
    
    tp1 = s850_rf_hot * dislum * dislum
    tp2 = (1.0+z) * kappa * bb2
    dust_mass_hot = con1 * (tp1/tp2)
    
    dust_mass_total = dust_mass_cold + dust_mass_hot
    
    dust_to_gas = dust_mass_total/Mstar

    return dust_to_gas


# Read in data
# ****************************************************************************
# ****************************************************************************

# Read in results of stacking
# ****************************************************************************

file='simstack-regular-results.dat'

data = open(file,'r')

ic = 0

for line in data.readlines():  
    info = line.split()
    zlow = np.append(zlow,float(info[0]))
    zup = np.append(zup,float(info[1]))
    Mlow = np.append(Mlow,float(info[2]))
    Mup = np.append(Mup,float(info[3]))
    S850 = np.append(S850,float(info[4]))
    Error = np.append(Error,float(info[5]))
    
    ic = ic + 1
data.close()

nbin = len(Mlow)
print(nbin)

# Read in estimates of the noise from the simulation
# ****************************************************************************

file='results_from_Matts_simulation.dat'

data = open(file,'r')

ic = 0

for line in data.readlines():  
    info = line.split()
    art_noise = np.append(art_noise,float(info[5]))
    
# Now read in the original JWST catalogue
# ****************************************************************************

file = 'COSMOS-Webb-stacking-catalogue.dat'

galaxy_data = []

galaxy_data = np.loadtxt(file)

print(galaxy_data.shape)
ngal = galaxy_data.shape[0]
print(ngal)

# Change the units from Jy to mJy
# ****************************************************************************

S850 = S850 * 1000.0
Error = Error * 1000.0
art_noise = art_noise*1000.0

# Correct the flux errors by adding the two sets of noise in quadrature
# and print out the results for all the bins
# ****************************************************************************
# ****************************************************************************

for i in range(0,nbin):
    tp1 = Error[i]
    tp2 = art_noise[i]
    tp3 = tp1*tp1 + tp2*tp2
    Error[i] = np.sqrt(tp3)

for i in range(0,nbin):
    print(zlow[i],zup[i],Mlow[i],Mup[i],S850[i],Error[i])
    
# Now calculate the mean stellar mass and mean dust-to-stellar mass in
# each bin using my new averaging method (see notebook)
# ****************************************************************************
# ****************************************************************************

Mass_mean = np.empty(nbin)
z_mean = np.empty(nbin)
dust_to_stellar_mass = np.empty(nbin)
error_mass_ratio = np.empty(nbin)

dust_to_stellar_mass_Tvar = np.empty(nbin)
error_mass_ratio_Tvar = np.empty(nbin)

dust_to_stellar_mass_bakx_model = np.empty(nbin)
error_mass_ratio_bakx_model = np.empty(nbin)

for i in range(0,nbin):
    sel3 = np.where((galaxy_data[:,3] >= Mlow[i]) &\
                    (galaxy_data[:,3] < Mup[i]) &\
                    (galaxy_data[:,2] >= zlow[i]) &\
                        (galaxy_data[:,2] < zup[i]))
    
    masses = galaxy_data[:,3][sel3]
    redshifts = galaxy_data[:,2][sel3]
    
    masses = np.power(10.0,masses)
    
    ntarget = len(masses)
    
    fluxes = np.empty(ntarget)
    fluxes[:] = S850[i]
    
    Mass_mean[i] = np.mean(masses)
    Mass_mean[i] = np.log10(Mass_mean[i])
    z_mean[i] = np.mean(redshifts)
    
# dust-to-stars with only cold dust and no evolution

    dtog = calc_ratio(redshifts,fluxes,masses)
    
    dust_to_stellar_mass[i] = np.mean(dtog)
    
# dust-to-stars with temperature evolution model
    
    dtog = calc_ratio_temp_var(redshifts,fluxes,masses)
    dust_to_stellar_mass_Tvar[i] = np.mean(dtog)

# dust-to-stars with two-component model

    dtog = calc_ratio_bakx(redshifts,fluxes,masses)
    dust_to_stellar_mass_bakx_model[i] = np.mean(dtog)

# Calculate error on dust-to-stellar mass

    error_mass_ratio[i] = (Error[i] / S850[i]) * dust_to_stellar_mass[i]
    error_mass_ratio_Tvar[i] = (Error[i] / S850[i]) * dust_to_stellar_mass_Tvar[i]    
    error_mass_ratio_bakx_model[i] = (Error[i] / S850[i]) *\
        dust_to_stellar_mass_bakx_model[i]

# Now arrange the results for each bin for plotting
# ****************************************************************************
# ****************************************************************************

# First redshift bin
# ****************************************************************************

x1 = np.empty(5)
y1 = np.empty(5)
error1 = np.empty(5)
xline = np.empty(2)
yline = np.empty(2)

x21 = np.empty(5)
y21 = np.empty(5)
error21 = np.empty(5)

x41 = np.empty(5)
y41 = np.empty(5)
error41 = np.empty(5)

for i in range(0,5):
    x1[i] = Mass_mean[i]
    y1[i] = dust_to_stellar_mass[i]
    error1[i] = error_mass_ratio[i]

    x21[i] = Mass_mean[i]
    y21[i] = dust_to_stellar_mass_Tvar[i]
    error21[i] = error_mass_ratio_Tvar[i]
    
    x41[i] = Mass_mean[i]
    y41[i] = dust_to_stellar_mass_bakx_model[i]
    error41[i] = error_mass_ratio_bakx_model[i]
        
# second redshift bin
# ****************************************************************************

x2 = np.empty(5)
y2 = np.empty(5)
error2 = np.empty(5)

x22 = np.empty(5)
y22 = np.empty(5)
error22 = np.empty(5)

x42 = np.empty(5)
y42 = np.empty(5)
error42 = np.empty(5)

for i in range(0,5):
    x2[i] = Mass_mean[5+i]
    y2[i] = dust_to_stellar_mass[5+i]
    error2[i] = error_mass_ratio[5+i]
    
    x22[i] = Mass_mean[5+i]
    y22[i] = dust_to_stellar_mass_Tvar[5+i]
    error22[i] = error_mass_ratio_Tvar[5+i]

    x42[i] = Mass_mean[5+i]
    y42[i] = dust_to_stellar_mass_bakx_model[5+i]
    error42[i] = error_mass_ratio_bakx_model[5+i]
          
# third redshift bin
# ****************************************************************************

x3 = np.empty(5)
y3 = np.empty(5)
error3 = np.empty(5)

x23 = np.empty(5)
y23 = np.empty(5)
error23 = np.empty(5)

x43 = np.empty(5)
y43 = np.empty(5)
error43 = np.empty(5)

for i in range(0,5):
    x3[i] = Mass_mean[10+i]
    y3[i] = dust_to_stellar_mass[10+i]
    error3[i] = error_mass_ratio[10+i]

    x23[i] = Mass_mean[10+i]
    y23[i] = dust_to_stellar_mass_Tvar[10+i]
    error23[i] = error_mass_ratio_Tvar[10+i]

    x43[i] = Mass_mean[10+i]
    y43[i] = dust_to_stellar_mass_bakx_model[10+i]
    error43[i] = error_mass_ratio_bakx_model[10+i]
        
# fourth redshift bin
# ****************************************************************************

x4 = np.empty(5)
y4 = np.empty(5)
error4 = np.empty(5)

x24 = np.empty(5)
y24 = np.empty(5)
error24 = np.empty(5)

x44 = np.empty(5)
y44 = np.empty(5)
error44 = np.empty(5)

for i in range(0,5):
    x4[i] = Mass_mean[15+i]
    y4[i] = dust_to_stellar_mass[15+i]
    error4[i] = error_mass_ratio[15+i]
    
    x24[i] = Mass_mean[15+i]
    y24[i] = dust_to_stellar_mass_Tvar[15+i]
    error24[i] = error_mass_ratio_Tvar[15+i]

    x44[i] = Mass_mean[15+i]
    y44[i] = dust_to_stellar_mass_bakx_model[15+i]
    error44[i] = error_mass_ratio_bakx_model[15+i]
        
# fifth redshift bin
# *****************************************************************************

x5 = np.empty(5)
y5 = np.empty(5)
error5 = np.empty(5)

x25 = np.empty(5)
y25 = np.empty(5)
error25 = np.empty(5)

x45 = np.empty(5)
y45 = np.empty(5)
error45 = np.empty(5)

for i in range(0,5):
    x5[i] = Mass_mean[20+i]
    y5[i] = dust_to_stellar_mass[20+i]
    error5[i] = error_mass_ratio[20+i]
    
    x25[i] = Mass_mean[20+i]
    y25[i] = dust_to_stellar_mass_Tvar[20+i]
    error25[i] = error_mass_ratio_Tvar[20+i]

    x45[i] = Mass_mean[20+i]
    y45[i] = dust_to_stellar_mass_bakx_model[20+i]
    error45[i] = error_mass_ratio_bakx_model[20+i]
    
        
# sixth redshift bin
# ****************************************************************************

x6 = np.empty(5)
y6 = np.empty(5)
error6 = np.empty(5)

x26 = np.empty(5)
y26 = np.empty(5)
error26 = np.empty(5)

x46 = np.empty(5)
y46 = np.empty(5)
error46 = np.empty(5)

for i in range(0,5):
    x6[i] = Mass_mean[25+i]
    y6[i] = dust_to_stellar_mass[25+i]
    error6[i] = error_mass_ratio[25+i]

    x26[i] = Mass_mean[25+i]
    y26[i] = dust_to_stellar_mass_Tvar[25+i]
    error26[i] = error_mass_ratio_Tvar[25+i]

    x46[i] = Mass_mean[25+i]
    y46[i] = dust_to_stellar_mass_bakx_model[25+i]
    error46[i] = error_mass_ratio_bakx_model[25+i]
    
        
# seventh redshift bin
# ****************************************************************************

x7 = np.empty(5)
y7 = np.empty(5)
error7 = np.empty(5)

x27 = np.empty(5)
y27 = np.empty(5)
error27 = np.empty(5)

x47 = np.empty(5)
y47 = np.empty(5)
error47 = np.empty(5)

for i in range(0,5):
    x7[i] = Mass_mean[30+i]
    y7[i] = dust_to_stellar_mass[30+i]
    error7[i] = error_mass_ratio[30+i]

    x27[i] = Mass_mean[30+i]
    y27[i] = dust_to_stellar_mass_Tvar[30+i]
    error27[i] = error_mass_ratio_Tvar[30+i]

    x47[i] = Mass_mean[30+i]
    y47[i] = dust_to_stellar_mass_bakx_model[30+i]
    error47[i] = error_mass_ratio_bakx_model[30+i]
    
        
# eigth redshift bin
# *****************************************************************************

x8 = np.empty(5)
y8 = np.empty(5)
error8 = np.empty(5)

x28 = np.empty(5)
y28 = np.empty(5)
error28 = np.empty(5)

x48 = np.empty(5)
y48 = np.empty(5)
error48 = np.empty(5)

print("8th bin")

for i in range(0,5):
    x8[i] = Mass_mean[35+i]
    y8[i] = dust_to_stellar_mass[35+i]
    error8[i] = error_mass_ratio[35+i]

    x28[i] = Mass_mean[35+i]
    y28[i] = dust_to_stellar_mass_Tvar[35+i]
    error28[i] = error_mass_ratio_Tvar[35+i]

    x48[i] = Mass_mean[35+i]
    y48[i] = dust_to_stellar_mass_bakx_model[35+i]
    error48[i] = error_mass_ratio_bakx_model[35+i]
    
# ninth redshift bin
# ****************************************************************************

x9 = np.empty(5)
y9 = np.empty(5)
error9 = np.empty(5)

x29 = np.empty(5)
y29 = np.empty(5)
error29 = np.empty(5)

x49 = np.empty(5)
y49 = np.empty(5)
error49 = np.empty(5)

for i in range(0,5):
    x9[i] = Mass_mean[40+i]
    y9[i] = dust_to_stellar_mass[40+i]
    error9[i] = error_mass_ratio[40+i]
    
    x29[i] = Mass_mean[40+i]
    y29[i] = dust_to_stellar_mass_Tvar[40+i]
    error29[i] = error_mass_ratio_Tvar[40+i]

    x49[i] = Mass_mean[40+i]
    y49[i] = dust_to_stellar_mass_bakx_model[40+i]
    error49[i] = error_mass_ratio_bakx_model[40+i]
    
        
# tenth redshift bin
# ****************************************************************************

x10 = np.empty(5)
y10 = np.empty(5)
error10 = np.empty(5)

x30 = np.empty(5)
y30 = np.empty(5)
error30 = np.empty(5)

x50 = np.empty(5)
y50 = np.empty(5)
error50 = np.empty(5)


for i in range(0,5):
    x10[i] = Mass_mean[45+i]
    y10[i] = dust_to_stellar_mass[45+i]   
    error10[i] = error_mass_ratio[45+i]

    x30[i] = Mass_mean[45+i]
    y30[i] = dust_to_stellar_mass_Tvar[45+i]
    error30[i] = error_mass_ratio_Tvar[45+i]

    x50[i] = Mass_mean[45+i]
    y50[i] = dust_to_stellar_mass_bakx_model[45+i]
    error50[i] = error_mass_ratio_bakx_model[45+i]
    
         
  # 11th redshift bin
# ****************************************************************************

x11 = np.empty(5)
y11 = np.empty(5)
error11 = np.empty(5)

x31 = np.empty(5)
y31 = np.empty(5)
error31 = np.empty(5)

x51 = np.empty(5)
y51 = np.empty(5)
error51 = np.empty(5)

for i in range(0,5):
    x11[i] = Mass_mean[50+i]
    y11[i] = dust_to_stellar_mass[50+i]   
    error11[i] = error_mass_ratio[50+i]

    x31[i] = Mass_mean[50+i]
    y31[i] = dust_to_stellar_mass_Tvar[50+i]
    error31[i] = error_mass_ratio_Tvar[50+i]

    x51[i] = Mass_mean[50+i]
    y51[i] = dust_to_stellar_mass_bakx_model[50+i]
    error51[i] = error_mass_ratio_bakx_model[50+i]
    
# 12th redshift bin
# ****************************************************************************

x12 = np.empty(5)
y12 = np.empty(5)
error12 = np.empty(5)

x32 = np.empty(5)
y32 = np.empty(5)
error32 = np.empty(5)

x52 = np.empty(5)
y52 = np.empty(5)
error52 = np.empty(5)


for i in range(0,5):
    x12[i] = Mass_mean[55+i]
    y12[i] = dust_to_stellar_mass[55+i]   
    error12[i] = error_mass_ratio[55+i]

    x32[i] = Mass_mean[55+i]
    y32[i] = dust_to_stellar_mass_Tvar[55+i]
    error32[i] = error_mass_ratio_Tvar[55+i]

    x52[i] = Mass_mean[55+i]
    y52[i] = dust_to_stellar_mass_bakx_model[55+i]
    error52[i] = error_mass_ratio_bakx_model[55+i]
    
        
# 13th redshift bin
# ****************************************************************************

x13 = np.empty(5)
y13 = np.empty(5)
error13 = np.empty(5)

x33 = np.empty(5)
y33 = np.empty(5)
error33 = np.empty(5)

x53 = np.empty(5)
y53 = np.empty(5)
error53 = np.empty(5)


for i in range(0,5):
    x13[i] = Mass_mean[60+i]
    y13[i] = dust_to_stellar_mass[60+i]   
    error13[i] = error_mass_ratio[60+i]

    x33[i] = Mass_mean[60+i]
    y33[i] = dust_to_stellar_mass_Tvar[60+i]
    error33[i] = error_mass_ratio_Tvar[60+i]
    
    x53[i] = Mass_mean[60+i]
    y53[i] = dust_to_stellar_mass_bakx_model[60+i]
    error53[i] = error_mass_ratio_bakx_model[60+i]
    
# 14th redshift bin
# ****************************************************************************

x14 = np.empty(5)
y14 = np.empty(5)
error14 = np.empty(5)

x34 = np.empty(5)
y34 = np.empty(5)
error34 = np.empty(5)

x54 = np.empty(5)
y54 = np.empty(5)
error54 = np.empty(5)

for i in range(0,5):
    x14[i] = Mass_mean[65+i]
    y14[i] = dust_to_stellar_mass[65+i]   
    error14[i] = error_mass_ratio[65+i]

    x34[i] = Mass_mean[65+i]
    y34[i] = dust_to_stellar_mass_Tvar[65+i]
    error34[i] = error_mass_ratio_Tvar[65+i]
 
    x54[i] = Mass_mean[65+i]
    y54[i] = dust_to_stellar_mass_bakx_model[65+i]
    error54[i] = error_mass_ratio_bakx_model[65+i]
    
# Set 3-sigma upper limits for negative dust-to-star ratios for panels 4,8,14
# *****************************************************************************

y4_uplims = np.array([0,0,0,0,0],dtype=bool)
y8_uplims = np.array([0,0,0,0,0],dtype=bool)
y14_uplims = np.array([0,0,0,0,0],dtype=bool)

y24_uplims = np.array([0,0,0,0,0],dtype=bool)
y28_uplims = np.array([0,0,0,0,0],dtype=bool)
y34_uplims = np.array([0,0,0,0,0],dtype=bool)

y44_uplims = np.array([0,0,0,0,0],dtype=bool)
y48_uplims = np.array([0,0,0,0,0],dtype=bool)
y54_uplims = np.array([0,0,0,0,0],dtype=bool)

for i in range(0,5):
    
    if y4[i]/error4[i] < 0.0:
        y4[i] = 3.0*error4[i]
        y4_uplims[i] = 1

    if y8[i]/error8[i] < 0.0:
        y8[i] = 3.0*error8[i]
        y8_uplims[i] = 1
        
    if y14[i]/error14[i] < 0.0:
        y14[i] = 3.0*error14[i]
        y14_uplims[i] = 1
        
    if y24[i]/error24[i] < 0.0:
        y24[i] = 3.0*error24[i]
        y24_uplims[i] = 1

    if y28[i]/error28[i] < 0.0:
        y28[i] = 3.0*error28[i]
        y28_uplims[i] = 1
        
    if y34[i]/error34[i] < 0.0:
        y34[i] = 3.0*error34[i]
        y34_uplims[i] = 1
        
    if y44[i]/error44[i] < 0.0:
        y44[i] = 3.0*error44[i]
        y44_uplims[i] = 1

    if y48[i]/error48[i] < 0.0:
        y48[i] = 3.0*error48[i]
        y48_uplims[i] = 1
        
    if y54[i]/error54[i] < 0.0:
        y54[i] = 3.0*error54[i]
        y54_uplims[i] = 1
        

for i in range(0,5):
    print(y14[i],y24[i],y54[i])
        


# Plot out results and do fits to relationships
# ****************************************************************************
# *****************************************************************************


fit_results = np.empty((14,4))
                 
fig = plt.figure(figsize=(10.0,10.0))

pos_ann = ([9.5,0.03])

# First panel
# *****************************************************************************

# fit straight line 

result = minimize(residuals,params,method='leastsq',args=(x41,y41,error41))
b_result = result.params['b'].value
a_result = result.params['a'].value
print("Panel 1: a and b: ",a_result,b_result)

# Calculate weighted mean and error

ndata = len(y1)
weights = np.zeros(ndata)
for i in range(0,ndata):
    weights[i] = 1.0/error1[i]**2.0

weights_norm= weights/np.sum(weights)

tp1 = weights_norm*y1
rmean = np.sum(tp1)
sd = np.sqrt(1.0/np.sum(weights))
print("mean and error: ",rmean,sd)

fit_results[0,0] = a_result
fit_results[0,1] = b_result
fit_results[0,2] = rmean
fit_results[0,3] = sd

xline = np.arange(9.0,11.5,0.01)
yline = a_result*xline + b_result
    
# Now do plot

f1 = plt.axes([0.15,0.1,0.2,0.2])
f1.set_xlim(9.0,11.5)
f1.set_ylim(0.0001,0.099)
f1.set_yscale('log')
f1.tick_params(axis='both',which='both',labelsize=18)

f1.plot(x1,y1,'bs')
f1.errorbar(x1,y1,yerr= error1,fmt='b',ls='none')
f1.plot(xline,yline,'k--')

x21 = x21+0.1
f1.plot(x21,y21,'co')
f1.errorbar(x21,y21,yerr= error21,fmt='c',ls='none')

x41 = x41-0.1
f1.plot(x41,y41,'ro')
f1.errorbar(x41,y41,yerr= error41,fmt='r',ls='none')


f1.annotate('0.2<z<0.5',pos_ann,fontweight='bold',size=15)

# Second panel
# ****************************************************************************

# fit straight line 

result = minimize(residuals,params,method='leastsq',args=(x42,y42,error42))
b_result = result.params['b'].value
a_result = result.params['a'].value
print("Panel a and b: ",a_result,b_result)

# Calculate weighted mean and error

ndata = len(y2)
weights = np.zeros(ndata)
for i in range(0,ndata):
    weights[i] = 1.0/error2[i]**2.0

weights_norm= weights/np.sum(weights)

tp1 = weights_norm*y2
rmean = np.sum(tp1)
sd = np.sqrt(1.0/np.sum(weights))
print("mean and error: ",rmean,sd)

fit_results[1,0] = a_result
fit_results[1,1] = b_result
fit_results[1,2] = rmean
fit_results[1,3] = sd

xline = np.arange(9.0,11.5,0.0099)
yline = a_result*xline + b_result

f1 = plt.axes([0.35,0.1,0.2,0.2])
f1.set_xlim(9.0,11.5)
f1.set_ylim(0.0001,0.1)

f1.tick_params(axis='x',labelsize=18)
f1.set_yscale('log')
f1.set_yticks([])
f1.plot(x2,y2,'bs')
f1.errorbar(x2,y2,yerr= error2,fmt='b',ls='none')

x22 = x22+0.1
f1.plot(x22,y22,'co')
f1.errorbar(x22,y22,yerr= error22,fmt='c',ls='none')

x42 = x42-0.1
f1.plot(x42,y42,'ro')
f1.errorbar(x42,y42,yerr= error42,fmt='r',ls='none')


f1.annotate('0.5<z<0.8',pos_ann,fontweight='bold',size=15)

f1.plot(xline,yline,'k--')

# Third panel
# ****************************************************************************

# fit straight line 

result = minimize(residuals,params,method='leastsq',args=(x43,y43,error43))
b_result = result.params['b'].value
a_result = result.params['a'].value
print("Panel 3: a and b: ",a_result,b_result)

# Calculate weighted mean and error

ndata = len(y3)
weights = np.zeros(ndata)
for i in range(0,ndata):
    weights[i] = 1.0/error3[i]**2.0

weights_norm= weights/np.sum(weights)

tp1 = weights_norm*y3
rmean = np.sum(tp1)
sd = np.sqrt(1.0/np.sum(weights))
print("mean and error: ",rmean,sd)

fit_results[2,0] = a_result
fit_results[2,1] = b_result
fit_results[2,2] = rmean
fit_results[2,3] = sd

xline = np.arange(9.0,11.5,0.01)
yline = a_result*xline + b_result

# do plot

f1 = plt.axes([0.55,0.1,0.2,0.2])
f1.set_xlim(9.0,11.5)
f1.set_ylim(0.0001,0.099)
f1.set_yscale('log')
f1.set_xlabel('$log(M_*)$',size=25)
f1.tick_params(axis='x',labelsize=18)
f1.set_yticks([])

f1.plot(x3,y3,'bs')
f1.errorbar(x3,y3,yerr= error3,fmt='b',ls='none')
f1.plot(xline,yline,'k--')

x23 = x23+0.1
f1.plot(x23,y23,'co')
f1.errorbar(x23,y23,yerr= error23,fmt='c',ls='none')

x43 = x43-0.1
f1.plot(x43,y43,'ro')
f1.errorbar(x43,y43,yerr= error43,fmt='r',ls='none')


f1.annotate('0.8<z<1.1',pos_ann,fontweight='bold',size=15)

# Fourth panel
# ****************************************************************************

# fit straight line 

result = minimize(residuals,params,method='leastsq',args=(x44,y44,error44))
b_result = result.params['b'].value
a_result = result.params['a'].value
print("Panel 4: a and b: ",a_result,b_result)

# Calculate weighted mean and error

ndata = len(y4)
weights = np.zeros(ndata)
for i in range(0,ndata):
    weights[i] = 1.0/error4[i]**2.0

weights_norm= weights/np.sum(weights)

tp1 = weights_norm*y4
rmean = np.sum(tp1)
sd = np.sqrt(1.0/np.sum(weights))
print("mean and error: ",rmean,sd)

fit_results[3,0] = a_result
fit_results[3,1] = b_result
fit_results[3,2] = rmean
fit_results[3,3] = sd

xline = np.arange(9.0,11.5,0.001)
yline = a_result*xline + b_result

# Now plot panel

f1 = plt.axes([0.75,0.1,0.2,0.2])
f1.set_xlim(9.0,11.5)
f1.set_ylim(0.0001,0.099)
f1.set_yscale('log')
f1.tick_params(axis='x',labelsize=18)
f1.set_yticks([])

f1.annotate('1.1<z<1.5',pos_ann,fontweight='bold',size=15)

f1.plot(x4,y4,'bs')
f1.errorbar(x4,y4,yerr= error4,fmt='b',ls='none',uplims=y4_uplims)
f1.plot(xline,yline,'k--')

x24 = x24+0.1
f1.plot(x24,y24,'co')
f1.errorbar(x24,y24,yerr= error24,fmt='c',ls='none',uplims=y24_uplims)

x44 = x44-0.1
f1.plot(x44,y44,'ro')
f1.errorbar(x44,y44,yerr= error44,fmt='r',ls='none',uplims=y44_uplims)


# Fifth panel
# ****************************************************************************

# fit straight line 

result = minimize(residuals,params,method='leastsq',args=(x45,y45,error45))
b_result = result.params['b'].value
a_result = result.params['a'].value
print("a and b: ",a_result,b_result)

# Calculate weighted mean and error

ndata = len(y5)
weights = np.zeros(ndata)
for i in range(0,ndata):
    weights[i] = 1.0/error5[i]**2.0

weights_norm= weights/np.sum(weights)

tp1 = weights_norm*y5
rmean = np.sum(tp1)
sd = np.sqrt(1.0/np.sum(weights))
print("mean and error: ",rmean,sd)

fit_results[4,0] = a_result
fit_results[4,1] = b_result
fit_results[4,2] = rmean
fit_results[4,3] = sd

xline = np.arange(9.0,11.5,0.01)
yline = a_result*xline + b_result

f1 = plt.axes([0.15,0.3,0.2,0.2])
f1.set_xlim(9.0,11.5)
f1.set_ylim(0.0001,0.099)

f1.set_ylabel('$M_d/M_*$',size=25)
f1.tick_params(axis='y',labelsize=18)
f1.set_xticks([])
f1.set_yscale('log')
f1.annotate('1.5<z<2.0',pos_ann,fontweight='bold',size=15)

f1.plot(x5,y5,'bs')
f1.errorbar(x5,y5,yerr= error5,fmt='b',ls='none')
f1.plot(xline,yline,'k--')

x45 = x45-0.1
f1.plot(x45,y45,'ro')
f1.errorbar(x45,y45,yerr= error45,fmt='r',ls='none')

x25 = x25+0.1
f1.plot(x25,y25,'co')
f1.errorbar(x25,y25,yerr= error25,fmt='c',ls='none')

# Sixth panel
# ****************************************************************************

# fit straight line 

result = minimize(residuals,params,method='leastsq',args=(x46,y46,error46))
b_result = result.params['b'].value
a_result = result.params['a'].value
print("Panel 6: a and b: ",a_result,b_result)

# Calculate weighted mean and error

ndata = len(y6)
weights = np.zeros(ndata)
for i in range(0,ndata):
    weights[i] = 1.0/error6[i]**2.0

weights_norm= weights/np.sum(weights)

tp1 = weights_norm*y6
rmean = np.sum(tp1)
sd = np.sqrt(1.0/np.sum(weights))
print("mean and error: ",rmean,sd)

fit_results[5,0] = a_result
fit_results[5,1] = b_result
fit_results[5,2] = rmean
fit_results[5,3] = sd

xline = np.arange(9.0,11.5,0.01)
yline = a_result*xline + b_result

# Plot panel

f1 = plt.axes([0.35,0.3,0.2,0.2])
f1.set_xlim(9.0,11.5)
f1.set_ylim(0.0001,0.099)
f1.set_yscale('log')
f1.set_xticks([])
f1.set_yticks([])

f1.annotate('2.0<z<2.5',pos_ann,fontweight='bold',size=15)

f1.plot(x6,y6,'bs')
f1.errorbar(x6,y6,yerr= error6,fmt='b',ls='none')
f1.plot(xline,yline,'k--')

x26 = x26+0.1
f1.plot(x26,y26,'co')
f1.errorbar(x26,y26,yerr= error26,fmt='c',ls='none')

x46 = x46-0.1
f1.plot(x46,y46,'ro')
f1.errorbar(x46,y46,yerr= error46,fmt='r',ls='none')

# Seventh panel
# ****************************************************************************

# fit straight line 

result = minimize(residuals,params,method='leastsq',args=(x47,y47,error47))
b_result = result.params['b'].value
a_result = result.params['a'].value
print("Panel 7: a and b: ",a_result,b_result)

# Calculate weighted mean and error

ndata = len(y7)
weights = np.zeros(ndata)
for i in range(0,ndata):
    weights[i] = 1.0/error7[i]**2.0

weights_norm= weights/np.sum(weights)

tp1 = weights_norm*y7
rmean = np.sum(tp1)
sd = np.sqrt(1.0/np.sum(weights))
print("mean and error: ",rmean,sd)

fit_results[6,0] = a_result
fit_results[6,1] = b_result
fit_results[6,2] = rmean
fit_results[6,3] = sd

xline = np.arange(9.0,11.5,0.01)
yline = a_result*xline + b_result

# Plot panel

f1 = plt.axes([0.55,0.3,0.2,0.2])
f1.set_xlim(9.0,11.5)
f1.set_ylim(0.0001,0.099)
f1.set_yscale('log')
f1.set_xticks([])
f1.set_yticks([])

f1.annotate('2.5<z<3.0',pos_ann,fontweight='bold',size=15)

f1.plot(x7,y7,'bs')
f1.errorbar(x7,y7,yerr= error7,fmt='b',ls='none')
f1.plot(xline,yline,'k--')

x27 = x27+0.1
f1.plot(x27,y27,'co')
f1.errorbar(x27,y27,yerr= error27,fmt='c',ls='none')

x47 = x47-0.1
f1.plot(x47,y47,'ro')
f1.errorbar(x47,y47,yerr= error47,fmt='r',ls='none')

# Eighth panel
# ****************************************************************************

# fit straight line 

result = minimize(residuals,params,method='leastsq',args=(x48,y48,error48))
b_result = result.params['b'].value
a_result = result.params['a'].value
print("Panel 8: a and b: ",a_result,b_result)

# Calculate weighted mean and error

ndata = len(y8)
weights = np.zeros(ndata)
for i in range(0,ndata):
    weights[i] = 1.0/error8[i]**2.0

weights_norm= weights/np.sum(weights)

tp1 = weights_norm*y8
rmean = np.sum(tp1)
sd = np.sqrt(1.0/np.sum(weights))
print("mean and error: ",rmean,sd)

fit_results[7,0] = a_result
fit_results[7,1] = b_result
fit_results[7,2] = rmean
fit_results[7,3] = sd

xline = np.arange(9.0,11.5,0.01)
yline = a_result*xline + b_result

# Plot panel

f1 = plt.axes([0.75,0.3,0.2,0.2])
f1.set_xlim(9.0,11.5)
f1.set_ylim(0.0001,0.099)
f1.set_yscale('log')
f1.set_xticks([])
f1.set_yticks([])

f1.annotate('3.0<z<3.5',pos_ann,fontweight='bold',size=15)

f1.plot(x8,y8,'bs')
f1.errorbar(x8,y8,yerr= error8,fmt='b',ls='none',uplims=y8_uplims)
f1.plot(xline,yline,'k--')

x28 = x28+0.1
f1.plot(x28,y28,'co')
f1.errorbar(x28,y28,yerr= error28,fmt='c',ls='none',uplims=y28_uplims)

x48 = x48-0.1
f1.plot(x48,y48,'ro')
f1.errorbar(x48,y48,yerr= error48,fmt='r',ls='none',uplims=y48_uplims)


# Ninth panel
# ****************************************************************************

# fit straight line 

result = minimize(residuals,params,method='leastsq',args=(x49,y49,error49))
b_result = result.params['b'].value
a_result = result.params['a'].value
print("Panel 9: a and b: ",a_result,b_result)

# Calculate weighted mean and error

ndata = len(y9)
weights = np.zeros(ndata)
for i in range(0,ndata):
    weights[i] = 1.0/error9[i]**2.0

weights_norm= weights/np.sum(weights)

tp1 = weights_norm*y9
rmean = np.sum(tp1)
sd = np.sqrt(1.0/np.sum(weights))
print("mean and error: ",rmean,sd)

fit_results[8,0] = a_result
fit_results[8,1] = b_result
fit_results[8,2] = rmean
fit_results[8,3] = sd

xline = np.arange(9.0,11.5,0.01)
yline = a_result*xline + b_result

# plot panel

pos_ann = ([9.5,0.0002])


f1 = plt.axes([0.15,0.5,0.2,0.2])
f1.set_xlim(9.0,11.5)
f1.set_ylim(0.0001,0.099)
f1.set_yscale('log')
f1.set_xticks([])
f1.tick_params(axis='y',labelsize=18)

f1.annotate('3.5<z<4.5',pos_ann,fontweight='bold',size=15)

f1.plot(x9,y9,'bs')
f1.errorbar(x9,y9,yerr= error9,fmt='b',ls='none')
f1.plot(xline,yline,'k--')

x29 = x29+0.1
f1.plot(x29,y29,'co')
f1.errorbar(x29,y29,yerr= error29,fmt='c',ls='none')

x49 = x49-0.1
f1.plot(x49,y49,'ro')
f1.errorbar(x49,y49,yerr= error49,fmt='r',ls='none')

# Tenth panel
# ****************************************************************************

# fit straight line 

result = minimize(residuals,params,method='leastsq',args=(x50,y50,error50))
b_result = result.params['b'].value
a_result = result.params['a'].value
print("Panel 10: a and b: ",a_result,b_result)

# Calculate weighted mean and error

ndata = len(y10)
weights = np.zeros(ndata)
for i in range(0,ndata):
    weights[i] = 1.0/error10[i]**2.0

weights_norm= weights/np.sum(weights)

tp1 = weights_norm*y10
rmean = np.sum(tp1)
sd = np.sqrt(1.0/np.sum(weights))
print("mean and error: ",rmean,sd)

fit_results[9,0] = a_result
fit_results[9,1] = b_result
fit_results[9,2] = rmean
fit_results[9,3] = sd

xline = np.arange(9.0,11.5,0.01)
yline = a_result*xline + b_result

# Plot panel

f1 = plt.axes([0.35,0.5,0.2,0.2])
f1.set_xlim(9.0,11.5)
f1.set_ylim(0.0001,0.099)
f1.set_yscale('log')

f1.set_xticks([])
f1.set_yticks([])
f1.tick_params(axis='y',labelsize=18)
f1.annotate('4.5<z<5.5',pos_ann,fontweight='bold',size=15)

f1.plot(x10,y10,'bs')
f1.errorbar(x10,y10,yerr= error10,fmt='b',ls='none')
f1.plot(xline,yline,'k--')

x30 = x30+0.1
f1.plot(x30,y30,'co')
f1.errorbar(x30,y30,yerr= error30,fmt='c',ls='none')

x50 = x50-0.1
f1.plot(x50,y50,'ro')
f1.errorbar(x50,y50,yerr= error50,fmt='r',ls='none')


# 11th panel
# ****************************************************************************

# fit straight line 

result = minimize(residuals,params,method='leastsq',args=(x51,y51,error51))
b_result = result.params['b'].value
a_result = result.params['a'].value
print("Panel 11: a and b: ",a_result,b_result)

# Calculate weighted mean and error

ndata = len(y11)
weights = np.zeros(ndata)
for i in range(0,ndata):
    weights[i] = 1.0/error11[i]**2.0

weights_norm= weights/np.sum(weights)

tp1 = weights_norm*y11
rmean = np.sum(tp1)
sd = np.sqrt(1.0/np.sum(weights))
print("mean and error: ",rmean,sd)

fit_results[10,0] = a_result
fit_results[10,1] = b_result
fit_results[10,2] = rmean
fit_results[10,3] = sd

xline = np.arange(9.0,11.5,0.01)
yline = a_result*xline + b_result

# Plot panel

pos_ann = ([9.5,0.0002])


f1 = plt.axes([0.55,0.5,0.2,0.2])
f1.set_xlim(9.0,11.5)
f1.set_ylim(0.0001,0.099)
f1.set_yscale('log')
f1.set_xticks([])
f1.set_yticks([])
f1.tick_params(axis='y',labelsize=18)

f1.annotate('5.5<z<6.5',pos_ann,fontweight='bold',size=15)

f1.plot(x11,y11,'bs')
f1.errorbar(x11,y11,yerr= error11,fmt='b',ls='none')
f1.plot(xline,yline,'k--')

x31 = x31+0.1
f1.plot(x31,y31,'co')
f1.errorbar(x31,y31,yerr= error31,fmt='c',ls='none')

x51 = x51-0.1
f1.plot(x51,y51,'ro')
f1.errorbar(x51,y51,yerr= error51,fmt='r',ls='none')

# 12th panel
# ****************************************************************************

# fit straight line 

result = minimize(residuals,params,method='leastsq',args=(x52,y52,error52))
b_result = result.params['b'].value
a_result = result.params['a'].value
print("Panel 12: a and b: ",a_result,b_result)

# Calculate weighted mean and error

ndata = len(y12)
weights = np.zeros(ndata)
for i in range(0,ndata):
    weights[i] = 1.0/error12[i]**2.0

weights_norm= weights/np.sum(weights)

tp1 = weights_norm*y12
rmean = np.sum(tp1)
sd = np.sqrt(1.0/np.sum(weights))
print("mean and error: ",rmean,sd)

fit_results[11,0] = a_result
fit_results[11,1] = b_result
fit_results[11,2] = rmean
fit_results[11,3] = sd

xline = np.arange(9.0,11.5,0.01)
yline = a_result*xline + b_result

# plot panel

pos_ann = ([9.8,0.05])

f1 = plt.axes([0.75,0.5,0.2,0.2])
f1.set_xlim(9.0,11.5)
f1.set_ylim(0.0001,0.099)
f1.set_yscale('log')
f1.set_xticks([])
f1.set_yticks([])
f1.tick_params(axis='y',labelsize=18)

f1.annotate('6.5<z<7.5',pos_ann,fontweight='bold',size=15)

f1.plot(x12,y12,'bs')
f1.errorbar(x12,y12,yerr= error12,fmt='b',ls='none')
f1.plot(xline,yline,'k--')

x32 = x32+0.1
f1.plot(x32,y32,'co')
f1.errorbar(x32,y32,yerr= error32,fmt='c',ls='none')

x52 = x52-0.1
f1.plot(x52,y52,'ro')
f1.errorbar(x52,y52,yerr= error52,fmt='r',ls='none')


# 13th panel
# ****************************************************************************

# fit straight line 

result = minimize(residuals,params,method='leastsq',args=(x53,y53,error53))
b_result = result.params['b'].value
a_result = result.params['a'].value
print("Panel 13: a and b: ",a_result,b_result)

# Calculate weighted mean and error

ndata = len(y13)
weights = np.zeros(ndata)
for i in range(0,ndata):
    weights[i] = 1.0/error12[i]**2.0

weights_norm= weights/np.sum(weights)

tp1 = weights_norm*y13
rmean = np.sum(tp1)
sd = np.sqrt(1.0/np.sum(weights))
print("mean and error: ",rmean,sd)

fit_results[12,0] = a_result
fit_results[12,1] = b_result
fit_results[12,2] = rmean
fit_results[12,3] = sd

xline = np.arange(9.0,11.5,0.01)
yline = a_result*xline + b_result

# Plot panel

pos_ann = ([9.9,0.0002])

f1 = plt.axes([0.15,0.7,0.2,0.2])
f1.set_xlim(9.0,11.5)
f1.set_ylim(0.0001,0.099)
f1.set_yscale('log')
f1.set_xticks([])
f1.tick_params(axis='y',labelsize=18)

f1.annotate('7.5<z<8.5',pos_ann,fontweight='bold',size=15)

f1.plot(x13,y13,'bs')
f1.errorbar(x13,y13,yerr= error13,fmt='b',ls='none')
f1.plot(xline,yline,'k--')

x33 = x33+0.1
f1.plot(x33,y33,'co')
f1.errorbar(x33,y33,yerr= error33,fmt='c',ls='none')

x53 = x53-0.1
f1.plot(x53,y53,'ro')
f1.errorbar(x53,y53,yerr= error53,fmt='r',ls='none')


# 14th panel
# *****************************************************************************

# fit straight line 

result = minimize(residuals,params,method='leastsq',args=(x54,y54,error54))
b_result = result.params['b'].value
a_result = result.params['a'].value
print("Panel 14: a and b: ",a_result,b_result)

# Calculate weighted mean and error

ndata = len(y14)
weights = np.zeros(ndata)
for i in range(0,ndata):
    weights[i] = 1.0/error12[i]**2.0

weights_norm= weights/np.sum(weights)

tp1 = weights_norm*y14
rmean = np.sum(tp1)
sd = np.sqrt(1.0/np.sum(weights))
print("mean and error: ",rmean,sd)

fit_results[13,0] = a_result
fit_results[13,1] = b_result
fit_results[13,2] = rmean
fit_results[13,3] = sd

xline = np.arange(9.0,11.5,0.01)
yline = a_result*xline + b_result

# plot panel

pos_ann = ([9.5,0.05])

f1 = plt.axes([0.35,0.7,0.2,0.2])
f1.set_xlim(9.0,11.5)
f1.set_ylim(0.0001,0.099)
f1.set_yscale('log')
f1.set_xticks([])
f1.set_yticks([])
f1.tick_params(axis='y',labelsize=18)

f1.annotate('8.5<z<12.0',pos_ann,fontweight='bold',size=15)

f1.plot(x14,y14,'bs')
f1.errorbar(x14,y14,yerr= error14,fmt='b',ls='none',uplims=y14_uplims)
#f1.plot(xline,yline,'k--')

x34 = x34+0.1
f1.plot(x34,y34,'co')
f1.errorbar(x34,y34,yerr= error34,fmt='c',ls='none',uplims=y34_uplims)

x54 = x54-0.1
f1.plot(x54,y54,'ro')
f1.errorbar(x54,y54,yerr= error54,fmt='r',ls='none',uplims=y54_uplims)


fig.savefig('Figure5.pdf')
plt.show()

# Save results to file
# ****************************************************************************
# ****************************************************************************

with open('dust_to_stellar_mass_results_SIMSTACK_averaging_method.dat', 'a') as file:
    
    for i in range(0,nbin):
        file.write(f"{zlow[i]} ")
        file.write(f"{zup[i]} ")
        file.write(f"{Mlow[i]} ")
        file.write(f"{Mup[i]} ")
        file.write(f"{dust_to_stellar_mass[i]} ")
        file.write(f"{error_mass_ratio[i]} ")
        file.write(f"{dust_to_stellar_mass_Tvar[i]} ")
        file.write(f"{error_mass_ratio_Tvar[i]} ")
        file.write(f"{dust_to_stellar_mass_bakx_model[i]} ")
        file.write(f"{error_mass_ratio_bakx_model[i]}" + "\n")
        
        





            
            
    
    

    
        
    
    



    





