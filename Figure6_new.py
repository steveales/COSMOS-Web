# This program calculates the ratio of dust-to-stellar mass versus stellar
# mass for the stacking data from COSMOS Webb using a new method that 
# averages over each bin of redshift and stellar mass. This program
# uses the results of the basic stacking method, including the estimates
# of the effect of the redshift and mass errors.

# This version includes a correction for the effect of the CMB

# This version also plots the values with the temperature-redshift relationship 

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
    fluxes = S850[i]
    
    Mass_mean[i] = np.mean(masses)
    Mass_mean[i] = np.log10(Mass_mean[i])
    z_mean[i] = np.mean(redshifts)
    
    dtog = calc_ratio(redshifts,fluxes,masses)
    
    dust_to_stellar_mass[i] = np.mean(dtog)

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

print("bin1")

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

    print(x21[i],y21[i],error21[i])

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

print("Second redshift bin")

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

    print(x22[i],y22[i],error22[i])

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


print("3rd bin")

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

    print(x23[i],y23[i],error23[i])

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

print("4th bin")

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

    print(x24[i],y24[i],error24[i])

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

print("5th bin")

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


    print(x25[i],y25[i],error25[i])

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

print("6th bin")

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

    print(x26[i],y26[i],error26[i])


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


print("7th bin")

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

    print(x27[i],y27[i],error27[i])


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


    print(x28[i],y28[i],error28[i])

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

print("Ninth bin")

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

    print(x9[i],y9[i],error9[i])
    print(x29[i],y29[i],error29[i])

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

print("Tenth bin")

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

    print(y10[i],y10[i],error10[i])
    print(x30[i],y30[i],error30[i])

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

print("11th bin")

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

    print(x11[i],y11[i],error11[i])
    print(x31[i],y31[i],error31[i])


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

print("12th bin")

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


    print(x12[i],y12[i],error12[i])
    print(x32[i],y32[i],error32[i])

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


print("13th bin")

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

    print(x13[i],y13[i],error13[i])
    print(x33[i],y33[i],error33[i])

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

print("14th bin")

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

    print(x14[i],y14[i],error14[i])
    print(x34[i],y34[i],error34[i])


# Gather the data to produce plots of dust-to-stellar mass versus redshift
# ****************************************************************************

redshift = np.zeros(14)

redshift[0] = 0.35
redshift[1] = 0.65
redshift[2] = 0.95
redshift[3] = 1.3
redshift[4] = 1.75
redshift[5] = 2.25
redshift[6] = 2.75
redshift[7] = 3.25
redshift[8] = 4.0
redshift[9] = 5.0
redshift[10] = 6.0
redshift[11] = 7.0
redshift[12] = 8.0
redshift[13] = 10.25

redshift_shift = redshift + 0.2
redshift_shift2 = redshift - 0.2

ratio1 = np.zeros(14)
errors1 = np.zeros(14)
ratio2 = np.zeros(14)
errors2 = np.zeros(14)
ratio3 = np.zeros(14)
errors3 = np.zeros(14)
ratio4 = np.zeros(14)
errors4 = np.zeros(14)
ratio5 = np.zeros(14)
errors5 = np.zeros(14)

ratio1_temp_var = np.zeros(14)
errors1_temp_var = np.zeros(14)
ratio2_temp_var = np.zeros(14)
errors2_temp_var = np.zeros(14)
ratio3_temp_var = np.zeros(14)
errors3_temp_var = np.zeros(14)
ratio4_temp_var = np.zeros(14)
errors4_temp_var = np.zeros(14)
ratio5_temp_var = np.zeros(14)
errors5_temp_var = np.zeros(14)

ratio1_Tmodel = np.zeros(14)
errors1_Tmodel = np.zeros(14)
ratio2_Tmodel = np.zeros(14)
errors2_Tmodel = np.zeros(14)
ratio3_Tmodel = np.zeros(14)
errors3_Tmodel = np.zeros(14)
ratio4_Tmodel = np.zeros(14)
errors4_Tmodel = np.zeros(14)
ratio5_Tmodel = np.zeros(14)
errors5_Tmodel = np.zeros(14)

# First relationship

for i in range(0,14):
    ratio1[i] = dust_to_stellar_mass[5*i]
    errors1[i] = error_mass_ratio[5*i]

    ratio2[i] = dust_to_stellar_mass[1+5*i]
    errors2[i] = error_mass_ratio[1+5*i]

    ratio3[i] = dust_to_stellar_mass[2+5*i]
    errors3[i] = error_mass_ratio[2+5*i]

    ratio4[i] = dust_to_stellar_mass[3+5*i]
    errors4[i] = error_mass_ratio[3+5*i]

    ratio5[i] = dust_to_stellar_mass[4+5*i]
    errors5[i] = error_mass_ratio[4+5*i]
    
    ratio1_temp_var[i] = dust_to_stellar_mass_Tvar[5*i]
    errors1_temp_var[i] = error_mass_ratio_Tvar[5*i]

    ratio2_temp_var[i] = dust_to_stellar_mass_Tvar[1+5*i]
    errors2_temp_var[i] = error_mass_ratio_Tvar[1+5*i]

    ratio3_temp_var[i] = dust_to_stellar_mass_Tvar[2+5*i]
    errors3_temp_var[i] = error_mass_ratio_Tvar[2+5*i]

    ratio4_temp_var[i] = dust_to_stellar_mass_Tvar[3+5*i]
    errors4_temp_var[i] = error_mass_ratio_Tvar[3+5*i]

    ratio5_temp_var[i] = dust_to_stellar_mass_Tvar[4+5*i]
    errors5_temp_var[i] = error_mass_ratio_Tvar[4+5*i]
    
    ratio1_Tmodel[i] = dust_to_stellar_mass_bakx_model[5*i]
    errors1_Tmodel[i] = error_mass_ratio_bakx_model[5*i]

    ratio2_Tmodel[i] = dust_to_stellar_mass_bakx_model[1+5*i]
    errors2_Tmodel[i] = error_mass_ratio_bakx_model[1+5*i]

    ratio3_Tmodel[i] = dust_to_stellar_mass_bakx_model[2+5*i]
    errors3_Tmodel[i] = error_mass_ratio_bakx_model[2+5*i]

    ratio4_Tmodel[i] = dust_to_stellar_mass_bakx_model[3+5*i]
    errors4_Tmodel[i] = error_mass_ratio_bakx_model[3+5*i]

    ratio5_Tmodel[i] = dust_to_stellar_mass_bakx_model[4+5*i]
    errors5_Tmodel[i] = error_mass_ratio_bakx_model[4+5*i]
    
# Print out results

print("9 < M < 9.5")

for i in range(0,14):
    print(redshift[i],ratio1[i],errors1[i])
                   
print("9.5 < M < 10.0")

for i in range(0,14):
    print(redshift[i],ratio2[i],errors2[i])
            
print("10.0 < M < 10.5")

for i in range(0,14):
    print(redshift[i],ratio3[i],errors3[i])
                   
print("10.5 < M < 11.0")

for i in range(0,14):
    print(redshift[i],ratio4[i],errors4[i])
                   
print("11.0 < M < 11.5")

for i in range(0,14):
    print(redshift[i],ratio5[i],errors5[i])
    
# Setting negatives to 3-sigma upper limits
# ****************************************************************************

ratio1_uplims = np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0],dtype=bool)
ratio1_temp_var_uplims = np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0],dtype=bool)
ratio1_Tmodel_uplims = np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0],dtype=bool)

ratio2_uplims = np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0],dtype=bool)
ratio2_temp_var_uplims = np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0],dtype=bool)
ratio2_Tmodel_uplims = np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0],dtype=bool)
                   
for i in range(0,14):
    
    if ratio1[i] < 0.0:
        ratio1[i] = 3.0*errors1[i]
        ratio1_uplims[i] = 1
        
    if ratio1_temp_var[i] < 0.0:
        ratio1_temp_var[i] = 3.0*errors1_temp_var[i]
        ratio1_temp_var_uplims[i] = 1

    if ratio1_Tmodel[i] < 0.0:
        ratio1_Tmodel[i] = 3.0*errors1_Tmodel[i]
        ratio1_Tmodel_uplims[i] = 1

    if ratio2[i] < 0.0:
        ratio2[i] = 3.0*errors2[i]
        ratio2_uplims[i] = 1
        
    if ratio2_temp_var[i] < 0.0:
        ratio2_temp_var[i] = 3.0*errors2_temp_var[i]
        ratio2_temp_var_uplims[i] = 1

    if ratio2_Tmodel[i] < 0.0:
        ratio2_Tmodel[i] = 3.0*errors2_Tmodel[i]
        ratio2_Tmodel_uplims[i] = 1

# Plot out results
# ****************************************************************************

fig = plt.figure(figsize=(10.0,10.0))

pos_ann = ([0.1,0.06])

# Panel 1

f1 = plt.axes([0.15,0.15,0.25,0.25])
f1.set_xlim(0.0,11.99)
f1.set_ylim(0.0001,0.099)
f1.set_yscale('log')
f1.tick_params(axis='both',which='both',labelsize=18)

f1.set_ylabel('$M_d/M_*$',size=25)

f1.plot(redshift,ratio1,'bs')
f1.errorbar(redshift,ratio1,yerr= errors1,fmt='b',ls='none',uplims=ratio1_uplims)

f1.plot(redshift_shift,ratio1_temp_var,'co')
f1.errorbar(redshift_shift,ratio1_temp_var,yerr= errors1_temp_var,fmt='c',ls='none',\
            uplims=ratio1_temp_var_uplims)

f1.plot(redshift_shift2,ratio1_Tmodel,'ro')
f1.errorbar(redshift_shift2,ratio1_Tmodel,yerr= errors1_Tmodel,fmt='r',ls='none',\
            uplims=ratio1_Tmodel_uplims)

f1.annotate('$9.0 < log_{10}M_*) < 9.5$',pos_ann,fontweight = 'bold',size=15)

# Panel 2

f1 = plt.axes([0.4,0.15,0.25,0.25])
f1.set_xlim(0.0,11.99)
f1.set_ylim(0.0001,0.099)
f1.set_yscale('log')
f1.tick_params(axis='x',labelsize=18)
f1.set_yticks([])

f1.set_xlabel('redshift',size=25)

f1.plot(redshift,ratio2,'bs')
f1.errorbar(redshift,ratio2,yerr= errors2,fmt='b',ls='none',uplims=ratio2_uplims)

f1.plot(redshift_shift,ratio2_temp_var,'co')
f1.errorbar(redshift_shift,ratio2_temp_var,yerr= errors2_temp_var,fmt='c',ls='none',\
            uplims=ratio2_temp_var_uplims)

f1.plot(redshift_shift2,ratio2_Tmodel,'ro')
f1.errorbar(redshift_shift2,ratio2_Tmodel,yerr= errors2_Tmodel,fmt='r',ls='none',\
            uplims=ratio2_Tmodel_uplims)

f1.annotate('$9.5 < log_{10}M_*) < 10.0$',pos_ann,fontweight = 'bold',size=15)

# Panel 3

f1 = plt.axes([0.65,0.15,0.25,0.25])
f1.set_xlim(0.0,11.99)
f1.set_ylim(0.0001,0.099)
f1.set_yscale('log')
f1.tick_params(axis='x',labelsize=18)
f1.set_yticks([])

f1.plot(redshift,ratio3,'bs')
f1.errorbar(redshift,ratio3,yerr= errors3,fmt='b',ls='none')

f1.plot(redshift_shift,ratio3_temp_var,'co')
f1.errorbar(redshift_shift,ratio3_temp_var,yerr= errors3_temp_var,fmt='c',ls='none')

f1.plot(redshift_shift2,ratio3_Tmodel,'ro')
f1.errorbar(redshift_shift2,ratio3_Tmodel,yerr= errors3_Tmodel,fmt='r',ls='none')

f1.annotate('$10.0 < log_{10}M_*) < 10.5$',pos_ann,fontweight = 'bold',size=15)

# Panel 4

f1 = plt.axes([0.15,0.4,0.25,0.25])
f1.set_xlim(0.0,11.99)
f1.set_ylim(0.0001,0.099)
f1.set_yscale('log')
f1.tick_params(axis='y',labelsize=18)
f1.set_xticks([])

f1.set_ylabel('$M_d/M_*$',size=25)
f1.plot(redshift,ratio4,'bs')
f1.errorbar(redshift,ratio4,yerr= errors4,fmt='b',ls='none')

f1.plot(redshift_shift,ratio4_temp_var,'co')
f1.errorbar(redshift_shift,ratio4_temp_var,yerr= errors4_temp_var,fmt='c',ls='none')

f1.plot(redshift_shift2,ratio4_Tmodel,'ro')
f1.errorbar(redshift_shift2,ratio4_Tmodel,yerr= errors4_Tmodel,fmt='r',ls='none')

f1.annotate('$10.5 < log_{10}M_*) < 11.0$',pos_ann,fontweight = 'bold',size=15)

# Panel 5

f1 = plt.axes([0.4,0.4,0.25,0.25])
f1.set_xlim(0.0,11.99)
f1.set_ylim(0.0001,0.099)
f1.set_yscale('log')
f1.set_xticks([])
f1.set_yticks([])

f1.plot(redshift,ratio5,'bs')
f1.errorbar(redshift,ratio5,yerr= errors5,fmt='b',ls='none')

f1.plot(redshift_shift,ratio5_temp_var,'co')
f1.errorbar(redshift_shift,ratio5_temp_var,yerr= errors5_temp_var,fmt='c',ls='none')

f1.plot(redshift_shift2,ratio5_Tmodel,'ro')
f1.errorbar(redshift_shift2,ratio5_Tmodel,yerr= errors5_Tmodel,fmt='r',ls='none')

f1.annotate('$11.0 < log_{10}M_*) < 11.5$',pos_ann,fontweight = 'bold',size=15)

fig.savefig('Figure6.pdf')
plt.show()

plt.show()