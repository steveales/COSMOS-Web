

# This program calculates the ratio of dust-to-stellar mass versus stellar
# mass for the stacking data from COSMOS Webb using a new method that 
# averages over each bin of redshift and stellar mass. This program
# uses the results from SIMSTACK, including the estimates
# of the effect of the redshift and mass errors.

# This version includes a correction for the effect of the CMB.
# This version also includes the T-z relation from Liang et al. (2019)
# and a two-temperature model from Bakx et al. (2018) 

# Last edited: 27nd October 2025

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
    Tdust = 22.0 * np.power(10.0,tp1)
       
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

# Read in dust model
# ****************************************************************************

file = 'results_dust_model.dat'

data = open(file,'r')

ic = 0

xmodel = np.array([])
ymodel1 = np.array([])
ymodel2 = np.array([])
ymodel3 = np.array([])

for line in data.readlines():
    info = line.split()
    xmodel = np.append(xmodel,float(info[0]))
    ymodel1 = np.append(ymodel1,float(info[4]))
    ymodel2 = np.append(ymodel2,float(info[5]))
    ymodel3 = np.append(ymodel3,float(info[6]))

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

# Now calculate mean dust-to-gas ratio in each redshift bin
# ****************************************************************************
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

ratio = np.zeros(14)
errors = np.zeros(14)
work = np.empty((5,4))

ratio_Tvar = np.zeros(14)
errors_Tvar = np.zeros(14)

ratio_Tmodel = np.zeros(14)
errors_Tmodel = np.zeros(14)

for i in range(0,14):

# First do this for constant temperature
# ****************************************************************************

    for j in range(0,5):
        work[j,0] = dust_to_stellar_mass[5*i+j]
        work[j,1] = error_mass_ratio[5*i+j]
        work[j,2] = 1.0/work[j,1]**2.0

    mean_dts = 0.0
    tot_weight = 0.0
    

    for j in range(0,5):
        mean_dts = mean_dts + work[j,2] * work[j,0]
        tot_weight = tot_weight + work[j,2]

    ratio[i] = mean_dts / tot_weight
    errors[i] = 1.0/np.sqrt(tot_weight)
    print(ratio[i],errors[i])

# Now do this for the T-z relationship
# *****************************************************************************

    for j in range(0,5):
        work[j,0] = dust_to_stellar_mass_Tvar[5*i+j]
        work[j,1] = error_mass_ratio_Tvar[5*i+j]
        work[j,2] = 1.0/work[j,1]**2.0

    mean_dts = 0.0
    tot_weight = 0.0

    for j in range(0,5):
        mean_dts = mean_dts + work[j,2] * work[j,0]
        tot_weight = tot_weight + work[j,2]

    ratio_Tvar[i] = mean_dts / tot_weight
    errors_Tvar[i] = 1.0/np.sqrt(tot_weight)

# Now do this for the two-temperature model
# ****************************************************************************

    for j in range(0,5):
        work[j,0] = dust_to_stellar_mass_bakx_model[5*i+j]
        work[j,1] = error_mass_ratio_bakx_model[5*i+j]
        work[j,2] = 1.0/work[j,1]**2.0

    mean_dts = 0.0
    tot_weight = 0.0

    for j in range(0,5):
        mean_dts = mean_dts + work[j,2] * work[j,0]
        tot_weight = tot_weight + work[j,2]

    ratio_Tmodel[i] = mean_dts / tot_weight
    errors_Tmodel[i] = 1.0/np.sqrt(tot_weight)

# Now plot results
# *****************************************************************************
# *****************************************************************************

fig = plt.figure(figsize=(10.0,10.0))

f1 = plt.axes([0.15,0.15,0.6,0.6])
f1.set_xlim(0.0,11.99)
f1.set_ylim(0.0001,0.1)
f1.set_yscale('log')
f1.tick_params(axis='both',which='both',labelsize=25)

f1.set_ylabel('$M_d/M_*$',size=25)
f1.set_xlabel('$Redshift$',size=25)

f1.plot(redshift,ratio,'bs')
f1.errorbar(redshift,ratio,yerr= errors,fmt='b',ls='none')

redshift = redshift + 0.1

f1.plot(redshift,ratio_Tvar,'co')
f1.errorbar(redshift,ratio_Tvar,yerr= errors_Tvar,fmt='c',ls='none')

redshift = redshift - 0.2

f1.plot(redshift,ratio_Tmodel,'ro')
f1.errorbar(redshift,ratio_Tmodel,yerr= errors_Tmodel,fmt='r',ls='none')

# Plot results of the chemical evolution model
# ****************************************************************************

f1.plot(xmodel,ymodel2,'g-')
f1.plot(xmodel,ymodel1,'m-')
f1.plot(xmodel,ymodel3,'y-')

pos_ann = ([2.0,5.0e-2])

f1.annotate('COLD',pos_ann,fontweight='bold',size=20,color='b')

pos_ann = ([2.0,3.0e-2])

f1.annotate('STAR-FORMER',pos_ann,fontweight='bold',size=20,color='r')

pos_ann = ([2.0,1.8e-2])

f1.annotate('EVOLVE',pos_ann,fontweight='bold',size=20,color='c')

xl = np.empty(2)
yl = np.empty(2)

xl[0] = 0.5
xl[1] = 1.5
yl[0] = 5.5e-2
yl[1] = 5.5e-2
f1.plot(xl,yl,'b-')

yl[0] = 3.3e-2
yl[1] = 3.3e-2
f1.plot(xl,yl,'r-')

yl[0] = 1.98e-2
yl[1] = 1.98e-2
f1.plot(xl,yl,'c-')


fig.savefig('Figure7.pdf')

plt.show()


