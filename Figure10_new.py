# This program plots mean star-formation rate versus redshift for
# each morphological class and range of galaxy stellar mass. The program
# adds the estimate of the star-formation rate from the JWST catalogue
# to an estimate from the 850-micron flux, based on the assumption of the
# standard Herschel SED 

# Important note: I've only used the bolometric luminosity as the SFR measure

# Note: 3sigma upper limits are shown as triangles

# Last edited: 10th November 2025

import numpy as np
import matplotlib.pyplot as plt
from astropy.cosmology import Planck15 as cosmo

# Constants
# ****************************************************************************
# *****************************************************************************

# General constants
# ****************************************************************************

# Constant in modified blackbody

con_bb = 48.01

# Constant from Kennicutt and Evans (2012)

con_KE = 43.41

# constant to calculate bolometric luminosity in ergs/s

con1 = 9.522e22

# value of beta

beta = 2.0

# SED model from Pearson et al. (2013)
# ****************************************************************************

Tcold = 23.9
Thot = 46.9
mass_ratio = 30.1

# Function to calculate flux densities from model SED
# ****************************************************************************
# ****************************************************************************

def create_model(wavelength,Tcold,Thot,mass_ratio):

# calculate frequency in units of 10^12 Hz

    freq = 3.0e2/wavelength
    tp1 = np.exp((con_bb * freq) / Tcold)

    tp2 = 1.0 / (tp1 - 1.0)
    
    tp3 = np.exp((con_bb * freq) / Thot)

    tp4 = 1.0 / (tp3 - 1.0)
    
    tp5 = tp4 + mass_ratio * tp2

    flux = tp5 * freq**(beta + 3.0)
    return flux

# Read in results from stacking
# ****************************************************************************
# ****************************************************************************

Mlow = np.array([])
Mup = np.array([])
zlow = np.array([])
zup = np.array([])
zmed = np.array([])
S850 = np.array([])
Error = np.array([])
Ngal = np.array([])
art_noise = np.array([])
        
file='SIMSTACK_results_collated_plus_number.dat'

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
    Ngal = np.append(Ngal,float(info[6]))
    
    ic = ic + 1
data.close()

nbin = len(Mlow)
print(nbin)

# Read in artificial data
# ****************************************************************************

file='SIMSTACK_artdata_collated.dat'

data = open(file,'r')

ic = 0

for line in data.readlines():
    info = line.split()
    art_noise = np.append(art_noise,float(info[5]))

    ic = ic + 1
data.close()

nbin = len(Mlow)
print(nbin)


# Convert fluxes into Jy
# ****************************************************************************

for i in range(0,nbin):
    S850[i] = S850[i] * 1000.0
    Error[i] = Error[i] * 1000.0
    art_noise[i] = art_noise[i] * 1000.0
    
# Add two noise estimates in quadrature
# ****************************************************************************

for i in range(0,nbin):
    tp1 = art_noise[i] * art_noise[i] + Error[i] * Error[i]
    Error[i] = np.sqrt(tp1)

# Read in the results for the SMGs
# ****************************************************************************
# ****************************************************************************

file = 'results_SMG.dat'

SMG_data = []

SMG_data = np.loadtxt(file)

n_smg = SMG_data.shape[0]

print(n_smg)
    
# Find the mean redshifts and stellar masses in each bin.
# ****************************************************************************

Mass_mean = np.empty(nbin)
z_mean = np.empty(nbin)

for i in range(0,nbin):
    Mass_mean[i] = (Mlow[i] + Mup[i]) / 2.0
    z_mean[i] = (zlow[i] + zup[i]) / 2.0
    
# Calculate hidden star-formation rate for each bin
# ****************************************************************************
# ****************************************************************************

# set up arrays

sfr_hidden = np.empty(nbin)
sfr_hidden[:] = -99.0
sfr_hidden_error = np.empty(nbin)
sfr_hidden_error[:] = -99.0

sfr_mean = np.empty(nbin)
sfr_mean[:] = 0.0

# wavelengths and frequencies

wave = np.arange(3.0,1100.0,1.0)

freq = 3.0e14 / wave

nlen = len(freq)

# calculate del(frequency)

del_freq = np.empty(nlen)

for i in range(0,nlen):

    if i==0:
        del_freq[0] = freq[0] - freq[1]
    elif i==(nlen-1):
        del_freq[nlen-1] = freq[nlen-2] - freq[nlen-1]
    else:
        tp1 = freq[i-1]
        tp2 = freq[i+1]
        del_freq[i] = (tp1 - tp2)/2.0
        
# Calculating bolometric luminosity
# ****************************************************************************

for j in range(0,nbin):
    z = z_mean[j]
    dislum = cosmo.luminosity_distance(z).value
    

# check whether there is a measurement of 850-micron flux

    if S850[j] < -98.0:
        continue
    
# Check whether there are enough galaxies for an estimate

    if Ngal[j] < 10.0:
        continue
    
# Find index coresponding to 850 microns

    wave_find = 850.0/(1.0 + z)

    for i in range(0,nlen):

        if wave[i] > wave_find:
            ind850 = i
            break
    
# find scaling for SED

    tp1 = create_model(wave[ind850],Tcold,Thot,mass_ratio)
    rnorm = S850[j] / tp1

    work1 = create_model(wave,Tcold,Thot,mass_ratio)
    work2 = work1 * del_freq

    tot = rnorm * np.sum(work2)

    lbol = con1 * 4.0 * 3.14159 * dislum * dislum * tot

# calculate star-formation rate
    if lbol > 0.0:
        tp1 = np.log10(lbol) - con_KE

        sfr_hidden[j] = np.power(10.0,tp1)
        sfr_hidden_error[j] = (Error[j]/S850[j]) * sfr_hidden[j]
    else:
        lbol = (Error[j]/S850[j]) * lbol
        tp1 = np.log10(lbol) - con_KE

        sfr_hidden_error[j] = np.power(10.0,tp1)
        sfr_hidden[j] = -1.0
        

# Change bins with -99.0s so they won't be plotted
# ****************************************************************************

for i in range(0,nbin):
    
    if sfr_hidden[i] < -98.0:
        sfr_hidden[i] = 1.0e-9
        sfr_hidden_error[i] = 1.0e-10
        
# Now arrange the results for each bin for plotting
# ****************************************************************************

# set the JWST estimates to zero

sfr_mean[:] = 0.0

# First mass range: 9.0 - 9.5

x1_m1 = np.empty(14)
y1_m1 = np.empty(14)
error1_m1 = np.empty(14)

x1_m2 = np.empty(14)
y1_m2 = np.empty(14)
error1_m2 = np.empty(14)

x1_m3 = np.empty(14)
y1_m3 = np.empty(14)
error1_m3 = np.empty(14)

x1_m4 = np.empty(14)
y1_m4 = np.empty(14)
error1_m4 = np.empty(14)

print("bin1")

for i in range(0,14):
    x1_m1[i] = z_mean[5*i]
    y1_m1[i] = sfr_mean[5*i] + sfr_hidden[5*i]
    error1_m1[i] = sfr_hidden_error[5*i]
    
for i in range(0,14):
    x1_m2[i] = z_mean[70+ 5*i]
    y1_m2[i] = sfr_mean[70 + 5*i] + sfr_hidden[70 + 5*i]
    error1_m2[i] = sfr_hidden_error[70 + 5*i]

for i in range(0,14):
    x1_m3[i] = z_mean[140+ 5*i]
    y1_m3[i] = sfr_mean[140 + 5*i] + sfr_hidden[140 + 5*i]
    error1_m3[i] = sfr_hidden_error[140 + 5*i]

for i in range(0,14):
    x1_m4[i] = z_mean[210+ 5*i]
    y1_m4[i] = sfr_mean[210 + 5*i] + sfr_hidden[210 + 5*i]
    error1_m4[i] = sfr_hidden_error[210 + 5*i]

# Second mass range: 9.5-10.0

x2_m1 = np.empty(14)
y2_m1 = np.empty(14)
error2_m1 = np.empty(14)

x2_m2 = np.empty(14)
y2_m2 = np.empty(14)
error2_m2 = np.empty(14)

x2_m3 = np.empty(14)
y2_m3 = np.empty(14)
error2_m3 = np.empty(14)

x2_m4 = np.empty(14)
y2_m4 = np.empty(14)
error2_m4 = np.empty(14)

print("bin2")

for i in range(0,14):
    x2_m1[i] = z_mean[5*i + 1]
    y2_m1[i] = sfr_mean[5*i + 1] + sfr_hidden[5*i + 1]
    error2_m1[i] = sfr_hidden_error[5*i + 1]
    
for i in range(0,14):
    x2_m2[i] = z_mean[71+ 5*i]
    y2_m2[i] = sfr_mean[71 + 5*i] + sfr_hidden[71 + 5*i]
    error2_m2[i] = sfr_hidden_error[71 + 5*i]

for i in range(0,14):
    x2_m3[i] = z_mean[141+ 5*i]
    y2_m3[i] = sfr_mean[141 + 5*i] + sfr_hidden[141 + 5*i]
    error2_m3[i] = sfr_hidden_error[141 + 5*i]

for i in range(0,14):
    x2_m4[i] = z_mean[211+ 5*i]
    y2_m4[i] = sfr_mean[211 + 5*i] + sfr_hidden[211 + 5*i]
    error2_m4[i] = sfr_hidden_error[211 + 5*i]

# Third mass range: 10.0-10.5

x3_m1 = np.empty(14)
y3_m1 = np.empty(14)
error3_m1 = np.empty(14)

x3_m2 = np.empty(14)
y3_m2 = np.empty(14)
error3_m2 = np.empty(14)

x3_m3 = np.empty(14)
y3_m3 = np.empty(14)
error3_m3 = np.empty(14)

x3_m4 = np.empty(14)
y3_m4 = np.empty(14)
error3_m4 = np.empty(14)

print("bin3")

for i in range(0,14):
    x3_m1[i] = z_mean[5*i + 2]
    y3_m1[i] = sfr_mean[5*i + 2] + sfr_hidden[5*i + 2]
    error3_m1[i] = sfr_hidden_error[5*i + 2]
    
for i in range(0,14):
    x3_m2[i] = z_mean[72+ 5*i]
    y3_m2[i] = sfr_mean[72 + 5*i] + sfr_hidden[72 + 5*i]
    error3_m2[i] = sfr_hidden_error[72 + 5*i]

for i in range(0,14):
    x3_m3[i] = z_mean[142+ 5*i]
    y3_m3[i] = sfr_mean[142 + 5*i] + sfr_hidden[142 + 5*i]
    error3_m3[i] = sfr_hidden_error[142 + 5*i]

for i in range(0,14):
    x3_m4[i] = z_mean[212+ 5*i]
    y3_m4[i] = sfr_mean[212 + 5*i] + sfr_hidden[212 + 5*i]
    error3_m4[i] = sfr_hidden_error[212 + 5*i]

# Fourth mass range: 10.5-11.0

x4_m1 = np.empty(14)
y4_m1 = np.empty(14)
error4_m1 = np.empty(14)

x4_m2 = np.empty(14)
y4_m2 = np.empty(14)
error4_m2 = np.empty(14)

x4_m3 = np.empty(14)
y4_m3 = np.empty(14)
error4_m3 = np.empty(14)

x4_m4 = np.empty(14)
y4_m4 = np.empty(14)
error4_m4 = np.empty(14)

print("bin4")

for i in range(0,14):
    x4_m1[i] = z_mean[5*i + 3]
    y4_m1[i] = sfr_mean[5*i + 3] + sfr_hidden[5*i + 3]
    error4_m1[i] = sfr_hidden_error[5*i + 3]
    
for i in range(0,14):
    x4_m2[i] = z_mean[73+ 5*i]
    y4_m2[i] = sfr_mean[73 + 5*i] + sfr_hidden[73 + 5*i]
    error4_m2[i] = sfr_hidden_error[73 + 5*i]

for i in range(0,14):
    x4_m3[i] = z_mean[143+ 5*i]
    y4_m3[i] = sfr_mean[143 + 5*i] + sfr_hidden[143 + 5*i]
    error4_m3[i] = sfr_hidden_error[143 + 5*i]

for i in range(0,14):
    x4_m4[i] = z_mean[213+ 5*i]
    y4_m4[i] = sfr_mean[213 + 5*i] + sfr_hidden[213 + 5*i]
    error4_m4[i] = sfr_hidden_error[213 + 5*i]

# Fifth mass range: 11.0-11.5

x5_m1 = np.empty(14)
y5_m1 = np.empty(14)
error5_m1 = np.empty(14)

x5_m2 = np.empty(14)
y5_m2 = np.empty(14)
error5_m2 = np.empty(14)

x5_m3 = np.empty(14)
y5_m3 = np.empty(14)
error5_m3 = np.empty(14)

x5_m4 = np.empty(14)
y5_m4 = np.empty(14)
error5_m4 = np.empty(14)

print("bin 5")

for i in range(0,14):
    x5_m1[i] = z_mean[5*i + 4]
    y5_m1[i] = sfr_mean[5*i + 4] + sfr_hidden[5*i + 4]
    error5_m1[i] = sfr_hidden_error[5*i + 4]
    
for i in range(0,14):
    x5_m2[i] = z_mean[74+ 5*i]
    y5_m2[i] = sfr_mean[74 + 5*i] + sfr_hidden[74 + 5*i]
    error5_m2[i] = sfr_hidden_error[74 + 5*i]

for i in range(0,14):
    x5_m3[i] = z_mean[144+ 5*i]
    y5_m3[i] = sfr_mean[144 + 5*i] + sfr_hidden[144 + 5*i]
    error5_m3[i] = sfr_hidden_error[144 + 5*i]

for i in range(0,14):
    x5_m4[i] = z_mean[214+ 5*i]
    y5_m4[i] = sfr_mean[214 + 5*i] + sfr_hidden[214 + 5*i]
    error5_m4[i] = sfr_hidden_error[214 + 5*i]
    
# Remover points that are not 2sigma detections and substitute limits
# *****************************************************************************

# Bin 1

y1_m1_uplims = np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0],dtype=bool)
y1_m2_uplims = np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0],dtype=bool)
y1_m3_uplims = np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0],dtype=bool)
y1_m4_uplims = np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0],dtype=bool)

for i in range(0,14):

    if y1_m1[i]/error1_m1[i] < 0.0:
        y1_m1[i] = 3.0 * error1_m1[i]
        y1_m1_uplims[i] = 1

    if y1_m2[i]/error1_m2[i] < 0.0:
        y1_m2[i] = 3.0 * error1_m2[i]
        y1_m2_uplims[i] = 1
        
    if y1_m3[i]/error1_m3[i] < 0.0:
        y1_m3[i] = 3.0 * error1_m3[i]
        y1_m3_uplims[i] = 1

    if y1_m4[i]/error1_m4[i] < 0.0:
        y1_m4[i] = 3.0 * error1_m4[i]
        y1_m4_uplims[i] = 1

for i in range(0,14):
    print(x1_m1[i],y1_m1[i],error1_m1[i],y1_m1_uplims[i])
        
# Bin 2

y2_m1_uplims = np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0],dtype=bool)
y2_m2_uplims = np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0],dtype=bool)
y2_m3_uplims = np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0],dtype=bool)
y2_m4_uplims = np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0],dtype=bool)

for i in range(0,14):

    if y2_m1[i]/error2_m1[i] < 0.0:
        y2_m1[i] = 3.0 * error2_m1[i]
        y2_m1_uplims[i] = 1

    if y2_m2[i]/error2_m2[i] < 0.0:
        y2_m2[i] = 3.0 * error2_m2[i]
        y2_m2_uplims[i] = 1
        
    if y2_m3[i]/error2_m3[i] < 0.0:
        y2_m3[i] = 3.0 * error2_m3[i]
        y2_m3_uplims[i] = 1

    if y2_m4[i]/error2_m4[i] < 0.0:
        y2_m4[i] = 3.0 * error2_m4[i]
        y2_m4_uplims[i] = 1
        
# Bin 3

y3_m1_uplims = np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0],dtype=bool)
y3_m2_uplims = np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0],dtype=bool)
y3_m3_uplims = np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0],dtype=bool)
y3_m4_uplims = np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0],dtype=bool)

for i in range(0,14):

    if y3_m1[i]/error3_m1[i] < 0.0:
        y3_m1[i] = 3.0 * error3_m1[i]
        y3_m1_uplims[i] = 1

    if y3_m2[i]/error3_m2[i] < 0.0:
        y3_m2[i] = 3.0 * error3_m2[i]
        y3_m2_uplims[i] = 1

    if y3_m3[i]/error3_m3[i] < 0.0:
        y3_m3[i] = 3.0 * error3_m3[i]
        y3_m3_uplims[i] = 1

    if y3_m4[i]/error3_m4[i] < 0.0:
        y3_m4[i] = 3.0 * error3_m4[i]
        y3_m4_uplims[i] = 1
        
# Bin 4

y4_m1_uplims = np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0],dtype=bool)
y4_m2_uplims = np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0],dtype=bool)
y4_m3_uplims = np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0],dtype=bool)
y4_m4_uplims = np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0],dtype=bool)

for i in range(0,14):

    if y4_m1[i]/error4_m1[i] < 0.0:
        y4_m1[i] = 3.0 * error4_m1[i]
        y4_m1_uplims[i] = 1

    if y4_m2[i]/error4_m2[i] < 0.0:
        y4_m2[i] = 3.0 * error4_m2[i]
        y4_m2_uplims[i] = 1
        
    if y4_m3[i]/error4_m3[i] < 0.0:
        y4_m3[i] = 3.0 * error4_m3[i]
        y4_m3_uplims[i] = 1

    if y4_m4[i]/error4_m4[i] < 0.0:
        y4_m4[i] = 3.0 * error4_m4[i]
        y4_m4_uplims[i] = 1
        
# Bin 5

y5_m1_uplims = np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0],dtype=bool)
y5_m2_uplims = np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0],dtype=bool)
y5_m3_uplims = np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0],dtype=bool)
y5_m4_uplims = np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0],dtype=bool)


for i in range(0,14):

    if y5_m1[i]/error5_m1[i] < 0.0:
        y5_m1[i] = 3.0 * error5_m1[i]
        y5_m1_uplims[i] = 1

    if y5_m2[i]/error5_m2[i] < 0.0:
        y5_m2[i] = 3.0 * error5_m2[i]
        y5_m2_uplims[i] = 1
        
    if y5_m3[i]/error5_m3[i] < 0.0:
        y5_m3[i] = 3.0 * error5_m3[i]
        y5_m3_uplims[i] = 1

    if y5_m4[i]/error5_m4[i] < 0.0:
        y5_m4[i] = 3.0 * error5_m4[i]
        y5_m4_uplims[i] = 1
        
# Plot out some results
# ************
        
    
# Plot out results
# ****************************************************************************
    
fig = plt.figure(figsize=(10.0,10.0))

pos_ann = ([8.0,100.0])

# First panel

f1 = plt.axes([0.15,0.1,0.75,0.15])
f1.set_xlim(0.0,12.0)
f1.set_ylim(1.0,9999.0)
f1.tick_params(axis='both',which='both',labelsize=18)
f1.set_xlabel('Redshift',size=25)
f1.set_yscale('log')
f1.annotate('$9.0<log_{10}M_* < 9.5$',pos_ann,fontweight='bold',size=17)

x1_m1 = x1_m1 + 0.033
x1_m2 = x1_m2 + 0.1
x1_m3 = x1_m3 - 0.033
x1_m4 = x1_m4 - 0.1

f1.plot(x1_m1,y1_m1,'rs')
f1.errorbar(x1_m1,y1_m1,yerr= error1_m1,fmt='r',ls='none',uplims=y1_m1_uplims)

f1.plot(x1_m2,y1_m2,'bo')
f1.errorbar(x1_m2,y1_m2,yerr= error1_m2,fmt='b',ls='none',uplims=y1_m2_uplims)

f1.plot(x1_m3,y1_m3,'g*')
f1.errorbar(x1_m3,y1_m3,yerr= error1_m3,fmt='g',ls='none',uplims=y1_m3_uplims)

f1.plot(x1_m4,y1_m4,'cp')
f1.errorbar(x1_m4,y1_m4,yerr= error1_m4,fmt='c',ls='none',uplims=y1_m4_uplims)

# Annotate with the number of spheroids in each bin

for i in range(0,10):
    pos_ann = ([z_mean[i*5],2.0])
    f1.annotate(str(int(Ngal[5*i])),pos_ann,fontweight='bold',size=7,color='r')


# Print the positions of SMGs

sel3 = np.where( (SMG_data[:,1] > 9.0) & (SMG_data[:,1] < 9.5))

xpoint = SMG_data[:,0][sel3]
ypoint = SMG_data[:,2][sel3]

f1.plot(xpoint,ypoint,'k.',alpha=0.2)


# second panel

f1 = plt.axes([0.15,0.25,0.75,0.15])
f1.set_xlim(0.001,11.999)
f1.set_ylim(1.0,9999.0)
f1.tick_params(axis='both',which='both',labelsize=18)
#f1.set_xticks([])
f1.set_yscale('log')
pos_ann = ([8.0,100.0])
f1.annotate('$9.5<log_{10}M_* < 10.0$',pos_ann,fontweight='bold',size=17)

x2_m1 = x2_m1 + 0.033
x2_m2 = x2_m2 + 0.1
x2_m3 = x2_m3 - 0.033
x2_m4 = x2_m4 - 0.1

f1.plot(x2_m1,y2_m1,'rs')
f1.errorbar(x2_m1,y2_m1,yerr= error2_m1,fmt='r',ls='none',uplims=y2_m1_uplims)

f1.plot(x2_m2,y2_m2,'bo')
f1.errorbar(x2_m2,y2_m2,yerr= error2_m2,fmt='b',ls='none',uplims=y2_m2_uplims)

f1.plot(x2_m3,y2_m3,'g*')
f1.errorbar(x2_m3,y2_m3,yerr= error2_m3,fmt='g',ls='none',uplims=y2_m3_uplims)

f1.plot(x2_m4,y2_m4,'cp')
f1.errorbar(x2_m4,y2_m4,yerr= error2_m4,fmt='c',ls='none',uplims=y2_m4_uplims)

# Annotate with the number of spheroids in each bin

for i in range(0,10):
    pos_ann = ([z_mean[1+i*5],2.0])
    f1.annotate(str(int(Ngal[1+5*i])),pos_ann,fontweight='bold',size=7,color='r')

for i in range(12,13):
    pos_ann = ([z_mean[1+i*5],2.0])
    f1.annotate(str(int(Ngal[1+5*i])),pos_ann,fontweight='bold',size=7,color='r')

# Print the positions of SMGs

sel3 = np.where( (SMG_data[:,1] > 9.5) & (SMG_data[:,1] < 10.0))

xpoint = SMG_data[:,0][sel3]
ypoint = SMG_data[:,2][sel3]

f1.plot(xpoint,ypoint,'k.',alpha=0.2)


# Third panel

f1 = plt.axes([0.15,0.4,0.75,0.15])
f1.set_xlim(0.001,11.999)
f1.set_ylim(1.0,9999.0)
f1.set_ylabel('$SFR/M_{\odot}\ year^{-1}$',size=25)
f1.tick_params(axis='both',which='both',labelsize=18)
f1.set_yscale('log')
pos_ann = ([8.0,100.0])
f1.annotate('$10.0<log_{10}M_* < 10.5$',pos_ann,fontweight='bold',size=17)

x3_m1 = x3_m1 + 0.033
x3_m2 = x3_m2 + 0.1
x3_m3 = x3_m3 - 0.033
x3_m4 = x3_m4 - 0.1

f1.plot(x3_m1,y3_m1,'rs')
f1.errorbar(x3_m1,y3_m1,yerr= error3_m1,fmt='r',ls='none',uplims=y3_m1_uplims)

f1.plot(x3_m2,y3_m2,'bo')
f1.errorbar(x3_m2,y3_m2,yerr= error3_m2,fmt='b',ls='none',uplims=y3_m2_uplims)

f1.plot(x3_m3,y3_m3,'g*')
f1.errorbar(x3_m3,y3_m3,yerr= error3_m3,fmt='g',ls='none',uplims=y3_m3_uplims)

f1.plot(x3_m4,y3_m4,'cp')
f1.errorbar(x3_m4,y3_m4,yerr= error3_m4,fmt='c',ls='none',uplims=y3_m4_uplims)

# Annotate with the number of spheroids in each bin

for i in range(0,10):
    pos_ann = ([z_mean[2+i*5],2.0])
    f1.annotate(str(int(Ngal[2+5*i])),pos_ann,fontweight='bold',size=7,color='r')


# Print the positions of SMGs

sel3 = np.where( (SMG_data[:,1] > 10.0) & (SMG_data[:,1] < 10.5))

xpoint = SMG_data[:,0][sel3]
ypoint = SMG_data[:,2][sel3]

f1.plot(xpoint,ypoint,'k.',alpha=0.2)

# Fourth panel

f1 = plt.axes([0.15,0.55,0.75,0.15])
f1.set_xlim(0.001,11.999)
f1.set_ylim(1.0,9999.0)
f1.tick_params(axis='both',which='both',labelsize=18)
f1.set_yscale('log')
pos_ann = ([8.0,100.0])
f1.annotate('$10.5<log_{10}M_* < 11.0$',pos_ann,fontweight='bold',size=17)

x4_m1 = x4_m1 + 0.033
x4_m2 = x4_m2 + 0.1
x4_m3 = x4_m3 - 0.033
x4_m4 = x4_m4 -0.1

f1.plot(x4_m1,y4_m1,'rs')
f1.errorbar(x4_m1,y4_m1,yerr= error4_m1,fmt='r',ls='none',uplims=y4_m1_uplims)

f1.plot(x4_m2,y4_m2,'bo')
f1.errorbar(x4_m2,y4_m2,yerr= error4_m2,fmt='b',ls='none',uplims=y4_m2_uplims)

f1.plot(x4_m3,y4_m3,'g*')
f1.errorbar(x4_m3,y4_m3,yerr= error4_m3,fmt='g',ls='none',uplims=y4_m3_uplims)

f1.plot(x4_m4,y4_m4,'cp')
f1.errorbar(x4_m4,y4_m4,yerr= error4_m4,fmt='c',ls='none',uplims=y4_m4_uplims)

# Annotate with the number of spheroids in each bin

for i in range(0,9):
    pos_ann = ([z_mean[3+i*5],2.0])
    f1.annotate(str(int(Ngal[3+5*i])),pos_ann,fontweight='bold',size=7,color='r')


# Print the positions of SMGs

sel3 = np.where( (SMG_data[:,1] > 10.5) & (SMG_data[:,1] < 11.0))

xpoint = SMG_data[:,0][sel3]
ypoint = SMG_data[:,2][sel3]

f1.plot(xpoint,ypoint,'k.',alpha=0.2)

# Fifth panel

f1 = plt.axes([0.15,0.7,0.75,0.15])
f1.set_xlim(0.001,11.999)
f1.set_ylim(1.0,9999.0)
f1.tick_params(axis='both',which='both',labelsize=18)
f1.set_yscale('log')
pos_ann = ([8,100.0])
f1.annotate('$11.0<log_{10}M_* < 11.5$',pos_ann,fontweight='bold',size=17)

x5_m1 = x5_m1 + 0.033
x5_m2 = x5_m2 + 0.1
x5_m3 = x5_m3 - 0.033
x5_m4 = x5_m4 - 0.1

f1.plot(x5_m1,y5_m1,'rs')
f1.errorbar(x5_m1,y5_m1,yerr= error5_m1,fmt='r',ls='none',uplims=y5_m1_uplims)

f1.plot(x5_m2,y5_m2,'bo')
f1.errorbar(x5_m2,y5_m2,yerr= error5_m2,fmt='b',ls='none',uplims=y5_m2_uplims)

f1.plot(x5_m3,y5_m3,'g*')
f1.errorbar(x5_m3,y5_m3,yerr= error5_m3,fmt='g',ls='none',uplims=y5_m3_uplims)

f1.plot(x5_m4,y5_m4,'cp')
f1.errorbar(x5_m4,y5_m4,yerr= error5_m4,fmt='c',ls='none',uplims=y5_m4_uplims)

# Annotate with the number of spheroids in each bin

for i in range(0,8):
    pos_ann = ([z_mean[4+i*5],3.0])
    f1.annotate(str(int(Ngal[4+5*i])),pos_ann,fontweight='bold',size=10,color='r')
    

# Print the positions of SMGs

sel3 = np.where( (SMG_data[:,1] > 11.0) & (SMG_data[:,1] < 11.5))

xpoint = SMG_data[:,0][sel3]
ypoint = SMG_data[:,2][sel3]

f1.plot(xpoint,ypoint,'k.',alpha=0.2)

fig.savefig('Figure12.pdf')
plt.show()





            
            
    
    

    
        
    
    



    





