# This program contains a model of the production of spheroidal galaxies
# from SMGs. It uses the far-infrared luminosity function from Duzeviciute
# et al. (2020), converting it into a star-formation function using the
# calibration from Dudzeviciute et al. (2020). The model assumes that the
# star-formation in the galaxy lasts long enough to make a galaxy with a mass
# of 10^11 solar masses.
#
# Last edited: 7th November 2025

import numpy as np
import matplotlib.pyplot as plt
from astropy.cosmology import Planck15 as cosmo

# Constants
# ****************************************************************************
# ****************************************************************************

# Lower and upper mass limits for the SMFs

mass_lower = 10.5
mass_upper = 11.5

# mass array

del_mass = 0.01
masses = np.arange(mass_lower,mass_upper,del_mass)
nmass = len(masses)

# Parameters for the far-infrared luminosity function

alpha = -1.5
thi_star = 6.6e-6
lstar = 1.35e13

# Redshift range for model

z_low = 2.0
z_up = 3.5

# Maximum mass for galaxy 

mass_max = 1.0e11

# depletion time for model

tdep = 10.0**8.5

# Function to calculate the luminosity function
# ****************************************************************************
# ****************************************************************************

def LF_calc(xdata,thistar,Lstar):
    xdata = np.power(10.0,xdata)

    tp1 = xdata/Lstar
    tp2 = np.exp(-tp1)

    tp3 = tp1**alpha

    lf_predict = thistar * tp3 * tp2

    return lf_predict

# Read in the parameters for the SMFs
# ****************************************************************************
# ****************************************************************************

# read in parameters for the spheroids

file = 'SMF_spheroid.dat'

SMF_sp = []

SMF_sp = np.loadtxt(file)

n_sp = SMF_sp.shape[0]

print(n_sp)

# read in parameters for the disks

file = 'SMF_disks.dat'

SMF_disks = []

SMF_d = np.loadtxt(file)

n_d = SMF_d.shape[0]

# read in parameters for the irregulars

file = 'SMF_pec.dat'

SMF_pec = []

SMF_pec = np.loadtxt(file)

n_pec = SMF_pec.shape[0]

# read in parameters for the bulge-dominated galaxies

file = 'SMF_BD.dat'

SMF_BD = []

SMF_BD = np.loadtxt(file)

n_BD = SMF_BD.shape[0]

# Model for elliptical production
# ****************************************************************************
# ****************************************************************************

# set up arrays
# ****************************************************************************

# luminosity array

delm = 0.01
lum = np.arange(12.3,13.5,delm)
nlum = len(lum)

# star-formation array, using KE calibration (see notebook)

sfr = lum - 9.827

# luminosity function

thi = LF_calc(lum,thi_star,lstar)

# Lifetimes of SMGs

lifetime = mass_max/np.power(10.0,sfr)

# convert into Gyr

lifetime = lifetime/1.0e9

# Redshift and time arrays

z = np.arange(z_low,z_up,0.1)
time = cosmo.age(z).value
nz = len(z)

for i in range(0,nlum):
    print(lum[i],lifetime[i])

# Calculate the increase in the number of elliptical galaxies during this
# redshift range
# ****************************************************************************

gen = np.empty(nz)

for i in range(0,nz):
    
    if i == 0:
        delt = time[0] - time[1]
    elif i==(nz-1):
        delt = time[nz-2] - time[nz-1]
    else:
        delt = time[i-1] - time[i+1]
        delt = delt / 2.0
        
    tot = 0.0
    
    for j in range(0,nlum):
        tot = tot + thi[j] * delm * (delt/lifetime[j])
        
    gen[i] = tot

# calculate the cumulative production rate

ellip_prod = np.empty(nz)

for i in range(0,nz):
    
    tot = 0.0
    
    for j in range(i,nz):
        tot = tot + gen[j]
        
    ellip_prod[i] = tot

for i in range(0,nz):
    print(z[i],ellip_prod[i])
    
# Calculate the increase in mass due to star formation during the quenching
# process
# ****************************************************************************

mass_increase = tdep/lifetime
mass_increase = mass_increase / 1.0e9

    
# Turn this into a probability distribution

mass_final = np.arange(11.0,12.0,0.1)
npoint = len(mass_final)

prob = np.zeros(npoint)

for i in range(0,nlum):
    print(lum[i],mass_increase[i])

for i in range(0,nlum):
    tp = 1.0e11 * (1.0+mass_increase[i])
    print(tp)
    tp2 = np.log10(tp)
    print(tp2)
    ind = int( (tp2-11.0)/0.1)
    print(ind)
    
    if ind < npoint:
        prob[ind] = thi[i] / lifetime[i]
    
# Normalise probability distribution

prob = prob / np.sum(prob)

# Calculate space-density of massive galaxies versus redshift
# ****************************************************************************
# ****************************************************************************

# First for the spheroids

x_sp = np.empty(n_sp)
y_sp = np.empty(n_sp)

for i in range(0,n_sp):
    x_sp[i] = ( SMF_sp[i,0] + SMF_sp[i,1]) / 2.0

    if SMF_sp[i,2] < 1.5:
        alpha = SMF_sp[i,3]
        mstar = np.power(10.0,SMF_sp[i,6])
        thi = np.power(10.0,SMF_sp[i,9])
        
        tot = 0.0
        
        for j in range(0,nmass):
            mass = np.power(10.0,masses[j])
            
            tp1 = mass/mstar
            tp2 = np.exp(-tp1)
            
            tp3 = thi * tp2 * tp1**(alpha+1.0)
            tot = tot + 2.30 * tp3 * del_mass
        
        y_sp[i] = tot
            
    if (SMF_sp[i,2] > 1.5) & (SMF_sp[i,2] < 2.5):
        alpha1 = SMF_sp[i,4]
        alpha2 = SMF_sp[i,5]
        mstar = np.power(10.0,SMF_sp[i,6])
        thi1 = np.power(10.0,SMF_sp[i,7])
        thi2 = np.power(10.0,SMF_sp[i,8])

        tot = 0.0
        
        for j in range(0,nmass):
            mass = np.power(10.0,masses[j])
            
            tp1 = mass/mstar
            tp2 = np.exp(-tp1)
            
            tp3 = thi1 * tp2 * tp1**(alpha1+1.0)
            tp4 = thi2 * tp2 * tp1**(alpha2+1.0)
            
            tot = tot + 2.30 * (tp3 + tp4) * del_mass
        
        y_sp[i] = tot
        
    if (SMF_sp[i,2] > 2.5) & (SMF_sp[i,2] < 3.5):
        alpha1 = SMF_sp[i,4]
        alpha2 = SMF_sp[i,5]
        mstar = np.power(10.0,SMF_sp[i,6])
        thi = np.power(10.0,SMF_sp[i,9])

        tot = 0.0
        
        for j in range(0,nmass):
            mass = np.power(10.0,masses[j])
            
            tp1 = mstar/mass
            tp2 = tp1**(alpha1+1.0)
            tp3 = tp1**(alpha2+1.0)
                        
            tot = tot + thi/ (tp2 + tp3) * del_mass
        
        y_sp[i] = tot
        
# Now for the disks

x_d = np.empty(n_d)
y_d = np.empty(n_d)

for i in range(0,n_d):
    x_d[i] = ( SMF_d[i,0] + SMF_d[i,1]) / 2.0

    if SMF_d[i,2] < 1.5:
        alpha = SMF_d[i,3]
        mstar = np.power(10.0,SMF_d[i,6])
        thi = np.power(10.0,SMF_d[i,9])
        
        tot = 0.0
        
        for j in range(0,nmass):
            mass = np.power(10.0,masses[j])
            
            tp1 = mass/mstar
            tp2 = np.exp(-tp1)
            
            tp3 = thi * tp2 * tp1**(alpha+1.0)
            tot = tot + 2.30 * tp3 * del_mass
        
        y_d[i] = tot
            
    if (SMF_d[i,2] > 1.5) & (SMF_d[i,2] < 2.5):
        alpha1 = SMF_d[i,4]
        alpha2 = SMF_d[i,5]
        mstar = np.power(10.0,SMF_d[i,6])
        thi1 = np.power(10.0,SMF_d[i,7])
        thi2 = np.power(10.0,SMF_d[i,8])
        
        tot = 0.0
        
        for j in range(0,nmass):
            mass = np.power(10.0,masses[j])
            
            tp1 = mass/mstar
            tp2 = np.exp(-tp1)
            
            tp3 = thi1 * tp2 * tp1**(alpha1+1.0)
            tp4 = thi2 * tp2 * tp1**(alpha2+1.0)
            
            tot = tot + 2.30 * (tp3 + tp4) * del_mass
            
        y_d[i] = tot
        
    if (SMF_d[i,2] > 2.5) & (SMF_d[i,2] < 3.5):
        alpha1 = SMF_d[i,4]
        alpha2 = SMF_d[i,5]
        mstar = np.power(10.0,SMF_d[i,6])
        thi = np.power(10.0,SMF_d[i,9])
        
        tot = 0.0
        
        for j in range(0,nmass):
            mass = np.power(10.0,masses[j])
            
            tp1 = mstar/mass
            tp2 = tp1**(alpha1+1.0)
            tp3 = tp1**(alpha2+1.0)
                        
            tot = tot + thi/ (tp2 + tp3) * del_mass

        y_d[i] = tot
        
# Now for the irregulars

x_pec = np.empty(n_pec)
y_pec = np.empty(n_pec)

for i in range(0,n_pec):
    x_pec[i] = ( SMF_pec[i,0] + SMF_pec[i,1]) / 2.0

    if SMF_pec[i,2] < 1.5:
        alpha = SMF_pec[i,3]
        mstar = np.power(10.0,SMF_pec[i,6])
        thi = np.power(10.0,SMF_pec[i,9])
        
        tot = 0.0
        
        for j in range(0,nmass):
            mass = np.power(10.0,masses[j])
            
            tp1 = mass/mstar
            tp2 = np.exp(-tp1)
            
            tp3 = thi * tp2 * tp1**(alpha+1.0)
            tot = tot + 2.30 * tp3 * del_mass
        
        y_pec[i] = tot
            
    if (SMF_pec[i,2] > 1.5) & (SMF_pec[i,2] < 2.5):
        alpha1 = SMF_pec[i,4]
        alpha2 = SMF_pec[i,5]
        mstar = np.power(10.0,SMF_pec[i,6])
        thi1 = np.power(10.0,SMF_pec[i,7])
        thi2 = np.power(10.0,SMF_pec[i,8])

        tot = 0.0
        
        for j in range(0,nmass):
            mass = np.power(10.0,masses[j])
            
            tp1 = mass/mstar
            tp2 = np.exp(-tp1)
            
            tp3 = thi1 * tp2 * tp1**(alpha1+1.0)
            tp4 = thi2 * tp2 * tp1**(alpha2+1.0)
            
            tot = tot + 2.30 * (tp3 + tp4) * del_mass
        
        y_pec[i] = tot
            
    if (SMF_pec[i,2] > 2.5) & (SMF_pec[i,2] < 3.5):
        alpha1 = SMF_pec[i,4]
        alpha2 = SMF_pec[i,5]
        mstar = np.power(10.0,SMF_pec[i,6])
        thi = np.power(10.0,SMF_pec[i,9])

        tot = 0.0
        
        for j in range(0,nmass):
            mass = np.power(10.0,masses[j])
            
            tp1 = mstar/mass
            tp2 = tp1**(alpha1+1.0)
            tp3 = tp1**(alpha2+1.0)
                        
            tot = tot + thi/ (tp2 + tp3) * del_mass
        
        y_pec[i] = tot
            
# Now for the bulge-dominated galaxies

x_BD = np.empty(n_BD)
y_BD = np.empty(n_BD)

for i in range(0,n_BD):
    x_BD[i] = ( SMF_BD[i,0] + SMF_BD[i,1]) / 2.0

    if SMF_BD[i,2] < 1.5:
        alpha = SMF_BD[i,3]
        mstar = np.power(10.0,SMF_BD[i,6])
        thi = np.power(10.0,SMF_BD[i,9])
        
        tot = 0.0
        
        for j in range(0,nmass):
            mass = np.power(10.0,masses[j])
            
            tp1 = mass/mstar
            tp2 = np.exp(-tp1)
            
            tp3 = thi * tp2 * tp1**(alpha+1.0)
            tot = tot + 2.30 * tp3 * del_mass
        
        y_BD[i] = tot
            
    if (SMF_BD[i,2] > 1.5) & (SMF_BD[i,2] < 2.5):
        alpha1 = SMF_BD[i,4]
        alpha2 = SMF_BD[i,5]
        mstar = np.power(10.0,SMF_BD[i,6])
        thi1 = np.power(10.0,SMF_BD[i,7])
        thi2 = np.power(10.0,SMF_BD[i,8])

        tot = 0.0
        
        for j in range(0,nmass):
            mass = np.power(10.0,masses[j])
            
            tp1 = mass/mstar
            tp2 = np.exp(-tp1)
            
            tp3 = thi1 * tp2 * tp1**(alpha1+1.0)
            tp4 = thi2 * tp2 * tp1**(alpha2+1.0)
            
            tot = tot + 2.30 * (tp3 + tp4) * del_mass
        
        y_BD[i] = tot
            
    if (SMF_BD[i,2] > 2.5) & (SMF_BD[i,2] < 3.5):
        alpha1 = SMF_BD[i,4]
        alpha2 = SMF_BD[i,5]
        mstar = np.power(10.0,SMF_BD[i,6])
        thi = np.power(10.0,SMF_BD[i,9])

        tot = 0.0
        
        for j in range(0,nmass):
            mass = np.power(10.0,masses[j])
            
            tp1 = mstar/mass
            tp2 = tp1**(alpha1+1.0)
            tp3 = tp1**(alpha2+1.0)
                        
            tot = tot + thi/ (tp2 + tp3) * del_mass
            
        y_BD[i] = tot
        
# Plot out results
# ****************************************************************************

fig = plt.figure(figsize=(10.0,10.0))

f1 = plt.axes([0.2,0.2,0.7,0.7])
f1.set_xlim(0.0,6.0)
f1.set_ylim(1.0e-6,1.0e-2)
f1.tick_params(axis='both',which='both',labelsize=25)
f1.set_xlabel('Redshift',size=25)
f1.set_ylabel('$\Phi/Mpc^{-3}$',size=25)
f1.set_yscale('log')
    
f1.plot(x_sp,y_sp,'r-',linewidth=3)
f1.plot(x_d,y_d,'b-',linewidth=3)
f1.plot(x_pec,y_pec,'g-',linewidth=3)
f1.plot(x_BD,y_BD,'c-',linewidth=3)

f1.plot(z,ellip_prod,'k--',linewidth=5)

fig.savefig('Figure15.pdf')
plt.show()

# Now plot the probability distribution for the final mass

# Only plot positive values

nprob = len(prob)

nl = 0

for i in range(0,nprob):
    
    if prob[i] > 0.0001:
        nl = nl + 1
        

mass_final2 = np.empty(nl)
prob2 = np.empty(nl)

nl = 0

for i in range(0,nprob):
    
    if prob[i] > 0.0001:
        prob2[nl] = prob[i]
        mass_final2[nl] = mass_final[i]
        nl = nl + 1

fig = plt.figure(figsize=(10.0,10.0))

f1 = plt.axes([0.2,0.2,0.7,0.7])
f1.set_xlim(11.0,12.0)
f1.set_ylim(0.0,0.5)
f1.tick_params(axis='both',which='both',labelsize=25)
f1.set_xlabel('$Final\ Mass/M_{\odot}$',size=25)
f1.set_ylabel('$Fraction$',size=25)

f1.plot(mass_final2,prob2,'k-',linewidth=3)
fig.savefig('Figure16.pdf')

    

