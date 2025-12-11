# This program calculates the evolution of the dust density using my estimates
# of the mean dust-to-gas ratio and the stellar mass functions from Shuntov
# et al. (2025)
# 
# Last edited: 11th September 2025

import numpy as np
import matplotlib.pyplot as plt
from astropy.cosmology import Planck15 as cosmo

# Constants and other key variables
# ****************************************************************************
# ****************************************************************************

# Number of redshift bins

nbin = 14

# Critical density in solar masses per cubic Mpc

rho_den = 1.273e11

### Redshift slices and the bin width

catalogInfo = [{"zmin":0.2, "zmax":0.5, "bin_width": 0.5, "nb":5},\
               {"zmin":0.5, "zmax":0.8, "bin_width": 0.5, "nb": 5},\
               {"zmin":0.8, "zmax":1.1, "bin_width": 0.5, "nb": 5},\
               {"zmin":1.1, "zmax":1.5, "bin_width": 0.5, "nb": 5},\
               {"zmin":1.5, "zmax":2.0, "bin_width": 0.5, "nb": 5},\
               {"zmin":2.0, "zmax":2.5, "bin_width": 0.5, "nb": 5},\
               {"zmin":2.5, "zmax":3.0, "bin_width": 0.5, "nb": 5},\
               {"zmin":3.0, "zmax":3.5, "bin_width": 0.5, "nb": 5},\
               {"zmin":3.5, "zmax":4.5, "bin_width": 0.5, "nb": 5},\
               {"zmin":4.5, "zmax":5.5, "bin_width": 0.5, "nb": 5},\
               {"zmin":5.5, "zmax":6.5, "bin_width": 0.5, "nb": 5},\
               {"zmin":6.5, "zmax":7.5, "bin_width": 0.5, "nb": 5},\
               {"zmin":7.5, "zmax":8.5, "bin_width": 0.5, "nb": 5},\
               {"zmin":8.5, "zmax":12.0, "bin_width": 0.5, "nb": 5},\
                   ]

# Read in data
# ****************************************************************************
# *****************************************************************************

# Read in stellar mass function
# *****************************************************************************

file = 'SMF_cutdown.dat'

mass_function_all = []

mass_function_all = np.loadtxt(file)

print(mass_function_all.shape)

# Read in mean dust-to-stellar mass ratios
# *****************************************************************************

file = 'dust_to_stellar_mass_results_SIMSTACK_averaging_method.dat'

dust_to_stars = []

dust_to_stars = np.loadtxt(file)

print(dust_to_stars.shape)

# Read in the data from the Madau and Dickinson ARAA paper
# ****************************************************************************

data_DM_col1 = np.array([])
data_DM_col2 = np.array([])
data_DM_col3 = np.array([])
data_DM_col4 = np.array([])
data_DM_col5 = np.array([])
data_DM_col6 = np.array([])
file = 'dickinson_madau.dat'

data = open(file,'r')

ic = 0

for line in data.readlines():
    info = line.split()
    data_DM_col1 = np.append(data_DM_col1,float(info[0]))
    data_DM_col2 = np.append(data_DM_col2,float(info[1]))
    data_DM_col3 = np.append(data_DM_col3,float(info[2]))
    data_DM_col4 = np.append(data_DM_col4,float(info[3]))
    data_DM_col5 = np.append(data_DM_col5,float(info[4]))
    data_DM_col6 = np.append(data_DM_col6,float(info[5]))
    ic = ic + 1
data.close()

n_DM = len(data_DM_col1)


# Read in the data from Peroux et al.
# ****************************************************************************

file = 'Peroux_data.dat'

data = open(file,'r')

ic = 0

data_Peroux_col1 = np.array([])
data_Peroux_col2 = np.array([])
data_Peroux_col3 = np.array([])
data_Peroux_col4 = np.array([])
data_Peroux_col5 = np.array([])
data_Peroux_col6 = np.array([])

for line in data.readlines():  
    info = line.split()
    data_Peroux_col1 = np.append(data_Peroux_col1,float(info[0]))
    data_Peroux_col2 = np.append(data_Peroux_col2,float(info[1]))
    data_Peroux_col3 = np.append(data_Peroux_col3,float(info[2]))
    data_Peroux_col4 = np.append(data_Peroux_col4,float(info[3]))
    data_Peroux_col5 = np.append(data_Peroux_col5,float(info[4]))
    data_Peroux_col6 = np.append(data_Peroux_col6,float(info[5]))
    
    ic = ic + 1
data.close()

# Read in the results from Eales and Ward (2024)

file = 'Eales_Ward_results.dat'

EW_results = []

EW_results = np.loadtxt(file)

print(EW_results.shape)

# Read in dust model
# ****************************************************************************

file = 'results_dust_model.dat'

data = open(file,'r')

ic = 0

xmodel = np.array([])
ymodel1 = np.array([])
ymodel2 = np.array([])
ymodel3 = np.array([])
ymodel4 = np.array([])

for line in data.readlines():
    info = line.split()
    xmodel = np.append(xmodel,float(info[0]))
    ymodel1 = np.append(ymodel1,float(info[1]))
    ymodel2 = np.append(ymodel2,float(info[2]))
    ymodel3 = np.append(ymodel3,float(info[3]))

# Set up variables and arrays
# ****************************************************************************
# ****************************************************************************

# start mass bin

lowest_mass = 9.0

# set up bin limits

mass_lower_limit = np.zeros(100)
mass_upper_limit = np.zeros(100)

z_lower_limit = np.zeros(100)

z_upper_limit = np.zeros(100)

nbin = 0

n_z_slice = len(catalogInfo)

for i in range(0,n_z_slice):
    width = catalogInfo[i]["bin_width"]
    nind = catalogInfo[i]["nb"]

    z1 = catalogInfo[i]["zmin"]
    z2 = catalogInfo[i]["zmax"]

    for j in range(0,nind):
        z_lower_limit[nbin+j] = z1
        z_upper_limit[nbin+j] = z2

        mass_lower_limit[nbin+j] = lowest_mass + float(j) * width
        mass_upper_limit[nbin+j] = mass_lower_limit[nbin+j] + width

    nbin = nbin + nind

# Separate data and stellar mass function into redshift bins
# ****************************************************************************

# Sort out the stellar mass function

n_thi = n_z_slice + 1

thi = np.empty((n_thi,100,4))

for i in range(0,n_thi):
    
    for j in range(0,100):
        thi[i,j,:] = mass_function_all[i*100+j,:] 
        
# average the stellar mass function in the final two redshfit bins

for j in range(0,100):
    thi[13,j,:] = (thi[13,j,:] + thi[14,j,:]) / 2.0
    
# Print out the stellar mass function for the 3rd and 4th bins

#for j in range(0,100):
#    print(thi[2,j,0],thi[2,j,1],thi[3,j,1])

# Calculate mean dust density
# ****************************************************************************
# ****************************************************************************

dust_mass_density = np.zeros(nbin)
error_mass_density = np.zeros(nbin)
redshift = np.empty(nbin)

dust_mass_density_Tvar = np.empty(nbin)
error_mass_density_Tvar = np.empty(nbin)

dust_mass_density_Tmodel = np.empty(nbin)
error_mass_density_Tmodel = np.empty(nbin)

for i in range(0,n_z_slice):

    work = np.zeros((5,10))
    
    for j in range(0,5):
        icount = i*5 + j
        
        z1 = z_lower_limit[icount]
        z2 = z_upper_limit[icount]
        m1 = mass_lower_limit[icount]
        m2 = mass_upper_limit[icount]
   
        redshift[i] = (z1 + z2)/2.0
        
        dts = dust_to_stars[icount,4]
        error_dts = dust_to_stars[icount,5]
        
        dts_Tvar = dust_to_stars[icount,6]
        error_dts_Tvar = dust_to_stars[icount,7]
        
        dts_Tmodel = dust_to_stars[icount,8]
        error_dts_Tmodel = dust_to_stars[icount,9]
        
# Integrate mass-weigted stellar mass function

        total = 0.0
        error = 0.0
        
        for k in range(0,100):
            mass = thi[i,k,0]
            if (mass > m1) and (mass < m2): 
                tp1 = np.power(10.0,mass)
                
                tp2 = mass - thi[i,k-1,0]
                tp3 = thi[i,k+1,0] - mass
                delm = (tp2 + tp3) / 2.0
                
                total = total + thi[i,k,1] * delm * tp1
                
                err1 = thi[i,k,1] - thi[i,k,2]
                err2 = thi[i,k,3] - thi[i,k,1]
                err = (err1 + err2)/2.0
                error = error + err*err
            
        work[j,0] = total
        work[j,1] = np.sqrt(error)
        work[j,2] = dts
        work[j,3] = error_dts
        
        work[j,4] = dts_Tvar
        work[j,5] = error_dts_Tvar
        
        work[j,6] = dts_Tmodel
        work[j,7] = error_dts_Tmodel
    
# Now combine these values to get the mean dust-to-stellar-mass ratio

    tot = 0.0
    err = 0.0
    
    tot_Tvar = 0.0
    err_Tvar = 0.0
    
    tot_Tmodel = 0.0
    err_Tmodel = 0.0
    
    for j in range(0,5):
        dust_mass = work[j,2] * work[j,0]
        
        if dust_mass < 0.0:
            dust_mass=0.0
        
#        print('here ',i,dust_mass,work[j,2],work[j,0])
        tot = tot + dust_mass
        
# calculate error

        tp1 = work[j,3] / work[j,2]
        tp2 = work[j,1] / work[j,0]
        tp3 = tp1 * tp1 + tp2*tp2
        tp3 = np.sqrt(tp3)
        tp4 = tp3 * dust_mass
        
        err = err + tp4*tp4
        
# do the same with the temperature variation

        dust_mass = work[j,4] * work[j,0]
        
        if dust_mass < 0.0:
            dust_mass=0.0
        
        tot_Tvar = tot_Tvar + dust_mass

# calculate error

        tp1 = work[j,5] / work[j,4]
        tp2 = work[j,1] / work[j,0]
        tp3 = tp1 * tp1 + tp2*tp2
        tp3 = np.sqrt(tp3)
        tp4 = tp3 * dust_mass
        
        err_Tvar = err_Tvar + tp4*tp4
        
# Now do the same with the other temperature model

        dust_mass = work[j,6] * work[j,0]
        
        if dust_mass < 0.0:
            dust_mass=0.0
            
        tot_Tmodel = tot_Tmodel + dust_mass

# calculate error

        tp1 = work[j,7] / work[j,6]
        tp2 = work[j,1] / work[j,0]
        tp3 = tp1 * tp1 + tp2*tp2
        tp3 = np.sqrt(tp3)
        tp4 = tp3 * dust_mass  
        
        err_Tmodel = err_Tmodel + tp4*tp4
    
    dust_mass_density[i] = tot
    error_mass_density[i] = np.sqrt(err)
    
    dust_mass_density_Tvar[i] = tot_Tvar
    error_mass_density_Tvar[i] = np.sqrt(err_Tvar)
    
    dust_mass_density_Tmodel[i] = tot_Tmodel
    error_mass_density_Tmodel[i] = np.sqrt(err_Tmodel)
     
    print(i,dust_mass_density_Tmodel[i],error_mass_density_Tmodel[i],\
          dust_mass_density_Tvar[i],error_mass_density_Tvar[i])
    

dust_mass_density = dust_mass_density / 1.0e5
error_mass_density = error_mass_density / 1.0e5

dust_mass_density_Tvar = dust_mass_density_Tvar / 1.0e5
error_mass_density_Tvar = error_mass_density_Tvar / 1.0e5

dust_mass_density_Tmodel = dust_mass_density_Tmodel / 1.0e5
error_mass_density_Tmodel = error_mass_density_Tmodel / 1.0e5

        
# Plot out results
# ****************************************************************************

redshift_Tvar = redshift - 0.1
redshift_Tmodel = redshift + 0.1

# set final measurement so it doesn't appear on graph

#dust_mass_density[13] = 0.0000001
#error_mass_density[13] = 0.0

#dust_mass_density_Tvar[13] = 0.0000001
#error_mass_density_Tvar[13] = 0.0

#dust_mass_density_Tmodel[13] = 0.000001
#error_mass_density_Tmodel[13] = 0.0

for i in range(0,14):
    print(dust_mass_density_Tvar[i],error_mass_density_Tvar[i])

fig = plt.figure(figsize=(10.0,10.0))
f1 = plt.axes([0.15,0.15,0.8,0.8])
f1.set_xlim(0.0,12.0)
f1.set_ylim(0.0001,10.0)
f1.set_yscale('log')
f1.set_xlabel('Redshift',size=25)
f1.set_ylabel('$Dust\ Density/10^5\ M_{\odot}\ Mpc^{-3}$',size=25)
f1.tick_params(axis='both',labelsize=20)

f1.plot(redshift,dust_mass_density,'bs',markersize=10)
f1.errorbar(redshift,dust_mass_density,yerr=error_mass_density,fmt='b',ls='none')

f1.plot(redshift_Tvar,dust_mass_density_Tvar,'co',alpha=1.0,markersize=10)
f1.errorbar(redshift_Tvar,dust_mass_density_Tvar,\
            yerr=error_mass_density_Tvar,fmt='c',alpha=1.0,ls='none')

f1.plot(redshift_Tmodel,dust_mass_density_Tmodel,'rs',alpha=1.0,markersize=10)
f1.errorbar(redshift_Tmodel,dust_mass_density_Tmodel,\
                yerr=error_mass_density_Tmodel,fmt='r',alpha=1.0,ls='none')

# plot redshift ranges

xp = np.empty(1)
yp = np.empty(1)
xl = np.empty(2)
yl = np.empty(2)

for i in range(0,n_z_slice):
    icount = i * 5
    xl[0] = z_lower_limit[icount]
    xl[1] = z_upper_limit[icount]
    yl[0] = dust_mass_density[i]
    yl[1] = dust_mass_density[i]
    f1.plot(xl,yl,'b-')
    
# Plot 3-sigma limit for highest redshift bin
# Note that this is the SIMSTACK value

xp[0] = 10.25
yp[0] = 0.02649
#f1.plot(xp,yp,'bv',markersize=10)

xl[0] = 8.5
xl[1] = 12.0
yl[0] = 0.02649
yl[1] = 0.02649
#f1.plot(xl,yl,'b-')

# Plot 3-sigma limit for highest redshift bin for the temperature variation 
# model
# Note that this is the SIMSTACK value

xp[0] = 10.25
yp[0] = 6.07e-4
#f1.plot(xp,yp,'cv',alpha=1.0,markersize=10)

xl[0] = 8.5
xl[1] = 12.0
yl[0] = 6.07e-4
yl[1] = 6.07e-4
#f1.plot(xl,yl,'c-',alpha=1.0)

# Plot 3-sigma limit for highest redshift bin with two-temperature model
# Note that this is the SIMSTACK value

xp[0] = 10.25
yp[0] = 0.01049
#f1.plot(xp,yp,'rv',alpha=1.0,markersize=10)

xl[0] = 8.5
xl[1] = 12.0
yl[0] = 0.01049
yl[1] = 0.01049
#f1.plot(xl,yl,'r-',alpha=1.0)

# Plot results of the chemical evolution model
# ****************************************************************************

f1.plot(xmodel,ymodel2,'g-')
f1.plot(xmodel,ymodel1,'m-')
f1.plot(xmodel,ymodel3,'y-')


# Plot results from Eales and Ward (2024)

x_ew = EW_results[:,2]
y_ew = EW_results[:,3]
errors_ew = EW_results[:,4]

f1.plot(x_ew,y_ew,'mo',markersize=10)
f1.errorbar(x_ew,y_ew,yerr=errors_ew,fmt='m',ls='none')

for i in range(0,10):
    xl[0] = EW_results[i,0]
    xl[1] = EW_results[i,1]
    yl[0] = y_ew[i]
    yl[1] = y_ew[i]
    f1.plot(xl,yl,'m-')
    

pos_ann = ([8.0,3.16])

f1.annotate('COLD',pos_ann,fontweight='bold',size=20,color='b')

pos_ann = ([8.0,1.77])

f1.annotate('STAR-FORMER',pos_ann,fontweight='bold',size=20,color='r')

pos_ann = ([8.0,1.0])

f1.annotate('EVOLVE',pos_ann,fontweight='bold',size=20,color='c')

xl = np.empty(2)
yl = np.empty(2)

xl[0] = 7.0
xl[1] = 7.5
yl[0] = 3.63
yl[1] = 3.63
f1.plot(xl,yl,'b-')

yl[0] = 2.04
yl[1] = 2.04
f1.plot(xl,yl,'r-')

yl[0] = 1.15
yl[1] = 1.15
f1.plot(xl,yl,'c-')


fig.savefig('Figure8.pdf')
plt.show()








    
    
    




            
            
    
    

    
        
    
    



    





