# This program plots flux versus redshift in a separate panel for
# each redshift slice. This version of the program uses the results
# of the basic method and adds (nin quadrature) the errors from
# the bootstrap versions of the JWST catalogue.
#
# Last edited: 2nd December 2025

import numpy as np
import matplotlib.pyplot as plt

# First collate the results from the basic stacking method
# ****************************************************************************
# ****************************************************************************

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

        
# Read in stacking data catalogue and catalogue of noise from the artificial
# JWST data
# ****************************************************************************

file='results_all_slices_monte_carlo.dat'

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

# Now read in the noise file

file='noise_values_artificial_data_basic_stacking.dat'

data = open(file,'r')

ic = 0

for line in data.readlines():  
    info = line.split()
    art_noise = np.append(art_noise,float(info[4]))
    
# Now add on the noise from the artificial JWST data in quadrature
# ****************************************************************************

for i in range(0,nbin):
    tp1 = Error[i]
    tp2 = art_noise[i]
    tp3 = tp1*tp1 + tp2*tp2
    Error[i] = np.sqrt(tp3)
    
# Print out the results for all the bins
# ****************************************************************************

for i in range(0,nbin):
    print(zlow[i],zup[i],Mlow[i],Mup[i],S850[i],Error[i])
    
# Find the mean redshifts and stellar masses in each bin. Note that this is
# ****************************************************************************

Mass_mean = np.empty(nbin)
z_mean = np.empty(nbin)

for i in range(0,nbin):
    Mass_mean[i] = (Mlow[i] + Mup[i]) / 2.0
    z_mean[i] = (zlow[i] + zup[i]) / 2.0


# Now arrange the results for each bin for plotting
# ****************************************************************************

# First redshift bin

x1 = np.empty(5)
y1 = np.empty(5)
error1 = np.empty(5)
xline = np.empty(2)
yline = np.empty(2)

print("bin1")

for i in range(0,5):
    x1[i] = Mass_mean[i]
    y1[i] = S850[i]
    error1[i] = Error[i]
    print(x1[i],y1[i],error1[i])
    
    
# second redshift bin

x2 = np.empty(5)
y2 = np.empty(5)
error2 = np.empty(5)

print("Second redshift bin")

for i in range(0,5):
    x2[i] = Mass_mean[5+i]
    y2[i] = S850[5+i]
    error2[i] = Error[5+i]
    print(x2[i],y2[i],error2[i])
    
# third redshift bin

x3 = np.empty(5)
y3 = np.empty(5)
error3 = np.empty(5)

print("3rd bin")


for i in range(0,5):
    x3[i] = Mass_mean[10+i]
    y3[i] = S850[10+i]
    error3[i] = Error[10+i]
    print(x3[i],y3[i],error3[i])

# fourth redshift bin

x4 = np.empty(5)
y4 = np.empty(5)
error4 = np.empty(5)

print("4th bin")


for i in range(0,5):
    x4[i] = Mass_mean[15+i]
    y4[i] = S850[15+i]
    error4[i] = Error[15+i]
    print(x4[i],y4[i],error4[i])

# fifth redshift bin

x5 = np.empty(5)
y5 = np.empty(5)
error5 = np.empty(5)

print("5th bin")


for i in range(0,5):
    x5[i] = Mass_mean[20+i]
    y5[i] = S850[20+i]
    error5[i] = Error[20+i]
    print(x5[i],y5[i],error5[i])

# sixth redshift bin

x6 = np.empty(5)
y6 = np.empty(5)
error6 = np.empty(5)

print("6th bin")


for i in range(0,5):
    x6[i] = Mass_mean[25+i]
    y6[i] = S850[25+i]
    error6[i] = Error[25+i]
    print(x6[i],y6[i],error6[i])

# seventh redshift bin

x7 = np.empty(5)
y7 = np.empty(5)
error7 = np.empty(5)

print("7th bin")

for i in range(0,5):
    x7[i] = Mass_mean[30+i]
    y7[i] = S850[30+i]
    error7[i] = Error[30+i]
    print(x7[i],y7[i],error7[i])
    
# eigth redshift bin

x8 = np.empty(5)
y8 = np.empty(5)
error8 = np.empty(5)

print("8th bin")

for i in range(0,5):
    x8[i] = Mass_mean[35+i]
    y8[i] = S850[35+i]
    error8[i] = Error[35+i]
    print(x8[i],y8[i],error8[i])
       
# ninth redshift bin

x9 = np.empty(5)
y9 = np.empty(5)
error9 = np.empty(5)

print("Ninth bin")

for i in range(0,5):
    x9[i] = Mass_mean[40+i]
    y9[i] = S850[40+i]
    error9[i] = Error[40+i]
    print(x9[i],y9[i],error9[i])
    
       
# tenth redshift bin

x10 = np.empty(5)
y10 = np.empty(5)
error10 = np.empty(5)

print("Tenth bin")

for i in range(0,5):
    x10[i] = Mass_mean[45+i]
    y10[i] = S850[45+i]   
    error10[i] = Error[45+i]
    print(x10[i],y10[i],error10[i])
    
  # 11th redshift bin

x11 = np.empty(5)
y11 = np.empty(5)
error11 = np.empty(5)

print("11th bin")

for i in range(0,5):
    x11[i] = Mass_mean[50+i]
    y11[i] = S850[50+i]   
    error11[i] = Error[50+i]
    print(x11[i],y11[i],error11[i])

# 12th redshift bin

x12 = np.empty(5)
y12 = np.empty(5)
error12 = np.empty(5)

print("12th bin")

for i in range(0,5):
    x12[i] = Mass_mean[55+i]
    y12[i] = S850[55+i]   
    error12[i] = Error[55+i]
    print(x12[i],y12[i],error12[i])

# 13th redshift bin

x13 = np.empty(5)
y13 = np.empty(5)
error13 = np.empty(5)

print("13th bin")

for i in range(0,5):
    x13[i] = Mass_mean[60+i]
    y13[i] = S850[60+i]   
    error13[i] = Error[60+i]
    print(x13[i],y13[i],error13[i])

# 14th redshift bin

x14 = np.empty(5)
y14 = np.empty(5)
error14 = np.empty(5)

print("14th bin")

for i in range(0,5):
    x14[i] = Mass_mean[65+i]
    y14[i] = S850[65+i]   
    error14[i] = Error[65+i]
    print(x14[i],y14[i],error14[i])

# Plot out results
# ****************************************************************************

# Put 3sigma limits for negative fluxes

y1_uplims = np.array([0,0,0,0,0],dtype=bool)
y2_uplims = np.array([0,0,0,0,0],dtype=bool)
y3_uplims = np.array([0,0,0,0,0],dtype=bool)
y4_uplims = np.array([0,0,0,0,0],dtype=bool)
y5_uplims = np.array([0,0,0,0,0],dtype=bool)
y6_uplims = np.array([0,0,0,0,0],dtype=bool)
y7_uplims = np.array([0,0,0,0,0],dtype=bool)
y8_uplims = np.array([0,0,0,0,0],dtype=bool)
y9_uplims = np.array([0,0,0,0,0],dtype=bool)
y10_uplims = np.array([0,0,0,0,0],dtype=bool)
y11_uplims = np.array([0,0,0,0,0],dtype=bool)
y12_uplims = np.array([0,0,0,0,0],dtype=bool)
y13_uplims = np.array([0,0,0,0,0],dtype=bool)
y14_uplims = np.array([0,0,0,0,0],dtype=bool)

for i in range(0,5):
    
    if y1[i] < 0.0:
        y1[i] = 3.0*error1[i]
        y1_uplims[i] = 1
        
    if y2[i] < 0.0:
        y2[i] = 3.0*error2[i]
        y2_uplims[i] = 1

    if y3[i] < 0.0:
        y3[i] = 3.0*error3[i]
        y3_uplims[i] = 1
        
    if y4[i] < 0.0:
        y4[i] = 3.0*error4[i]
        y4_uplims[i] = 1
        
    if y5[i] < 0.0:
        y5[i] = 3.0*error5[i]
        y5_uplims[i] = 1
        
    if y6[i] < 0.0:
        y6[i] = 3.0*error6[i]
        y6_uplims[i] = 1
        
    if y7[i] < 0.0:
        y7[i] = 3.0*error7[i]
        y7_uplims[i] = 1
        
    if y8[i] < 0.0:
        y8[i] = 3.0*error8[i]
        y8_uplims[i] = 1

    if y9[i] < 0.0:
        y9[i] = 3.0*error9[i]
        y9_uplims[i] = 1

    if y10[i] < 0.0:
        y10[i] = 3.0*error10[i]
        y10_uplims[i] = 1
        
    if y11[i] < 0.0:
        y11[i] = 3.0*error11[i]
        y11_uplims[i] = 1

    if y12[i] < 0.0:
        y12[i] = 3.0*error12[i]
        y12_uplims[i] = 1

    if y13[i] < 0.0:
        y13[i] = 3.0*error13[i]
        y13_uplims[i] = 1
        
    if y14[i] < 0.0:
        y14[i] = 3.0*error14[i]
        y14_uplims[i] = 1


# Now collate the results for SIMSTACK
# ****************************************************************************
# ****************************************************************************
        
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

        
# Read in stacking data catalogue and catalogue of noise from the artificial
# JWST data
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

# Now read in the noise file

file='results_from_Matts_simulation.dat'

data = open(file,'r')

ic = 0

for line in data.readlines():  
    info = line.split()
    art_noise = np.append(art_noise,float(info[5]))
    
    ic = ic + 1
data.close()

# Change the units to mJy
# ****************************************************************************

S850 = S850 * 1000.0
Error = Error * 1000.0
art_noise = art_noise * 1000.0
    
# Now add on the noise from the artificial JWST data in quadrature
# ****************************************************************************

for i in range(0,nbin):
    tp1 = Error[i]
    tp2 = art_noise[i]
    tp3 = tp1*tp1 + tp2*tp2
    Error[i] = np.sqrt(tp3)
    
# Print out the results for all the bins
# ****************************************************************************

for i in range(0,nbin):
    print(zlow[i],zup[i],Mlow[i],Mup[i],S850[i],Error[i])
    
# Find the mean redshifts and stellar masses in each bin. Note that this is
# ****************************************************************************

Mass_mean = np.empty(nbin)
z_mean = np.empty(nbin)

for i in range(0,nbin):
    Mass_mean[i] = (Mlow[i] + Mup[i]) / 2.0
    z_mean[i] = (zlow[i] + zup[i]) / 2.0


# Now arrange the results for each bin for plotting
# ****************************************************************************

# First redshift bin

x21 = np.empty(5)
y21 = np.empty(5)
error21 = np.empty(5)
xline = np.empty(2)
yline = np.empty(2)

for i in range(0,5):
    x21[i] = Mass_mean[i]
    y21[i] = S850[i]
    error21[i] = Error[i]
        
# second redshift bin

x22 = np.empty(5)
y22 = np.empty(5)
error22 = np.empty(5)

for i in range(0,5):
    x22[i] = Mass_mean[5+i]
    y22[i] = S850[5+i]
    error22[i] = Error[5+i]
    
# third redshift bin

x23 = np.empty(5)
y23 = np.empty(5)
error23 = np.empty(5)


for i in range(0,5):
    x23[i] = Mass_mean[10+i]
    y23[i] = S850[10+i]
    error23[i] = Error[10+i]

# fourth redshift bin

x24 = np.empty(5)
y24 = np.empty(5)
error24 = np.empty(5)

print("4th bin")


for i in range(0,5):
    x24[i] = Mass_mean[15+i]
    y24[i] = S850[15+i]
    error24[i] = Error[15+i]

# fifth redshift bin

x25 = np.empty(5)
y25 = np.empty(5)
error25 = np.empty(5)



for i in range(0,5):
    x25[i] = Mass_mean[20+i]
    y25[i] = S850[20+i]
    error25[i] = Error[20+i]

# sixth redshift bin

x26 = np.empty(5)
y26 = np.empty(5)
error26 = np.empty(5)


for i in range(0,5):
    x26[i] = Mass_mean[25+i]
    y26[i] = S850[25+i]
    error26[i] = Error[25+i]

# seventh redshift bin

x27 = np.empty(5)
y27 = np.empty(5)
error27 = np.empty(5)

for i in range(0,5):
    x27[i] = Mass_mean[30+i]
    y27[i] = S850[30+i]
    error27[i] = Error[30+i]
    
# eigth redshift bin

x28 = np.empty(5)
y28 = np.empty(5)
error28 = np.empty(5)

for i in range(0,5):
    x28[i] = Mass_mean[35+i]
    y28[i] = S850[35+i]
    error28[i] = Error[35+i]
       
# ninth redshift bin

x29 = np.empty(5)
y29 = np.empty(5)
error29 = np.empty(5)


for i in range(0,5):
    x29[i] = Mass_mean[40+i]
    y29[i] = S850[40+i]
    error29[i] = Error[40+i]
    
       
# tenth redshift bin

x30 = np.empty(5)
y30 = np.empty(5)
error30 = np.empty(5)


for i in range(0,5):
    x30[i] = Mass_mean[45+i]
    y30[i] = S850[45+i]   
    error30[i] = Error[45+i]
    
  # 11th redshift bin

x31 = np.empty(5)
y31 = np.empty(5)
error31 = np.empty(5)

for i in range(0,5):
    x31[i] = Mass_mean[50+i]
    y31[i] = S850[50+i]   
    error31[i] = Error[50+i]

# 12th redshift bin

x32 = np.empty(5)
y32 = np.empty(5)
error32 = np.empty(5)


for i in range(0,5):
    x32[i] = Mass_mean[55+i]
    y32[i] = S850[55+i]   
    error32[i] = Error[55+i]

# 13th redshift bin

print('13th bin')

x33 = np.empty(5)
y33 = np.empty(5)
error33 = np.empty(5)

for i in range(0,5):
    x33[i] = Mass_mean[60+i]
    y33[i] = S850[60+i]   
    error33[i] = Error[60+i]
    print(x33[i],y33[i],error33[i])

# 14th redshift bin

x34 = np.empty(5)
y34 = np.empty(5)
error34 = np.empty(5)

print('14th bin')

for i in range(0,5):
    x34[i] = Mass_mean[65+i]
    y34[i] = S850[65+i]   
    error34[i] = Error[65+i]
    print(x34[i],y34[i],error34[i])

y21_uplims = np.array([0,0,0,0,0],dtype=bool)
y22_uplims = np.array([0,0,0,0,0],dtype=bool)
y23_uplims = np.array([0,0,0,0,0],dtype=bool)
y24_uplims = np.array([0,0,0,0,0],dtype=bool)
y25_uplims = np.array([0,0,0,0,0],dtype=bool)
y26_uplims = np.array([0,0,0,0,0],dtype=bool)
y27_uplims = np.array([0,0,0,0,0],dtype=bool)
y28_uplims = np.array([0,0,0,0,0],dtype=bool)
y29_uplims = np.array([0,0,0,0,0],dtype=bool)
y30_uplims = np.array([0,0,0,0,0],dtype=bool)
y31_uplims = np.array([0,0,0,0,0],dtype=bool)
y32_uplims = np.array([0,0,0,0,0],dtype=bool)
y33_uplims = np.array([0,0,0,0,0],dtype=bool)
y34_uplims = np.array([0,0,0,0,0],dtype=bool)

for i in range(0,5):
    
    if y21[i] < 0.0:
        y21[i] = 3.0*error21[i]
        y21_uplims[i] = 1
        
    if y22[i] < 0.0:
        y22[i] = 3.0*error22[i]
        y22_uplims[i] = 1

    if y23[i] < 0.0:
        y23[i] = 3.0*error23[i]
        y23_uplims[i] = 1
        
    if y24[i] < 0.0:
        y24[i] = 3.0*error24[i]
        y24_uplims[i] = 1
        
    if y25[i] < 0.0:
        y25[i] = 3.0*error25[i]
        y25_uplims[i] = 1
        
    if y26[i] < 0.0:
        y26[i] = 3.0*error26[i]
        y26_uplims[i] = 1
        
    if y27[i] < 0.0:
        y27[i] = 3.0*error27[i]
        y27_uplims[i] = 1
        
    if y28[i] < 0.0:
        y28[i] = 3.0*error28[i]
        y28_uplims[i] = 1

    if y29[i] < 0.0:
        y29[i] = 3.0*error29[i]
        y29_uplims[i] = 1

    if y30[i] < 0.0:
        y30[i] = 3.0*error30[i]
        y30_uplims[i] = 1
        
    if y31[i] < 0.0:
        y31[i] = 3.0*error31[i]
        y31_uplims[i] = 1

    if y32[i] < 0.0:
        y32[i] = 3.0*error32[i]
        y32_uplims[i] = 1

    if y33[i] < 0.0:
        y33[i] = 3.0*error33[i]
        y33_uplims[i] = 1
        
    if y34[i] < 0.0:
        y34[i] = 3.0*error34[i]
        y34_uplims[i] = 1

        
fig = plt.figure(figsize=(10.0,10.0))

pos_ann = ([9.5,4.0])

xl = np.empty(2)
yl = np.empty(2)
xl[0] = 9.0
xl[1] = 11.5
yl[0] = 0.1
yl[1] = 0.1

# First panel

f1 = plt.axes([0.15,0.1,0.2,0.2])
f1.set_xlim(9.0,11.5)
f1.set_ylim(0.01,9.0)
f1.set_yscale('log')
f1.tick_params(axis='both',which='both',labelsize=18)

f1.plot(x1,y1,'ms')
f1.errorbar(x1,y1,yerr= error1,fmt='m',ls='none',uplims=y1_uplims)

x21 = x21 + 0.1
f1.plot(x21,y21,'co')
f1.errorbar(x21,y21,yerr= error21,fmt='c',ls='none',uplims=y21_uplims)

f1.plot(xl,yl,'k--')

f1.annotate('0.2<z<0.5',pos_ann,fontweight='bold',size=15)

# second panel

f1 = plt.axes([0.35,0.1,0.2,0.2])
f1.set_xlim(9.0,11.5)
f1.set_ylim(0.01,9.0)
f1.set_yscale('log')
f1.tick_params(axis='x',labelsize=18)
f1.set_yticks([])

f1.plot(x2,y2,'ms')
f1.errorbar(x2,y2,yerr= error2,fmt='m',ls='none',uplims=y2_uplims)

x22 = x22 + 0.1
f1.plot(x22,y22,'co')
f1.errorbar(x22,y22,yerr= error22,fmt='c',ls='none',uplims=y22_uplims)

f1.plot(xl,yl,'k--')

f1.annotate('0.5<z<0.8',pos_ann,fontweight='bold',size=15)

# Third panel

f1 = plt.axes([0.55,0.1,0.2,0.2])
f1.set_xlim(9.0,11.5)
f1.set_ylim(0.01,9.0)
f1.set_yscale('log')
f1.set_xlabel('$log(M_*)$',size=25)
f1.tick_params(axis='x',labelsize=18)
f1.set_yticks([])

f1.plot(x3,y3,'ms')
f1.errorbar(x3,y3,yerr= error3,fmt='m',ls='none',uplims=y3_uplims)

x23 = x23 + 0.1
f1.plot(x23,y23,'co')
f1.errorbar(x23,y23,yerr= error23,fmt='c',ls='none',uplims=y23_uplims)

f1.plot(xl,yl,'k--')

f1.annotate('0.8<z<1.1',pos_ann,fontweight='bold',size=15)


# Fourth panel

f1 = plt.axes([0.75,0.1,0.2,0.2])
f1.set_xlim(9.0,11.5)
f1.set_ylim(0.01,9.0)
f1.set_yscale('log')
f1.tick_params(axis='x',labelsize=18)
f1.set_yticks([])

f1.annotate('1.1<z<1.5',pos_ann,fontweight='bold',size=15)

f1.plot(x4,y4,'ms')
f1.errorbar(x4,y4,yerr= error4,fmt='m',ls='none',uplims=y4_uplims)

x24 = x24 + 0.1
f1.plot(x24,y24,'co')
f1.errorbar(x24,y24,yerr= error24,fmt='c',ls='none',uplims=y24_uplims)

f1.plot(xl,yl,'k--')

# Fifth panel

f1 = plt.axes([0.15,0.3,0.2,0.2])
f1.set_xlim(9.0,11.5)
f1.set_ylim(0.01,9.0)
f1.set_yscale('log')
f1.set_ylabel('$Flux/mJy$',size=25)
f1.tick_params(axis='y',labelsize=18)
f1.set_xticks([])

f1.annotate('1.5<z<2.0',pos_ann,fontweight='bold',size=15)

f1.plot(x5,y5,'ms')
f1.errorbar(x5,y5,yerr= error5,fmt='m',ls='none',uplims=y5_uplims)

x25 = x25 + 0.1
f1.plot(x25,y25,'co')
f1.errorbar(x25,y25,yerr= error25,fmt='c',ls='none',uplims=y25_uplims)

f1.plot(xl,yl,'k--')

# Sixth panel

f1 = plt.axes([0.35,0.3,0.2,0.2])
f1.set_xlim(9.0,11.5)
f1.set_ylim(0.01,9.0)
f1.set_yscale('log')
f1.set_xticks([])
f1.set_yticks([])

f1.annotate('2.0<z<2.5',pos_ann,fontweight='bold',size=15)

f1.plot(x6,y6,'ms')
f1.errorbar(x6,y6,yerr= error6,fmt='m',ls='none',uplims=y6_uplims)

x26 = x26 + 0.1
f1.plot(x26,y26,'co')
f1.errorbar(x26,y26,yerr= error26,fmt='c',ls='none',uplims=y26_uplims)

f1.plot(xl,yl,'k--')

# Seventh panel

f1 = plt.axes([0.55,0.3,0.2,0.2])
f1.set_xlim(9.0,11.5)
f1.set_ylim(0.01,9.0)
f1.set_yscale('log')
f1.set_xticks([])
f1.set_yticks([])

f1.annotate('2.5<z<3.0',pos_ann,fontweight='bold',size=15)

f1.plot(x7,y7,'ms')
f1.errorbar(x7,y7,yerr= error7,fmt='m',ls='none',uplims=y7_uplims)

x27 = x27 + 0.1
f1.plot(x27,y27,'co')
f1.errorbar(x27,y27,yerr= error27,fmt='c',ls='none',uplims=y27_uplims)

f1.plot(xl,yl,'k--')

# Eighth panel

f1 = plt.axes([0.75,0.3,0.2,0.2])
f1.set_xlim(9.0,11.5)
f1.set_ylim(0.01,9.0)
f1.set_yscale('log')
f1.set_xticks([])
f1.set_yticks([])

f1.annotate('3.0<z<3.5',pos_ann,fontweight='bold',size=15)

f1.plot(x8,y8,'ms')
f1.errorbar(x8,y8,yerr= error8,fmt='m',ls='none',uplims=y8_uplims)

x28 = x28 + 0.1
f1.plot(x28,y28,'co')
f1.errorbar(x28,y28,yerr= error28,fmt='c',ls='none',uplims=y28_uplims)

f1.plot(xl,yl,'k--')


# Ninth panel

f1 = plt.axes([0.15,0.5,0.2,0.2])
f1.set_xlim(9.0,11.5)
f1.set_ylim(0.01,9.0)
f1.set_yscale('log')
f1.set_xticks([])
f1.tick_params(axis='y',labelsize=18)

f1.annotate('3.5<z<4.5',pos_ann,fontweight='bold',size=15)

f1.plot(x9,y9,'ms')
f1.errorbar(x9,y9,yerr= error9,fmt='m',ls='none',uplims=y9_uplims)

x29 = x29 + 0.1
f1.plot(x29,y29,'co')
f1.errorbar(x29,y29,yerr= error29,fmt='c',ls='none',uplims=y29_uplims)

f1.plot(xl,yl,'k--')

# Tenth panel

f1 = plt.axes([0.35,0.5,0.2,0.2])
f1.set_xlim(9.0,11.5)
f1.set_ylim(0.01,9.0)
f1.set_yscale('log')
f1.set_xticks([])
f1.set_yticks([])

f1.annotate('4.5<z<5.5',pos_ann,fontweight='bold',size=15)

f1.plot(x10,y10,'ms')
f1.errorbar(x10,y10,yerr= error10,fmt='m',ls='none',uplims=y10_uplims)

x30 = x30 + 0.1
f1.plot(x30,y30,'co')
f1.errorbar(x30,y30,yerr= error30,fmt='c',ls='none',uplims=y30_uplims)

f1.plot(xl,yl,'k--')

# 11th panel

f1 = plt.axes([0.55,0.5,0.2,0.2])
f1.set_xlim(9.0,11.5)
f1.set_ylim(0.01,9.0)
f1.set_yscale('log')
f1.set_xticks([])
f1.set_yticks([])
f1.tick_params(axis='y',labelsize=12)

f1.annotate('5.5<z<6.5',pos_ann,fontweight='bold',size=15)

f1.plot(x11,y11,'ms')
f1.errorbar(x11,y11,yerr= error11,fmt='m',ls='none',uplims=y11_uplims)

x31 = x31 + 0.1
f1.plot(x31,y31,'co')
f1.errorbar(x31,y31,yerr= error31,fmt='c',ls='none',uplims=y31_uplims)

f1.plot(xl,yl,'k--')

# 12th panel

f1 = plt.axes([0.75,0.5,0.2,0.2])
f1.set_xlim(9.0,11.5)
f1.set_ylim(0.01,9.0)
f1.set_yscale('log')
f1.set_xticks([])
f1.set_yticks([])
f1.tick_params(axis='y',labelsize=12)

f1.annotate('6.5<z<7.5',pos_ann,fontweight='bold',size=15)

f1.plot(x12,y12,'ms')
f1.errorbar(x12,y12,yerr= error12,fmt='m',ls='none',uplims=y12_uplims)

x32 = x32 + 0.1
f1.plot(x32,y32,'co')
f1.errorbar(x32,y32,yerr= error32,fmt='c',ls='none',uplims=y32_uplims)

f1.plot(xl,yl,'k--')


# 13th panel

f1 = plt.axes([0.15,0.7,0.2,0.2])
f1.set_xlim(9.0,11.5)
f1.set_ylim(0.01,9.0)
f1.set_yscale('log')
f1.set_xticks([])
f1.tick_params(axis='y',labelsize=18)

f1.annotate('7.5<z<8.5',pos_ann,fontweight='bold',size=15)

f1.plot(x13,y13,'ms')
f1.errorbar(x13,y13,yerr= error13,fmt='m',ls='none',uplims=y13_uplims)

x33 = x33 + 0.1
f1.plot(x33,y33,'co')
f1.errorbar(x33,y33,yerr= error33,fmt='c',ls='none',uplims=y33_uplims)

f1.plot(xl,yl,'k--')


# 14th panel

f1 = plt.axes([0.35,0.7,0.2,0.2])
f1.set_xlim(9.0,11.5)
f1.set_ylim(0.01,9.0)
f1.set_yscale('log')
f1.set_xticks([])
f1.set_yticks([])
f1.tick_params(axis='y',labelsize=12)

f1.annotate('8.5<z<12.0',pos_ann,fontweight='bold',size=15)

f1.plot(x14,y14,'ms')
f1.errorbar(x14,y14,yerr= error14,fmt='m',ls='none',uplims=y14_uplims)

x34 = x34 + 0.1
f1.plot(x34,y34,'co')
f1.errorbar(x34,y34,yerr= error34,fmt='c',ls='none',uplims=y34_uplims)

f1.plot(xl,yl,'k--')

fig.savefig('Figure2.pdf')
plt.show()

with open('TableA1.dat', 'a') as file:
    for i in range(0,nbin):
        file.write(f"{zlow[i]} ")
        file.write(f"{zup[i]} ")
        file.write(f"{Mlow[i]} ")
        file.write(f"{Mup[i]} ")
        file.write(f"{S850[i]} ")
        file.write(f"{Error[i]}" + "\n")



            
            
    
    

    
        
    
    



    





