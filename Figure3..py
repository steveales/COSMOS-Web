# This program plots flux versus redshift in a separate panel for
# each redshift slice. This version of the program uses the results
# of the basic method and plots the results for star-forming and
# quiescent galaxies in the same diagram. It does not yet include the
# errors from the bootstrap of the JWST catalogue.

# This version of the program plots measurements detected at <3sigma as
# upper limits

# Last edited: 28th November 2025

import numpy as np
import matplotlib.pyplot as plt

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

        
# Read in stacking data catalogue
# ****************************************************************************

file='SIMSTACK_results_collated.dat'

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

Error1 = np.empty(nbin)
Error1 = Error

# Read in artificial data
# *****************************************************************************

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


# Convert fluxes to mJy
# ****************************************************************************

for i in range (0,nbin):
    S850[i] = S850[i] * 1000.0
    Error[i] = Error[i] * 1000.0
    art_noise[i] = art_noise[i] * 1000.0
    
for i in range(0,nbin):
    print(zlow[i],zup[i],S850[i],Error1[i],Error[i],art_noise[i])

Error1 = Error
    
# Add the noise in quadrature

for i in range(0,nbin):
    tp1 = Error[i] * Error[i] + art_noise[i]*art_noise[i]
    Error[i] = np.sqrt(tp1)

    
# Print out the results for all the bins
# ****************************************************************************

#for i in range(0,nbin):
#    print(zlow[i],zup[i],Mlow[i],Mup[i],S850[i],Error[i],Error1[i],art_noise[i])
    
# Find the mean redshifts and stellar masses in each bin. Note that this is
# ****************************************************************************

Mass_mean = np.empty(nbin)
z_mean = np.empty(nbin)

for i in range(0,nbin):
    Mass_mean[i] = (Mlow[i] + Mup[i]) / 2.0
    z_mean[i] = (zlow[i] + zup[i]) / 2.0

# Change -99s
# *****************************************************************************

for i in range(0,nbin):
    if S850[i] < -98.0:
        Error[i] = 0.0

# Now arrange the results for each bin for plotting
# ****************************************************************************

# First redshift bin

x1 = np.empty(5)
y1 = np.empty(5)
error1 = np.empty(5)
xline = np.empty(2)
yline = np.empty(2)


for i in range(0,5):
    x1[i] = Mass_mean[i]
    y1[i] = S850[i]
    error1[i] = Error[i]

x1_q = np.empty(5)
y1_q = np.empty(5)
error1_q = np.empty(5)

for i in range(0,5):
    x1_q[i] = Mass_mean[70+i]
    y1_q[i] = S850[70+i]
    error1_q[i] = Error[70+i]
    
  
# second redshift bin

x2 = np.empty(5)
y2 = np.empty(5)
error2 = np.empty(5)


for i in range(0,5):
    x2[i] = Mass_mean[5+i]
    y2[i] = S850[5+i]
    error2[i] = Error[5+i]
    
x2_q = np.empty(5)
y2_q = np.empty(5)
error2_q = np.empty(5)

for i in range(0,5):
    x2_q[i] = Mass_mean[75+i]
    y2_q[i] = S850[75+i]
    error2_q[i] = Error[75+i]
    
# third redshift bin

x3 = np.empty(5)
y3 = np.empty(5)
error3 = np.empty(5)



for i in range(0,5):
    x3[i] = Mass_mean[10+i]
    y3[i] = S850[10+i]
    error3[i] = Error[10+i]
    
x3_q = np.empty(5)
y3_q = np.empty(5)
error3_q = np.empty(5)

for i in range(0,5):
    x3_q[i] = Mass_mean[80+i]
    y3_q[i] = S850[80+i]
    error3_q[i] = Error[80+i]
    
# fourth redshift bin

x4 = np.empty(5)
y4 = np.empty(5)
error4 = np.empty(5)



for i in range(0,5):
    x4[i] = Mass_mean[15+i]
    y4[i] = S850[15+i]
    error4[i] = Error[15+i]
    
x4_q = np.empty(5)
y4_q = np.empty(5)
error4_q = np.empty(5)

for i in range(0,5):
    x4_q[i] = Mass_mean[85+i]
    y4_q[i] = S850[85+i]
    error4_q[i] = Error[85+i]
    
# fifth redshift bin

x5 = np.empty(5)
y5 = np.empty(5)
error5 = np.empty(5)



for i in range(0,5):
    x5[i] = Mass_mean[20+i]
    y5[i] = S850[20+i]
    error5[i] = Error[20+i]
    
x5_q = np.empty(5)
y5_q = np.empty(5)
error5_q = np.empty(5)

for i in range(0,5):
    x5_q[i] = Mass_mean[90+i]
    y5_q[i] = S850[90+i]
    error5_q[i] = Error[90+i]

# sixth redshift bin

x6 = np.empty(5)
y6 = np.empty(5)
error6 = np.empty(5)



for i in range(0,5):
    x6[i] = Mass_mean[25+i]
    y6[i] = S850[25+i]
    error6[i] = Error[25+i]
    
x6_q = np.empty(5)
y6_q = np.empty(5)
error6_q = np.empty(5)

for i in range(0,5):
    x6_q[i] = Mass_mean[95+i]
    y6_q[i] = S850[95+i]
    error6_q[i] = Error[95+i]

# seventh redshift bin

x7 = np.empty(5)
y7 = np.empty(5)
error7 = np.empty(5)


for i in range(0,5):
    x7[i] = Mass_mean[30+i]
    y7[i] = S850[30+i]
    error7[i] = Error[30+i]
    
x7_q = np.empty(5)
y7_q = np.empty(5)
error7_q = np.empty(5)

for i in range(0,5):
    x7_q[i] = Mass_mean[100+i]
    y7_q[i] = S850[100+i]
    error7_q[i] = Error[100+i]

# eigth redshift bin

x8 = np.empty(5)
y8 = np.empty(5)
error8 = np.empty(5)


for i in range(0,5):
    x8[i] = Mass_mean[35+i]
    y8[i] = S850[35+i]
    error8[i] = Error[35+i]
    
x8_q = np.empty(5)
y8_q = np.empty(5)
error8_q = np.empty(5)

for i in range(0,5):
    x8_q[i] = Mass_mean[105+i]
    y8_q[i] = S850[105+i]
    error8_q[i] = Error[105+i]
 
# ninth redshift bin

x9 = np.empty(5)
y9 = np.empty(5)
error9 = np.empty(5)


for i in range(0,5):
    x9[i] = Mass_mean[40+i]
    y9[i] = S850[40+i]
    error9[i] = Error[40+i]
    
x9_q = np.empty(5)
y9_q = np.empty(5)
error9_q = np.empty(5)

for i in range(0,5):
    x9_q[i] = Mass_mean[110+i]
    y9_q[i] = S850[110+i]
    error9_q[i] = Error[110+i]
 
# tenth redshift bin

x10 = np.empty(5)
y10 = np.empty(5)
error10 = np.empty(5)


for i in range(0,5):
    x10[i] = Mass_mean[45+i]
    y10[i] = S850[45+i]   
    error10[i] = Error[45+i]
    
x10_q = np.empty(5)
y10_q = np.empty(5)
error10_q = np.empty(5)

for i in range(0,5):
    x10_q[i] = Mass_mean[115+i]
    y10_q[i] = S850[115+i]
    error10_q[i] = Error[115+i]
    
  # 11th redshift bin

x11 = np.empty(5)
y11 = np.empty(5)
error11 = np.empty(5)


for i in range(0,5):
    x11[i] = Mass_mean[50+i]
    y11[i] = S850[50+i]   
    error11[i] = Error[50+i]

x11_q = np.empty(5)
y11_q = np.empty(5)
error11_q = np.empty(5)

for i in range(0,5):
    x11_q[i] = Mass_mean[120+i]
    y11_q[i] = S850[120+i]
    error11_q[i] = Error[120+i]

# 12th redshift bin

x12 = np.empty(5)
y12 = np.empty(5)
error12 = np.empty(5)


for i in range(0,5):
    x12[i] = Mass_mean[55+i]
    y12[i] = S850[55+i]   
    error12[i] = Error[55+i]

x12_q = np.empty(5)
y12_q = np.empty(5)
error12_q = np.empty(5)

for i in range(0,5):
    x12_q[i] = Mass_mean[125+i]
    y12_q[i] = S850[125+i]
    error12_q[i] = Error[125+i]

# 13th redshift bin

x13 = np.empty(5)
y13 = np.empty(5)
error13 = np.empty(5)


for i in range(0,5):
    x13[i] = Mass_mean[60+i]
    y13[i] = S850[60+i]   
    error13[i] = Error[60+i]

    
x13_q = np.empty(5)
y13_q = np.empty(5)
error13_q = np.empty(5)

for i in range(0,5):
    x13_q[i] = Mass_mean[130+i]
    y13_q[i] = S850[130+i]
    error13_q[i] = Error[130+i]

# 14th redshift bin

x14 = np.empty(5)
y14 = np.empty(5)
error14 = np.empty(5)


for i in range(0,5):
    x14[i] = Mass_mean[65+i]
    y14[i] = S850[65+i]   
    error14[i] = Error[65+i]
 
    
x14_q = np.empty(5)
y14_q = np.empty(5)
error14_q = np.empty(5)

for i in range(0,5):
    x14_q[i] = Mass_mean[135+i]
    y14_q[i] = S850[135+i]
    error14_q[i] = Error[135+i]


# Replace negative measurements <3sigma with upper limits
# ****************************************************************************

# Bin 1


y1_uplims = np.array([0,0,0,0,0],dtype=bool)
y1_q_uplims = np.array([0,0,0,0,0],dtype=bool)


for i in range(0,5):
    
    if y1[i]/error1[i] < 0.0:
        y1[i] = 3.0*error1[i]
        y1_uplims[i] = 1

    if y1_q[i]/error1_q[i] < 0.0:
        y1_q[i] = 3.0*error1_q[i]
        y1_q_uplims[i] = 1

# Bin 2

y2_uplims = np.array([0,0,0,0,0],dtype=bool)
y2_q_uplims = np.array([0,0,0,0,0],dtype=bool)


for i in range(0,5):
    
    if y2[i]/error2[i] < 0.0:
        y2[i] = 3.0*error2[i]
        y2_uplims[i] = 1

    if y2_q[i]/error2_q[i] < 0.0:
        y2_q[i] = 3.0*error2_q[i]
        y2_q_uplims[i] = 1

# Bin 3

y3_uplims = np.array([0,0,0,0,0],dtype=bool)
y3_q_uplims = np.array([0,0,0,0,0],dtype=bool)


for i in range(0,5):
    
    if y3[i]/error3[i] < 0.0:
        y3[i] = 3.0*error3[i]
        y3_uplims[i] = 1

    if y3_q[i]/error3_q[i] < 0.0:
        y3_q[i] = 3.0*error3_q[i]
        y3_q_uplims[i] = 1
        
# Bin 4

y4_uplims = np.array([0,0,0,0,0],dtype=bool)
y4_q_uplims = np.array([0,0,0,0,0],dtype=bool)


for i in range(0,5):
    
    if y4[i]/error4[i] < 0.0:
        y4[i] = 3.0*error4[i]
        y4_uplims[i] = 1

    if y4_q[i]/error4_q[i] < 0.0:
        y4_q[i] = 3.0*error4_q[i]
        y4_q_uplims[i] = 1

# Bin 5

y5_uplims = np.array([0,0,0,0,0],dtype=bool)
y5_q_uplims = np.array([0,0,0,0,0],dtype=bool)


for i in range(0,5):
    
    if y5[i]/error5[i] < 0.0:
        y5[i] = 3.0*error5[i]
        y5_uplims[i] = 1

    if y5_q[i]/error5_q[i] < 0.0:
        y5_q[i] = 3.0*error5_q[i]
        y5_q_uplims[i] = 1

# Bin 6

y6_uplims = np.array([0,0,0,0,0],dtype=bool)
y6_q_uplims = np.array([0,0,0,0,0],dtype=bool)


for i in range(0,5):
    
    if y6[i]/error6[i] < 0.0:
        y6[i] = 3.0*error6[i]
        y6_uplims[i] = 1

    if y6_q[i]/error6_q[i] < 0.0:
        y6_q[i] = 3.0*error6_q[i]
        y6_q_uplims[i] = 1

# Bin 7

y7_uplims = np.array([0,0,0,0,0],dtype=bool)
y7_q_uplims = np.array([0,0,0,0,0],dtype=bool)


for i in range(0,5):
    
    if y7[i]/error7[i] < 0.0:
        y7[i] = 3.0*error7[i]
        y7_uplims[i] = 1

    if y7_q[i]/error7_q[i] < 0.0:
        y7_q[i] = 3.0*error7_q[i]
        y7_q_uplims[i] = 1

# Bin 8

y8_uplims = np.array([0,0,0,0,0],dtype=bool)
y8_q_uplims = np.array([0,0,0,0,0],dtype=bool)


for i in range(0,5):
    
    if y8[i]/error8[i] < 0.0:
        y8[i] = 3.0*error8[i]
        y8_uplims[i] = 1

    if y8_q[i]/error8_q[i] < 0.0:
        y8_q[i] = 3.0*error8_q[i]
        y8_q_uplims[i] = 1

# Bin 9

y9_uplims = np.array([0,0,0,0,0],dtype=bool)
y9_q_uplims = np.array([0,0,0,0,0],dtype=bool)


for i in range(0,5):
    
    if y9[i]/error9[i] < 0.0:
        y9[i] = 3.0*error9[i]
        y9_uplims[i] = 1

    if y9_q[i]/error9_q[i] < 0.0:
        y9_q[i] = 3.0*error9_q[i]
        y9_q_uplims[i] = 1

# Bin 10

y10_uplims = np.array([0,0,0,0,0],dtype=bool)
y10_q_uplims = np.array([0,0,0,0,0],dtype=bool)


for i in range(0,5):
    
    if y10[i]/error10[i] < 0.0:
        y10[i] = 3.0*error10[i]
        y10_uplims[i] = 1

    if y10_q[i]/error10_q[i] < 0.0:
        y10_q[i] = 3.0*error10_q[i]
        y10_q_uplims[i] = 1


# Bin 11

y11_uplims = np.array([0,0,0,0,0],dtype=bool)
y11_q_uplims = np.array([0,0,0,0,0],dtype=bool)


for i in range(0,5):
    
    if y11[i]/error11[i] < 0.0:
        y11[i] = 3.0*error11[i]
        y11_uplims[i] = 1

    if y11_q[i]/error11_q[i] < 0.0:
        y11_q[i] = 3.0*error11_q[i]
        y11_q_uplims[i] = 1
        
# Bin 12

y12_uplims = np.array([0,0,0,0,0],dtype=bool)
y12_q_uplims = np.array([0,0,0,0,0],dtype=bool)


for i in range(0,5):
    
    if y12[i]/error12[i] < 0.0:
        y12[i] = 3.0*error12[i]
        y12_uplims[i] = 1

    if y12_q[i]/error12_q[i] < 0.0:
        y12_q[i] = 3.0*error12_q[i]
        y12_q_uplims[i] = 1
        

# Bin 13

y13_uplims = np.array([0,0,0,0,0],dtype=bool)
y13_q_uplims = np.array([0,0,0,0,0],dtype=bool)


for i in range(0,5):
    
    if y13[i]/error13[i] < 0.0:
        y13[i] = 3.0*error13[i]
        y13_uplims[i] = 1

    if y13_q[i]/error13_q[i] < 0.0:
        y13_q[i] = 3.0*error13_q[i]
        y13_q_uplims[i] = 1

# Bin 14

y14_uplims = np.array([0,0,0,0,0],dtype=bool)
y14_q_uplims = np.array([0,0,0,0,0],dtype=bool)


for i in range(0,5):
    
    if y14[i]/error14[i] < 0.0:
        y14[i] = 3.0*error14[i]
        y14_uplims[i] = 1

    if y14_q[i]/error14_q[i] < 0.0:
        y14_q[i] = 3.0*error14_q[i]
        y14_q_uplims[i] = 1
        

# Plot out results
# ****************************************************************************
    
fig = plt.figure(figsize=(10.0,10.0))

pos_ann = ([9.5,4.0])

xl = np.empty(2)
yl = np.empty(2)
xl[0] = 9.0
xl[1] = 11.5
yl[0] = 0.1
yl[1]= 0.1

# First panel

f1 = plt.axes([0.15,0.1,0.2,0.2])
f1.set_xlim(9.0,11.5)
f1.set_ylim(0.01,9.0)
f1.set_yscale('log')
f1.tick_params(axis='both',which='both',labelsize=18)

f1.plot(x1,y1,'ms')
f1.errorbar(x1,y1,yerr= error1,fmt='m',ls='none',uplims=y1_uplims)

f1.annotate('0.2<z<0.5',pos_ann,fontweight='bold',size=15)

x1_q = x1_q + 0.1
f1.plot(x1_q,y1_q,'co')
f1.errorbar(x1_q,y1_q,yerr= error1_q,fmt='c',ls='none',uplims=y1_q_uplims)


f1.plot(xl,yl,'k--')


# second panel

f1 = plt.axes([0.35,0.1,0.2,0.2])
f1.set_xlim(9.0,11.5)
f1.set_ylim(0.01,9.0)
f1.set_yscale('log')
f1.tick_params(axis='x',labelsize=18)
f1.set_yticks([])

f1.plot(x2,y2,'ms')
f1.errorbar(x2,y2,yerr= error2,fmt='m',ls='none',uplims=y2_uplims)

f1.annotate('0.5<z<0.8',pos_ann,fontweight='bold',size=15)

x2_q = x2_q + 0.1
f1.plot(x2_q,y2_q,'co')
f1.errorbar(x2_q,y2_q,yerr= error2_q,fmt='c',ls='none',uplims=y2_q_uplims)

f1.plot(xl,yl,'k--')


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

f1.annotate('0.8<z<1.1',pos_ann,fontweight='bold',size=15)

x3_q = x3_q + 0.1
f1.plot(x3_q,y3_q,'co')
f1.errorbar(x3_q,y3_q,yerr= error3_q,fmt='c',ls='none',uplims=y3_q_uplims)

f1.plot(xl,yl,'k--')

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

x4_q = x4_q + 0.1
f1.plot(x4_q,y4_q,'co')
f1.errorbar(x4_q,y4_q,yerr= error4_q,fmt='c',ls='none',uplims=y4_q_uplims)


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

x5_q = x5_q + 0.1
f1.plot(x5_q,y5_q,'co')
f1.errorbar(x5_q,y5_q,yerr= error5_q,fmt='c',ls='none',uplims=y5_q_uplims)

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

x6_q = x6_q + 0.1
f1.plot(x6_q,y6_q,'co')
f1.errorbar(x6_q,y6_q,yerr= error6_q,fmt='c',ls='none',uplims=y6_q_uplims)

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

x7_q = x7_q + 0.1
f1.plot(x7_q,y7_q,'co')
f1.errorbar(x7_q,y7_q,yerr= error7_q,fmt='c',ls='none',uplims=y7_q_uplims)

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

x8_q = x8_q + 0.1
f1.plot(x8_q,y8_q,'co')
f1.errorbar(x8_q,y8_q,yerr= error8_q,fmt='c',ls='none',uplims=y8_q_uplims)


f1.plot(xl,yl,'k--')

# Ninth panel

pos_ann = ([9.5,4.0e-2])

f1 = plt.axes([0.15,0.5,0.2,0.2])
f1.set_xlim(9.0,11.5)
f1.set_ylim(0.01,9.0)
f1.set_yscale('log')
f1.set_xticks([])
f1.tick_params(axis='y',labelsize=18)

f1.annotate('3.5<z<4.5',pos_ann,fontweight='bold',size=15)

f1.plot(x9,y9,'ms')
f1.errorbar(x9,y9,yerr= error9,fmt='m',ls='none',uplims=y9_uplims)

x9_q = x9_q + 0.1
f1.plot(x9_q,y9_q,'co')
f1.errorbar(x9_q,y9_q,yerr= error9_q,fmt='c',ls='none',uplims=y9_q_uplims)

f1.plot(xl,yl,'k--')


# Tenth panel

f1 = plt.axes([0.35,0.5,0.2,0.2])
f1.set_xlim(9.0,11.5)
f1.set_ylim(0.01,9.0)
f1.set_yscale('log')
f1.set_xticks([])
f1.set_yticks([])
f1.tick_params(axis='y',labelsize=12)

f1.annotate('4.5<z<5.5',pos_ann,fontweight='bold',size=15)

f1.plot(x10,y10,'ms')
f1.errorbar(x10,y10,yerr= error10,fmt='m',ls='none',uplims=y10_uplims)

x10_q = x10_q + 0.1
f1.plot(x10_q,y10_q,'co')
f1.errorbar(x10_q,y10_q,yerr= error10_q,fmt='c',ls='none',uplims=y10_q_uplims)


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

x11_q = x11_q + 0.1
f1.plot(x11_q,y11_q,'co')
f1.errorbar(x11_q,y11_q,yerr= error11_q,fmt='c',ls='none',uplims=y11_q_uplims)


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

x12_q = x12_q + 0.1
f1.plot(x12_q,y12_q,'co')
f1.errorbar(x12_q,y12_q,yerr= error12_q,fmt='c',ls='none',uplims=y12_q_uplims)

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

x13_q = x13_q + 0.1
f1.plot(x13_q,y13_q,'co')
f1.errorbar(x13_q,y13_q,yerr= error13_q,fmt='c',ls='none',uplims=y13_q_uplims)

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

x14_q = x14_q + 0.1
f1.plot(x14_q,y14_q,'co')
f1.errorbar(x14_q,y14_q,yerr= error14_q,fmt='c',ls='none',uplims=y14_q_uplims)

f1.plot(xl,yl,'k--')

fig.savefig('Figure3.pdf')
plt.show()





            
            
    
    

    
        
    
    



    





