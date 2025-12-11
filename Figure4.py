# This program plots flux versus redshift in a separate panel for
# each redshift slice. This version of the program uses the results
# of the basic method and plots the results for four morphological
# classes in the same diagram. It does not yet include the
# errors from the bootstrap of the JWST catalogue.

# This version only plots points if they are detected at >2 sigma. The
# upper limits are 3sigma upper limits. No values are plotted if there are
# less than 10 objects in a bin

# Last edited: 25th November 2025

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
Ngal = np.array([])
        
# Read in stacking data catalogue
# ****************************************************************************

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

# Convert fluxes to mJy
# ****************************************************************************

for i in range (0,nbin):
    S850[i] = S850[i] * 1000.0
    Error[i] = Error[i] * 1000.0
    art_noise[i] = art_noise[i] * 1000.0
    
# Add two noise estimates in quadrature
# ****************************************************************************

for i in range(0,nbin):
    tp1 = art_noise[i] * art_noise[i] + Error[i] * Error[i]
    Error[i] = np.sqrt(tp1)
    
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

# Change bins with less than 10 objects and ones with -99s so they
# won't be plotted
# *****************************************************************************

for i in range(0,nbin):
    
    if Ngal[i] < 10:
        S850[i] = 1.0e-9
        Error[i] = 1.0e-10
    
for i in range(0,nbin):
    if S850[i] < -98.0:
        S850[i] = 1.0e-9
        Error[i] = 1.0e-10

# Now arrange the results for each bin for plotting
# ****************************************************************************

# First redshift bin

x1_m1 = np.empty(5)
y1_m1 = np.empty(5)
error1_m1 = np.empty(5)

x1_m2 = np.empty(5)
y1_m2 = np.empty(5)
error1_m2 = np.empty(5)

x1_m3 = np.empty(5)
y1_m3 = np.empty(5)
error1_m3 = np.empty(5)

x1_m4 = np.empty(5)
y1_m4 = np.empty(5)
error1_m4 = np.empty(5)


print("bin1")

for i in range(0,5):
    x1_m1[i] = Mass_mean[i]
    y1_m1[i] = S850[i]
    error1_m1[i] = Error[i]

for i in range(0,5):
    x1_m2[i] = Mass_mean[70+i]
    y1_m2[i] = S850[70+i]
    error1_m2[i] = Error[70+i]
    
for i in range(0,5):
    x1_m3[i] = Mass_mean[140+i]
    y1_m3[i] = S850[140+i]
    error1_m3[i] = Error[140+i]
    
for i in range(0,5):
    x1_m4[i] = Mass_mean[210+i]
    y1_m4[i] = S850[210+i]
    error1_m4[i] = Error[210+i] 
    
# second redshift bin

x2_m1 = np.empty(5)
y2_m1 = np.empty(5)
error2_m1 = np.empty(5)

x2_m2 = np.empty(5)
y2_m2 = np.empty(5)
error2_m2 = np.empty(5)

x2_m3 = np.empty(5)
y2_m3 = np.empty(5)
error2_m3 = np.empty(5)

x2_m4 = np.empty(5)
y2_m4 = np.empty(5)
error2_m4 = np.empty(5)


print("bin2")

for i in range(0,5):
    x2_m1[i] = Mass_mean[5+i]
    y2_m1[i] = S850[5+i]
    error2_m1[i] = Error[5+i]

for i in range(0,5):
    x2_m2[i] = Mass_mean[75+i]
    y2_m2[i] = S850[75+i]
    error2_m2[i] = Error[75+i]
    
for i in range(0,5):
    x2_m3[i] = Mass_mean[145+i]
    y2_m3[i] = S850[145+i]
    error2_m3[i] = Error[145+i]
    
for i in range(0,5):
    x2_m4[i] = Mass_mean[215+i]
    y2_m4[i] = S850[215+i]
    error2_m4[i] = Error[215+i] 
    
    
# third redshift bin

x3_m1 = np.empty(5)
y3_m1 = np.empty(5)
error3_m1 = np.empty(5)

x3_m2 = np.empty(5)
y3_m2 = np.empty(5)
error3_m2 = np.empty(5)

x3_m3 = np.empty(5)
y3_m3 = np.empty(5)
error3_m3 = np.empty(5)

x3_m4 = np.empty(5)
y3_m4 = np.empty(5)
error3_m4 = np.empty(5)

print("bin3")

for i in range(0,5):
    x3_m1[i] = Mass_mean[10+i]
    y3_m1[i] = S850[10+i]
    error3_m1[i] = Error[10+i]

for i in range(0,5):
    x3_m2[i] = Mass_mean[80+i]
    y3_m2[i] = S850[80+i]
    error3_m2[i] = Error[80+i]
    
for i in range(0,5):
    x3_m3[i] = Mass_mean[150+i]
    y3_m3[i] = S850[150+i]
    error3_m3[i] = Error[150+i]
    
for i in range(0,5):
    x3_m4[i] = Mass_mean[220+i]
    y3_m4[i] = S850[220+i]
    error3_m4[i] = Error[220+i] 
    
    
# fourth redshift bin

x4_m1 = np.empty(5)
y4_m1 = np.empty(5)
error4_m1 = np.empty(5)

x4_m2 = np.empty(5)
y4_m2 = np.empty(5)
error4_m2 = np.empty(5)

x4_m3 = np.empty(5)
y4_m3 = np.empty(5)
error4_m3 = np.empty(5)

x4_m4 = np.empty(5)
y4_m4 = np.empty(5)
error4_m4 = np.empty(5)

print("bin 4")

for i in range(0,5):
    x4_m1[i] = Mass_mean[15+i]
    y4_m1[i] = S850[15+i]
    error4_m1[i] = Error[15+i]

for i in range(0,5):
    x4_m2[i] = Mass_mean[85+i]
    y4_m2[i] = S850[85+i]
    error4_m2[i] = Error[85+i]
    
for i in range(0,5):
    x4_m3[i] = Mass_mean[155+i]
    y4_m3[i] = S850[155+i]
    error4_m3[i] = Error[155+i]
    
for i in range(0,5):
    x4_m4[i] = Mass_mean[225+i]
    y4_m4[i] = S850[225+i]
    error4_m4[i] = Error[225+i] 
    
# fifth redshift bin

x5_m1 = np.empty(5)
y5_m1 = np.empty(5)
error5_m1 = np.empty(5)

x5_m2 = np.empty(5)
y5_m2 = np.empty(5)
error5_m2 = np.empty(5)

x5_m3 = np.empty(5)
y5_m3 = np.empty(5)
error5_m3 = np.empty(5)

x5_m4 = np.empty(5)
y5_m4 = np.empty(5)
error5_m4 = np.empty(5)

print("bin 5")

for i in range(0,5):
    x5_m1[i] = Mass_mean[20+i]
    y5_m1[i] = S850[20+i]
    error5_m1[i] = Error[20+i]

for i in range(0,5):
    x5_m2[i] = Mass_mean[90+i]
    y5_m2[i] = S850[90+i]
    error5_m2[i] = Error[90+i]
    
for i in range(0,5):
    x5_m3[i] = Mass_mean[160+i]
    y5_m3[i] = S850[160+i]
    error5_m3[i] = Error[160+i]
    
for i in range(0,5):
    x5_m4[i] = Mass_mean[230+i]
    y5_m4[i] = S850[230+i]
    error5_m4[i] = Error[230+i] 

# sixth redshift bin

x6_m1 = np.empty(5)
y6_m1 = np.empty(5)
error6_m1 = np.empty(5)

x6_m2 = np.empty(5)
y6_m2 = np.empty(5)
error6_m2 = np.empty(5)

x6_m3 = np.empty(5)
y6_m3 = np.empty(5)
error6_m3 = np.empty(5)

x6_m4 = np.empty(5)
y6_m4 = np.empty(5)
error6_m4 = np.empty(5)

print("bin 6")

for i in range(0,5):
    x6_m1[i] = Mass_mean[25+i]
    y6_m1[i] = S850[25+i]
    error6_m1[i] = Error[25+i]

for i in range(0,5):
    x6_m2[i] = Mass_mean[95+i]
    y6_m2[i] = S850[95+i]
    error6_m2[i] = Error[95+i]
    
for i in range(0,5):
    x6_m3[i] = Mass_mean[165+i]
    y6_m3[i] = S850[165+i]
    error6_m3[i] = Error[165+i]
    
for i in range(0,5):
    x6_m4[i] = Mass_mean[235+i]
    y6_m4[i] = S850[235+i]
    error6_m4[i] = Error[235+i] 

# seventh redshift bin

x7_m1 = np.empty(5)
y7_m1 = np.empty(5)
error7_m1 = np.empty(5)

x7_m2 = np.empty(5)
y7_m2 = np.empty(5)
error7_m2 = np.empty(5)

x7_m3 = np.empty(5)
y7_m3 = np.empty(5)
error7_m3 = np.empty(5)

x7_m4 = np.empty(5)
y7_m4 = np.empty(5)
error7_m4 = np.empty(5)

print("bin 7")

for i in range(0,5):
    x7_m1[i] = Mass_mean[30+i]
    y7_m1[i] = S850[30+i]
    error7_m1[i] = Error[30+i]

for i in range(0,5):
    x7_m2[i] = Mass_mean[100+i]
    y7_m2[i] = S850[100+i]
    error7_m2[i] = Error[100+i]
    
for i in range(0,5):
    x7_m3[i] = Mass_mean[170+i]
    y7_m3[i] = S850[170+i]
    error7_m3[i] = Error[170+i]
    
for i in range(0,5):
    x7_m4[i] = Mass_mean[240+i]
    y7_m4[i] = S850[240+i]
    error7_m4[i] = Error[240+i] 

# eighth redshift bin

x8_m1 = np.empty(5)
y8_m1 = np.empty(5)
error8_m1 = np.empty(5)

x8_m2 = np.empty(5)
y8_m2 = np.empty(5)
error8_m2 = np.empty(5)

x8_m3 = np.empty(5)
y8_m3 = np.empty(5)
error8_m3 = np.empty(5)

x8_m4 = np.empty(5)
y8_m4 = np.empty(5)
error8_m4 = np.empty(5)

print("bin 8")

for i in range(0,5):
    x8_m1[i] = Mass_mean[35+i]
    y8_m1[i] = S850[35+i]
    error8_m1[i] = Error[35+i]

for i in range(0,5):
    x8_m2[i] = Mass_mean[105+i]
    y8_m2[i] = S850[105+i]
    error8_m2[i] = Error[105+i]
    
for i in range(0,5):
    x8_m3[i] = Mass_mean[175+i]
    y8_m3[i] = S850[175+i]
    error8_m3[i] = Error[175+i]
    
for i in range(0,5):
    x8_m4[i] = Mass_mean[245+i]
    y8_m4[i] = S850[245+i]
    error8_m4[i] = Error[245+i] 

# ninth redshift bin

x9_m1 = np.empty(5)
y9_m1 = np.empty(5)
error9_m1 = np.empty(5)

x9_m2 = np.empty(5)
y9_m2 = np.empty(5)
error9_m2 = np.empty(5)

x9_m3 = np.empty(5)
y9_m3 = np.empty(5)
error9_m3 = np.empty(5)

x9_m4 = np.empty(5)
y9_m4 = np.empty(5)
error9_m4 = np.empty(5)

print("bin 9")

for i in range(0,5):
    x9_m1[i] = Mass_mean[40+i]
    y9_m1[i] = S850[40+i]
    error9_m1[i] = Error[40+i]

for i in range(0,5):
    x9_m2[i] = Mass_mean[110+i]
    y9_m2[i] = S850[110+i]
    error9_m2[i] = Error[110+i]
    
for i in range(0,5):
    x9_m3[i] = Mass_mean[180+i]
    y9_m3[i] = S850[180+i]
    error9_m3[i] = Error[180+i]
    
for i in range(0,5):
    x9_m4[i] = Mass_mean[250+i]
    y9_m4[i] = S850[250+i]
    error9_m4[i] = Error[250+i] 

# tenth redshift bin

x10_m1 = np.empty(5)
y10_m1 = np.empty(5)
error10_m1 = np.empty(5)

x10_m2 = np.empty(5)
y10_m2 = np.empty(5)
error10_m2 = np.empty(5)

x10_m3 = np.empty(5)
y10_m3 = np.empty(5)
error10_m3 = np.empty(5)

x10_m4 = np.empty(5)
y10_m4 = np.empty(5)
error10_m4 = np.empty(5)

print("bin 10")

for i in range(0,5):
    x10_m1[i] = Mass_mean[45+i]
    y10_m1[i] = S850[45+i]
    error10_m1[i] = Error[45+i]

for i in range(0,5):
    x10_m2[i] = Mass_mean[115+i]
    y10_m2[i] = S850[115+i]
    error10_m2[i] = Error[115+i]
    
for i in range(0,5):
    x10_m3[i] = Mass_mean[185+i]
    y10_m3[i] = S850[185+i]
    error10_m3[i] = Error[185+i]
    
for i in range(0,5):
    x10_m4[i] = Mass_mean[255+i]
    y10_m4[i] = S850[255+i]
    error10_m4[i] = Error[255+i] 
  
# 11th redshift bin

x11_m1 = np.empty(5)
y11_m1 = np.empty(5)
error11_m1 = np.empty(5)

x11_m2 = np.empty(5)
y11_m2 = np.empty(5)
error11_m2 = np.empty(5)

x11_m3 = np.empty(5)
y11_m3 = np.empty(5)
error11_m3 = np.empty(5)

x11_m4 = np.empty(5)
y11_m4 = np.empty(5)
error11_m4 = np.empty(5)

print("bin 11")

for i in range(0,5):
    x11_m1[i] = Mass_mean[50+i]
    y11_m1[i] = S850[50+i]
    error11_m1[i] = Error[50+i]

for i in range(0,5):
    x11_m2[i] = Mass_mean[120+i]
    y11_m2[i] = S850[120+i]
    error11_m2[i] = Error[120+i]
    
for i in range(0,5):
    x11_m3[i] = Mass_mean[190+i]
    y11_m3[i] = S850[190+i]
    error11_m3[i] = Error[190+i]
    
for i in range(0,5):
    x11_m4[i] = Mass_mean[260+i]
    y11_m4[i] = S850[260+i]
    error11_m4[i] = Error[260+i] 

# 12th redshift bin

x12_m1 = np.empty(5)
y12_m1 = np.empty(5)
error12_m1 = np.empty(5)

x12_m2 = np.empty(5)
y12_m2 = np.empty(5)
error12_m2 = np.empty(5)

x12_m3 = np.empty(5)
y12_m3 = np.empty(5)
error12_m3 = np.empty(5)

x12_m4 = np.empty(5)
y12_m4 = np.empty(5)
error12_m4 = np.empty(5)

print("bin 12")

for i in range(0,5):
    x12_m1[i] = Mass_mean[55+i]
    y12_m1[i] = S850[55+i]
    error12_m1[i] = Error[55+i]

for i in range(0,5):
    x12_m2[i] = Mass_mean[125+i]
    y12_m2[i] = S850[125+i]
    error12_m2[i] = Error[125+i]
    
for i in range(0,5):
    x12_m3[i] = Mass_mean[195+i]
    y12_m3[i] = S850[195+i]
    error12_m3[i] = Error[195+i]
    
for i in range(0,5):
    x12_m4[i] = Mass_mean[265+i]
    y12_m4[i] = S850[265+i]
    error12_m4[i] = Error[265+i] 

# 13th redshift bin

x13_m1 = np.empty(5)
y13_m1 = np.empty(5)
error13_m1 = np.empty(5)

x13_m2 = np.empty(5)
y13_m2 = np.empty(5)
error13_m2 = np.empty(5)

x13_m3 = np.empty(5)
y13_m3 = np.empty(5)
error13_m3 = np.empty(5)

x13_m4 = np.empty(5)
y13_m4 = np.empty(5)
error13_m4 = np.empty(5)

print("bin 13")

for i in range(0,5):
    x13_m1[i] = Mass_mean[60+i]
    y13_m1[i] = S850[60+i]
    error13_m1[i] = Error[60+i]

for i in range(0,5):
    x13_m2[i] = Mass_mean[130+i]
    y13_m2[i] = S850[130+i]
    error13_m2[i] = Error[130+i]
    
for i in range(0,5):
    x13_m3[i] = Mass_mean[200+i]
    y13_m3[i] = S850[200+i]
    error13_m3[i] = Error[200+i]
    
for i in range(0,5):
    x13_m4[i] = Mass_mean[270+i]
    y13_m4[i] = S850[270+i]
    error13_m4[i] = Error[270+i] 

# 14th redshift bin

x14_m1 = np.empty(5)
y14_m1 = np.empty(5)
error14_m1 = np.empty(5)

x14_m2 = np.empty(5)
y14_m2 = np.empty(5)
error14_m2 = np.empty(5)

x14_m3 = np.empty(5)
y14_m3 = np.empty(5)
error14_m3 = np.empty(5)

x14_m4 = np.empty(5)
y14_m4 = np.empty(5)
error14_m4 = np.empty(5)

print("bin 14")

for i in range(0,5):
    x14_m1[i] = Mass_mean[65+i]
    y14_m1[i] = S850[65+i]
    error14_m1[i] = Error[65+i]

for i in range(0,5):
    x14_m2[i] = Mass_mean[135+i]
    y14_m2[i] = S850[135+i]
    error14_m2[i] = Error[135+i]
    
for i in range(0,5):
    x14_m3[i] = Mass_mean[205+i]
    y14_m3[i] = S850[205+i]
    error14_m3[i] = Error[205+i]
    
for i in range(0,5):
    x14_m4[i] = Mass_mean[275+i]
    y14_m4[i] = S850[275+i]
    error14_m4[i] = Error[275+i] 
    
# Transform negative fluxes into 3sigma upper limits
# *****************************************************************************

# Bin 1

y1_m1_uplims = np.array([0,0,0,0,0],dtype=bool)
y1_m2_uplims = np.array([0,0,0,0,0],dtype=bool)
y1_m3_uplims = np.array([0,0,0,0,0],dtype=bool)
y1_m4_uplims = np.array([0,0,0,0,0],dtype=bool)

for i in range(0,5):
    
    if y1_m1[i]/error1_m1[i] < 0.0:
        y1_m1[i] = 3.0*error1_m1[i]
        y1_m1_uplims[i] = 1
        
    if y1_m2[i]/error1_m2[i] < 0.0:
        y1_m2[i] = 3.0*error1_m2[i]
        y1_m2_uplims[i] = 1
        
    if y1_m3[i]/error1_m3[i] < 0.0:
        y1_m3[i] = 3.0*error1_m3[i]
        y1_m3_uplims[i] = 1
        
        
    if y1_m4[i]/error1_m4[i] < 0.0:
        y1_m4[i] = 3.0 * error1_m4[i]
        y1_m4_uplims[i] = 1
         
        
# Bin 2
     
y2_m1_uplims = np.array([0,0,0,0,0],dtype=bool)
y2_m2_uplims = np.array([0,0,0,0,0],dtype=bool)
y2_m3_uplims = np.array([0,0,0,0,0],dtype=bool)
y2_m4_uplims = np.array([0,0,0,0,0],dtype=bool)

for i in range(0,5):
    
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
     
y3_m1_uplims = np.array([0,0,0,0,0],dtype=bool)
y3_m2_uplims = np.array([0,0,0,0,0],dtype=bool)
y3_m3_uplims = np.array([0,0,0,0,0],dtype=bool)
y3_m4_uplims = np.array([0,0,0,0,0],dtype=bool)

for i in range(0,5):
    
    if y3_m1[i]/error3_m1[i] < 0.0:
        y3_m1[i] = 3.0 * error3_m1[i]
        y3_m1_uplims[i] = 1
        
    if y3_m2[i]/error2_m3[i] < 0.0:
        y3_m2[i] = 3.0 * error3_m2[i]
        y3_m2_uplims[i] = 1  

    if y3_m3[i]/error3_m3[i] < 0.0:
        y3_m3[i] = 3.0 * error3_m3[i]
        y3_m3_uplims[i] = 1
        
    if y3_m4[i]/error3_m4[i] < 0.0:
        y3_m4[i] = 3.0 * error3_m4[i]  
        y3_m4_uplims[i] = 1


# Bin 4
     
y4_m1_uplims = np.array([0,0,0,0,0],dtype=bool)
y4_m2_uplims = np.array([0,0,0,0,0],dtype=bool)
y4_m3_uplims = np.array([0,0,0,0,0],dtype=bool)
y4_m4_uplims = np.array([0,0,0,0,0],dtype=bool)

for i in range(0,5):
    
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
     
y5_m1_uplims = np.array([0,0,0,0,0],dtype=bool)
y5_m2_uplims = np.array([0,0,0,0,0],dtype=bool)
y5_m3_uplims = np.array([0,0,0,0,0],dtype=bool)
y5_m4_uplims = np.array([0,0,0,0,0],dtype=bool)

for i in range(0,5):
    
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
        
# Bin 6
     
y6_m1_uplims = np.array([0,0,0,0,0],dtype=bool)
y6_m2_uplims = np.array([0,0,0,0,0],dtype=bool)
y6_m3_uplims = np.array([0,0,0,0,0],dtype=bool)
y6_m4_uplims = np.array([0,0,0,0,0],dtype=bool)

for i in range(0,5):
    
    if y6_m1[i]/error6_m1[i] < 0.0:
        y6_m1[i] = 3.0 * error6_m1[i]
        y6_m1_uplims[i] = 1
        
    if y6_m2[i]/error6_m2[i] < 0.0:
        y6_m2[i] = 3.0 * error6_m2[i]
        y6_m2_uplims[i] = 1

    if y6_m3[i]/error6_m3[i] < 0.0:
        y6_m3[i] = 3.0 * error6_m3[i]
        y6_m3_uplims[i] = 1

        
    if y6_m4[i]/error6_m4[i] < 0.0:
        y6_m4[i] = 3.0 * error6_m4[i] 
        y6_m4_uplims[i] = 1

        
# Bin 7
     
y7_m1_uplims = np.array([0,0,0,0,0],dtype=bool)
y7_m2_uplims = np.array([0,0,0,0,0],dtype=bool)
y7_m3_uplims = np.array([0,0,0,0,0],dtype=bool)
y7_m4_uplims = np.array([0,0,0,0,0],dtype=bool)

for i in range(0,5):
    
    if y7_m1[i]/error7_m1[i] < 0.0:
        y7_m1[i] = 3.0 * error7_m1[i]
        y7_m1_uplims[i] = 1
        
    if y7_m2[i]/error7_m2[i] < 0.0:
        y7_m2[i] = 3.0 * error7_m2[i]
        y7_m2_uplims[i] = 1


    if y7_m3[i]/error7_m3[i] < 0.0:
        y7_m3[i] = 3.0 * error7_m3[i]
        y7_m3_uplims[i] = 1
        
    if y7_m4[i]/error7_m4[i] < 0.0:
        y7_m4[i] = 3.0 * error7_m4[i] 
        y7_m4_uplims[i] = 1

# Bin 8
     
y8_m1_uplims = np.array([0,0,0,0,0],dtype=bool)
y8_m2_uplims = np.array([0,0,0,0,0],dtype=bool)
y8_m3_uplims = np.array([0,0,0,0,0],dtype=bool)
y8_m4_uplims = np.array([0,0,0,0,0],dtype=bool)

for i in range(0,5):
    
    if y8_m1[i]/error8_m1[i] < 0.0:
        y8_m1[i] = 3.0 * error8_m1[i]
        y8_m1_uplims[i] = 1
        
    if y8_m2[i]/error8_m2[i] < 0.0:
        y8_m2[i] = 3.0 * error8_m2[i]
        y8_m2_uplims[i] = 1

    if y8_m3[i]/error8_m3[i] < 0.0:
        y8_m3[i] = 3.0 * error8_m3[i]
        y8_m3_uplims[i] = 1
        
    if y8_m4[i]/error8_m4[i] < 0.0:
        y8_m4[i] = 3.0 * error8_m4[i]
        y8_m4_uplims[i] = 1

# Bin 9

y9_m1_uplims = np.array([0,0,0,0,0],dtype=bool)
y9_m2_uplims = np.array([0,0,0,0,0],dtype=bool)
y9_m3_uplims = np.array([0,0,0,0,0],dtype=bool)
y9_m4_uplims = np.array([0,0,0,0,0],dtype=bool)

for i in range(0,5):
    
    if y9_m1[i]/error9_m1[i] < 0.0:
        y9_m1[i] = 3.0 * error9_m1[i]
        y9_m1_uplims[i] = 1
        
    if y9_m2[i]/error9_m2[i] < 0.0:
        y9_m2[i] = 3.0 * error9_m2[i]
        y9_m2_uplims[i] = 1

    if y9_m3[i]/error9_m3[i] < 0.0:
        y9_m3[i] = 3.0 * error9_m3[i]
        y9_m3_uplims[i] = 1
        
        
    if y9_m4[i]/error9_m4[i] < 0.0:
        y9_m4[i] = 3.0 * error9_m4[i]
        y9_m4_uplims[i] = 1
        
# Bin 10
   
y10_m1_uplims = np.array([0,0,0,0,0],dtype=bool)
y10_m2_uplims = np.array([0,0,0,0,0],dtype=bool)
y10_m3_uplims = np.array([0,0,0,0,0],dtype=bool)
y10_m4_uplims = np.array([0,0,0,0,0],dtype=bool)  

for i in range(0,5):
    
    if y10_m1[i]/error10_m1[i] < 0.0:
        y10_m1[i] = 3.0 * error10_m1[i]
        y10_m1_uplims[i] = 1
        
    if y10_m2[i]/error10_m2[i] < 0.0:
        y10_m2[i] = 3.0 * error10_m2[i]
        y10_m2_uplims[i] = 1

    if y10_m3[i]/error10_m3[i] < 0.0:
        y10_m3[i] = 3.0 * error10_m3[i]
        y10_m3_uplims[i] = 1
        
    if y10_m4[i]/error10_m4[i] < 0.0:
        y10_m4[i] = 3.0 * error10_m4[i]
        y10_m4_uplims[i] = 1
        
# Bin 11
 
y11_m1_uplims = np.array([0,0,0,0,0],dtype=bool)
y11_m2_uplims = np.array([0,0,0,0,0],dtype=bool)
y11_m3_uplims = np.array([0,0,0,0,0],dtype=bool)
y11_m4_uplims = np.array([0,0,0,0,0],dtype=bool)      

for i in range(0,5):
    
    if y11_m1[i]/error11_m1[i] < 0.0:
        y11_m1[i] = 3.0 * error11_m1[i]
        y11_m1_uplims[i] = 1
             
    if y11_m2[i]/error11_m2[i] < 0.0:
        y11_m2[i] = 3.0 * error11_m2[i]
        y11_m2_uplims[i] = 1
        
    if y11_m3[i]/error11_m3[i] < 0.0:
        y11_m3[i] = 3.0 * error11_m3[i]
        y11_m3_uplims[i] = 1
        
    if y11_m4[i]/error11_m4[i] < 0.0:
        y11_m4[i] = 3.0 * error11_m4[i]
        y11_m4_uplims[i] = 1
        
# Bin 12
     
y12_m1_uplims = np.array([0,0,0,0,0],dtype=bool)
y12_m2_uplims = np.array([0,0,0,0,0],dtype=bool)
y12_m3_uplims = np.array([0,0,0,0,0],dtype=bool)
y12_m4_uplims = np.array([0,0,0,0,0],dtype=bool)  

for i in range(0,5):
    
    if y12_m1[i]/error12_m1[i] < 0.0:
        y12_m1[i] = 3.0 * error12_m1[i]
        y12_m1_uplims[i] = 1
        
    if y12_m2[i]/error12_m2[i] < 0.0:
        y12_m2[i] = 3.0 * error12_m2[i]
        y12_m2_uplims[i] = 1

    if y12_m3[i]/error12_m3[i] < 0.0:
        y12_m3[i] = 3.0 * error12_m3[i]
        y12_m3_uplims[i] = 1
        
    if y12_m4[i]/error12_m4[i] < 0.0:
        y12_m4[i] = 3.0 * error12_m4[i]   
        y12_m4_uplims[i] = 1

# Bin 13

y13_m1_uplims = np.array([0,0,0,0,0],dtype=bool)
y13_m2_uplims = np.array([0,0,0,0,0],dtype=bool)
y13_m3_uplims = np.array([0,0,0,0,0],dtype=bool)
y13_m4_uplims = np.array([0,0,0,0,0],dtype=bool)  
     

for i in range(0,5):
    
    if y13_m1[i]/error13_m1[i] < 0.0:
        y13_m1[i] = 3.0 * error13_m1[i]
        y13_m1_uplims[i] = 1
        
    if y13_m2[i]/error13_m2[i] < 0.0:
        y13_m2[i] = 3.0 * error13_m2[i]
        y13_m2_uplims[i] = 1

    if y13_m3[i]/error13_m3[i] < 0.0:
        y13_m3[i] = 3.0 * error13_m3[i]
        y13_m3_uplims[i] = 1
        
    if y13_m4[i]/error13_m4[i] < 0.0:
        y13_m4[i] = 3.0 * error13_m4[i] 
        y13_m4_uplims[i] = 1

# Bin 14
     
y14_m1_uplims = np.array([0,0,0,0,0],dtype=bool)
y14_m2_uplims = np.array([0,0,0,0,0],dtype=bool)
y14_m3_uplims = np.array([0,0,0,0,0],dtype=bool)
y14_m4_uplims = np.array([0,0,0,0,0],dtype=bool)  

for i in range(0,5):
    
    if y14_m1[i]/error14_m1[i] < 0.0:
        y14_m1[i] = 3.0 * error14_m1[i]
        y14_m1_uplims[i] = 1
        
    if y14_m2[i]/error14_m2[i] < 0.0:
        y14_m2[i] = 3.0 * error14_m2[i]
        y14_m2_uplims[i] = 1

    if y14_m3[i]/error14_m3[i] < 0.0:
        y14_m3[i] = 3.0 * error14_m3[i]
        y14_m3_uplims[i] = 1
        
    if y14_m4[i]/error14_m4[i] < 0.0:
        y14_m4[i] = 3.0 * error14_m4[i]
        y14_m4_uplims[i] = 1

# Plot out results
# ****************************************************************************
    
fig = plt.figure(figsize=(10.0,10.0))

pos_ann = ([9.5,4.0])

xl = np.empty(2)
yl = np.empty(2)
xl[0] = 9.0
xl[1] = 11.5
yl[:] = 0.1


# First panel

f1 = plt.axes([0.15,0.1,0.2,0.2])
f1.set_xlim(9.0,11.5)
f1.set_ylim(0.01,9.0)
f1.set_yscale('log')
f1.tick_params(axis='both',which='both',labelsize=18)

x1_m1 = x1_m1 + 0.033
x1_m2 = x1_m2 + 0.1
x1_m3 = x1_m3 - 0.033
x1_m4 = x1_m4 - 0.1

f1.annotate('0.2<z<0.5',pos_ann,fontweight='bold',size=15)

f1.plot(x1_m1,y1_m1,'rs')
f1.errorbar(x1_m1,y1_m1,yerr= error1_m1,fmt='r',ls='none',uplims=y1_m1_uplims)

f1.plot(x1_m2,y1_m2,'bo')
f1.errorbar(x1_m2,y1_m2,yerr= error1_m2,fmt='b',ls='none',uplims=y1_m2_uplims)

f1.plot(x1_m3,y1_m3,'g*')
f1.errorbar(x1_m3,y1_m3,yerr= error1_m3,fmt='g',ls='none',uplims=y1_m3_uplims)

f1.plot(x1_m4,y1_m4,'cp')
f1.errorbar(x1_m4,y1_m4,yerr= error1_m4,fmt='c',ls='none',uplims=y1_m4_uplims)

f1.plot(xl,yl,'k--')

# second panel

f1 = plt.axes([0.35,0.1,0.2,0.2])
f1.set_xlim(9.0,11.5)
f1.set_ylim(0.01,9.0)
f1.set_yscale('log')
f1.tick_params(axis='x',labelsize=18)
f1.set_yticks([])

f1.annotate('0.5<z<0.8',pos_ann,fontweight='bold',size=15)

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

f1.plot(xl,yl,'k--')

# Third panel

f1 = plt.axes([0.55,0.1,0.2,0.2])
f1.set_xlim(9.0,11.5)
f1.set_ylim(0.01,9.0)
f1.set_xlabel('$log(M_*)$',size=25)
f1.set_yscale('log')
f1.tick_params(axis='x',labelsize=18)
f1.set_yticks([])


f1.annotate('0.8<z<1.1',pos_ann,fontweight='bold',size=15)

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

f1.plot(xl,yl,'k--')


# Fourth panel

f1 = plt.axes([0.75,0.1,0.2,0.2])
f1.set_xlim(9.0,11.5)
f1.set_ylim(0.01,9.0)
f1.tick_params(axis='x',labelsize=18)
f1.set_yscale('log')
f1.set_yticks([])

f1.annotate('1.1<z<1.5',pos_ann,fontweight='bold',size=15)

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

f1.plot(xl,yl,'k--')


# Fifth panel

f1 = plt.axes([0.15,0.3,0.2,0.2])
f1.set_xlim(9.0,11.5)
f1.set_ylim(0.01,9.0)
f1.set_ylabel('$Flux/mJy$',size=25)
f1.tick_params(axis='y',labelsize=18)
f1.set_yscale('log')
f1.set_xticks([])

f1.annotate('1.5<z<2.0',pos_ann,fontweight='bold',size=15)

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


f1.plot(xl,yl,'k--')


# Sixth panel

f1 = plt.axes([0.35,0.3,0.2,0.2])
f1.set_xlim(9.0,11.5)
f1.set_ylim(0.01,9.0)
f1.set_yscale('log')
f1.set_xticks([])
f1.set_yticks([])

f1.annotate('2.0<z<2.5',pos_ann,fontweight='bold',size=15)

x6_m1 = x6_m1 + 0.033
x6_m2 = x6_m2 + 0.1
x6_m3 = x6_m3 - 0.033
x6_m4 = x6_m4 - 0.1

f1.plot(x6_m1,y6_m1,'rs')
f1.errorbar(x6_m1,y6_m1,yerr= error6_m1,fmt='r',ls='none',uplims=y6_m1_uplims)

f1.plot(x6_m2,y6_m2,'bo')
f1.errorbar(x6_m2,y6_m2,yerr= error6_m2,fmt='b',ls='none',uplims=y6_m2_uplims)

f1.plot(x6_m3,y6_m3,'g*')
f1.errorbar(x6_m3,y6_m3,yerr= error6_m3,fmt='g',ls='none',uplims=y6_m3_uplims)

f1.plot(x6_m4,y6_m4,'cp')
f1.errorbar(x6_m4,y6_m4,yerr= error6_m4,fmt='c',ls='none',uplims=y6_m4_uplims)


f1.plot(xl,yl,'k--')

# Seventh panel

f1 = plt.axes([0.55,0.3,0.2,0.2])
f1.set_xlim(9.0,11.5)
f1.set_ylim(0.01,9.0)
f1.set_yscale('log')
f1.set_xticks([])
f1.set_yticks([])

f1.annotate('2.5<z<3.0',pos_ann,fontweight='bold',size=15)

x7_m1 = x7_m1 + 0.033
x7_m2 = x7_m2 + 0.1
x7_m3 = x7_m3 - 0.033
x7_m4 = x7_m4 - 0.1

f1.plot(x7_m1,y7_m1,'rs')
f1.errorbar(x7_m1,y7_m1,yerr= error7_m1,fmt='r',ls='none',uplims=y7_m1_uplims)

f1.plot(x7_m2,y7_m2,'bo')
f1.errorbar(x7_m2,y7_m2,yerr= error7_m2,fmt='b',ls='none',uplims=y7_m2_uplims)

f1.plot(x7_m3,y7_m3,'g*')
f1.errorbar(x7_m3,y7_m3,yerr= error7_m3,fmt='g',ls='none',uplims=y7_m3_uplims)

f1.plot(x7_m4,y7_m4,'cp')
f1.errorbar(x7_m4,y7_m4,yerr= error7_m4,fmt='c',ls='none',uplims=y7_m4_uplims)

f1.plot(xl,yl,'k--')


# Eighth panel

f1 = plt.axes([0.75,0.3,0.2,0.2])
f1.set_xlim(9.0,11.5)
f1.set_ylim(0.01,9.0)
f1.set_yscale('log')
f1.set_xticks([])
f1.set_yticks([])

f1.annotate('3.0<z<3.5',pos_ann,fontweight='bold',size=15)

x8_m1 = x8_m1 + 0.033
x8_m2 = x8_m2 + 0.1
x8_m3 = x8_m3 - 0.033
x8_m4 = x8_m4 - 0.1

f1.plot(x8_m1,y8_m1,'rs')
f1.errorbar(x8_m1,y8_m1,yerr= error8_m1,fmt='r',ls='none',uplims=y8_m1_uplims)

f1.plot(x8_m2,y8_m2,'bo')
f1.errorbar(x8_m2,y8_m2,yerr= error8_m2,fmt='b',ls='none',uplims=y8_m2_uplims)

f1.plot(x8_m3,y8_m3,'g*')
f1.errorbar(x8_m3,y8_m3,yerr= error8_m3,fmt='g',ls='none',uplims=y8_m3_uplims)

f1.plot(x8_m4,y8_m4,'cp')
f1.errorbar(x8_m4,y8_m4,yerr= error8_m4,fmt='c',ls='none',uplims=y8_m4_uplims)


f1.plot(xl,yl,'k--')

# Ninth panel

f1 = plt.axes([0.15,0.5,0.2,0.2])
f1.set_xlim(9.0,11.5)
f1.set_ylim(0.01,9.0)
f1.set_yscale('log')
f1.set_xticks([])
f1.tick_params(axis='y',labelsize=18)

f1.annotate('3.5<z<4.5',pos_ann,fontweight='bold',size=15)

x9_m1 = x9_m1 + 0.033
x9_m2 = x9_m2 + 0.1
x9_m3 = x9_m3 - 0.033
x9_m4 = x9_m4 - 0.1

f1.plot(x9_m1,y9_m1,'rs')
f1.errorbar(x9_m1,y9_m1,yerr= error9_m1,fmt='r',ls='none',uplims=y9_m1_uplims)

f1.plot(x9_m2,y9_m2,'bo')
f1.errorbar(x9_m2,y9_m2,yerr= error9_m2,fmt='b',ls='none',uplims=y9_m2_uplims)

f1.plot(x9_m3,y9_m3,'g*')
f1.errorbar(x9_m3,y9_m3,yerr= error9_m3,fmt='g',ls='none',uplims=y9_m3_uplims)

f1.plot(x9_m4,y9_m4,'cp')
f1.errorbar(x9_m4,y9_m4,yerr= error9_m4,fmt='c',ls='none',uplims=y9_m4_uplims)

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

x10_m1 = x10_m1 + 0.033
x10_m2 = x10_m2 + 0.1
x10_m3 = x10_m3 - 0.033
x10_m4 = x10_m4 - 0.1

f1.plot(x10_m1,y10_m1,'rs')
f1.errorbar(x10_m1,y10_m1,yerr= error10_m1,fmt='r',ls='none',uplims=y10_m1_uplims)

f1.plot(x10_m2,y10_m2,'bo')
f1.errorbar(x10_m2,y10_m2,yerr= error10_m2,fmt='b',ls='none',uplims=y10_m2_uplims)

f1.plot(x10_m3,y10_m3,'g*')
f1.errorbar(x10_m3,y10_m3,yerr= error10_m3,fmt='g',ls='none',uplims=y10_m3_uplims)

f1.plot(x10_m4,y10_m4,'cp')
f1.errorbar(x10_m4,y10_m4,yerr= error10_m4,fmt='c',ls='none',uplims=y10_m4_uplims)

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

x11_m1 = x11_m1 + 0.033
x11_m2 = x11_m2 + 0.1
x11_m3 = x11_m3 - 0.033
x11_m4 = x11_m4 - 0.1

f1.plot(x11_m1,y11_m1,'rs')
f1.errorbar(x11_m1,y11_m1,yerr= error11_m1,fmt='r',ls='none',uplims=y11_m1_uplims)

f1.plot(x11_m2,y11_m2,'bo')
f1.errorbar(x11_m2,y11_m2,yerr= error11_m2,fmt='b',ls='none',uplims=y11_m2_uplims)

f1.plot(x11_m3,y11_m3,'g*')
f1.errorbar(x11_m3,y11_m3,yerr= error11_m3,fmt='g',ls='none',uplims=y11_m3_uplims)

f1.plot(x11_m4,y11_m4,'cp')
f1.errorbar(x11_m4,y11_m4,yerr= error11_m4,fmt='c',ls='none',uplims=y11_m4_uplims)

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

x12_m1 = x12_m1 + 0.033
x12_m2 = x12_m2 + 0.1
x12_m3 = x12_m3 - 0.033
x12_m4 = x12_m4 - 0.1

f1.plot(x12_m1,y12_m1,'rs')
f1.errorbar(x12_m1,y12_m1,yerr= error12_m1,fmt='r',ls='none',uplims=y12_m1_uplims)

f1.plot(x12_m2,y12_m2,'bo')
f1.errorbar(x12_m2,y12_m2,yerr= error12_m2,fmt='b',ls='none',uplims=y12_m2_uplims)

f1.plot(x12_m3,y12_m3,'g*')
f1.errorbar(x12_m3,y12_m3,yerr= error12_m3,fmt='g',ls='none',uplims=y12_m3_uplims)

f1.plot(x12_m4,y12_m4,'cp')
f1.errorbar(x12_m4,y12_m4,yerr= error12_m4,fmt='c',ls='none',uplims=y12_m4_uplims)


f1.plot(xl,yl,'k--')

# 13th panel

f1 = plt.axes([0.15,0.7,0.2,0.2])
f1.set_xlim(9.0,11.5)
f1.set_ylim(0.01,9.0)
f1.set_yscale('log')
f1.set_xticks([])
f1.tick_params(axis='y',labelsize=18)

f1.annotate('7.5<z<8.5',pos_ann,fontweight='bold',size=15)

x13_m1 = x13_m1 + 0.033
x13_m2 = x13_m2 + 0.1
x13_m3 = x13_m3 - 0.033
x13_m4 = x13_m4 - 0.1

f1.plot(x13_m1,y13_m1,'rs')
f1.errorbar(x13_m1,y13_m1,yerr= error13_m1,fmt='r',ls='none',uplims=y13_m1_uplims)

f1.plot(x13_m2,y13_m2,'bo')
f1.errorbar(x13_m2,y13_m2,yerr= error13_m2,fmt='b',ls='none',uplims=y13_m2_uplims)

f1.plot(x13_m3,y13_m3,'g*')
f1.errorbar(x13_m3,y13_m3,yerr= error13_m3,fmt='g',ls='none',uplims=y13_m3_uplims)

f1.plot(x13_m4,y13_m4,'cp')
f1.errorbar(x13_m4,y13_m4,yerr= error13_m4,fmt='c',ls='none',uplims=y13_m4_uplims)


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

x14_m1 = x14_m1 + 0.033
x14_m2 = x14_m2 + 0.1
x14_m3 = x14_m3 - 0.033
x14_m4 = x14_m4 - 0.1

f1.plot(x14_m1,y14_m1,'rs')
f1.errorbar(x14_m1,y14_m1,yerr= error14_m1,fmt='r',ls='none',uplims=y14_m1_uplims)

f1.plot(x13_m2,y13_m2,'bo')
f1.errorbar(x14_m2,y14_m2,yerr= error14_m2,fmt='b',ls='none',uplims=y14_m2_uplims)

f1.plot(x14_m3,y14_m3,'g*')
f1.errorbar(x14_m3,y14_m3,yerr= error14_m3,fmt='g',ls='none',uplims=y14_m3_uplims)

f1.plot(x14_m4,y14_m4,'cp')
f1.errorbar(x14_m4,y14_m4,yerr= error14_m4,fmt='c',ls='none',uplims=y14_m4_uplims)

f1.plot(xl,yl,'k--')

fig.savefig('Figure4.pdf')
plt.show()





            
            
    
    

    
        
    
    



    





