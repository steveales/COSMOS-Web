# Program to do stacking analysis on S2COSMOS image
# using the COSMOS-Web catalogue. The program uses a Monte-Carlo
# simulation to generate the errors.
#
# Last edited: 29th August 2025

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename
from astropy import wcs

# Constants
# ****************************************************************************
# ****************************************************************************

# No. of trials in the Monte-Carlo simulations

nmonte = 1000

# size of array

nx = 3307
ny = 3307


### choose the redshift slices and the bin width

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

# Set up variables and arrays
# ****************************************************************************

# start mass bin
lowest_mass = 9.0

# Set FWHM 850u,
fwhm850 = 14.0

# set pixel size
pix850 = 2.0

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

# Check the bin limits

for i in range(0,nbin):
    print(z_lower_limit[i],z_upper_limit[i],mass_lower_limit[i],mass_upper_limit[i])

print(nbin)

# Read in galaxy data
# ****************************************************************************

file = 'COSMOS-web 0.021-stacking_catalogue_flag_star.dat'

galaxy_data = []

galaxy_data = np.loadtxt(file)

print(galaxy_data.shape)
ngal = galaxy_data.shape[0]
print(ngal)

# Print out the first ten values as a check

for i in range(0,10):
    print(galaxy_data[i,0],galaxy_data[i,1],galaxy_data[i,2],\
          galaxy_data[i,3],galaxy_data[i,4])

# Read in images
# ****************************************************************************

hdu1 = fits.open(get_pkg_data_filename('S2COSMOS_20180927_850_fcf_mf_crop.fits'))
hdu2 = fits.open(get_pkg_data_filename('S2COSMOS_20180927_850_err_mf_crop.fits'))

# Extract data, noise and header
# ****************************************************************************

image850_header = hdu1[0].header

image850_3D = hdu1[0].data
print(image850_3D.shape)
noise850_3D = hdu2[0].data
print(noise850_3D.shape)

# reduce 3D SCUBA images to 2D images

image850 = np.empty((3307,3307))
noise850 = np.empty((3307,3307))

for i in range(0,3307):
    
    for j in range(0,3307):
        image850[i,j] = image850_3D[0,i,j]
        noise850[i,j] = noise850_3D[0,i,j]
        
# Do Matt's fixes to sort out the headers

image850_header['NAXIS'] = 2
image850_header["i_naxis"] = 2
del(image850_header['NAXIS3'])
del(image850_header["CRPIX3"])

try:
    del(image850_header["CDELT3"])
except:
    del(image850_header['CD3_3'])

del(image850_header["CRVAL3"])
del(image850_header["CTYPE3"])
del(image850_header["LBOUND3"])
del(image850_header["CUNIT3"])

# Extract the astrometric information

c850 = wcs.WCS(image850_header)

# Calculate the statistics of the image excluding the 'nans'
# ****************************************************************************

# How many nans are there?

oned = np.ravel(image850)
num_pix = len(oned)

print("Number of pixels in image: ",num_pix)

isnan_array = np.isnan(oned)

num_nan = sum(isnan_array)

num_good = num_pix - num_nan

print("Number of data values: ",num_good)

# Do the statistics

rmean = np.nanmean(image850)

image850 = image850 - rmean

rstd = np.nanstd(image850)

rstd = rstd / np.sqrt(float(num_good))

print("Mean of image is: ",rmean,' +/- ',rstd)

# Calculate stacking signal for each bin
# ****************************************************************************

binInfo = {}
work = np.zeros((ngal,5))

results = np.zeros((nbin,5))

print(results.shape)

for k in range(0,nbin):
    binID = k
    
# Find all the galaxies in this range of stellar mass

    num_target = 0
    
    mass_low = mass_lower_limit[k]
    mass_high = mass_upper_limit[k]
    z_low = z_lower_limit[k]
    z_high = z_upper_limit[k]

# save info to binInfo

    binInfo[binID] = {"mass_low": mass_low,\
                      "mass_high": mass_high, "z_low": z_low, "z_high": z_high}

    sel3 = np.where((galaxy_data[:,3] >= mass_low) & (galaxy_data[:,3] < mass_high) &\
                    (galaxy_data[:,2] >= z_low) & (galaxy_data[:,2] < z_high))

    num_target = len(sel3[0])

        # save info to binInfo

    binInfo[binID] = {"mass_low": mass_low, "mass_high": mass_high, "z_low": z_low,\
                      "z_high": z_high, "num_in_bin": num_target}
    print(f"Mass: {mass_low}-{mass_high}, Redshift: {z_low}-{z_high}, No. galaxies: {num_target}")
    
    ra = galaxy_data[:,0][sel3]
    dec = galaxy_data[:,1][sel3]

        # convert the positions of the targets in RA and dec into x and y positions

    xpos,ypos = c850.wcs_world2pix(ra,dec,0)
    xgal = np.round(xpos).astype(int)
    ygal = np.round(ypos).astype(int)
    
# measure the signals and noises at the position of each galaxy

    work[0:ngal,0:3] = -99.0

    num_good = 0

    for i in range(0,num_target):
    
        if ygal[i] > -1 and ygal[i] < ny and xgal[i] > -1 and xgal[i] < nx:
    
            if np.isnan(image850[ygal[i],xgal[i]]) != True:
        
                if np.isnan(noise850[ygal[i],xgal[i]]) != True:
                    num_good = num_good + 1
                    work[i,0] = image850[ygal[i],xgal[i]]
                    work[i,1] = noise850[ygal[i],xgal[i]]
                    work[i,2] = 1.0/(noise850[ygal[i],xgal[i]]**2.0)

# combining the signals

    totweight = 0.0
    totmean = 0.0

    for i in range(0,num_target):
    
        if work[i,0] > -98.0:
            totmean = totmean + work[i,0]*work[i,2]
            totweight = totweight + work[i,2]
        
    totmean = totmean/totweight

    err = np.sqrt(1.0/totweight)
    
    results[k,0] = totmean
    results[k,1] = err
    
    print("mean is: ",totmean)
#   print("error is: ",error)

# Monte Carlo simulation
# *****************************************************************************

    nsample = num_good

    data_monte = np.empty((2,nmonte))
    work_monte = np.empty((3,nsample))

    for i in range(0,nmonte):
        work_monte[0:3,0:nsample] = -99.0

# set up random positions

        xpos = np.random.random(nsample)
        xpos = xpos * (nx-1)

        ypos = np.random.random(nsample)
        ypos = ypos * (ny-1)

# Measure signals and noises at the position of each galaxy

        num_good = 0

        for j in range(0,nsample):

            xgal = int(xpos[j]+0.5)
            ygal = int(ypos[j]+0.5)

            if np.isnan(image850[ygal,xgal]) != True:

                if np.isnan(noise850[ygal,xgal]) != True:
                    num_good = num_good + 1
                    work_monte[0,j] = image850[ygal,xgal]
                    work_monte[1,j] = noise850[ygal,xgal]
                    work_monte[2,j] = 1.0/(noise850[ygal,xgal]**2.0)

# Now combine the signals

        totweight = 0.0
        totmean = 0.0

        for j in range(0,nsample):

            if work_monte[0,j] > -98.0:
                totmean = totmean + work_monte[0,j] * work_monte[2,j]
                totweight = totweight + work_monte[2,j]

        totmean = totmean/totweight

        error = np.sqrt(1.0/totweight)

        data_monte[0,i] = totmean
        data_monte[1,i] = error

# Combine results to calculate the mean of the means and the error

    totmean = 0.0
    totvar = 0.0

    for i in range(0,nmonte):
        totmean = totmean + data_monte[0,i]

    totmean = totmean/float(nmonte)

    for i in range(0,nmonte):
        tp = data_monte[0,i]
        totvar = totvar + tp*tp

    totvar = totvar / float(nmonte)
    sd = np.sqrt(totvar)

    error_mean = sd /np.sqrt(float(nmonte))
    
    results[k,2] = totmean
    results[k,3] = error_mean
    results[k,4] = sd
    
#   print("Monte Carlo mean: ",totmean," +/- ",error_mean)
    
    print("Monte Carlo error: ",sd)
    
    mass_low = binInfo[k]['mass_low']
    mass_high = binInfo[k]['mass_high']
    z_low = binInfo[k]['z_low']
    z_high = binInfo[k]['z_high']
    print(mass_low,mass_high,z_low,z_high,results[k,0],results[k,4])

print(results.shape)

with open('results_flag_star.dat', 'a') as file:
    print(nbin)

    for i in range(0,nbin):
        mass_low = binInfo[i]['mass_low']
        mass_high = binInfo[i]['mass_high']
        z_low = binInfo[i]['z_low']
        z_high = binInfo[i]['z_high']
        num_in_bin = binInfo[i]['num_in_bin']
               
        file.write(f"{z_low} ")
        file.write(f"{z_high} ")
        file.write(f"{mass_low} ")
        file.write(f"{mass_high} ")
        file.write(f"{results[i,0]} ")
        file.write(f"{results[i,4]}" + "\n")
        
        print(mass_low,mass_high,z_low,z_high,results[i,0],results[i,4])
        
        
        
        
        
        
        