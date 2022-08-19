#################################################################################
#
# This script will walk you through the steps on how to generate value-added
#	products from a spectral cube, including moment maps, P-V diagrams,
#	line profiles and fluxes, and radial profiles.
#
# We will use the MAPS images of HC3N 11-10 and 29-28 emission.
#
# We will also estimate the disk-averaged column density and rotational
#	temperature using the disk-integrated fluxes of the two lines.
#
# There are blanks purposefully left in the script for you to make decisions
# 	or that require some research from the literature.
#
# We recommend that you execute each block of the script in python independently 
#	- the script provided here is not meant to be executed in its entirety
#	- copy and paste into a python command prompt
#	- alternatively you can copy and paste into an executible python script
#
# John and Catherine and other tutors/lecturers are on hand to answer any
#	questions, or help to troubleshoot any issues!  Have fun!
#

#################################################################################
# Preparations
#################################################################################

# Before you begin make sure you have a recent version of CASA installed, and that
#	you have Analysis Utilities and Reduction Utilities also downloaded and
#	accessible from within CASA
# Also make sure you have installed gofish and casatasks as python modules

#################################################################################
# Data inspection and S/N
#################################################################################

# As done in the imaging step, we can use casa/imview or casaviewer to inspect the 
#	image cubes. This is your first step.

# From the visual inspection of both cubes, consider the following:

#	- what is roughly the spatial extent of the emission (in au)?
#	- over how many channels is emission clearly detected?
#	- what is the likely spatial morphology of the emission?

# We can estimate the S/N in both images from within casa using a function in
#	reduction, start-up casa if needed

% casa

# We first need to determine the line free channels from visual inspection

# WARNING - can only select one range combined with a region!

# HD 163296 29-28 line-free channels: 
# HD 163296 11-10line-free channels: 

# or

# MWC 480 29-28 line-free channels: 
# MWC 480 11-10 line-free channels: 

# And we need to set a region to only cover the spatial extent of the source

# This is because the noise varies across the image

# HD_163296 coordinates => 

# or

# MWC 480 coordinates => 

CASA <X>: hd163296_ra  = ' '
CASA <X>: hd163296_dec = ' '
CASA <X>: hd163296_extent = 

# or

CASA <X>: mwc480_ra  = ' '
CASA <X>: mwc480_dec = ' '
CASA <X>: mwc480_extent = 

CASA <X>: hd163296_region  = 'circle[[%s, %s], %.1farcsec]' % (hd163296_ra,hd163296_dec,hd163296_extent)

# or

CASA <X>: mwc480_region    = 'circle[[%s, %s], %.1farcsec]' % (mwc480_ra,mwc480_dec,mwc480_extent)

# Now we load the function and modules from reduction_utils.py

CASA <X>: execfile('PATH-TO/reduction_utils.py') 

CASA <X>: importfits(fitsimage=' ',
		imagename=' ')

CASA <X>: estimate_MAPSSNR(imagename=' ',
		region = ,
		chan_range=' ')

# COPY FIRST TRANSITION IMAGE STATS HERE

CASA <X>: importfits(fitsimage='',
		imagename='')

CASA <X>: estimate_MAPSSNR(imagename=' ',
		region = ,
		chan_range='')

# COPY SECOND TRANSITION IMAGE STATS HERE

# Compare your numbers with those reported in the MAPS overview paper.

#################################################################################
# Moment maps
#################################################################################

# To make moment maps of the emission we use the casa task, immoments:
#	https://casadocs.readthedocs.io/en/latest/api/tt/casatasks.analysis.immoments.html?casatasks.analysis.immoments

# We can do this step outside of casa, using the python module, casatasks

# You can either execute the following within interactive python, or create
#	a script to be executed on the command line

import os
import casatools
from casatasks import immoments
from casatasks import exportfits
from casatasks import specflux

my_image_cube = ' '

# Using the in-built CASA routines
###################################

# Measure a noise level in the image cube to set cut levels
rms =     # Jy/beam

# Define a range of channels over which to calculate the moments
mychans = ' '

# Moment 0 (integrated intensity) - use no intensity cuts

immoments(imagename=my_image_cube, 
	moments=[0], 
	chans=mychans, 
	includepix=[-1000,1000], 
	outfile=my_image_cube+'.mom0')

# Moment 1 (intensity weighted velocity) - currently set to only use pixels >3*rms

immoments(imagename=my_image_cube, 
	moments=[1], 
	chans=mychans, 
	includepix=[3*rms,1000], 
	outfile=my_image_cube+'.mom1')

# Moment 2 (intensity weighted dispersion) - currently set to only use pixels >3*rms

immoments(imagename=my_image_cube, 
	moments=[2], 
	chans=mychans, 
	includepix=[3*rms,1000], 
	outfile=my_image_cube+'.mom2')

# Moment 8 (peak values) - use no intensity cuts

immoments(imagename=my_image_cube, 
	moments=[8], 
	chans=mychans, 
	includepix=[-1000,1000], 
	outfile=my_image_cube+'.mom8')

# As done in the imaging step, we can use casa/imview or casaviewer to inspect the 
#	moment maps.

# When you are happy with your images, you can export these to fits format for further analysis:

for infile in [my_image_cube+'.mom0', my_image_cube+'.mom1', my_image_cube+'.mom2', my_image_cube+'.mom8']:
    exportfits(imagename=infile, fitsimage=infile+".fits", dropdeg=True, overwrite=True) 

# Using the bettermoments package
##################################

# We can also use bettermoments to make moment maps

# Need to use bettermoments on the command line via os

# Need a fits file for better moments
my_image_fits = ' '

# Moment 0 (integrated intensity) - use no intensity cuts

# With a mask 
# mom0_cmd = "bettermoments " + my_image_fits + " -clip 0 -mask " + my_mask_fits + " -method zeroth"

# Without a mask 
mom0_cmd = "bettermoments " + my_image_fits + " -method zeroth"
os.system(mom0_cmd)

# Now look at the other methods available in the bettermoments documentation and 
# see how they compare to the CASA implementations

# Bettermoments makes fits images by default so no need for export

#################################################################################
# Disk-integrated line profile and flux
#################################################################################

# Using the in-built CASA routines
###################################

# We can do this also two ways, first using the casa task, specflux:
#	https://casadocs.readthedocs.io/en/stable/api/tt/casatasks.analysis.specflux.html

specflux(imagename=' ',
	region =  ,
	logfile='')

specflux(imagename='',
	region = ,
	logfile='')

# This will print the total flux and line profile to the logfile

# Run this for both lines and make a note of the total flux

# HD 163296 HC3N 11-10 => 171 mJy km/s
# HD 163296 HC3N 29-28 => 187 mJy km/s

# or

# MWC 480 HC3N 11-10 => 168 mJy km/s
# MWC 480 HC3N 11-10 => 99.4 mJy km/s

# Compare these with the values from Ilee et al. (2021)

# How would one estimate the error on these fluxes?

# Have a think then consult a tutor (hint - think about the mathematics 
#	used to compute the total flux and error propagation ...)

# Now make a plot of the disk-integrated line profiles using matplotlib
#	or your favourite plotting package

# Compare your line profiles from those presented in Ilee et al. 2021

# Using the gofish package
##########################

# gofish also enables extraction of an averaged line profile

import numpy as np
import matplotlib.pyplot as plt
from gofish import imagecube

my_image_cube = ' '

inclination =        # degrees, 90 = edge on 
PA =                 # degrees
Mstar =              # solar masses
distance =           # parsecs
radius =             # distance out to which to stack in disk, arcsec

cube = imagecube(my_image_cube, FOV=5.0)

x, y, dy = cube.average_spectrum(r_min=0.0, 
				 r_max=radius, 
				 inc=inclination,
                                   PA=PA, 
				 mstar=Mstar, 
				 dist=distance)

# Now we can plot the data

fig, ax = plt.subplots()
ax.errorbar(x, y, dy, fmt=' ', capsize=1.25, capthick=1.25, color='k', lw=1.0)
ax.step(x, y, where='mid', color='k', lw=1.0)
ax.set_xlabel('Velocity (m/s)')
ax.set_ylabel('Line Flux (Jy/beam)')
fig.show()

# Save the line profile data as .txt files
avelineprofDATA = np.array([x, y, dy]).T

np.savetxt('ave_line_profile.txt', 
	avelineprofDATA, header='v [velocity], F [Jy/beam], F_err [Jy/beam]', 
	fmt=['%f','%f','%f'], 
	delimiter=', ')

# We can also use gofish to calculate an integrated spectrum

x, y, dy = cube.integrated_spectrum(r_min=0.0, 
                                    r_max=radius, 
                                    inc=inclination, 
                                    PA=PA,
                                    mstar=Mstar, 
                                    dist=distance,
                                    resample=1,
                                    include_spectral_decorrelation=True)

# Now we can plot the data

fig, ax = plt.subplots()
ax.errorbar(x, y, dy, fmt=' ', capsize=1.25, capthick=1.25, color='k', lw=1.0)
ax.step(x, y, where='mid', color='k', lw=1.0)
ax.set_xlabel('Velocity (m/s)')
ax.set_ylabel('Integrated Flux (Jy)')
fig.show()

# Save the profile data as .txt files
intlineprofDATA = np.array([x, y, dy]).T

np.savetxt('int_line_profile.txt', 
	intlineprofDATA, header='v [velocity], F [Jy], F_err [Jy]', 
	fmt=['%f','%f','%f'], 
	delimiter=', ')

# You will need to integrate the spectrum (e.g. x vs. y) numerically to calculate
# an integrated line flux in mJy km/s for the rotational analysis

# We leave this as an exercise for you to do in python and to compare the numbers 
#	obtained with the specflux method

#################################################################################
# Radial profiles
#################################################################################

# Go fish also allows us to export radial profiles of the emission

my_image_cube = ' '

cube = imagecube(my_image_cube, FOV=5.0)

inclination =       # degrees, 90 = edge on 
PA =                # degrees

# Calculate the radial profile (currently assumes a flat emission surface, full azimuth)

x, y, dy = cube.radial_profile(inc=inclination, PA=PA)

# Check for and remove nans
iR = ~np.isnan(y)

# Save the profile data as .txt files
radprofDATA = np.array([x[iR], y[iR], dy[iR]]).T

np.savetxt('radial_profile.txt', 
	radprofDATA, header='R [arcsec], I [mJy/beam km/s], I_err [mJy/beam km/s]', 
	fmt=['%f','%f','%f'], 
	delimiter=', ')
	
# Compare the radial profiles obtained with those presented in Ilee et al.

# What do the radial profiles tell us about the morphology of the emission?	

#################################################################################
# END OF ANALYSIS EXERCISE
#################################################################################

# Once you have measured the disk-integrated flux of both lines, you can now 
#	estimate the column density and rotational temperature
#
# Now move over to to rotation_diagram_script.py
#
# Any concerns or questions at this point, contact a tutor

#################################################################################





















