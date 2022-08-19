#################################################################################
#
# This script will walk you through the steps on how to generate value-added
#	products from a spectral cube, including moment maps, P-V diagrams,
#	line profiles and fluxes, and radial profiles.
#
# We will use your simulated images of HC3N 11-10 and 29-28 emission.
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

# Start up python an import the following modules

import os
import casatools
from casatasks import imsmooth
from casatasks import immoments
from casatasks import exportfits
from casatasks import specflux

#################################################################################
# Smooth images
#################################################################################

# We first need to smooth the simulated cubes to the same resolution as the data

# We can do this using the casa task, imsmooth:
#	https://casadocs.readthedocs.io/en/stable/api/tt/casatasks.analysis.imsmooth.html

my_model_fits = ' '

imsmooth(imagename=my_model_fits,
	major = ' ',
	minor = ' ',
	pa = ' ',
	outfile = my_model_fits+'.im.smooth')
	
# Export the smoothed maps to fits to use with bettermoments

exportfits(imagename=my_model_fits+'.im.smooth',
	fitsimage = my_model_fits+'.im.smooth.fits')

#################################################################################
# Data inspection and S/N
#################################################################################

# As done in the imaging step, we can use casa/imview or casaviewer to inspect the 
#	image cubes. This is your next step.

# If the emission looks very blobby or pixelated and not smooth, it is recommended
#	to go back to the lime step and increase the resolution of the lime grid

# If running on your laptop, there will be a limit on the size of model you can run
#	so we recommend running lime on the cluster if available

# From the visual inspection of both cubes, consider the following:

#	- what is roughly the spatial extent of the emission (in au)?
#	- over how many channels is emission clearly detected?
#	- what is the likely spatial morphology of the emission?
#	- how does the emission compare with that presented in Ilee et al. (2021)

#################################################################################
# Moment maps
#################################################################################

# To make moment maps of the emission we use the casa task, immoments:
#	https://casadocs.readthedocs.io/en/latest/api/tt/casatasks.analysis.immoments.html?casatasks.analysis.immoments

# We can do this step outside of casa, using the python module, casatasks

# You can either execute the following within interactive python, or create
#	a script to be executed on the command line

my_model_image = ' '

# Using the in-built CASA routines
###################################

# Define a range of channels over which to calculate the moments
#	selected by visual inspection of the smoothed model cubes

mychans = ' '

# Moment 0 (integrated intensity) - use no intensity cuts

immoments(imagename=my_model_image, 
	moments=[0], 
	chans=mychans, 
	outfile=my_model_cube+'.mom0')

# Moment 1 (intensity weighted velocity) - use no intensity cuts for models

immoments(imagename=my_model_image, 
	moments=[1], 
	chans=mychans, 
	outfile=my_model_cube+'.mom1')

# Moment 2 (intensity weighted dispersion) - use no intensity cuts for models

immoments(imagename=my_model_image, 
	moments=[2], 
	chans=mychans, 
	outfile=my_model_cube+'.mom2')

# Moment 8 (peak values) - use no intensity cuts

immoments(imagename=my_model_image, 
	moments=[8], 
	chans=mychans, 
	outfile=my_model_cube+'.mom8')

# We can use casa/imview or casaviewer to inspect the moment maps.

# When you are happy with your images, you can export these to fits format for further analysis:

for infile in [my_model_image+'.mom0', my_model_image+'.mom1', my_model_image+'.mom2', my_model_imagee+'.mom8']:
    exportfits(imagename=infile, fitsimage=infile+".fits", dropdeg=True, overwrite=True) 

# Using the bettermoments package
##################################

# We can also use bettermoments to make moment maps

# Need to use bettermoments on the command line via os

# Set a fits file to work with
my_model_fits = ' '

# Moment 0 (integrated intensity) - use no intensity cuts

# With a mask 
# mom0_cmd = "bettermoments " + my_image_fits + " -clip 0 -mask " + my_mask_fits + " -method zeroth"

# Without a mask 
mom0_cmd = "bettermoments " + my_model_fits + " -method zeroth"
os.system(mom0_cmd)

# Now look at the other methods available in the bettermoments documentation and 
# see how they compare to the CASA implementations

# Bettermoments makes fits images by default so no need for export

#################################################################################
# Disk-integrated line profile and flux
#################################################################################

# Decide on the regions size you will extract line profiles and fluxes

flux_region    = 'circle[[%sdeg, %sdeg], %.1farcsec]' % ( , ,)

# Using the in-built CASA routines
###################################

my_model_image = ' '

# We can do this also two ways, first using the casa task, specflux:
#	https://casadocs.readthedocs.io/en/stable/api/tt/casatasks.analysis.specflux.html

specflux(imagename=my_model_image,
	region = flux_region,
	logfile=' ')

# This will print the total flux and line profile to the logfile

# Run this for both lines and make a note of the total flux

# HD 163296 HC3N 29-28 => 
# HD 163296 HC3N 11-10 => 

# Compare these with the values from Ilee et al. (2021)

# Now make a plot of the disk-integrated line profiles using matplotlib
#	or your favourite plotting package

# Compare your line profiles from those presented in Ilee et al. 2021 and
#	those produced by the other team

# Using the gofish package
##########################

# gofish also enables extraction of an averaged line profile

import numpy as np
import matplotlib.pyplot as plt
from gofish import imagecube

my_model_fits = ' '

inclination =        # degrees, 90 = edge on 
PA =                 # degrees
Mstar =              # solar masses
distance =           # parsecs
radius =             # distance out to which to stack in disk, arcsec

cube = imagecube(my_model_fits, FOV=5.0)

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

my_model_fits = ' '

cube = imagecube(my_model_fits, FOV=5.0)

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
	
# Compare the radial profiles obtained with those presented in Ilee et al. and
#	those produced by the other group

# What do the radial profiles tell us about the morphology of the simulated 
#	emission? How might you adapt the models to better reproduce the 
#	observations?

#################################################################################
# END OF ANALYSIS EXERCISE
#################################################################################

# Once you have extracted the disk-integrated flux of both lines, you can now 
#	estimate the column density and rotational temperature
#
# Now move over to to rotation_diagram_script.py
#
# Any concerns or questions at this point, contact a tutor

#################################################################################





















