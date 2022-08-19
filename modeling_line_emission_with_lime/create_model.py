#!/usr/bin/python3

import numpy as np

# Script to convert a protoplanetary disk model to a c header file for reading into lime

# Model file
infile = 'HD_163296_disk_structure.npz'

# Read in model data into master data array
data = np.load('HD_163296_disk_structure.npz')

# These data are in 2D spherical coordinates, need to 
#	convert to cylindrical coordinated

# Assign columns to model variables
#	rmod = radius (AU)
#	theta = angle(rad)
#	hmod = height (AU)
# 	dmod = density (cm-3)
#	tgmod = gas temperature (K)
#	tdmod = dust temperature (K)

r = data['r']
a = data['theta']

rsph,theta = np.meshgrid(r,a,indexing='ij')

rmod = rsph*np.sin(theta)
hmod = rsph*np.cos(theta)

dmod  = data['n_H2']
tgmod = data['tgas']
tdmod = data['td2']

# Create an empty array with the same shape for the abundance

amod = np.full_like(rmod,0.0)

# Set nans in density array to zero

dmod[np.isnan(dmod)] = 0

# Set total number of grid points

nmod = len(r)*len(a)

# Set here the rules for the abundance structure, you can set
#	- a temperature range
#	- a radial range
#	- a vertical range, or
#	- a combination of the above

# In this sample script, we will place the molecule in a warm molecular layer
#	with a constant fractional abundance

tgmin = 30.0
tgmax = 70.0

amin = 1e-16
amax = 1e-09

for i in range(len(r)):
	for j in range(len(a)):
		amod[i][j] = amin	
		if ((tgmod[i][j] >= tgmin) and (tgmod[i][j] <= tgmax)):			
			amod[i][j] = amax
		if(dmod[i][j] == 'nan'):
			dmod[i][j] = 0.0
	
# Lime uses SI unit: convert radius and height from AU to m 
#	and density from cm-3 to m-3

# Conversion factors

AUtoM	= 1.4959787e+11
CM3toM3	= 1.0e-06

rmod = rmod*AUtoM
hmod = hmod*AUtoM

dmod = dmod/CM3toM3

# Set here the molecular weight (amu) and the star mass (Msol)

mmol  = 51.0
mstar = 2.40

# Set here the disk inclination and position angle

dinc = '46.7*PI/180.0'
dpa  = '133.3*PI/180.0'

# Open c head file for writing

f = open('HC3N_model.h','w')

# Write preamble for model.h

f.write('#include <stdio.h>\n')
f.write('#ifndef _model_h\n')
f.write('#define _model_h\n')
f.write('\n')

# Set model variables here

f.write('#define mmol  %4.1f \t\t /* Molecular weight (amu) */ \n'     % (mmol))
f.write('#define mstar %4.2f \t\t /* Stellar mass (Msol) */ \n'        % (mstar))
f.write('#define dinc  %s \t /* Disk inclination (radians) */ \n'    % (dinc))
f.write('#define dpa   %s \t /* Disk position angle (radians) */ \n' % (dpa))

f.write('#define nmod  %i   \n' % (nmod))
f.write('\n')

# Set physical structure arrays

f.write('double rmod[nmod] = {\n')
for i in range(len(r)):
	for j in range(len(a)): 
		f.write('%12.6e,\n' % (rmod[i][j]))	
f.write('};\n')

f.write('\n')

f.write('double hmod[nmod] = {\n')
for i in range(len(r)):
	for j in range(len(a)): 
		f.write('%12.6e,\n' % (hmod[i][j]))	
f.write('};\n')

f.write('\n')

f.write('double dmod[nmod] = {\n')
for i in range(len(r)):
	for j in range(len(a)): 
		f.write('%12.6e,\n' % (dmod[i][j]))	
f.write('};\n')

f.write('\n')

f.write('double tgmod[nmod] = {\n')
for i in range(len(r)):
	for j in range(len(a)): 
		f.write('%12.6e,\n' % (tgmod[i][j]))	
f.write('};\n')

f.write('\n')

f.write('double tdmod[nmod] = {\n')
for i in range(len(r)):
	for j in range(len(a)): 
		f.write('%12.6e,\n' % (tdmod[i][j]))	
f.write('};\n')

f.write('\n')

f.write('double amod[nmod] = {\n')
for i in range(len(r)):
	for j in range(len(a)):  
		f.write('%12.6e,\n' % (amod[i][j]))	
f.write('};\n')

f.write('\n')

f.write('#endif\n')
