#################################################################################
#
# This script will walk you through the steps on how to inspect, self-calibrate, 
# 	and image molecular line emission from a protoplanetary disk.
#
# We will use the ALMA science verification data for TW Hya as an example.
#
# We will follow the self-calibration procedure as outlined by the MAPS large
# 	program team (Molecules with ALMA on Planet-forming Scales).
#
# There are blanks purposefully left in the script for you to make decisions
# 	or that require some research from the literature.
#
# We recommend that you execute each block of the script in CASA independently 
#	- the script provided here is not meant to be executed in its entirety
#	- copy and paste into the casa command prompt
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

# If it does not already exist make a new file "startup.py" in ~/.casa/.
# Add the following code

import sys, os
sys.path.append("/PATH-TO-THE/analysis_scripts/")
import analysisUtils as au

# Make sure you give the full path with no short-cuts

# Start up casa in a new terminal window

% casa

# This will open a casa prompt and the casa logger
# Check that Analysis Utilities has loaded correctly

CASA <X>: au.help()    

# You should see a big list of the available commands if successful
# Contact a tutor if you are having difficulties at this point

#################################################################################
# Inspection of the data
#################################################################################

# First set the path to the TW Hya measurement set
# Make sure to change the path if the data directory is not located where you have
#	started up casa

CASA <X>: inpvis = ' '

# Before proceeding we will first take a look at the calibrated data
# We can do this in casa using the command, listobs: 
#	https://casadocs.readthedocs.io/en/stable/api/tt/casatasks.information.listobs.html#casatasks.information.listobs

# You can check  what input a casa task needs using the command, inp:
	
CASA <X>: inp listobs

# Run listobs on the TW Hya measurement set

CASA <X>: listobs(vis=inpvis,verbose=True)

# This will output a huge amount of information to the casa logger
# Let's save it to a file so that we can inspect it properly

CASA <X>: listobs(vis=inpvis,verbose=True,listfile=' ')

# From the output of listobs ,determine the following:

#	- how many executions/scans were conducted?
#	- how many spectral windows were set?
#	- the spectral resolution of the data?
#	- how many antennas were used?

# We can use some of our supplementary scripts to get additional information

# Let's determine the number and range of baselines for these data

CASA <X>: au.getBaselineLengths(inpvis,sort=True)

# From the output of this task, determine the following:

#	- the number of independent baselines? => 36
#	- the minumum baseline? => 16.5m
#	- the maximum baseline? => 91.0m

# Given the maximum baseline and observing frequency, can you estimate the 
#	spatial resolution of the data?

# Let's now determine the time spent on source (total observing time on source)

CASA <X>: au.timeOnSource(inpvis,verbose=False)

# From the output of this task, determine the following:

#	- the time spent on source? => 139 min
#	- the percentage of observing time spent on source? => 49.5%

# Now we'll take a look at the actual data: we do this using the casa task, plotms:
#	https://casadocs.readthedocs.io/en/stable/api/tt/casaplotms.plotms.html

CASA <X>: inp plotms

# This has a huge number of options, first we will look at the data in frequency space
# 	with no averaging

# This might take a little while depending on your machine ...

CASA <X>: plotms(vis=inpvis,xaxis='frequency',yaxis='amplitude')

# You should see clearly the number of spectral windows and how they are arranged 
#	in frequency

# Let's now export a plot after applying some averaging to reduce the datasize
# We will average the data over scan and time (in seconds)

CASA <X>: plotms(vis=inpvis,xaxis='frequency',yaxis='amplitude',avgscan=True,avgtime='60',plotfile=' ')

# This should be much faster ...

# From this plot, determine the following:

#	- how many spectral lines are clearly detected?

# We can also take a look at the uv coverage of the data in plotms

CASA <X>: plotms(vis=inpvis,xaxis='uwave',yaxis='vwave',avgscan=True,avgtime='60',plotfile=' ')

# From this plot, consider the following:

#	- how does the uv coverage compare to an idealised case?
#	- what impact might that have on the imaging?

#################################################################################
# Flagging of the data
#################################################################################

# The first step to conduct self-calibration is to extract the continuum from the
#	data, i.e., merging all line-free channels

# We can create a new measurement set in casa using the casa task, split:
#	https://casadocs.readthedocs.io/en/stable/api/tt/casatasks.manipulation.split.html#casatasks.manipulation.split

CASA <X>: inp split

# We will average over channels to reduce the file size, cut off the noisy
#	edges in the spectral windows, and only retain the 'data' column

CASA <X>: split(vis=inpvis,outputvis=' ', spw=' ', width=' ', datacolumn='data')  

# This might take a little while ...

# Let's now inspect this new measurement set using plotms and averaging over the
#	total time on source

CASA <X>: contvis = ' '

CASA <X>: plotms(vis=contvis,spw='0~3',
		xaxis='channel',
		yaxis='amp',
		avgtime='1e8',
		avgscan=True,
		coloraxis='spw',
		iteraxis='spw',
		xselfscale=True,
		showgui=True)

# Cycle through the individual spws using the green arrows in the gui

# We now need to flag the data before further processing:
#	- we flag/mask any channels containing lines
#	- we flag/mask any erroneous or particularly noisy data

# We can flag the data automatically using the casa task, flagdata:
#	https://casadocs.readthedocs.io/en/stable/api/tt/casatasks.flagging.flagdata.html#casatasks.flagging.flagdata

# Consult with a tutor on which data you should flag

CASA <X>: flagdata(vis=contvis, mode='manual',spw=' ')

# Re-run plotms again to check that the correct data have been flagged

CASA <X>: plotms(vis=contvis,spw='0~3',
		xaxis='channel',
		yaxis='amp',
		avgtime='1e8',
		avgscan=True,
		coloraxis='spw',
		iteraxis='spw',
		xselfscale=True,
		showgui=True)

# If happy, we now have a continuum only measurement set!

# Next we will examine the data and averaging over channels

CASA <X>: plotms(vis=contvis,
		xaxis='time',
		yaxis='amp',
		ydatacolumn='data',
		avgchannel='38',
		coloraxis='spw',
		plotfile='TW-Hya_B7cont_time-v-amp.png',
		showgui=True)

CASA <X>: plotms(vis=contvis,
		xaxis='uvdist',
		yaxis='amp',
		ydatacolumn='data',
		avgchannel='38',
		coloraxis='spw',
		plotfile='TW-Hya_B7cont_uvdist-v-amp.png',
		showgui=True)

# From this plot, consider the following:

#	- how does the ammplitude vary with uv distance (or baseline)?
#	- what does this tell you about the source emission?
#	- how does the amplitude vary with spw/frequency?
#	- is this what is expected?

# Finally run listobs over these data to extract information needed later

CASA <X>: listobs(vis=contvis, verbose=True, listfile=' ')

#################################################################################
# Initial imaging of the data
#################################################################################

# Need to load the reduction_utils.py script

CASA <X>: execfile('PATH-TO/reduction_utils.py') 	

# Now that we have a continuum only measurement set, we will proceed with 
#	the self-calibration step

# Let's first estimate the theoretical sensitivity reached by the data, using the 
#	ALMA sensitivity calculator
#	https://almascience.nrao.edu/proposing/sensitivity-calculator

# Theoretical sensitivity continuum =>  mJy/beam for   GHz bandwidth
# Theoretical sensitivity lines     =>  mJy/beam for   km/s channel resolution

# We will now make our first image of the continuum data using the casa task, tclean:
#	https://casadocs.readthedocs.io/en/stable/api/tt/casatasks.imaging.tclean.html#casatasks.imaging.tclean

# First we define a mask, within which we expect to find emission, to be used
#	later for imaging

# Casa regions: https://casadocs.readthedocs.io/en/latest/notebooks/image_analysis.html#Region-Files

# TW_Hya coordinates => 

CASA <X>: mask_pa  =    #position angle of mask in degrees                                                          
CASA <X>: mask_maj =    #semimajor axis of mask in arcsec                                                      
CASA <X>: mask_min =    #semiminor axis of mask in arcsec                                                      
CASA <X>: mask_ra  = ' '
CASA <X>: mask_dec = ' '
CASA <X>: tclean_mask  = 'ellipse[[%s, %s], [%.1farcsec, %.1farcsec], %.1fdeg]' % (mask_ra, mask_dec, mask_maj, mask_min, mask_pa)

# Now we make a so-called dirty image with no "cleaning", simply the Fourier transform of the visibilities to an image
# To do this we set niter = 0

# We also need to set some parameters for our images

# Cell size: ~ expected beam divided by 6-7
# Expected beam =>  ~  rad ~  arcsec

CASA <X>: cellsize_im =  ' '

# Image size: ~ primary beam 
# Primary beam @ 345 GHz => 

CASA <X>: imsize_im = 

CASA <X>: tclean(vis=contvis,
		imagename=' ',
      		specmode='mfs', 
      		imsize=imsize_im,
		cell=cellsize_im,
      		weighting='briggs',
		robust=0.5, 
		niter=0)

# We can inspect this image using the casa task, imview

CASA <X>: imview

# This opens a gui where you can open up and overlay images

# Have a play with imview for a short while, can you figure out how to:

# 	- change the colour map range, scaling, and palette?
#	- change the x and y axis from absolute to relative coordinates?
#	- turn on the colour wedge?

# WARNING - the casa viewer is very buggy and crashes randomly, you may have to 
#	quit and restart your casa session ...

# From this image, consider the following:

#	- does the data look how you might expect?
#	- does the noise off source appear as you might expect?

# Next we estimate the rms noise and signal-to-noise in this image; first we define a noise annulus

CASA <X>: noise_annulus = "annulus[[%s, %s],['%.2farcsec', '8.0arcsec']]" % \
                	(mask_ra, mask_dec, 1.1*mask_maj) 

CASA <X>: estimate_SNR(' ', disk_mask=tclean_mask, noise_mask=noise_annulus)

# LIST IMAGE STATISTICS HERE

# How does the rms from the dirty image compare with the theoretical sensivity?

# Now we do our first proper "clean"!
# Set threshold to ~ 3 sigma (theoretical noise for the continuum) for self-cal

CASA <X>: noise_im = ' '

# We have already defined a mask and a threshold, so we will run a non-interactive clean

CASA <X>: tclean(vis=contvis,
		imagename=' ',
      		specmode='mfs', 
      		imsize=imsize_im,
		cell=cellsize_im,
      		weighting='briggs',
		robust=0.5, 
		mask = tclean_mask,
		niter=10000,
		threshold=noise_im,
		interactive=False)
		
# Use imview to visually inspect this image, what has noticeably changed?

# Next we estimate the rms noise and signal-to-noise in this cleaned image
		
CASA <X>: estimate_SNR(' ', disk_mask=tclean_mask, noise_mask=noise_annulus)

# LIST IMAGE STATISTICS HERE

# Compare these numbers with those from the dirty image:

#	- what has changed and why?
#	- how does the image rms compare with the theoretical sensitivity?

#################################################################################
# Self-calibration of the data
#################################################################################

# Now we use the clean model generated in the initial imaging to self-calibrate 
#	the data

# The first task we run is gaincal:
#	https://casadocs.readthedocs.io/en/latest/api/tt/casatasks.calibration.gaincal.html#casatasks.calibration.gaincal

# Set solint initially to ~ 3 times the integration time step (check listobs)
# gaintype = 'T' averages the polarisations

CASA <X>: inp gaincal

CASA <X>: gaincal(vis=contvis,
		caltable=' ',
        		solint=' ',
		gaintype='T',
        		refant='DV06',
		spw='',
		minblperant=4,
        		calmode='p',
		minsnr=2)

# Now we inspect the solutions using the casa task, plotms:

CASA <X>: plotms(vis=' ',
		xaxis='time',
		yaxis='phase',
        		spw='',
		field='',
		antenna='1~8',
		showmajorgrid=True,
		showminorgrid=True,
		coloraxis = 'antenna1',
		iteraxis='antenna',
		gridrows=3,
		gridcols=3,
		plotrange=[0,0,-80,80],
		plotfile=' ',
		overwrite=True,
		showgui=True)

# This creates an image with an overly long filename, let's change it back using mv:

CASA <X>: !mv long_filename.png short_filename.png

# Let's inspect the calibration tables

CASA <X>: !open short_filename.png

# You will see the magnitude of corrections needed for each antenna with respect 
#	to the reference antenna, DV06, as a function of observing time

# Let's now apply these corrections to the data using the casa task, applycal:
#	https://casadocs.readthedocs.io/en/latest/api/tt/casatasks.calibration.applycal.html?casatasks.calibration.applycal

CASA <X>: applycal(vis=contvis,
		gaintable=[' '],
		applymode='calonly',
		interp='linear',
		calwt=True)

# We run a tclean again to see whether or not this has improved the image

CASA <X>: tclean(vis=contvis,
		imagename=' ',
      		specmode='mfs', 
      		imsize=imsize_im,
		cell=cellsize_im,
      		weighting='briggs',
		robust=0.5, 
		mask = tclean_mask,
		niter=10000,
		threshold=noise_im,
		interactive=False,
		savemodel='modelcolumn')

CASA <X>: estimate_SNR(' ', disk_mask=tclean_mask, noise_mask=noise_annulus)

# LIST IMAGE STATISTICS HERE

# Inspect this new image using imview, and consider the following:
#	https://casadocs.readthedocs.io/en/v6.2.1/api/tt/casatasks.visualization.imview.html

CASA <X>: imview(raster=[{'file':' ','range':[-0.01,1.0]},
			{'file':' ','range':[-0.01,1.0]},
       			{'file':' ','range':[-0.01,1.0]}])

#	- has the image quality improved?
#	- how have the rms and S/N changed after one round of phase calibration?

# You might also want to inspect the residual images using imview ...

# Now that we have a better model, we do another round of phase calibration, this
#	time using the integration timestep

CASA <X>: gaincal(vis=contvis,
		caltable=' ',
        		solint='int',
		gaintype='T',
        		refant='DV06',
		spw='',
		minblperant=4,
        		calmode='p',
		minsnr=2)
		
# As before, we plot the new solutions

CASA <X>: plotms(vis=' ',
		xaxis='time',
		yaxis='phase',
        		spw='',
		field='',
		antenna='1~8',
		showmajorgrid=True,
		showminorgrid=True,
		coloraxis = 'antenna1',
		iteraxis='antenna',
		gridrows=3,
		gridcols=3,
		plotrange=[0,0,-80,80],
		plotfile=' ',
		overwrite=True,
		showgui=True)l

# This creates an image with an overly long filename, let's change it back:

CASA <X>: !mv long_filename.png short_filename.png

# Let's inspect the calibration tables

CASA <X>: !open short_filename.png

# And if happy, we apply the solutions again		
		
CASA <X>: applycal(vis=contvis,
		gaintable=[' '],
		applymode='calonly',
		interp='linear',
		calwt=True)

# We run a tclean again to see whether or not this has improved the image

CASA <X>: tclean(vis=contvis,
		imagename=' ',
      		specmode='mfs', 
      		imsize=imsize_im,
		cell=cellsize_im,
      		weighting='briggs',
		robust=0.5, 
		mask = tclean_mask,
		niter=10000,
		threshold=noise_im,
		interactive=False,
		savemodel='modelcolumn')

CASA <X>: estimate_SNR(' ', disk_mask=tclean_mask, noise_mask=noise_annulus)

# LIST IMAGE STATS HERE

# From the image statistics, consider the following:

#	- has the second round of phase calibration improved the signal-to-noise ratio?
#	- if so by how much?

# Also visually inspect the two recent images for an improvement in image quality
# We use a negative scaling to highlight small-scale features ...

CASA <X>: imview(raster=[{'file':' ','range':[-0.01,1.0],'scaling':-2},
       			{'file':' ','range':[-0.01,1.0],'scaling':-2}])

# Should we continue with phase calibration?  

# If not, we move on to amplitude calibration
# Amplitude changes on longer timescales than phase, so we use solint='inf'

CASA <X>: gaincal(vis=contvis,
		caltable=' ',
       	 	solint='inf',
		gaintype='T',
        		refant='DV06',
		minblperant=4,
        		gaintable=[' '],
       		calmode='ap',
		minsnr=2)

# Now we inspect the phase solutions using the casa task, plotms:

CASA <X>: plotms(vis=' ',
		xaxis='time',
		yaxis='phase',
        		spw='',
		field='',
		antenna='1~8',
		showmajorgrid=True,
		showminorgrid=True,
		coloraxis = 'antenna1',
		iteraxis='antenna',
		gridrows=3,
		gridcols=3,
		plotrange=[0,0,-1,1],
		plotfile=' ',
		overwrite=True,
		showgui=True)

# ... and also amplitude

CASA <X>: plotms(vis=' ',
		xaxis='time',
		yaxis='amp',
        		spw='',
		field='',
		antenna='1~8',
		showmajorgrid=True,
		showminorgrid=True,
		coloraxis = 'antenna1',
		iteraxis='antenna',
		gridrows=3,
		gridcols=3,
		plotrange=[0,0,0.6,1.4],
		plotfile=' ',
		overwrite=True,
		showgui=True)

# And finally apply all solutions to the data

CASA <X>: applycal(vis=contvis,
		gaintable=[' ',' '],
		interp='linear',
		calwt=True)

# Let's do a final inspection of the self-calibrated data in plotms

CASA <X>: plotms(vis=contvis,
		xaxis='time',
		yaxis='amp',
		ydatacolumn='corrected',
		avgchannel='38',
		coloraxis='spw',
		plotfile=' ',
		showgui=True)

CASA <X>: plotms(vis=contvis,
		xaxis='uvdist',
		yaxis='amp',
		ydatacolumn='corrected',
		avgchannel='38',
		coloraxis='spw',
		plotfile=' ',
		showgui=True)

# Compare these images with the ones from earlier, what has changed?

#################################################################################
# Final continuum image
#################################################################################

# Before we make the final continuum image, we will first fix the phase centre of 
#	the data (this will also centre the image)

# We first make a test image to determine the phase shift needed

CASA <X>: tclean(vis=contvis,
		imagename=' ',
		datacolumn='corrected',
      		specmode='mfs', 
      		imsize=imsize_im,
		cell=cellsize_im,
      		weighting='briggs',
		robust=0.5, 
		mask=tclean_mask,
		niter=10000,
		threshold=noise_im,
		interactive=False)

# We first fit a gaussian to this test continuum image

CASA <X>: fit_gaussian(' ', region=tclean_mask)

# LIST RESULTS OF IMFIT HERE

# We now apply this phase shift to the measurement set

CASA <X>: peak_pos   = ' '
CASA <X>: newvis     = ' '

CASA <X>: fixvis(vis=contvis,outputvis=newvis,phasecenter=peak_pos)
CASA <X>: fixplanets(vis=newvis,field='0',direction=peak_pos)

# Now we make our final image 

CASA <X>: tclean(vis=newvis,
		imagename=' ',
		datacolumn='corrected',
      		specmode='mfs', 
      		imsize=imsize_im,
		cell=cellsize_im,
      		weighting='briggs',
		robust=0.5, 
		mask=tclean_mask,
		niter=10000,
		threshold=noise_im,
		interactive=False)

CASA <X>: estimate_SNR(' ', disk_mask=tclean_mask, noise_mask=noise_annulus)

# LIST IMAGE STATS HERE

# From these image statistics, consider the following:

#	- by how much has the rms improved from i) the dirty image, and ii) the first clean image
# 	- by how much has the SNR improved from i) the dirty image, and ii) the first clean image

# - rms:  mJy/beam /  mJy/beam  = 
# - rms:   mJy/beam / mJy/beam  = 

# - SNR:  /  = 
# - SNR:  /  = 

#################################################################################
# Apply self-calibration solutions to full unaveraged dataset
#################################################################################

# First make a copy in the current directory of the original data using split

CASA <X>: allvis = ' '

CASA <X>: split(vis=inpvis,outputvis=allvis,datacolumn = 'all')

# Now apply the solutions of the phase and amplitude calibration

CASA <X>: applycal(vis=allvis,
		gaintable=[' ',' '],
		interp=['linear','linear'],
		applymode='calonly',
		calwt=True)

# This might take a wee while ...

# And now apply the phase shift identified in the continuum fitting

CASA <X>: finalvis = ' '

CASA <X>: fixvis(vis=allvis,outputvis=finalvis,phasecenter=peak_pos)
CASA <X>: fixplanets(vis=finalvis,field='0',direction=peak_pos)

# Now we split out the line-containing spws for further processing
# Use the output of listobs to determine the number of the spw 

# CO 3-2   @ 345.796 GHz => spw =  
# HCO+ 4-3 @ 356.734 GHz => spw =  

CASA <X>: CO32vis   = ' '
CASA <X>: HCOp43vis = ' '

CASA <X>: split(vis=finalvis, outputvis=CO32vis, spw=' ' ,datacolumn = 'corrected')
CASA <X>: split(vis=finalvis, outputvis=HCOp43vis, spw=' ' ,datacolumn = 'corrected')

#################################################################################
# Continuum subtraction of molecular line data
#################################################################################

# To image the molecular line emission we need to subtract the continuum baseline
#	in uv space

# First we identify the line-free (i.e., continuum-only) channels in the data, and
#	averaging over time

CASA <X>: plotms(vis=CO32vis,
		xaxis='channel',
		yaxis='amp',
      		avgtime='1e8',
		avgscan=True,
		plotfile=' ',
		showgui=True)

# From this plot consider the following:

#	- what do you see?
#	- what channels should be excluded if considering continuum only?

# Let's exclude channels from   -> , no need to be precise with plenty of bandwidth

# Next we subtract the continuum using the casa task, uvcontsub:
#	https://casadocs.readthedocs.io/en/latest/api/tt/casatasks.manipulation.uvcontsub.html?casatasks.manipulation.uvcontsub

CASA <X>: uvcontsub(vis=CO32vis,
		fitorder=1,
		excludechans=True,
		fitspw=' ')
		
# This might take a wee while ... it will make a new measurement set with the suffix, ".contsub"

# Now we do the same for HCO+ ...

CASA <X>: plotms(vis=HCOp43vis,
		xaxis='channel',
		yaxis='amp',
      		avgtime='1e8',
		avgscan=True,
		plotfile=' ',
		showgui=True)

# From this plot consider the following:

#	- what channels should be excluded if considering continuum only?

# Let's exclude channels from    ->  , again no need to be precise with plenty of bandwidth

CASA <X>: uvcontsub(vis=HCOp43vis,
		fitorder=1,
		excludechans=True,
		fitspw=' ')
		
# Again, this might take a wee while ...

#################################################################################
# Imaging of molecular line data
#################################################################################

# Finally we make spectral cubes containing the molecular line emission

# We use the casa task, tclean, but now with the mode set to "cube"

# We will also have a few more parameters to set
# We will want to centre the cube on the source velocity of TW Hya =>  km/s

CASA <X>: velstart = ' '
CASA <X>: velwidth = ' '
CASA <X>: nvelchan =  
CASA <X>: freqCO32 = ' '
CASA <X>: freqHCOp43 = ' '

# Set threshold to ~ sigma (theoretical noise for the lines)
# We often need to clean lines a bit deeper than the continuum ...

CASA <X>: noise_line = ' '

# We also need to make a mask to clean the lines, we will use an already generated keplerian mask

# I have provided a python script that will create a FITS image from which we will make the mask
# Ouside of casa, open up keplerian_mask.py in a text editor, scroll to the bottom and
# 	edit the parameters so that they match those you have chosen for your image

# When you are satisfied you have a good keplerian mask, run the following in casa

# First import the fits image using the casa task, importfits:
#	https://casadocs.readthedocs.io/en/latest/api/tt/casatasks.data.importfits.html?casatasks.data.importfits

CASA <X>: importfits(fitsimage=' ',
		imagename=' ',
		defaultaxes=True,
		defaultaxesvalues=['','','','I'])

# Then use the image to make a mask for cleaning, using the casa task, makemask:
#	https://casadocs.readthedocs.io/en/latest/api/tt/casatasks.imaging.makemask.html?casatasks.imaging.makemask

CASA <X>: makemask(mode='copy',
		inpimage=' ',
		inpmask=' ',
		output=' ')

# WARNING - you may get some errors here, but the mask is still created!

# Now we image the CO 3-2 line, again running an non-interative clean with a mask and threshold

CASA <X>: tclean(vis=CO32vis+'.contsub',
		imagename=' ',
      		imsize=imsize_im,
		cell=cellsize_im,
      		specmode='cube',
		start=velstart,
		width=velwidth,
		nchan=nvelchan,
      		restfreq=freqCO32,
		outframe='LSRK',
      		weighting='briggs',
		robust=0.5,
      		interactive=False,
		mask=' ',
		threshold=noise_line,
		niter=1000000)

# Using casaviewer/imview, overlay a contour map of the residuals onto a raster image of the emission

# Examine the residuals, if there remains significant emission therein, try the following:

#	- clean to a deeper threshold (not too deep, no more than ~ sigma)
#	- use velocity smoothing to create a larger keplerian mask

# Using casaviewer/imview, inspect the final image and consider the following:
#
#	- how does the emission pattern change with velocity?
#	- is this what you would expect?

# Next let's estimate the rms noise and SNR for the CO 3-2 line
# Inspect the cude and determine the line-free channels


CASA <X>: estimate_lineSNR(imagename=' ',chan_range=' ')

# LIST IMAGE STATS HERE
 
# Now we image the HCO+ 4-3 line, again running an non-interative clean with a mask and threshold

# You will need to generate a new keplerian mask at the right frequency/velocity range

CASA <X>: importfits(fitsimage=' ',
		imagename=' ',
		defaultaxes=True,
		defaultaxesvalues=['','','','I'])

CASA <X>: makemask(mode='copy',
		inpimage=' ',
		inpmask=' ',
		output=' ')

CASA <X>: tclean(vis=HCOp43vis+'.contsub',
		imagename=' ',
      		imsize=imsize_im,
		cell=cellsize_im,
      		specmode='cube',
		start=velstart,
		width=velwidth,
		nchan=nvelchan,
      		restfreq=freqHCOp43,
		outframe='LSRK',
      		weighting='briggs',
		robust=0.5,
      		interactive=False,
		mask=' ',
		threshold=noise_line,
		niter=1000000)

# Next let's estimate the rms noise and SNR for the HCO+ 4-3 line using the line-free channels

CASA <X>: estimate_lineSNR(imagename=' ',chan_range=' ')

# LIST IMAGE STATS HERE

# For these image statistics, consider the following:
#
#	- how does the HCO+ line compare with the CO line?
#	- is this what you would expect?

#################################################################################
# END OF IMAGING EXERCISE
#################################################################################

# Now you have two spectral line cubes with which to do further analysis!

# For the analysis step we will use higher resolution data from the MAPS project

# Now switch to P1_analysis_script_student.py and the directory "Analysis_Step"

# Any concerns or questions at this point, contact a tutor

#################################################################################





















