# To use as a mask in CASA, first set imcen to RA and Dec, in decimal degrees, 
#    of the phase center of the image and then execute the following commands within CASA

# importfits(fitsimage='mask.fits',imagename='mask.im',defaultaxes=True,defaultaxesvalues=['','','','I'])
# makemask(mode='copy',inpimage='mask.im',inpmask='mask.im',output='mask.im:mask0')

# This should then be used in clean with the mask keyword (mask='mask.im') in theory, 
#      but in practice it doesn't always seem to work...

def butterfly(mstar,Rin,Rout,incl,PA,dist,dv,npix,imres,offs,beam,nchans,chanstep,voff,outfile,freq0,imcen):
    
    import numpy as np
    from astropy import constants as const
    from astropy.io import fits
    from scipy import ndimage
    from astropy.convolution import Gaussian2DKernel,convolve

    AU = const.au.cgs.value
    pc = const.pc.cgs.value
    Msun = const.M_sun.cgs.value
    G = const.G.cgs.value
    incl = np.radians(incl) #convert inclination from degrees to radians
    dist *= pc
    Rin *= AU
    Rout *= AU
    mstar *= Msun

# Create image
    ra_image = np.linspace(-npix/2*imres,npix/2*imres,npix)
    dec_image = np.linspace(-npix/2*imres,npix/2*imres,npix)
    ra_image,dec_image = np.meshgrid(ra_image,dec_image)

# Define radius and phi for each pixel
    radius = np.sqrt((ra_image)**2+((dec_image)/np.cos(incl))**2)
    radius = np.radians(radius/3600.)*dist
    phi = np.arctan2((dec_image),(ra_image))
    
    vkep = np.sqrt(G*mstar/radius)/1e5 
    vlos = vkep*np.cos(phi)*np.sin(incl)
    dv = dv*(radius/(100*AU))**(-.5)
    
# Remove points outside of disk
    wremove = (radius<Rin) | (radius>Rout)
    vlos[wremove] = -100
    vmax = np.max(vlos[vlos>-90])
    vmin = np.min(vlos[vlos>-90])

    chanmin = (-nchans/2.+.5)*chanstep
    velo_steps = np.arange(nchans)*chanstep+chanmin
    image = np.zeros((npix,npix,nchans))
        
    for i in range(nchans):
        w = (vlos+dv/2.>velo_steps[i]-chanstep/2.) & (vlos-dv/2.<(velo_steps[i]+chanstep/2.))
        image[:,:,i][w] = 1

# Convolve with finite spatial resolution
        gauss = Gaussian2DKernel(beam/(imres*2*np.sqrt(2*np.log(2))))
    for i in range(nchans):
        image[:,:,i] = convolve(image[:,:,i],gauss,boundary=None)
    
# Shift and rotate image to correct orientation
    image = ndimage.rotate(image,PA+180,reshape=False)
    pixshift = np.array([-1.,1.])*offs/(np.array([imres,imres]))
    image = ndimage.shift(image,(pixshift[0],pixshift[1],0),mode='nearest')
    image[image<.01]=0. #clean up artifacts from interpolation 
    image[image>.01]=1.

# Make fits file
    hdr = fits.Header()
    hdr['TELESCOP'] = 'ALMA'
    hdr['SIMPLE'] = 'T'
    hdr['BITPIX'] = 32
    hdr['NAXIS'] = 4
    hdr['NAXIS1'] = npix
    hdr['NAXIS2'] = npix
    hdr['NAXIS3'] = nchans
    hdr['NAXIS4'] = 1
    hdr['CDELT1'] = -1*imres/3600.
    hdr['CRPIX1'] = npix/2.+1
    hdr['CRVAL1'] = imcen[0]
    hdr['CTYPE1'] = 'RA---SIN'
    hdr['CUNIT1'] = 'deg'
    hdr['CDELT2'] = imres/3600.
    hdr['CRPIX2'] = npix/2.+1
    hdr['CRVAL2'] = imcen[1]
    hdr['CTYPE2'] = 'DEC--SIN'
    hdr['CUNIT2'] = 'deg'
    hdr['CTYPE3'] = 'VELO-LSR'
    hdr['CDELT3'] = chanstep*1e3
    hdr['CRPIX3'] = 1.
    hdr['CRVAL3'] = velo_steps[0]*1.0e3 + voff*1.0e3
    hdr['CTYPE4'] = 'STOKES'
    hdr['CDELT4'] = 1.
    hdr['CRPIX4'] = 1.
    hdr['CRVAL4'] = 1.
    hdr['RESTFREQ'] = freq0*1.0e9
    hdu = fits.PrimaryHDU(image.T,hdr)
    hdu.writeto(outfile,overwrite=True,output_verify='fix')

# Call the butterfly function

butterfly(mstar=0.8,
	Rin=4.0,
	Rout=400.0,
	incl=7.0,
	PA=-20.0,
	dist=54.0,
	dv=1.0,
	npix=64,
	imres=0.3,
	offs=[0.0,0.0],
	beam=1.6,
	nchans=99,
	chanstep=0.12,
	voff=2.82,
	outfile='HCOp_4-3_keplerian_mask.fits',
	freq0=356.734288,
	imcen=[165.4660208,-34.7047667]) 
	
	

