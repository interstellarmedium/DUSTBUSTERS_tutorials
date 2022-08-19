# Script to calculate a disk integrated rotation diagram using MCMC
# John Ilee, University of Leeds, April 2022
#
# Based on Ilee et al. 2021, ApJS, 257, 9
#

# Import the required packages (pip install if you don't have these)
import matplotlib
import sys
import emcee
import corner
import numpy as np
import pylab as pl
from multiprocessing import Pool
from scipy.constants import c as c
from scipy.constants import h as planck
from scipy.constants import k as boltz
from matplotlib.font_manager import FontProperties
from scipy import interpolate

# This file contains lots of the molecular data needed for the calculations
import moleculardictionary as md

# Set font size and type
font = {'size'   : 16}
pl.rc('font', **font)
pl.rc('text', usetex=False)
fontsize = 16

#################################################################
# Define some functions for later

def error_gradient(E_u, y):
    percentiles = np.percentile(y, np.arange(101), axis=0)

    for i in np.arange(49):
        pl.fill_between(E_u, percentiles[i], percentiles[-(i+1)], facecolor='lightblue', edgecolor='none', alpha=0.04*i**0.05)

    return

# Interpolate partition function from CDMS entries
def Q_rot(T):
    # This will take the Q information and linearly interpolate in log space
    Ti = np.log10(md.molec_dict['Band6'][molecule]['T_Q'])  # Band here doesn't matter
    Qi = md.molec_dict['Band6'][molecule]['Q']
    fQi = interpolate.interp1d(Ti,Qi)
    return 10**fQi(np.log10(T))

def N_u_div_g_u_calc(E_u, N_t, T_r, sym):
    ln_N_u = np.log(N_t/Q_rot(T_r))-E_u/T_r
    return np.exp(ln_N_u)

def tau_calc(E_u, N_t, T_r, A_ul, nu_ul, g_u, sym, lw):
    N_u =  N_u_div_g_u_calc(E_u, N_t, T_r, sym)*g_u
    tau = ((A_ul*c**3)/(8.*np.pi*nu_ul**3*lw))*N_u*1e4*(np.exp((planck*nu_ul)/(boltz*T_r))-1)
    return tau

# Set priors 
def lnprior(theta):
    N_t, T_r = theta

    # Make these priors are reasonable for your source
    if 1e6 < N_t < 1e18 and 20. < T_r < 150.:
        return 0.
    return -np.inf

def lnlike(theta, y, yerr, A_ul, E_u, nu_ul, g_u, sym, lw):
    y = y/g_u
    yerr = yerr/g_u
    N_t, T_r = theta
    tau = tau_calc(E_u, N_t, T_r, A_ul, nu_ul, g_u, sym, lw)

    # if tau is less than 1e-4, treat as optically thin
    if np.any(tau<0.0001):
        C_t = np.ones(tau.shape)
    else:
        C_t = tau/(1-np.exp(-tau))
    #print('C_t in ln like = ' + str(C_t))

    y = y*C_t
    model = np.exp(np.log(N_t/Q_rot(T_r))-E_u/T_r)
    inv_sigma2 = 1/(yerr**2)
    return -0.5*(np.sum((y-model)**2*inv_sigma2 - np.log(inv_sigma2)))

def lnprob(theta, y, yerr, A_ul, E_u, nu_ul, g_u, sym, lw):
    lp = lnprior(theta)
    if not np.isfinite(lp):
        return -np.inf
    #added to by-pass lnprob returned NaN.
    if np.isnan(lp + lnlike(theta, y, yerr, A_ul, E_u, nu_ul, g_u, sym, lw)):
        print('ValueError (by-passed) : lnprob returned NaN.')
        return -np.inf
    #print('Ln prior = ' + str(lp))
    #print('Ln likelyhood = ' + str(lnlike(theta, y, yerr, A_ul, E_u, nu_ul, g_u, sym, lw)))
    #print('')
    return lp + lnlike(theta, y, yerr, A_ul, E_u, nu_ul, g_u, sym, lw)

# Choose number of walkers
def fit_rot_temp(ys, yerr, init, A_ul, E_u, nu_ul, g_u, sym, lw, nruns):
    ndim, nwalkers = 2, 150
    initial = init
    pos = [initial + 1e-4*np.array(initial)*np.random.randn(ndim) for i in range(nwalkers)]

    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(ys, yerr, A_ul, E_u, nu_ul, g_u, sym, lw), pool=Pool())
    # ~ sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(ys, yerr, A_ul, E_u, nu_ul, g_u, sym, lw), threads=6)

    for i in range(nruns):
        sampler.run_mcmc(pos, 1500)
        #file_name = "multi_gauss_samples_chain_" + str(i) + ".npy"
        #np.save(file_name, np.array(sampler.chain))
        pos = sampler.chain[:,-1,:]    
    return sampler.chain


def rot_diag_emcee(fluxes, flux_uncs, int_area, init, A_ul, E_u, nu_ul, g_u, sym, lw, Emin, Emax, Nmin, Nmax, axis_loc, colorisomer, radius):
    #ax = pl.axes([axis_loc[0], axis_loc[1], axis_loc[2], axis_loc[3]])

    ax = pl.axes()

    #~ N_u = fluxes/1000.*1e-23*4*np.pi/(planck*1e4)/c/A_ul/int_area
    #~ N_uncs = flux_uncs/1000.*1e-23*4*np.pi/(planck*1e4)/c/A_ul/int_area
    
    N_u = fluxes/1000.*1e-23*4*np.pi/(planck*1e4)/c/A_ul/int_area
    N_uncs = (fluxes + flux_uncs)/1000.*1e-23*4*np.pi/(planck*1e4)/c/A_ul/int_area - N_u

    chain = fit_rot_temp(N_u, N_uncs, init, A_ul, E_u, nu_ul, g_u, sym, lw, 1)
        
    np.save('rot_diagram_disk_avg_chains.npy', chain)

#     np.save(savefolder + "rot_diagram_" + str(sym) + '_' + str(radius) + "_chains.npy", chain)


    samples = chain[:,1000:,:].reshape((-1,2))

    ax.set_yscale('log')

    y = []
    for N_t, T_r in samples:
        y.append(N_u_div_g_u_calc(np.arange(Emax+1), N_t, T_r, sym))
    y = np.array(y)
    error_gradient(np.arange(Emax+1), y)

    #for N_t, T_r in samples[np.random.randint(len(samples), size=500)]:
    #    pl.plot(np.arange(Emax+1), N_u_div_g_u_calc(np.arange(Emax+1), N_t, T_r, sym), color='darkslateblue', alpha=0.01)

    N_t = np.percentile(samples[:,0], 50)
    T_r = np.percentile(samples[:,1], 50)
    T_err_p = np.percentile(samples[:,1], 84) - np.percentile(samples[:,1], 50)
    T_err_n = np.percentile(samples[:,1], 50) - np.percentile(samples[:,1], 16)
    N_err_p = np.percentile(samples[:,0], 84) - np.percentile(samples[:,0], 50)
    N_err_n = np.percentile(samples[:,0], 50) - np.percentile(samples[:,0], 16)
    print()
    print("Column density = " + '{:.2e}'.format(N_t) + " +/- " + '{:.2e}'.format(N_err_p) + ", " + '{:.2e}'.format(N_err_n) + " cm^-2")
    print("Rotational temperature = " + '{:.1f}'.format(T_r) + " +/- " + '{:.1f}'.format(T_err_p) + ", " + '{:.1f}'.format(T_err_n) + " K")
   
    tau = tau_calc(E_u, N_t, T_r, A_ul, nu_ul, g_u, sym, lw)
    #print('Tau values')
    print("Optical depths = ", tau)
    print()

    # if tau is less than 1e-4, treat as optically thin
    if np.any(tau<0.0001):
        C_t = np.ones(tau.shape)
    else:
        C_t = tau/(1-np.exp(-tau))
        #continue
        

    #~ print('N_u, C_t values')
    #~ print(N_u, C_t)

    # could have multiple types of points / colors if there are multiple k ladders
    img = pl.errorbar(E_u, (N_u*C_t)/g_u, yerr=(N_uncs*C_t)/g_u, fmt='o', markersize=8, markeredgecolor='none', markerfacecolor='k', ecolor='k', capsize=5, elinewidth=1, zorder=10000)

    
    pl.xlim([Emin, Emax])
    pl.ylim([Nmin, Nmax])

    return ax, img, [N_t, N_err_p, N_err_n, T_r, T_err_p, T_err_n]


def annotate(ipan, ax, vals_o):
    N_t_o, N_err_p_o, N_err_n_o, T_r_o, T_err_p_o, T_err_n_o = vals_o

    ax.minorticks_on()

    pl.ylabel(r"N$_{\rm u}$ / g$_{\rm u}$ (cm$^{{\rm -2}}$)")

    pl.xlabel(r"E$_{\rm u}$ (K)")

    pl.title(obj_lab)
    
    ax.text(0.95, 0.95, mol_lab, ha='right', va='top', transform = ax.transAxes, size=fontsize)

    oom_o = int(np.log10(N_t_o))
#     oom_p = int(np.log10(N_t_p))

#     ax.text(0.05,0.375,r'N$_{{\rm T}}$~=~{:.2f}$_{{-{:.2f}}}^{{+{:.2f}}}$$\times$10$^{{{:d}}}$~cm$^{{-2}}$'.format(N_t_o/(10**oom_o), N_err_n_o/(10**oom_o), N_err_p_o/(10**oom_o), oom_o), transform = ax.transAxes, size='7')
#     ax.text(0.05,0.175,r'T$_{{\rm rot}}$~=~{:.1f}$_{{-{:.1f}}}^{{+{:.1f}}}$~K'.format(T_r_o, T_err_n_o, T_err_p_o), transform = ax.transAxes, size='7')

    ax.text(0.05,0.25,r'N$_{{\rm T}}$  = {:.1f}$_{{-{:.1f}}}^{{+{:.1f}}}$$ \times$10$^{{{:d}}}$ cm$^{{\rm -2}}$'.format(N_t_o/(10**oom_o), N_err_n_o/(10**oom_o), N_err_p_o/(10**oom_o), oom_o), transform = ax.transAxes, size=fontsize)
    ax.text(0.05,0.1,r'T$_{{\rm rot}}$ = {:.1f}$_{{-{:.1f}}}^{{+{:.1f}}}$ K'.format(T_r_o, T_err_n_o, T_err_p_o), transform = ax.transAxes, size=fontsize)

    ax.tick_params(axis='both', which='major', pad=4)
    
    pl.setp(ax.get_xticklabels(), size=fontsize)
    pl.setp(ax.get_yticklabels(), size=fontsize)
    return()

 
def cornerplot(chain, corner_TF, ndim, sym, radius):
    # Create walker path and corner plots
    if corner_TF == True:
        # split the sampler chain into the different parameters 
        # and remove the first 1000 steps as a burn-in 
        samples = chain[:,1000:,:].reshape((-1,ndim))

        # plot the path of each walker
        # can be used to see if the distribution is multimodal and whether the walkers have converged
        fig, axes = pl.plt.subplots(ndim, 1, figsize = (9, ndim*2+1))
        labels = ['N$_T$', 'T$_{rot}$']
        if ndim == 2:
            labels += [r'$\Delta$V$_{\rm int}$']

        for walker in chain:
            for i,(ax,lab) in enumerate(zip(axes,labels)):
                ax.plot(walker[:,i])
                ax.set_xlabel('step number',fontsize=14)
                ax.set_ylabel(lab,fontsize=14)

        pl.tight_layout()
        fig.savefig(savefolder + 'walkers_' + sym + '_' + str(radius) + '.png')
        pl.close()

        # Making the "famous" corner plot
        # NOTE: the corner plots shows the 0.16, 0.50 and 0.84 quantile in both the title and the figures itself      
        fig = corner.corner(samples, labels=labels, fontsize=12, show_titles=True, quantiles = [0.16,0.5,0.84])
        fig.savefig(savefolder + 'corner_' + sym + '_' + str(radius) + '.png')
        pl.close()


# Calculate a linewidth based on a temperature (Pegues+ 2020)
def tot_lw(Tex):
    amu = 1.6605390666e-27 # kg   
    mu = 2.37
    t0 = 0.01
    mH = 1.00784*amu
    mX = md.molec_dict['Band6'][molecule]['m']*amu
                           
    first  = (np.sqrt((2.0*boltz*Tex)/mX))**2
    second = (t0*np.sqrt((boltz*Tex)/(mu*mH)))**2

    dV = np.sqrt(first+second) * 2.0*np.sqrt(np.log(2))

    return dV
 
#################################################################
# Main Program
#################################################################

# What object are you looking at?
my_obj = 'HD_163296'
obj_lab = 'HD 163296' # no underscores for labels

# Assumed distance to object
distance = 101 # parsecs

# To what radial extent did you measure fluxes?
mask_radius = 500.0     # au

# Which molecule are you looking at?
molecule = 'HC3N'                          
mol_lab = 'HC$_{\mathregular{3}}$N'  # latex'd for labels

# Set the bounds for the plot (make sure the priors above are consistent)
Emin, Emax = 0, 200
Nmin, Nmax = 1e8, 1e12

print()
print("Rotation diagram analysis for ", molecule, " in ", obj_lab)
print()

# Set a linewidth (m/s) based on a temperature.  Iterate this with Trot after first run.
lw = tot_lw(30.0)       # 30 K

print('Assuming an intrisic linewidth of ', lw, ' m/s') 
print()

# Disk integrated flux measurements
#####################################

# Integrated flux densities (mJy km/s)

fluxes_o    = np.array([171, 187])   # order must match molecular data for other arrays
noises      = np.array([5.8, 3.1])   # measured uncertainties on above

# Observational parameters
###########################

# Beam major axis and minor axis (need this for each transition)

bmaj = 0.3      # arcsec
bmin = 0.3
beam_area_1 = (bmaj*bmin*np.pi)/(4.*np.log(2.))

bmaj = 0.3      # arcsec
bmin = 0.3
beam_area_2 = (bmaj*bmin*np.pi)/(4.*np.log(2.))

beam_areas_o    = np.array([beam_area_1,beam_area_2])


# Molecular data parameters (grabbed from molecular dictionary)
################################################################

# Symmetry parameter for molecule
sym = 1

# Frequencies of the transitions (Hz)
nu_ul_o  = (10**9)*np.array([md.molec_dict['Band3'][molecule]['freq'],
                             md.molec_dict['Band6'][molecule]['freq']])

# Einstein A co-efficients of the transitions ()
A_ul_o   = 10**np.array([md.molec_dict['Band3'][molecule]['logAul'],
                         md.molec_dict['Band6'][molecule]['logAul']])
                         
# Upper energy of the transitions (K)
E_u_o    = np.array([md.molec_dict['Band3'][molecule]['Eu'],
                     md.molec_dict['Band6'][molecule]['Eu']])

# Upper state degeneracies of the transitions
g_u_o    = np.array([md.molec_dict['Band3'][molecule]['gu'],
                     md.molec_dict['Band6'][molecule]['gu']])


# Fitting parameters
#####################

# Initial guess of column density (cm^-2) and rotational temp (K)
init = [4.e13, 20.]

# Some switches
corner_TF = True
crash = False
mode = 'disk'

# Specify a save folder if you want
savefolder = ''

######################################
# Disk-averaged calculation and plot #
######################################

if mode == 'disk':
   
    fig = pl.figure(figsize=(5,4))

    # Calculate integrated area
    int_area = np.pi*((mask_radius/distance)/206265)**2
    
    # Fold in a 10% inter-band flux calibration uncertainty
    flux_uncs_o = np.sqrt((noises**2) + (0.1*fluxes_o)**2)

    # Make the rotation diagram
    ax, img_o, vals_o = rot_diag_emcee(fluxes_o, flux_uncs_o, int_area, init, A_ul_o, E_u_o, nu_ul_o, g_u_o, 'ortho', lw, Emin, Emax, Nmin, Nmax, [0.145, 0.14, 0.83, 0.8], 'darkorange', 'disk_avg')
    prop=FontProperties(size=10)
    annotate(1, ax, vals_o)
    
    pl.tight_layout()

    pl.savefig('Rot_diag_'+my_obj+'_'+molecule+'.pdf')

    pl.show()
    
    # Finally, make the corner plot
    cornerplot(np.load("rot_diagram_disk_avg_chains.npy"), corner_TF, 2, 'ortho', 'disk-avg')
    
    print("Done.")
