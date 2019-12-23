#
# Import NumPy for array handling
#
from trapezoidal import trapezoidal
import numpy as np
#import pandas as pd
#
# Import plotting libraries (start Python with ipython --matplotlib)
#
from mpl_toolkits.mplot3d import axes3d
from matplotlib import pyplot as plt
#
# Some natural constants
#
#Monte Carlo parameters
nphot    = 1000000

#
# Some natural constants
#
au  = 1.49598e13     # Astronomical Unit       [cm]
pc  = 3.08572e18     # Parsec                  [cm]
ms  = 1.98892e33     # Solar mass              [g]
ts  = 5.78e3         # Solar temperature       [K]
ls  = 3.8525e33      # Solar luminosity        [erg/s]
rs  = 6.96e10        # Solar radius            [cm]
cc  = 3.0e10         # speed of light
pi  = 3.1415
#
# Model parameters
#
radius   = 0.5*pc  #change
rho0     = 1e-10  #change from -16

#grid parameters for modelling; after adjusting should be removed

nx       = 100
ny       = 100
nz       = 100
sizex    = 1*pc
sizey    = 1*pc
sizez    = 1*pc




#################################################################
#


#
# Star parameters
#
mstar    = ms*1e+6  #change
rstar    = 0.001*pc #change from 0.001*pc
tstar    = ts*1   #change from ts
pstar    = np.array([0.,0.,0.])



# Write the wavelength_micron.inp file
#
lam1     = 0.00001e0
lam2     = 0.01e0
lam3     = 1.e1
lam4     = 10.e1
lam5     = 100.e1
n12      = 20
n23      = 40
n34      = 20
n45      = 20
lam12    = np.logspace(np.log10(lam1),np.log10(lam2),n12,endpoint=False)
lam23    = np.logspace(np.log10(lam2),np.log10(lam3),n23,endpoint=False)
lam34    = np.logspace(np.log10(lam3),np.log10(lam4),n34,endpoint=False)
lam45    = np.logspace(np.log10(lam4),np.log10(lam5),n45,endpoint=True)
lam      = np.concatenate([lam12,lam23,lam34,lam45])
nlam     = lam.size

#
# Write the wavelength file (in microns)
#
with open('wavelength_micron.inp','w+') as f:
    f.write('%d\n'%(nlam))
    for value in lam:
        f.write('%13.6e\n'%(value))

#
#
# Write the stars.inp file
#
F0 = 25.7*1e-14
with open('stars.inp','w+') as f:
    f.write('2\n')
    f.write('1 %d\n\n'%(nlam))
    f.write('%13.6e %13.6e %13.6e %13.6e %13.6e\n\n'%(rstar,mstar,pstar[0],pstar[1],pstar[2]))
    for value in lam:
        f.write('%13.6e\n'%(value))
    for value in lam:
        if value <= 0.01:
            F = F0*(value/0.01)**(1.2)
        f.write('%13.6e\n'%(F))
        if 0.01<value<0.1:  
          F = F0*1   
          f.write('%13.6e\n'%(F))
        if 0.1<value<1:
            F = F0*(value/0.1)**(-0.5)   
            f.write('%13.6e\n'%(F))
        if 1<value:
            F = F0*(0.1/1)**(0.5)*(value/1)**(-3)   
            f.write('%13.6e\n'%(F))


#Write the input luminosity
#d = 4*pi*(0.0001*pc)**2*1e7
#with open('inputSpectrum2.inp', 'w+') as f:
#      for value in lam:
#          if lam1<= value <=0.01:
#            v = lambda value: d*(value)**0.2*(1/0.01)**(1.2)
#            n = 20
#            I =trapezoidal(v, lam1, 0.01, n)   
#            f.write('%13.6e, %13.6e \n'%(value, I))
#          if 0.01< value <=0.1:
#            v = lambda value: d*(value)**(-1)
#            n = 20
#            I = trapezoidal(v, 0.01, 0.1, n)   
#          f.write('%13.6e, %13.6e \n'%(value, I))
#          if 0.1< value <=1:
#            v = lambda value: d*(value)**(-1.5)
#            n = 20
#            I = trapezoidal(v, 0.1, 1, n)   
#          f.write('%13.6e, %13.6e \n'%(value, I))
#          if 1<= value<=lam4:
#            v = lambda value: d*(0.1)**0.5*(value)**(-3)
#           n =20
#           I = trapezoidal(v, value, lam5, n)   
#          f.write('%13.6e, %13.6e \n'%(value, I))


# Write the radmc3d.inp control file
with open('radmc3d.inp','w+') as f:
    f.write('nphot = %d\n'%(nphot))
    f.write('scattering_mode_max = 0\n')   # Put this to 1 for isotropic scattering
    f.write('iranfreqmode = 1\n')
    f.write('modified_random_walk = 1\n')
    f.write('rto_style = 3')  # For the binary output
# Dust opacity control file
with open('dustopac.inp','w+') as f:
    f.write('2               Format number of this file\n')
    f.write('1               Nr of dust species\n')
    f.write('============================================================================\n')
    f.write('1               Way in which this dust species is read\n')
    f.write('0               0=Thermal grain\n')
    f.write('silicate        Extension of name of dustkappa_***.inp file\n')
    f.write('----------------------------------------------------------------------------\n')


print("==============================================")
print("               SUCCESS Part 2. :)             ")
print("==============================================")
