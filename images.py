#!/home/maria/anaconda2/bin/python

import matplotlib as plb
import numpy as np
import pandas as pd
import os
from mpl_toolkits.mplot3d import axes3d
from matplotlib import pyplot as plt
from matplotlib import cm
from radmc3dPy.image import *    # Make sure that the shell variable PYTHONPATH points to the RADMC-3D python directory
from radmc3dPy.analyze import *  # Make sure that the shell variable PYTHONPATH points to the RADMC-3D python directory
from radmc3dPy import *
pc  = 3.08572e18
#dust density contours
fig1=plb.figure()
q = analyze.readGrid(sgrid = True)
xx   = q.x[:]
yy   = q.y[:]
zz   = q.z[:]
#qq = np.meshgrid(xx, yy, zz, indexing='ij')
#xc = qq[0]
#yc = qq[1]
#zc = qq[2]
#x1 = xc[:,:,0]
#y1 = yc[:,:,0]
#z1 = zc[:,:,0]
#z2 = np.pi/2 - y1 
data = analyze.readData(ddens=True, binary=False)
v = data.rhodust[0,:,:,0]
c = plb.contourf(zz/natconst.pc,  xx/natconst.pc, (data.rhodust[0,:,:,0]), 290)
#plb.xlabel('r [au]')
#plb.axis([-500000, 500000, -250000, 250000])
#plb.ylabel(r' [au]')
#plb.xscale('log')
#plb.yscale('log')
cb = plb.colorbar(c)
cb.set_label(r'$\log_{10}{\rho}$', rotation=270.)


# View a 2-D slice of the 3-D array of the setup
#
#fig2 = plt.figure()
#q = analyze.readGrid(sgrid = True)
#v = analyze.readData(ddens=True,  binary=False)
#ax   = fig2.add_subplot(111, projection='3d')
#ax.plot_wireframe((y1), (x1), (v), rstride=1, cstride=1)
#surf = ax.plot_surface((y1),  (x1) ,  (v), cmap=cm.coolwarm,linewidth=0, antialiased=False)# change from data1 to zz
#fig2.colorbar(surf, shrink=0.5, aspect=5)
# Make and plot an example image
#
#fig2  = plt.figure()
#makeImage(npix=200,incl=60.,phi=30.,wav=30.,sizeau=300)   # This calls radmc3d 
#a=readImage()
#plotImage(a,log=True,au=True,maxlog=6,cmap='hot')


# dust temperature contours
#opac = analyze.readOpac(ext=['silicate'])
#fig2=plb.figure()
#data = analyze.readData(dtemp=True)
#c= plb.contourf(data.grid.x/natconst.pc, data.grid.y/natconst.pc, (data.dusttemp[:,:,0,0].T), 298)
#plb.xlabel('r [pc]')
#plb.ylabel(r'$\pi$')
#plb.xscale('log')
#plb.yscale('log')
#plb.axis([-2.5, 2.5, -2, 4])
#cb = plb.colorbar(c)
#cb.set_label('T [K]', rotation=270.)
#c= plb.contour(data.grid.x/natconst.pc, np.pi/2.-data.grid.y/natconst.pc, (data.dusttemp[:,:,0,0].T), 100,  colors='k', linestyles='solid')
#plb.clabel(c, inline=1, fontsize=10)


#Plotting it "by hand", the SED as seen at 1 pc distance
#
#fig2  = plt.figure()
#s     = readSpectrum()
#lam   = s[:,0]
#nu    = 1e4*cc/lam
#fnu   = s[:,1]
#nufnu = nu*fnu
#plt.plot(lam,nufnu)
#plt.xscale('log')
#plt.yscale('log')
#plt.axis([1e-3, 1e6, 1e+1, 5e17])
#plt.xlabel('$\lambda\; [\mu \mathrm{m}$]')
#plt.ylabel('$\\nu F_\\nu \; [\mathrm{erg}\,\mathrm{cm}^{-2}\,\mathrm{s}^{-1}]$')
# Use the radmc3dPy.analyze tool set for plotting the SED,
# this time let's plot nuLnu in units of Lsun
#
#fig3  = plt.figure()
#plotSpectrum(s,nulnu=True,lsun=True,xlg=True,ylg=False,micron=True)
#plotSpectrum(s,nulnu=True,lsun=True,xlg=True,ylg=True,kev=True)
#plt.axis([1e-8,1e+2,1e+6,3.0e16])

#plot of GIVEN spectrum
#fig4  = plt.figure()
#table = pd.read_csv('inputSpectrum.inp')
#x = table.values[:,0]
#y = table.values[:,1]
#plt.figure(figsize=(15, 7))
#plt.plot(x,y)
#plt.xscale('log')
#plt.yscale('log')
#plt.axis([1e-3, 1e6, 1e+39, 5e+45])
#plt.xlabel('$\lambda')
#plt.ylabel('$\\ \lambda L_\\lambda')

#plot L_lambda dlambda/ L_bol
#fig5  = plt.figure()
#table1 = pd.read_csv('inputSpectrum1.inp')
#x = table1.values[:,0]
#y = table1.values[:,1]
#L_l = y*4*3.1415*(0.001*pc)**2
#L_bol = 10**42
#delta_L = L_l/L_bol
#plt.plot(x,y)
#plt.xscale('log')
#plt.yscale('linear')
#plt.axis([1e-4, 1e6, 9.0e+34, 15.1e+35])
#plt.xlabel('$\lambda')
#plt.ylabel('$\\ L_\\nu')
plb.show()
