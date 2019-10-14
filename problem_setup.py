#!/Users/dora/Library/Enthought/Canopy_32bit/User/bin/python2.7

import numpy as np
from scipy import *
from numpy import ndarray, zeros,array,size,meshgrid,flipud,floor,where,amin,argmin,int
import os
import hdfFileToModel
import math

from SaveModelToFile import SaveModelToFile, ReadModelFromFile,ConvertAnsSaveCylindricalToCart,ConvertAnsSaveCylindricalToCart_1 





def formatFileName(timeNum):
    time = str(timeNum)
    strLen= len(time)

    if strLen==1:
        return('00'+time)
    elif strLen==2:
        return('0'+time)
    else:
        return(time)

#from mpl_toolkits.mplot3d import Axes3D
#from pylab import *
#from matplotlib import pyplot as plt
#from mpl_toolkits.axes_grid1 import make_axes_locatable
#from mpl_toolkits.axes_grid1 import AxesGrid,ImageGrid
#from docutils.parsers.rst.directives.misc import Replace
# from apptools.appscripting.action import stop_recording_action
# from coverage import stop


put0= '/home/maria/DATA/'

locdir = 'runN300x300_L001fuv05fx05n1e8/'

put0 +=locdir

# creating a data model from datafile
dat=hdfFileToModel.HdfRadHydroModel()

dat.loadAtorusParam(put0 + '/atorus_inp', printDetail=True)

file_Name = "put0 + '/atorus_inp',"
#hdfaa000000
fileNamePrefix ="hdfaa000000."
dat.loadZmp_inp(put0 + 'zmp_inp')
dthdf=dat.dthdf
#print("dthdf=", dthdf)
 
timeNumeric =[2, 40, 60, 90] 
id0 = []
arrayTimeTeX = []

for time, i in zip(timeNumeric,xrange(len(timeNumeric))):    
    # timePhys = timeNumeric[i]*t_0/YR
    # arrayTimeTeX.append(roundThenStringToLatexFormat(timePhys))              
    # print('time in yrs=',  timePhys, arrayTimeTeX)
    timeStr= formatFileName(timeNumeric[i])
    id0.append(timeStr)
#    print (timeStr)
fileNamePrefix_1 = "dust_density"

hdfFileName_pure = fileNamePrefix+id0[0]
hdfFileName_pure_1 = fileNamePrefix_1

hdfFileName = put0 + hdfFileName_pure
hdfFileName1 = put0 + hdfFileName_pure_1

# copying data from file to a data model
try:
    dat.loadDataIn(hdfFileName, dat,radiation =True)
   


#exit()


except Exception as inst:
    print(Exception)
    





#file(1) = hdfFileName1+'.ascii'






Nt = 5  #change from 5 to 1
Nr = Nt
Rout = 5
eps = 1e-5

# r = linspace(0.1, Rout, Nr)
# t = linspace(0, 2.*pi - eps, Nt)

N_cart1 =dat.Nx                  #dimensions
N_cart2 =max(dat.Nz, dat.Ny)
N_cart3 = dat.Ny = max(dat.Nx, dat.Nz)

# Nx_cart = 10
# Ny_cart = 10
# Nz_cart = 10

x = zeros(N_cart1)    #axes

y = zeros(N_cart1)

z = zeros(N_cart3)

DD = zeros(shape=(N_cart1, N_cart2, N_cart3))   #unfilled density array
#print("DD=",DD)
x = linspace(-abs(dat.x[-1]), abs(dat.x[-1]), N_cart1)
#print("x=",x)
         #maximal radius#
y = x

z = x  #dat.z - 0.5*dat.z[-1] #allowed to change

# print (x)
# print (y)
# exit()
# print(z, 'z=')



# writing coordinates
file = hdfFileName+'.ascii'
file = os.getcwd() + '/' + hdfFileName_pure + '.inp'

ConvertAnsSaveCylindricalToCart(dat,  x,y,z, file)
#writing density
file = hdfFileName1+'.ascii'
file = os.getcwd() + '/' + hdfFileName_pure_1 + '.inp'
ConvertAnsSaveCylindricalToCart_1(dat,  DD, x,y,z, file)
print("==============================================")
print("               SUCCESS Part 1. :)             ")
print("==============================================")

exit()





