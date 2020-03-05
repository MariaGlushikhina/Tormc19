import numpy as np
import math
import pickle 
import struct
def SaveModelToFile(dat, fileToWrite, code = 'zeus'):
   
    if code=='zeus':
        f = open(fileToWrite, 'w')

        # f.write("Fortran-style indexing: F(i,j,k) \n:\
        # do z=>k =slowest, x=>i=fastest \n\n " )

        strToWrite ="nr:" +str(dat.Nx-1) + ' nz:'+str(dat.Nz-1)+"\n" # nr and nz are not declared
        f.write(strToWrite)

        print("writing x-coord")        
        for i in range(0,dat.Nx): #change from -1
            f.write(str(dat.x[i])+"\n"+"\n")       
        
        print("writing z-coord")
        f.write("\n")
        f.write("k-coord: z \n")
        for k in range(0, dat.Ny):
            f.write(str(dat.z[k])+"\n"+"\n")       

        print("writing number density")
        f.write("\n")
        f.write("density\n")
        for k in range(0,dat.Nz): #dat.Nz-1
                for i in range(0,dat.Nx): #dat.Nx-1
                    strToWrite = str(dat.dd[i,k])+ "\n"+"\n"
                    f.write(strToWrite)       

       
    
    f.close()



#  ====================================================================

  

def ReadModelFromFile(dat, file, code = 'zeus'):
    from re import findall

    if code=='zeus':
        f = open(file, 'r').read()
        splitData = f.split('\n')
        for line in splitData:
            if 'nr' in line:
                Nx = int( findall('nr:(\d*)', line)[0])  #depends on nr??
                Nz = int( findall('nz:(\d*)', line)[0])  #depends on nz???
               #print(Nx, Nz)

def ConvertAnsSaveCylindricalToCart(dat, x,y,z, fileToWrite): #writing coordinates file
    
    f = open(fileToWrite, 'w')
    
    strToWrite = str(dat.Nx-2)+ " "  +str(dat.Ny-2)+ " "+ str(dat.Nz-2)+ '\n' #change from 2
          
    
    f.write("1 \n")  #iformat
    f.write("0 \n")  # grid style regular = 0
    f.write("99 \n")  #coord system cartesian
    f.write("0 \n")   # gridinfo
    f.write("1  1  1 \n") #number of active coordinate axis, inactive axis = 0
    f.write(strToWrite)  # number of grid cells for each dimension
# print "writing x-coord"
    #f.write("i-coord: x \n\n")
    for i in range(1, dat.Nx): #dat.Nx
        f.write(str(x[i]*3e+18)+"\n")   #change from 0 to 2    
       #print("xi", x[i])
#  "writing y-coord"
    #f.write("\n")
    #f.write("j-coord: y \n")    
    for j in range(1, dat.Ny): #dat.Ny

        f.write(str(y[j]*3e+18)+"\n" )  #change from 0 to 2 #3e+18 - r_scale make a change
        
    
# "writing z-coord"
    #f.write("\n")
    #f.write( "k-coord: z \n" )    
    for k in range (1, dat.Nz):  
        if k< dat.Nz/4: 
         f.write(str(z[k]*3e18)+"\n")  
        if  dat.Nz/4<= k <= 3*dat.Nz/4:
         f.write(str(0*z[k]*3e18)+"\n")
        if k> 3*dat.Nz/4:
         f.write(str(z[k]*3e18)+"\n")
        
def ConvertAnsSaveCylindricalToCart_1(dat,  DD, x,y,z,  fileToWrite):# writing density file
    
    f = open(fileToWrite, 'w')
    strToWrite = str((dat.Nx-2)*(dat.Ny-2)*(dat.Nz-2)) #change from 2
    f.write("1 \n")      #iformat
    f.write(strToWrite+"\n") #number of cells
    f.write("1 \n")      #nrspec - number of independent dust species
    print 'dat.x', (dat.x+1) 
    for i,xi in zip(range(dat.Nx-2), x):  #change from 2
   
      for j,yj in zip(range(dat.Ny-2), y): 
   
         for k,zk in zip(range(dat.Nz-2), z):
               
                rc = np.sqrt(1.*(xi)**2 + 1.*(yj)**2)  # 
               # print rc, 'rc'                          
                r_index = np.where((dat.x+1.1) >=0.999*rc)[0][0]  # no more than 2
               # print 'r_index',r_index        
                z_index = np.where(1.0*z >= zk)[0][0] #change from z>=zk 
                #print 'z_index', z_index          
                DD[i,j,k] = dat.dd[r_index, z_index]
                f.write(str(DD[i,j,k]*1.67266e-16)+"\n") #print into file
                

    f.close()
               
               





            




         



