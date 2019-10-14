import numpy as np
import math
def SaveModelToFile(dat, fileToWrite, code = 'zeus'):
   
    if code=='zeus':
        f = open(fileToWrite, 'w')

        # f.write("Fortran-style indexing: F(i,j,k) \n:\
        # do z=>k =slowest, x=>i=fastest \n\n " )

        strToWrite ="nr:" +str(dat.Nx-1) + ' nz:'+str(dat.Nz-1)+"\n" # nr and nz are not declared
        f.write(strToWrite)

        print("writing x-coord")        
        for i in range(0,dat.Nx-1):
            f.write(str(dat.x[i])+"\n"+"\n")       
        
        print("writing z-coord")
        f.write("\n")
        f.write("k-coord: z \n")
        for k in range(0, dat.Ny-1):
            f.write(str(dat.z[k])+"\n"+"\n")       

        print("writing number density")
        f.write("\n")
        f.write("density\n")
        for k in range(0,dat.Nz-1): #dat.Nz-1
                for i in range(0,dat.Nx-1): #dat.Nx-1
                    strToWrite = str(dat.dd[i,k])+ "\n"+"\n"
                    f.write(strToWrite)       

        # print("writing Vr_cyl velocity [cm s^-1]")
        # f.write("\n")
        # f.write("Vr_cyl velocity [cm s^-1]\n")
        # for k in range(0,dat.nz-1):
        #     for j in range(0,dat.nt-1):
        #         for i in range(0,dat.nr-1):
        #             strToWrite = str(dat.Usc*Mx[i,j,k])+ "\n"
        #             f.write(strToWrite)   

        # print("writing Vt velocity [cm s^-1]")
        # f.write("\n")
        # f.write("Vt_cyl velocity [cm s^-1]\n")
        # for k in range(0,dat.nz-1):
        #     for j in range(0,dat.nt-1):
        #         for i in range(0,dat.nr-1):
        #             strToWrite = str(dat.Usc*Mt[i,j,k])+ "\n"
        #             f.write(strToWrite)    

        # print("writing Vz velocity [cm s^-1]")
        # f.write("\n")
        # f.write("Vr_cyl velocity [cm s^-1]\n")
        # for k in range(0,dat.nz-1):
        #     for j in range(0,dat.nt-1):
        #         for i in range(0,dat.nr-1):
        #             # strToWrite = str(i)+str(j)+str(k)+' '+str(dat.Usc*Mz[i,j,k])+ "\n"
        #             strToWrite = str(dat.Usc*Mz[i,j,k])+ "\n"
        #             f.write(strToWrite)    


        # print("SaveModelToFile ... done")
    
    
    f.close()



#  ====================================================================

    #  f = open(fileToWrite, 'w')

    #     f.write("Fortran-style indexing: F(i,j,k) \n:\
    #     do z=>k =slowest, phi=>j, x=>i=fastest \n\n " )

    #     strToWrite ="nr:" +str(dat.nr-1)+' nph:'+str(dat.nt-1)+' nz:'+str(dat.nz-1)+"\n"
    #     f.write(strToWrite)

    #     print("writing x-coord")
    #     f.write("i-coord: cyl. radius [cm]\n")
    #     for i in range(0,dat.nr):
    #         f.write(str(dat.Rsc*dat.r[i,0,0])+"\n")       

    #     print("writing phi-coord")
    #     f.write("\n")
    #     f.write("j-coord: theta [rad] \n")
    #     for j in range(0,dat.nt):
    #         f.write(str(dat.t[0,j,0])+"\n")       

    #     print("writing z-coord, [cm]")
    #     f.write("\n")
    #     f.write("k-coord: z \n")
    #     for k in range(0,dat.nz):
    #         f.write(str(dat.Rsc*dat.z[0,0,k])+"\n")       

    #     print("writing number density [cm^-3]")
    #     f.write("\n")
    #     f.write("number density [cm^-3] \n")
    #     for k in range(0,dat.nz-1):
    #         for j in range(0,dat.nt-1):
    #             for i in range(0,dat.nr-1):
    #                 strToWrite = str(dat.n0*dat.d[i,j,k])+ "\n"
    #                 f.write(strToWrite)       

    #     print("writing Vr_cyl velocity [cm s^-1]")
    #     f.write("\n")
    #     f.write("Vr_cyl velocity [cm s^-1]\n")
    #     for k in range(0,dat.nz-1):
    #         for j in range(0,dat.nt-1):
    #             for i in range(0,dat.nr-1):
    #                 strToWrite = str(dat.Usc*Mx[i,j,k])+ "\n"
    #                 f.write(strToWrite)   

    #     print("writing Vt velocity [cm s^-1]")
    #     f.write("\n")
    #     f.write("Vt_cyl velocity [cm s^-1]\n")
    #     for k in range(0,dat.nz-1):
    #         for j in range(0,dat.nt-1):
    #             for i in range(0,dat.nr-1):
    #                 strToWrite = str(dat.Usc*Mt[i,j,k])+ "\n"
    #                 f.write(strToWrite)    

    #     print("writing Vz velocity [cm s^-1]")
    #     f.write("\n")
    #     f.write("Vr_cyl velocity [cm s^-1]\n")
    #     for k in range(0,dat.nz-1):
    #         for j in range(0,dat.nt-1):
    #             for i in range(0,dat.nr-1):
    #                 # strToWrite = str(i)+str(j)+str(k)+' '+str(dat.Usc*Mz[i,j,k])+ "\n"
    #                 strToWrite = str(dat.Usc*Mz[i,j,k])+ "\n"
    #                 f.write(strToWrite)    


def ReadModelFromFile(dat, file, code = 'zeus'):
    from re import findall

    if code=='zeus':
        f = open(file, 'r').read()
        splitData = f.split('\n')
        for line in splitData:
            if 'nr' in line:
                Nx = int( findall('nr:(\d*)', line)[0])  #depends on nr??
                Nz = int( findall('nz:(\d*)', line)[0])  #depends on nz???
                print(Nx, Nz)



        # nr:299 nz:299

    # strToWrite ="nr:" +str(dat.Nx-1) + ' nz:'+str(dat.Nz-1)+"\n"
    # f.write(strToWrite)

def ConvertAnsSaveCylindricalToCart(dat, x,y,z, fileToWrite):
    
    f = open(fileToWrite, 'w')
    
    strToWrite = str(dat.Nx-1)+ " "  +str(dat.Ny-1)+ " "+ str(dat.Nz-1)+ '\n'
          
    
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

#  "writing y-coord"
    #f.write("\n")
    #f.write("j-coord: y \n")    
    for j in range(1, dat.Ny): #dat.Ny
        f.write(str(y[j]*3e+18)+"\n" )  #change from 0 to 2 
        
    
# "writing z-coord"
    #f.write("\n")
    #f.write( "k-coord: z \n" )    
    for k in range (1, dat.Nz):  #dat.Nz
        f.write(str(z[k]*3e+18)+"\n")  #change from 0 to 2 

        
def ConvertAnsSaveCylindricalToCart_1(dat,  DD, x,y,z,  fileToWrite):
    
    f = open(fileToWrite, 'w')
    strToWrite = str((dat.Nx-2)**3)
    f.write("1 \n")      #iformat
    f.write(strToWrite+"\n") #number of cells
    f.write("1 \n")      #nrspec - number of independent dust species

    for i,xi in zip(range(dat.Nx-1), x):
   
      for j,yj in zip(range(dat.Ny-1), y): 
                 
            for k,zk in zip(range(dat.Nz-1), z):
               
                rc = 1.0*np.sqrt(xi**2 + yj**2)                                 
                r_index = np.where(dat.x+2 >=rc)[0][0]            
                z_index = np.where(z >= zk)[0][0]               
                DD[i,j,k] = dat.dd[r_index, z_index]
                
                f.write(str(DD[i,j,k]*1.67266e-16)+"\n") #print into file
                

    f.close()
               
               




            




         


# for all x,y,z do 

#     rc2 = x**2 + y**2 + z**2    
#     rc = sqrt(rc2)
    
#     Rcyl = sqrt( x**2 + y**2 )    
#     phi = atan2(y,x)
#     phi= math.atan2(yj,xi) 
    
    
#     interpolate phi -> t find j_t
    
#     interpolate Rcyl -> r find i_r    
    
#     interpolate z -> z' find k_z
                
# end

# r = interp1d(Rcyl, r)

# itemindex = numpy.where(array==item)
