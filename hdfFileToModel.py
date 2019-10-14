from numpy import linspace,sqrt, pi, arccos, zeros,searchsorted
from pyhdf import SD
from numpy import size
import numpy as nm
from re import findall
import re
from physics1 import *
import pprint


class HdfRadHydroModel:
  """ reads data from hdf file
        fields: r, rth, dd, ee, u1,u2,u3"""        
  def __init__(self):        
    self.r=None
    self.z=None
    self.th=None
    self.dd=None
    self.ee=None
    self.erad=None
    self.u1=None
    self.u2=None
    self.u3=None
    
    self.kpd=None   #dust opacity
    self.tau_e=None #optical depth
    self.Td=None    #dust temperature
    self.geff_z=None    #rad pressure -grav
    self.geff_r=None    #rad pressure -grav
    
    self.i_s = None
    self.ie = None    
    self.js = None
    self.je = None
    
    self.Nz=None
    self.Nx=None
    self.Ny=None
    self.Mbh=None
    self.F2Fedd=None
    self.Td_sblm=None
    self.kpa=None
    self.n0 =None   #number density
    self.R0 = None
    self.Dsc =None   #density
    self.tsc =None
    self.Eesc = None
    self.Rsc =None
    self.Usc =None
    self.zmax = None
    self.zmin = None
    self.dthdf = None

  def loadDataIn(self, hdfFileName,dat,radiation=False, printDetail = True):
    try:         
        dat.load_from_ZMP_CYL(hdfFileName, radiation)        
        if printDetail: print('loaded ', hdfFileName)        
    except Exception as e:
        print('failed to read', hdfFileName)
        print(str(e))                
  
  def load_from_ZMP_CYL(self, file_name, radiation=True, printDetail = True):
        """for ZEUS-2D hdf """   
        if printDetail: print(file_name)
        hdf = SD.SD(file_name)        
        sds=hdf.select("fakeDim1")    
        self.x=sds.get()
        
        sds=hdf.select("fakeDim2")
        self.z=sds.get()
        
        sds=hdf.select("Data-Set-2")    #u1
        x=(sds.get())
        self.u1 = x[0,:,:].T
        
        sds=hdf.select("Data-Set-3")    #u2
        x=(sds.get())
        self.u2 = x[0,:,:].T

        sds=hdf.select("Data-Set-4")    #u3
        x=(sds.get())
        self.u3 = x[0,:,:].T

        sds=hdf.select("Data-Set-5")    #density
        x=(sds.get())
        self.dd = x[0,:,:].T
   
        sds=hdf.select("Data-Set-6")     #gas energy
        x=(sds.get())
        self.ee = x[0,:,:].T
        
        self.Nz=size(self.z) # how the massives are declared??
        self.Nx=size(self.x)
        
        self.i_s = 1
        self.ie = self.Nz-1
        
        self.js = 1
        self.je = self.Nx-1
                        
        sds=hdf.select("Data-Set-7")     #rad energy
        x=(sds.get())
        self.erad = x[0,:,:].T
        
        if radiation:
            sds=hdf.select("Data-Set-8")     #dust opac
            x=(sds.get())
            self.kpd = x[0,:,:].T
    
            sds=hdf.select("Data-Set-9")     #tau_e
            x=(sds.get())        
            self.tau_e =x[0,:,:].T
    
    
            sds=hdf.select("Data-Set-10")     #Tdust
            x=(sds.get())
            self.Td = x[0,:,:].T
    
            sds=hdf.select("Data-Set-11")     #geff_z
            x=(sds.get())
            self.geff_z = x[0,:,:].T
    
            sds=hdf.select("Data-Set-12")     #geff_r
            x=(sds.get())
            self.geff_r = x[0,:,:].T
         
        
  def load_from_fileZ2D(self, file_name):        
    print( file_name )
    hdf = SD.SD(file_name)		
    sds=hdf.select("fakeDim0")	
    self.th=sds.get()
    
    sds=hdf.select("fakeDim1")
    self.r=sds.get()
                    
    sds=hdf.select("Data-Set-2")	#density
    self.dd = (sds.get()).T

    sds=hdf.select("Data-Set-3")	 #total energy
    self.ee = (sds.get()).T
    
    sds=hdf.select("Data-Set-4")	 #ur
    self.ur = (sds.get()).T
    
    sds=hdf.select("Data-Set-5")	 #uth
    self.uth = (sds.get()).T
    
    sds=hdf.select("Data-Set-6")	 #uph
    self.uph = (sds.get()).T
    
        

  def projectPolar2Cart(self, X,Y, f):
      Nr = self.r.size
      Nth = self.th.size
      Res = zeros((self.Nx,self.Ny))
      for i in range(self.Nx):
        for k in range(self.Ny):
          
          x=X[0,:]          
          y=Y[:,0]
                  
          ri = sqrt(x[i]**2 + y[k]**2)
          if ri==0: ri =1e-10
          
          thj = arccos(y[k]/ri)
                  
          i0= searchsorted(self.r, ri, side='left')
          j0= searchsorted(self.th, thj, side='right')
        
          if j0 == 0: t=0; j0=1
          if j0 == Nth: t=1; j0=j0-1
        
          if i0 == 0: u=0; i0=1
          if i0 == Nr: t=1; i0=i0-1
        
          t = ( thj - self.th[j0] ) / ( self.th[j0] - self.th[j0-1] )
          u = ( ri - self.r[i0] ) / (self.r[i0] - self.r[i0-1] )
        
          if j0 == 0: t=0; j0=1
          if j0 == Nth: t=1; j0=j0-1
        
          if i0 == 0: u=0; i0=1
          if i0 == Nr: t=1; i0=i0-1
            
          Res[i,k] = (1.-t) * (1.-u) * f[i0-1, j0-1] + \
             t * (1.-u) * f[i0, j0-1] + \
             t * u * f[i0,j0] + \
             (1.-t) * u * f[i0-1,j0]
    
      return(Res)

  def loadAtorusParam(self, fileName, printDetail = True):
                        
        try:
            atorData = open(fileName, 'r').read()
            splitData = atorData.split('\n')
                    
            for dataLine in splitData:                
                if 'MBH' in dataLine:            
                     self.Mbh = float( findall('MBH\s*=(\s*\d*.\d*e\d*)', dataLine)[0] )                     
                                   
                if 'nc0' in dataLine:                                               
                     self.n0 = float( findall('nc0\s*=(\s*\d*.\d*e\d*)', dataLine)[0] )                                      
                                    
                if 'R0' in dataLine:                                 
                     self.Rsc =  float(  findall('R0\s*=\s*(\d*.*\d)', dataLine)[0] )        
                     self.Rsc *= PARSEC         
                                                        
                if 'F2Fedd' in dataLine:         
                     self. F2Fedd = float(  findall('F2Fedd\s*=\s*(\d*.\d*)', dataLine)[0] )                 
                       
                if 'Td_sblm' in dataLine:         
                     self. Td_sblm = float(  findall('Td_sblm\s*=\s*(\d*.\d*)', dataLine)[0] )                 
                                                 
                if printDetail:
 #                  print('Mbh =', self.Mbh);
 #                  print( 'nc0 = ', self.n0)
 #                  print(' Rsc =',self.Rsc)      
 #                  print('F2Fedd=', self.F2Fedd)
                   print('Td_sblm=', self.Td_sblm)          
        except Exception as e:
          print(str(e))
        
#        self.Nz=size(self.z)
#        self.Nx=size(self.x)     
        
        self.Dsc = MP* self.n0
        self.tsc=sqrt( self.Rsc**3/G/self.Mbh/MSOL)
        self.Usc = sqrt(G*self.Mbh*MSOL/self.Rsc)
        mass = self.Mbh*MSOL
        self.Esc=self.Dsc*G*mass/self.Rsc
        self.Usc = sqrt( G*mass / self.Rsc )
        
        self.par= {'Dsc': self.Dsc, 'tsc':self.tsc,  'Usc':self.Usc, 'Mbh':self.Mbh, 'Esc':self.Esc, 'Usc' :self.Usc,
                    'nc0':self.n0,' Rsc': self.Rsc, 'F2Fedd': self.F2Fedd, 'Td_sblm': self.Td_sblm
                    }

#        print('Nx = ', self.Nx)
#         print('Nz = ', self.Nz)
#         print('Rsc = ', self.Rsc)
#         print('Usc = ', self.Usc)
#         print('Dsc = ', self.Dsc)
    
  def loadZmp_inp(self, fileName):
        try:
            atorData = open(fileName, 'r').read()
            splitData = atorData.split('\n')
            
            for dataLine in splitData:
                if '&iocon' in dataLine:
                    self. dthdf = float(  findall('dthdf\s*=\s*(\d*.\d*)', dataLine)[0] )
#                     print('dthdf=', self.dthdf)                    
                if '&ggen1' in dataLine:
                    self.zmin = float(  findall('x1min\s*=\s*(\d*.\d*)', dataLine)[0] )
                    self.zmax = float(  findall('x1max\s*=\s*(\d*.\d*)', dataLine)[0] )
        except e:
            print(str(e))
   
  def read_atorus_inp(self, file1):
        f1 = open(file1)    
        
        lines= f1.readlines()

        for line in lines:
            if 'MBH' in line:                                
                tmp=re.findall(r'\s*(MBH)\s*=\s*([\d\.\w\d*]+)',line)            
                MBH = float(tmp[0][1])
            if 'R0' in line:
                tmp=re.findall(r'(R0)\s*=\s*([\d* \.\w* \d*]+)',line)            
                R0 = float(tmp[0][1])
            if 'nc0' in line:
                tmp=re.findall(r'(nc0)\s*=\s*([\d\.\w\d*]+)',line)        
                n0 = float(tmp[0][1])
            if 'F2Fedd' in line:    
                tmp=re.findall(r'(F2Fedd)\s*=\s*([\d\.\w*\d*]+)',line)
                F2Fedd = float(tmp[0][1])
        
        f1.close()
        return(MBH,R0,n0,F2Fedd)  
     
  def velocAvr(self, mtor):
            
        i_eq = nm.argmin(abs(self.z-1.))
    #     print(i_eq, '  ', mtor);exit()    
        vAvr=0.
            
        for i in range(i_eq, self.Nz):        
            for j in range(1, self.Nx):
                
                r =  self.x[i]            
                dr = (self.x[i] - self.x[i-1])            
                dh = (self.z[i] - self.z[i-1])             
                dvol = 2*pi*r*dr*dh                        
                
                dvzAvr1 = self.u1[i, j]*self.dd[i,j] *dvol            
                dvzAvr1 *= self.Rsc**3* self.Dsc*self.Usc
                
                #print(dvzAvr, dvol, self.u1[i, j])
                
                vAvr += dvzAvr1
         
        res = vAvr/mtor
        
        print('v1Avr=', res)
        
        return(res)
        
  def torusMass(self):
    from math import pi
    
    y = zeros(self.Nz)                        
    for i in range(1, self.Nz):            
      y[i] = nm.trapz( 2*pi*self.x*self.dd[i,:], self.x )                  
      m=self.Rsc**3* self.Dsc * nm.trapz(y, self.z)        
      return (m)    
    
    
    
    
    
