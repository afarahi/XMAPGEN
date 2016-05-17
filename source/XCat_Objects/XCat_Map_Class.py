from XCat_Objects import FDIR_SB_MAP_SAVE
import numpy as np


class Map_Parameters_Class: 

   def __init__(self):
       # This is well define by options number
       self.XrayBand = 1
       self.XrayBandstr = '[0.5-2.0]'
       # This is defining which class of surface brightness we are using
       self.ModelOption = 1
       # Map size in degree (it can not be lager than 5 degree!!)
       self.DECMapSize = 0.5 # in degree
       self.RAMapSize = 0.5 # in degree
       # Number of pixels in each side
       self.PixSide = 3600  #7200/120
       # Effective area of detector in cm^2
       self.effArea = 1.0
       # Exposure time in second
       self.exposureTime = 1000.0 
       # Map Intensity Unit
       self.mapUnit = 'ergs / s / cm^2 / pixel'


class Map_Class:
 
   def __init__(self):
      self.MPa = Map_Parameters_Class()
      # Map is nxn array which firs element is the RA and second one is DEC
      self.MAP = np.zeros([self.MPa.PixSide,self.MPa.PixSide])
      # Map file name
      self.fname = 'sample.fits'
      from Surface_Brightness_Models import Surface_Brightness_Class 
      SBC = Surface_Brightness_Class()
      self.flux2SurfaceBrightness = SBC.F2SB_Class.flux2SurfaceBrightness
      del SBC

      # Center of the map
      self.DECc   = 0.0
      self.RAc    = 0.0

      self.DECmin = 0.0
      self.RAmin  = 0.0
      self.DECmax = 0.0
      self.RAmax  = 0.0

      self.dDEC   = 0.0
      self.dRA    = 0.0
      # Number of pixels = nside x nside
      # ( DECmax - DECmin ) / dDEC = nside 
      # ( RAmax  - RAmin  ) / dRA  = nside 

   def update(self):#,DECmin,RAmin,DECmax,RAmax,dDEC,dRA,nside):
      self.i,self.j=np.mgrid[0:self.MPa.PixSide,0:self.MPa.PixSide]
      self.DECmin = - self.MPa.DECMapSize / 2.
      self.RAmin  = - self.MPa.RAMapSize / 2.
      self.DECmax = self.MPa.DECMapSize / 2.
      self.RAmax  = self.MPa.RAMapSize / 2.
      self.nside  = self.MPa.PixSide
      self.dDEC   = self.MPa.DECMapSize / float(self.nside)
      self.dRA    = self.MPa.RAMapSize / float(self.nside) 
      self.DECc   = 0. 
      self.RAc    = 0. 
      from XCat_Objects import StertoArcmin2,pi
      #area = (self.DECmax-self.DECmin)*(self.RAmax-self.RAmin)*4.0*pi*pi**2/(180.0*180.0*self.nside*self.nside)          
      #self.area = StertoArcmin2*area
      self.fname = 'sample'
   
   def ang2pix(self,RA,DEC):
      # RA  component
      i = int((RA-self.RAmin)/self.dRA)
      # DEC component
      j = int((DEC-self.DECmin)/self.dDEC)
      return (i,j)

   def pix2ang(self,i,j):
      RA  = self.RAmin +  self.dRA/2.0 + float(i)*self.dRA
      DEC = self.DECmin+ self.dDEC/2.0 + float(j)*self.dDEC
      return (RA,DEC)

   # Exact formula for angular distance
   #(Note that RA and DEC should be in radian unit)
   def angularDistance(RA1,DEC1,RA2,DEC2):
     return np.arccos( np.sin(DEC1)*np.sin(DEC2)
                      + np.cos(DEC1)*np.cos(DEC2)*np.cos(RA1-RA2) )
      
   def pixInsideArea(self,RA,DEC,angularDis):
      pixi = [];      pixj = []
      angularDis2 = angularDis*angularDis
      (imin,jmin) = self.ang2pix(RA-angularDis,DEC-angularDis)
      (imax,jmax) = self.ang2pix(RA+angularDis,DEC+angularDis)
      for i in range(imin-1,imax+1):
         for j in range(jmin-1,jmax+1):
           (newRA,newDEC) = self.pix2ang(i,j)
           if ( ((RA-newRA)**2+(DEC-newDEC)**2) < angularDis2 ): 
              pixi.append(i)
              pixj.append(j)
      return (pixi,pixj)
 
   # Working in real space (not log space)
   def fluxCalculator(self,RA,DEC,angularDis):
      (pixi,pixj) = self.pix_inside_halo(RA,DEC,angularDis)
      flux = 0.0
      for i in range(len(pixi)):
         if (pixi[i] < 0 or pixj[i] < 0 
             or pixi[i] > self.nside-1 
             or pixj[i] > self.nside-1): pass
         else: flux += self.Map[pixi[i],pixj[i]]
      return flux*self.MPa.effArea

   def showMap(self):
      import matplotlib.pyplot as plt 
      plt.figure()
      plt.clf()
      plt.imshow(self.MAP.T,aspect='auto',origin='lower'
                ,extent=(self.RAmin,self.RAmax,self.DECmin,self.DECmax))
      plt.colorbar()
      plt.xlabel('RA',{'fontsize':20})
      plt.ylabel('DEC',{'fontsize':20})
      plt.title('DEC =%.2f , RA =%.2f '%(self.DECc,self.RAc),{'fontsize':20})
      plt.grid()
      plt.show()
      print 'done'

   def saveMapPic(self):
      import matplotlib.pyplot as plt
      import os
      plt.figure()
      plt.clf()
      plt.imshow(np.log10(self.MAP/6.24e8+1e-20).T,\
                aspect='auto',origin='lower',\
                vmin=-20,vmax=-16,\
                #vmin=-18.,vmax=-15.,\
                extent=(self.RAmin,self.RAmax,self.DECmin,self.DECmax))
      plt.colorbar()
      plt.xlabel('RA',{'fontsize':20})
      plt.ylabel('DEC',{'fontsize':20})
      plt.title('DEC =%.2f , RA =%.2f '%(self.DECc,self.RAc),{'fontsize':20})
      plt.grid()
      path = FDIR_SB_MAP_SAVE
      if not os.path.exists(path): os.makedirs(path)
      #plt.savefig(path+self.fname+'_'+self.MPa.XrayBandstr+'.pdf')
      plt.savefig(path+self.fname+'.pdf')
      plt.close()

   def saveMapFits(self):
      import os
      import pyfits
      path = FDIR_SB_MAP_SAVE
      if not os.path.exists(path): os.makedirs(path)
      #fdir = path + self.fname + r'_' + self.MPa.XrayBandstr + r'.fit'
      fdir = path + self.fname + r'.fit'
      print "Removing old halo/map catalog (if exist) :", fdir
      os.system('rm -r -f %s'%fdir)
      print 'Saving data', self.fname

      #hduListObj = pyfits.open('./Catalog/Output_File/XCAT_DMONLY_L400N1024_zle3p0_02_[0.5-2.0].fit')#,
      #mode='update', save_backup='True')
      hduListObj = pyfits.open('./Catalog/Output_File/'+self.fname+'.fit')[1:]
      hdr = pyfits.open('./Catalog/Output_File/'+self.fname+'.fit')[0].header
      hdr.add_comment('------- MAP INFO -------- ')
      hdr.set('XRAY_MAP','TRUE')     
      hdr.set('DEC_SIZE',str(self.MPa.DECMapSize),'unit:degree')
      hdr.set('RA_SIZE',str(self.MPa.RAMapSize),'unit:degree')
      hdr.set('PIX_NUM',str(self.MPa.PixSide))
      hdr.set('PIX_UNIT',self.MPa.mapUnit)
      hdr.set('EFF_AREA',str(self.MPa.effArea),'unit:cm^2') 
      hdr.add_comment(' --------------------------- ')
      hdr.add_comment('This map is created by XCAT version 0.0.2')

      hdu = pyfits.PrimaryHDU(data=self.MAP/6.24e8+1e-20,header=hdr)
      hdu = [hdu]
      for i in range(len(hduListObj)): hdu.append(hduListObj[i])
      thdulist = pyfits.HDUList(hdu) 
      thdulist.writeto(fdir)
      thdulist.close()
      #hduListObj.insert(0,hdu)
      print "Map file is created."

   def addHalos2Map(self,Halos):
      self.MPa.XrayBandstr = Halos.XrayBandstr 
      for iHalo in range(len(Halos.RA[::])):
         #print iHalo
         if (Halos.F[iHalo] == 0.): pass
         else: self.flux2SurfaceBrightness(self,Halos,iHalo,self.MPa)
         #else: self.flux2SurfaceBrightness(self,Halos,iHalo,self.i,self.j)
      self.fname = Halos.fname

   #def makeEventMap(self):
   #   self.eventClass.test()
   #   print "OK"


