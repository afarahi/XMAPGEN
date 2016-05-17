import numpy as np
from XCat_Objects import DtoR, RtoD

class Surface_Brightness_Model():

   def __init__(self):
      self.clusterEffRadCOEFF = 20.0
      pass


   def flux2SurfaceBrightness(self,Map,Halos,samplei,MPa):
 
      pixelSide = MPa.PixSide 
      mapCenterShift = MPa.DECMapSize/2. # in degree
      mapSize = MPa.DECMapSize # in degree

      RApix   = (Halos.RA[samplei]+mapCenterShift)*pixelSide/mapSize
      DECpix  = (Halos.DEC[samplei]+mapCenterShift)*pixelSide/mapSize
      rad     = Halos.Rc[samplei]*pixelSide/mapSize
      effRad  = self.clusterEffRadCOEFF*rad
      thetac2 = rad**2

      bpower = -3.0*Halos.beta[samplei]+0.5

      iMax = int(RApix+effRad);   iMin = int(RApix-effRad)
      jMax = int(DECpix+effRad);  jMin = int(DECpix-effRad)
      if (iMin>=pixelSide-1 or jMin>=pixelSide-1): return 0
      if (iMax<=0 or jMax<=0): return 0

      i,j = np.mgrid[iMin:iMax,jMin:jMax]
      mask = (i-RApix)**2 + (j-DECpix)**2 < effRad**2
      mapA = mask*( 1. + ((RApix - i)**2 + (DECpix - j)**2)/thetac2 )**bpower
      mapA *= Halos.F[samplei] / sum(sum(mapA))

      iMaxA = min(iMax,pixelSide)
      iMinA = max(iMin,0)
      jMaxA = min(jMax,pixelSide)
      jMinA = max(jMin,0)

      #print (iMinA-iMin),(iMaxA-iMin)-(iMinA-iMin),\
      #      (jMinA-jMin),(jMaxA-jMin)-(jMinA-jMin),\
      #       "***" ,iMaxA-iMinA, jMaxA-jMinA

      Map.MAP[iMinA:iMaxA,jMinA:jMaxA] += \
            mapA[(iMinA-iMin):(iMaxA-iMin),(jMinA-jMin):(jMaxA-jMin)] 
        
   '''
   def flux2SurfaceBrightness(self,Map,Halos,samplei):
  
      RA    = Halos.RA[samplei]
      DEC   = Halos.DEC[samplei]
 
      bpower = -3.0*Halos.beta[samplei]+0.5

      thetac2 = Halos.Rc[samplei]**2
      #thetac2 = ( RtoD * Halos.Rc[samplei] * (1.0+Halos.Z[samplei]) / Halos.DIS[samplei] )**2

      (pixi,pixj) = Map.pixInsideArea(RA,DEC,10.0*np.sqrt(thetac2))
      pixValue = np.zeros(len(pixi))

      summ = 0.0
      sigmaRA = Map.dRA/2.0
      sigmaDEC= Map.dDEC/2.0

      #print "Sample Number : ", len(pixi)
      for i in range(len(pixi)):
         (sRA,sDEC) = Map.pix2ang(pixi[i],pixj[i]) 
         stheta2 = ((sRA - RA)**2 + (sDEC-DEC)**2)
         pixValue[i] = ( 1.0 + stheta2/thetac2 )**bpower
      summ = np.sum(pixValue[:])
      pixValue = Halos.F[samplei] * pixValue / summ
      for i in range(len(pixValue[:])): 
         if (pixi[i] < 0 or pixj[i] < 0 
            or pixi[i] > Map.nside-1 
            or pixj[i] > Map.nside-1): pass
         else: Map.MAP[pixi[i],pixj[i]] += pixValue[i] 
   '''
