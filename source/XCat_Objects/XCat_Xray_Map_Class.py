import numpy as np
from random       import uniform
from XCat_Objects import DtoR, RtoD

class Xray_Surface_Brightness_Class():

   def __init__(self,option=1):

      import sys 
      if (option==1): 
         self.F2SB = self.flux2SurfaceBrightnessSingleBetaModel
      elif (option==2): pass
      else:
         print "ERROR: Input option does not exist!"
         print "Please change XML file and try again."
         raw_input("Press enter to exit! ")
         sys.exit(2)

   def flux2SurfaceBrightnessSingleBetaModel(self,Map,Halos,samplei):
 
      power = -3.0*Halos.beta[samplei]+0.5
      RA    = Halos.RA[samplei]
      DEC   = Halos.DEC[samplei]

      # Rc = R200 * 0.075 * (M200/10^12)^0.6
      # prescrition astro-ph/9703027
      # self.Rc = Halo_data.R200[samplei]*0.075 * (Halo_data.M200[samplei]/10**12)**0.6
      Rc          = Halos.Rc[samplei]
      thetac2     = (RtoD*Rc*(1.0+Halos.Z[samplei])/Halos.DIS[samplei])**2
      (pixi,pixj) = Map.pixInsideArea(RA,DEC,10.0*np.sqrt(thetac2))
      pixValue    = np.zeros(len(pixi))
      summ        = 0.0
      #sigmaRA = Map.dRA/2.0
      #sigmaDEC= Map.dDEC/2.0

      sample_num = 1
      print "Sample Number : ", sample_num*len(pixi)
      for i in range(len(pixi)):
         (sRA,sDEC) = Map.pix2ang(pixi[i],pixj[i]) 
         #for sSBi in range(sample_num):
         RAi      = sRA  #+ uniform(-sigmaRA,sigmaRA)
         DECi     = sDEC #+ uniform(-sigmaDEC,sigmaDEC)
         stheta2  = ((RAi - RA)**2 + (DECi-DEC)**2)
         pixValue[i] += ( 1.0 + stheta2/thetac2 )**power
      summ = np.sum(pixValue[:])
      pixValue = Halos.F[samplei] * pixValue / summ / Map.MPa.effArea
      for i in range(len(pixValue[:])): 
         if (pixi[i] < 0 or pixj[i] < 0 
            or pixi[i] > Map.nside-1 
            or pixj[i] > Map.nside-1): pass
         else: Map.MAP[pixi[i],pixj[i]] += pixValue[i] 

