from XCat_Objects          import fluxfac, ln10, ergs2keV

def fluxUnitCanvFac(FluxUnitStr = 'ergs'):
    from XCat_Objects import ergs2keV
    if    (   FluxUnitStr == 'ergs' or FluxUnitStr == 'erg'
           or FluxUnitStr == 'ERGS' or FluxUnitStr == 'ERG'):
       return 1, 'ergs/s/cm^2'
    elif (    FluxUnitStr == 'keV' or FluxUnitStr == 'kev'
           or FluxUnitStr == 'KeV' or FluxUnitStr == 'KEV'):
       return ergs2keV, 'keV/s/cm^2'
    else:
       print 'WARNING: Unit %s was not implemented!'%FluxUnitStr
       print 'This code using default energy unit for saving flux.'
       print 'The default energy unit is is ergs.'
       raw_input("Press enter to continue. ")
       return 1, 'ergs/s/cm^2'


class Xray_Band_Class:
 
   def __init__(self,option=1):

      from XCat_Utilities  import read_data_string
      import sys 

      FluxUnitStr    = read_data_string(tag_name = 'Flux_Unit',file_name = 'parameters/Input_Parameters.xml')
      self.fluxFac, self.fluxFacStr = fluxUnitCanvFac(FluxUnitStr=FluxUnitStr)

      self.Conv_fac  = 1.0 #0.62
      self.Str       = '[0.5-2.0]'
      self.opt_num   = 1
      self.Descrption= 'Soft band'

      if (option==1): 
         #K-Correction function
         self.K_corr= self.K_corr_1
      elif (option==2): 
         self.Conv_fac  = 1.
         self.Str       = '[0.1-2.4]'
         self.opt_num   = 2
         self.Descrption= 'Soft band'
         #K-Correction function
         self.K_corr= self.K_corr_2
      elif (option==3):
         self.K_corr= self.K_corr_1
         import pyfits; import numpy as np
         fdir = './parameters/Models/K_Correction/'
         #x = pyfits.open(fdir+'L_bol_cube_CR_1.fits')[1].data
         x = pyfits.open(fdir+'CRb05-2_overdl2.fits')[1].data
         #x = pyfits.open(fdir+'alldet_b05-2_crtab.fits')[1].data
         self.Str = '[0.5-2.0]'
         self.lgTCUBE = np.log10(x[0]['T'])
         self.ZCUBE = x[0]['Z']
         #self.NHCUBE= x[0]['NH']
         #self.L_CUBE = x[0]['L_CUBE'][:][:][:]
         self.L_CUBE = x[0]['CR_PERL_TIMESD2'][:][:]
         del x
      else:
         print "ERROR: Input XRay band option does not exist!"
         print "Please change XML file and try again."
         raw_input("Press enter to exit! ")
         sys.exit(2)


   #!--- K-Correction --> K_corr = Lobs / Lrest
   # For option = 1 (Nord et al. 2008 --- arXiv:0706.2189)
   # This fit is accurate to a few percent within z=2. (T in keV)
   def K_corr_1(self,Z,lgT):
      K1 = -0.209 + lgT * (  1.18  - 0.39*lgT  )
      K2 = -0.098 + lgT * ( -0.092 + 0.085*lgT )
      return 1. + Z*(K1 + K2*Z)


   # For option = 2 (Stanek et al. 2005 --- arXiv:astro-ph/0602324)
   # This is for band [0.1-2.4] keV (equation 4 in paper). (T in keV)
   #! --- set K_corr = sqrt(1+(1+log10(Tx/5))*z) 
   def K_corr_2(self,Z,lgT):
      from numpy import sqrt
      return sqrt(1. + (1.+lgT-0.69897)*Z)


   def fluxCalculator(self,r2,z_redshift,lglobs): 
      # convert comoving distance to Mpc for assumed h0
      dlMpc2 = r2 * (1.0+z_redshift)**2
      fx = fluxfac * (10**lglobs) / dlMpc2
      return fx


   def Luminosity2FluxWithKcorrection(self,Z,lgT,lgLx,comdis2):
      # Lx = Lb * fac 
      import pyfits; import numpy as np
      fdir = './parameters/Models/K_Correction/'
      x = pyfits.open(fdir+'L_0.5-2_L_bol_factor.fits')[1].data
     
      lgTCUBE = np.log10(x[0]['T'])
      fac = x[0]['LBAND_LBOL']; del x
      iT=int(float(len(lgTCUBE)-1) * (lgT-lgTCUBE[0])/(lgTCUBE[-1]-lgTCUBE[0]))
      T2 = lgTCUBE[iT+1]; T1 = lgTCUBE[iT]
      fac1 = fac[iT];      fac2 = fac[iT+1]
      fac = np.log10( ( fac1*(T2-lgT) + fac2*(lgT-T1) ) / (T2-T1) )
      lgLx= fac + lgLx 
      
      from numpy import log10
      #! --- set Lx[0.5-2.0] = 0.62 Lx[0.1-2.4]
      #! --- set Lx[0.5-2.0] = Lx[0.1-2.4] * conv_fac
      #! --- set K_corr = Lx_obs / Lx_rest
      LobsLrest = self.K_corr(Z,lgT)
      lglobs = lgLx + log10(self.Conv_fac * LobsLrest)
      return self.fluxCalculator(comdis2,Z,lglobs)

   def LCUBE_INDEX(self,Z=0.1,lgT=0.0,iNH=0):
      iZ = int( float(len(self.ZCUBE)-1) \
                * (Z-self.ZCUBE[0])/(self.ZCUBE[-1]-self.ZCUBE[0]) )
      ilgT = int( float(len(self.lgTCUBE)-1) \
                  * (lgT-self.lgTCUBE[0])/(self.lgTCUBE[-1]-self.lgTCUBE[0]) )
      return iZ, ilgT, iNH


   def Luminosity2FluxWithCube(self,Z,lgT,lgLx):
      iZ, ilgT, iNH = self.LCUBE_INDEX(Z=Z,lgT=lgT)
      Z2  = self.ZCUBE[iZ+1];     Z1 = self.ZCUBE[iZ]
      T2  = self.lgTCUBE[ilgT+1]; T1 = self.lgTCUBE[ilgT]
      L11 = self.L_CUBE[iZ][ilgT]
      L21 = self.L_CUBE[iZ+1][ilgT]
      L12 = self.L_CUBE[iZ][ilgT+1]
      L22 = self.L_CUBE[iZ+1][ilgT+1]

      #LxCR1 = ( L11*(Z2-Z)*(T2-lgT) + L21*(Z-Z1)*(T2-lgT) + \
      #          L12*(Z2-Z)*(lgT-T1) + L22*(Z-Z1)*(lgT-T1) \
      #        ) / ( (Z2-Z1)*(T2-T1) )
      CRPERLx = ( L11*(Z2-Z)*(T2-lgT) + L21*(Z-Z1)*(T2-lgT) + \
                L12*(Z2-Z)*(lgT-T1) + L22*(Z-Z1)*(lgT-T1) \
                ) / ( (Z2-Z1)*(T2-T1) )

      return CRPERLx * 10.**lgLx
      #return 10.**lgLx/LxCR1



''' COMMENTS
   #! --- Coversion factor Lx_[0.1-2.4] to Lx_[0.5-2.4]
   #from numpy import loadtxt
   #kt, c_fact = loadtxt("./Output/tabulated_data/Lx_0.1-2.4_to_Lx_0.5-2.0.txt"
   #                      , comments='#', unpack=True)
   #len_kt = len(kt)

      #! --- set Lx[0.5-2.0] = Lx[0.1-2.4]/conv_fac  
      #! data extracted from : ArXiv:astro-ph/0405546
      #if (lgT[i] <= kt[0]):
      #   conv_fac = c_fact[0]
      #elif (lgT[i] > kt[len_kt-1]):
      #   conv_fac = c_fact[len_kt-1]
      #else:
      #   for kt_n in range(len_kt-1):
      #      if (lgT[i] > kt[kt_n] and lgT[i] <= kt[kt_n+1]):
      #         conv_fac = c_fact[kt_n] + (lgT[i]-kt[kt_n])*(c_fact[kt_n+1]-c_fact[kt_n])/(kt[kt_n+1]-kt[kt_n])
      #         break
      #lnlobs[i] = lnl[i] + log(1.0/conv_fac) + InputParameters.K_corr*log(LobsLrest)


'''


