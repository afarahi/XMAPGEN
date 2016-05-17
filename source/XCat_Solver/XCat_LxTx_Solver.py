from numpy                 import log, log10, exp, sqrt, zeros, loadtxt
from XCat_Evolution_factor import Evolution_factors
from XCat_Normal_Deviates  import normal_dev, correlatedNormalDistribution
from XCat_Objects          import ln10

# Evolutionary Factor
def Evolution_factors(Cosmology,zred):
   Ez = sqrt(Cosmology.Omega_M*(1.0+zred)**3 + Cosmology.Omega_DE)
   # Ez = sqrt(0.25*(1.0+zred)**3 + 0.75)
   return Ez

#h(z)
def Hubz(Cosmology,zred):
   Ez = Evolution_factors(Cosmology,zred)
   return Ez #Cosmology.h_0*Ez


import numpy as np


class LxTx_Solver:


    def __init__(self):
        pass

    def solve(self,Halos):
        import numpy.random as npr; npr.seed(0)
        from XCat_Objects import RtoD
        import random; random.seed(0)

        n_data = Halos.number_of_halos
        IP = Halos.InputParam
        LxRcCF = 0.0
        RcBetaCF = 0.0
        betaBar = 0.6667
        betaSig = 0.000001
        sigRc = 0.000001 


        #! --- set up Lx , Tx scaling relations
        LxS = IP.LxScaling #Lxm_Parameters
        TxS = IP.TxScaling #Txm_Parameters

        # I am using Andrea notes (reference)
        #Efac = Evolution_factors(IP,Halos.Z_red)
        #Halos.lgT  = np.log10(Halos.M500*Efac/1.0e14/IP.h_0) / 1.49 + 0.2877 \
        #       + 0.0 #+ npr.normal(0.0,0.1,len(Halos.M500))/np.log(10.0)
        #Halos.lgLx = 0.4 + np.log10(Efac) + 2.89*(Halos.lgT - 0.602) \
        #       + 0.0 #+ npr.normal(0.0,0.267,len(Halos.M500))/np.log(10.0)

        #Told  =  Halos.lgT[:]
        #Lxold = Halos.lgLx[:]

        # calcualte evolution factor
        E  = Evolution_factors(IP,Halos.Z_red)

        # calculate mean-log-Lx + scatter
        Ep = Evolution_factors(IP,LxS.z_p)
        Mp = LxS.M_p* IP.h_0
           
        Halos.lgLx = np.log(LxS.Norm) \
                      + LxS.M_slope*np.log(Halos.M500/Mp) \
                      + LxS.E_slope*np.log(E/Ep) \
                      + npr.normal(0.0, LxS.sig,n_data)
        Halos.lgLx /=  np.log(10.0)
        Halos.lgLx /=  - 44.

        # calculate mean-log-Tx + scatter
        Ep = Evolution_factors(IP,TxS.z_p)
        Mp = TxS.M_p* IP.h_0

        Halos.lgT = np.log(TxS.Norm) \
                     + TxS.M_slope*np.log(Halos.M500/Mp) \
                     + TxS.E_slope*np.log(E/Ep) \
                     + npr.normal(0.0,TxS.sig,n_data)
        Halos.lgT /=  np.log(10.0)

        #print np.mean(Halos.lgT[:] - Told),
        #print np.std(Halos.lgT[:] - Told)
        #print np.mean(Halos.lgLx[:] - Lxold),
        #print np.std(Halos.lgLx[:] - Lxold)

        x = [];  y = []


        # calculate core radius, flux, and beta
        dis2 = (Halos.pd/IP.h_0)**2
        for i in range(n_data):
            #fx  = IP.xray_band.Luminosity2FluxWithKcorrection(\
            #                  Halos.Z_red[i],Halos.lgT[i],Halos.lgLx[i],dis2[i])\
            #     * 6.2415 * 1e8 * 1387.71
            fx = IP.xray_band.Luminosity2FluxWithCube(\
                         Halos.Z_red[i],Halos.lgT[i],Halos.lgLx[i])
            Halos.lgFx[i] = log10(fx+1e-40)                                                 
            Rcbar = RtoD * (0.14*Halos.R500[i]) *\
                      (1.0 + Halos.Z_red[i]) / Halos.pd[i]
            Halos.Rc[i] = npr.lognormal(np.log(Rcbar),sigRc)
 
            sigBeta = abs(betaSig*np.log(betaBar))
            betaBarC = np.log(betaBar) + RcBetaCF * sigBeta / sigRc *\
                     (np.log(Halos.Rc[i]) - np.log(Rcbar)) 
            sigBeta = np.sqrt(1.0-RcBetaCF**2) * sigBeta
            Halos.beta[i] = npr.lognormal(betaBarC,sigBeta)

        print "Fluxes, Tempratures, and Luminosities are assigned successfully!" 




def LxTx_Solver_ahan(Halos):

   
   IP = Halos.InputParam
   n_data = Halos.number_of_halos

   LxRcCF = 0.0
   RcBetaCF = 0.0
 
   betaBar = 0.6667
   betaSig = 0.000001
   sigRc = 0.000001 #abs(4.49-.134*np.log(Halos.M500[i])) 

   #! --- set up Lx at nominal log-mean relation at z=0.23
   LxS = IP.LxScaling() #Lxm_Parameters
   TxS = IP.TxScaling() #Txm_Parameters

   # I am using Andrea notes

   Efac = Evolution_factors(IP,Halos.Z_red)
   Halos.lgT  = np.log10(Halos.M500*Efac/1.0e14/IP.h_0) / 1.49 + 0.2877 \
               + 0.0 + npr.normal(0.0,0.1,len(Halos.M500))/np.log(10.0)
   Halos.lgLx = 0.4 + np.log10(Efac) + 2.89*(Halos.lgT - 0.602) \
               + 0.0 + npr.normal(0.0,0.267,len(Halos.M500))/np.log(10.0)
   #Efac23= Evolution_factors(IP,0.23)
   #Halos.lgLx = 1.16 + 1.61*np.log10(Halos.M500*Efac**2/4.8e14/Efac23**2) \
   #            + 0.0 #+ npr.normal(0.0,0.4,len(Halos.M500))/np.log(10.0)

   x = []
   y = []

   dis2 = (Halos.pd/IP.h_0)**2
   for i in range(n_data):
      #fx  = IP.xray_band.Luminosity2FluxWithKcorrection(\
      #                  Halos.Z_red[i],Halos.lgT[i],Halos.lgLx[i],dis2[i])\
      #     * 6.2415 * 1e8 * 1387.71
      fx = IP.xray_band.Luminosity2FluxWithCube(\
                        Halos.Z_red[i],Halos.lgT[i],Halos.lgLx[i])
      Halos.lgFx[i] = log10(fx+1e-40)                                                 
      Rcbar = RtoD * (0.14*Halos.R500[i]) *\
                     (1.0 + Halos.Z_red[i]) / Halos.pd[i]
      Halos.Rc[i] = npr.lognormal(np.log(Rcbar),sigRc)

      sigBeta = abs(betaSig*np.log(betaBar))
      betaBarC = np.log(betaBar) + RcBetaCF * sigBeta / sigRc *\
                     (np.log(Halos.Rc[i]) - np.log(Rcbar)) 
      sigBeta = np.sqrt(1.0-RcBetaCF**2) * sigBeta
      Halos.beta[i] = npr.lognormal(betaBarC,sigBeta)


   #                     np.log10(Halos.Rc) \
   #                     + npr.normal(0.,4.24-.134*np.log10.0(Halos.M500))
   #Halos.beta = betaBar*np.ones(n_data)
   #thetac2 = ( RtoD * Halos.Rc[samplei] * (1.0+Halos.Z[samplei]) / Halos.DIS[samplei] )**2

   '''
   for i in range(n_data):

      rcBar = 0.2 * Halos.R500[i]
      rcSig = 0.000001 * rcBar
      
      Halos.Rc[i] = correlatedNormalDistribution(rcBar,lgLxbar,\
                                                 rcSig,lglxSig,\
                                                 LxRcCF,lgLx)
      Halos.beta[i] = correlatedNormalDistribution(betaBar,rcBar,\
                                                   betaSig,rcSig,\
                                                   RcBetaCF,Halos.Rc[i])
      #correlatedNormalDistribution(xBar,yBar,xSig,ySig,corr,y)
   '''
   print "Fluxes, Tempratures, and Luminosities are assigned successfully!" 



def LxTx_Solver_old(Halos):

   IP = Halos.InputParam
   n_data = Halos.number_of_halos

   LxRcCF = -0.7
   RcBetaCF = 0.6
  
   betaBar = 0.667
   betaSig = 0.000001 * betaBar

   #! --- set up Lx at nominal log-mean relation at z=0.23
   Lxm_Para = IP.Lxm_para_obj() #Lxm_Parameters
   Txm_Para = IP.Txm_para_obj() #Txm_Parameters
   (gdev,gdevX) = normal_dev(n_data)
   lglxSig = Lxm_Para.siglnL/ln10

   for i in range(n_data):

      (Ez,Ez_n) = Evolution_factors(IP,Halos.Z_red[i])
      lgLxbar = (   Lxm_Para.a 
                    + Lxm_Para.alpha*log(Halos.nM500[i]/Lxm_Para.mass)
                    + Lxm_Para.LxEzPow*log(Ez/Ez_n)
                ) / ln10
      lghM15 = log10(IP.h_0 * Ez * Halos.nM500[i] / 10.)
      lgTbar = log10(Txm_Para.T_15) + Txm_Para.alpha_T*lghM15

      #! ---- include correlation coefficient here
      gdevT = IP.rTL*gdev[i] + sqrt(1.0-IP.rTL**2)*gdevX[i]
      #! --- assign random temps w/ covariance (thru gdevT)
      lgT = lgTbar #+ Txm_Para.siglnT*gdevT/ln10
      #! ---- random scatter in luminosities 
      lgLx = lgLxbar #+ Lxm_Para.siglnL*gdev[i]/ln10

      dis2 = (Halos.pd[i]/IP.h_0)**2
      fx = IP.xray_band.Luminosity2FluxWithKcorrection(\
                      Halos.Z_red[i],lgT,lgLx,dis2)
      #fx = IP.xray_band.Luminosity2FluxWithCube(\
      #                Halos.Z_red[i],lgT,lgLx,dis2)

      Halos.lgLx[i] = lgLx
      Halos.lgT[i]  = lgT                                        
      Halos.lgFx[i] = log10(fx)                                                 

      rcBar = 0.2 * Halos.R500[i]
      rcSig = 0.000001 * rcBar
      
      Halos.Rc[i] = correlatedNormalDistribution(rcBar,lgLxbar,\
                                                 rcSig,lglxSig,\
                                                 LxRcCF,lgLx)
      Halos.beta[i] = correlatedNormalDistribution(betaBar,rcBar,\
                                                   betaSig,rcSig,\
                                                   RcBetaCF,Halos.Rc[i])

      #correlatedNormalDistribution(xBar,yBar,xSig,ySig,corr,y)

   del gdev, gdevX

   print "Fluxes, Tempratures, and Luminosities are assigned successfully!" 


