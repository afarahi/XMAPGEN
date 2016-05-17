import numpy as np
import pyfits
import sys 
Omega_M = 0.25
Omega_DE= 0.75


def LxErgsPerStoLxPhotonCounts(x):
    ergs2kev = 624.15*1000.
    return ergs2kev*x

def flux_calculator(r2,z_redshift,lobs): 
   #! --- convert comoving distance to Mpc for assumed h0
   Mpc2cm = 3.08567758e24
   fluxfac= ((1./Mpc2cm)*(1./Mpc2cm)/(4.0*np.pi))  #for cgs flux (ergs)
   dlMpc2 = r2 * (1.0+z_redshift)**2
   fx = fluxfac * lobs / dlMpc2
   return fx

def distance_ie(z):
    dl = np.sqrt( Omega_M*(1.0+z)**3 + Omega_DE )
    return 1.0/dl

def Proper_Distance(z_max): 
   #! --- convert redshift to proper distance
   integral_n = 200
   if (z_max < 0.5): integral_n = 50
   z = np.linspace(0.0,z_max,integral_n+1)
   dz= z[1]-z[0]
   integral_dl = 0.0
   for i in range(1,integral_n+1):   
      integral_dl += dz * ( distance_ie(z[i])+distance_ie(z[i-1]) ) / 2.0
   return (integral_dl * 2997.92458)

def Lx2Flux(Lx,lgT,Z):
   K1 = -0.209 + lgT * (  1.18  - 0.39 *lgT )
   K2 = -0.098 + lgT * ( -0.092 + 0.085*lgT )
   Kcorr = 1. + Z*(K1 + K2*Z)  #LobsLrest = Kcorr
   LxObs = Kcorr * Lx
   dis2 = ( Proper_Distance(Z) / 0.724 )**2
   fx   = flux_calculator(dis2,Z,LxObs)
   return fx


fdir = '../../../parameters/Models/K_Correction/'
print pyfits.open(fdir+'L_bol_cube_CR_1.fits')[1].header
x = pyfits.open(fdir+'L_bol_cube_CR_1.fits')[1].data
lgTCUBE = np.log10(x[0]['T'])
ZCUBE = x[0]['Z']
NHCUBE= x[0]['NH']
L_CUBE = x[0]['L_CUBE'][:][:][:]

#TLCUBE = np.log10(x[0]['T'][15:30])
#ZLCUBE = x[0]['Z'][5:150:5]
#L_CUBE = x[0]['L_CUBE'][:][:][:]

def LCUBE_INDEX(Z=0.1,lgT=0.0,iNH=0):
  iZ = int( float(len(ZCUBE)) * (Z-ZCUBE[0])/(ZCUBE[-1]-ZCUBE[0]) )
  ilgT = int( float(len(lgTCUBE)) * (lgT-lgTCUBE[0])/(lgTCUBE[-1]-lgTCUBE[0]) )
  return iZ, ilgT, iNH

def LCUBE_FLUX(Lx,Z=0.1,lgT=0.0,iNH=0):
   iZ, ilgT, iNH = LCUBE_INDEX(Z=Z,lgT=lgT,iNH=iNH)
   Z2 = ZCUBE[iZ+1]; Z1 = ZCUBE[iZ]
   T2 = lgTCUBE[ilgT+1]; T1 = lgTCUBE[ilgT]
   L11 = L_CUBE[iNH][iZ][ilgT]
   L21 = L_CUBE[iNH][iZ+1][ilgT]
   L12 = L_CUBE[iNH][iZ][ilgT+1]
   L22 = L_CUBE[iNH][iZ+1][ilgT+1]

   LxCR1 = ( L11*(Z2-Z)*(T2-lgT) + L21*(Z-Z1)*(T2-lgT) + \
             L12*(Z2-Z)*(lgT-T1) + L22*(Z-Z1)*(lgT-T1) \
           ) / ( (Z2-Z1)*(T2-T1) )

   return Lx/LxCR1


def ARYA_FLUX(Lx,Z=0.1,lgT=0.0):
   flux = Lx2Flux(Lx,lgT,Z)
   effArea  = 1387.71
   ergs2kev = 6.2415 * 10**8
   return flux * effArea * ergs2kev   

print LCUBE_FLUX(5e44,Z=0.45,lgT=-0.1,iNH=0)
print 0.62*ARYA_FLUX(5e44,Z=0.45,lgT=-0.1)
















