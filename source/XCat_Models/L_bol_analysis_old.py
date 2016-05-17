import numpy as np
Omega_M = 0.25
Omega_DE= 0.75

def Read_Halos_Catalog_fit_file(fname = None):
   import pyfits
   import sys 
   try:
      fitsdata = pyfits.open(fname)
      fitsdata.info()
      hdr = fitsdata[1].header
      #print hdr
      data = fitsdata[1].data[:]
      print "Reading Halos Catalog is done."
   except IOError:
      print "Error: can\'t find file or read data"
      raw_input("Press enter to continue. ")
      return (False, 0)
   return data

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
   #lgT = np.log10(T)
   K1 = -0.209 + lgT * (  1.18  - 0.39 *lgT )
   K2 = -0.098 + lgT * ( -0.092 + 0.085*lgT )
   Kcorr = 1. + Z*(K1 + K2*Z)
   #LobsLrest = Kcorr
   LxObs = Kcorr * Lx
   dis2 = ( Proper_Distance(Z) / 0.724 )**2
   #print np.sqrt(dis2), Z
   fx   = flux_calculator(dis2,Z,LxObs)
   return fx

x = Read_Halos_Catalog_fit_file('./L_0.5-2_L_bol_factor.fits')
Lband_Lbol = x[0]['LBAND_LBOL'][:]
del x

#print Lband_Lbol

x = Read_Halos_Catalog_fit_file('./L_bol_cube_CR_1.fits')

TLCUBE = np.log10(x[0]['T'][15:30])
ZLCUBE = x[0]['Z'][5:150:5]
L_CUBE = x[0]['L_CUBE'][:][:][:]

del x

nCR = np.zeros([len(TLCUBE),len(ZLCUBE)])
nAR = np.zeros([len(TLCUBE),len(ZLCUBE)])
for i in range(len(TLCUBE)):
  for j in range(len(ZLCUBE)):
    nCR[i,j] = L_CUBE[0][j*5+5][15+i] * Lband_Lbol[15+i]
    flux = Lx2Flux(1e44,TLCUBE[i],ZLCUBE[j])
    #print flux
    nAR[i,j] = 1e44 / flux
  
effArea  = 1387.71
ergs2kev = 6.2415 * 10**8.
SD2      = 41253.

for i in range(len(TLCUBE)):
   for j in range(len(ZLCUBE)):
      nAR[i,j] = nAR[i,j] / effArea / ergs2kev
      print nAR[i,j], nCR[i,j]

#print nAR

Ratio = np.zeros([len(TLCUBE),len(ZLCUBE)])

for i in range(len(TLCUBE)):
  for j in range(len(ZLCUBE)):
    Ratio[i,j] = nAR[i,j] / nCR[i,j]
    #print i, j, Ratio[i,j]

import matplotlib
matplotlib.use('Agg')
import matplotlib.pylab as plt

plt.imshow(Ratio[:,:].T,interpolation='nearest',
             vmin = 1., vmax=1.25,
             extent=[min(TLCUBE), max(TLCUBE), max(ZLCUBE), min(ZLCUBE)]
          )
plt.ylabel('Z')
plt.xlabel('log10(T)[keV]')
plt.colorbar()
#plt.show()
plt.savefig('result.png')










































'''
x = Read_Halos_Catalog_fit_file('./L_0.5-2_L_bol_factor.fits')
Lband_Lbol = x[0]['LBAND_LBOL'][:]
#T          = np.log10(x[0]['T'][18:30])
#print T
del x

print Lband_Lbol

x = Read_Halos_Catalog_fit_file('./L_bol_cube_CR_1.fits')

TLCUBE = np.log10(x[0]['T'][15:30])
ZLCUBE = x[0]['Z'][5:150:5]
L_CUBE = x[0]['L_CUBE'][:][:][:]

del x

nCR = np.zeros([len(TLCUBE),len(ZLCUBE)])
for i in range(len(TLCUBE)):
  for j in range(len(ZLCUBE)):
    nCR[i,j] = L_CUBE[0][j*5+5][15+i] * Lband_Lbol[15+i]

print nCR
#print ZLCUBE

x = Read_Halos_Catalog_fit_file('./XCAT_DMONLY_L400N1024_zle3p0_02_[0.5-2.0].fit')

TARYA = np.log10(x.field('T')[0:7500])
ZARYA = x.field('Z')[0:7500]
L_A   = x.field('Lx')[0:7500]
Flux_A= x.field('FLUX')[0:7500]
nCR_A = np.zeros(len(Flux_A))

#print len(TARYA)

effArea  = 1387.71
ergs2kev = 6.2415 * 10**8.
SD2      = 41253.

for i in range(len(Flux_A)):
  nCR_A[i] = L_A[i] / effArea / ergs2kev / Flux_A[i]


#ZLCUBE = np.linspace(0.0,1.9,20)
Ratio = np.zeros([len(TLCUBE),len(ZLCUBE),2])
dT = TLCUBE[1] - TLCUBE[0]
dZ = ZLCUBE[1] - ZLCUBE[0]

for i in range(len(nCR_A)):
  Tindex = int( (TARYA[i]-TLCUBE[0]-dT/2.) / dT )
  Zindex = int( (ZARYA[i]-ZARYA[0]) / dZ )
  if (Tindex < 0 or Zindex < 0): continue
  try: 
    Ratio[Tindex,Zindex,0] += nCR_A[i]
    Ratio[Tindex,Zindex,1] += 1.
  except IndexError: continue

for i in range(len(TLCUBE)):
  for j in range(len(ZLCUBE)):
    Ratio[i,j,0] = Ratio[i,j,0] / Ratio[i,j,1] / nCR[i,j]
    print i, j, Ratio[i,j,0], Ratio[i,j,1]

import matplotlib.pylab as plt
#plt.subplot(121)
#plt.imshow(Ratio[:,:,0].T,interpolation='nearest',
#              extent=[min(TLCUBE), max(TLCUBE), max(ZLCUBE), min(ZLCUBE)]
#          )

plt.imshow(Ratio[:,:,1].T,interpolation='nearest',
              extent=[min(TLCUBE), max(TLCUBE), max(ZLCUBE), min(ZLCUBE)]
          )
plt.ylabel('Z')
plt.xlabel('log10(T)[keV]')
plt.colorbar()
#plt.show()
plt.savefig('NumberOfObjects.png')


x = Read_Halos_Catalog_fit_file('./CRtab_b205-2_L05-2.fits')
TLCUBE = np.log10(x[0]['T'][13:26])
ZLCUBE = x[0]['Z'][0:200]
CR_CUBE = x[0]['CR_PERL_TIMESD2'][:][:][:]

del x

print np.size(CR_CUBE)
print len(CR_CUBE)
#print TLCUBE
#print len(ZLCUBE)
#print CR_CUBE

nCR = np.zeros([len(TLCUBE),len(ZLCUBE)])
for i in range(len(TLCUBE)):
  for j in range(len(ZLCUBE)):
    nCR[i,j] = CR_CUBE[0][j][13+i]

print nCR
#print ZLCUBE

x = Read_Halos_Catalog_fit_file('./XCAT_DMONLY_L400N1024_zle3p0_02_[0.5-2.0].fit')

TARYA = np.log10(x.field('T')[0:7500])
ZARYA = x.field('Z')[0:7500]
L_A   = x.field('Lx')[0:7500]
Flux_A= x.field('FLUX')[0:7500]
nCR_A = np.zeros(len(Flux_A))

#print len(TARYA)

effArea  = 1387.71
ergs2kev = 6.2415 * 10**8.
SD2      = 1.0 #41253.

for i in range(len(Flux_A)):
  nCR_A[i] = effArea * ergs2kev * 10**44 * Flux_A[i] / L_A[i] 


ZLCUBE = np.linspace(0.0,1.9,20)
Ratio = np.zeros([len(TLCUBE),len(ZLCUBE),2])
dT = TLCUBE[1] - TLCUBE[0]
dZ = ZLCUBE[1] - ZLCUBE[0]

for i in range(len(nCR_A)):
  Tindex = int( (TARYA[i]-TLCUBE[0]-dT/2.) / dT )
  Zindex = int( (ZARYA[i]-ZARYA[0]) / dZ )
  try: 
    Ratio[Tindex,Zindex,0] += nCR_A[i]
    Ratio[Tindex,Zindex,1] += 1.
  except IndexError: continue

for i in range(len(TLCUBE)):
  for j in range(len(ZLCUBE)):
    Ratio[i,j,0] = Ratio[i,j,0] / Ratio[i,j,1] / nCR[i,10*j+5]
    print i, j, Ratio[i,j,0], Ratio[i,j,1]

#print nCR_A
#print min(TARYA)
#print max(TARYA)



x = Read_Halos_Catalog_fit_file('./L_bol_cube_CR_1.fits')
TLCUBE = x[0]['T']
ZLCUBE = x[0]['Z']
L_C = x[0]['L_CUBE']
del x

L_CUBE = np.zeros([len(TLCUBE),len(ZLCUBE)])
for i in range(len(TLCUBE)):
   for j in range(len(ZLCUBE)): L_CUBE[i,j] = np.log10(L_C[i+len(TLCUBE)*j])


x = Read_Halos_Catalog_fit_file('./XCAT_DMONLY_L400N1024_zle3p0_02_[0.5-2.0].fit')
TARYA = x.field('T')
ZARYA = x.field('Z')
L_A   = x.field('Lx')
Flux_A= x.field('FLUX')




import matplotlib.pylab as plt


#print len(TARYA)
Flux_A = LxErgsPerStoLxPhotonCounts(Flux_A)

L_A_CUBE = L_A/Flux_A 

Zi = 20
Ti = 15
# Redshift Bin
Zmin = ZLCUBE[Zi]; Zmax = ZLCUBE[Zi+1]
# Temprature Bin
Tmin = TLCUBE[Ti]; Tmax = TLCUBE[Ti+5]
print Tmin, Tmax, Zmin, Zmax

LSampleArya = []
for i in range(len(Flux_A)):
  if ( TARYA[i] < Tmax and TARYA[i] > Tmin and  
       ZARYA[i] < Zmax and ZARYA[i] > Zmin ):
     LSampleArya.append(L_A_CUBE[i])
LSampleCUBE = ( L_C[Zi*len(TLCUBE)+Ti] + L_C[(Zi+1)*len(TLCUBE)+Ti] +
                L_C[Zi*len(TLCUBE)+Ti+5] + L_C[(Zi+1)*len(TLCUBE)+Ti+5]
              ) / 4.



print LSampleArya, LSampleCUBE


plt.figure(figsize=(10,8))
ax = plt.subplot(211)
ax.imshow(L_CUBE)
#plt.colorbar()
plt.show()
'''
