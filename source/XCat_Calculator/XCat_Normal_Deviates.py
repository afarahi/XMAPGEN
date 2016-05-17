# -- make normal deviates (g)
def normal_dev(n):
   from random import random
   from numpy  import zeros, sqrt, log
   gdev   = zeros(n)
   gdevx  = zeros(n)
   i      = 0
   while (i<n):
      v1 = random()
      v2 = random()
      r  = v1**2 + v2**2
      if (r < 1.0 and r > 0.0):
        fac      = sqrt(-2.0*log(r)/r)
        gdev[i]  = v1*fac
        gdevx[i] = v2*fac
        i       += 1
   return (gdev,gdevx)


def correlatedNormalDistribution(xBar,yBar,xSig,ySig,corr,y,ymin=0.0):
   import random as rd
   from numpy import sqrt
   mean = xBar + xSig * corr * (y-yBar) / ySig
   sig  = xSig * sqrt(1. - corr*corr)
   samp = rd.normalvariate(mean,sig)
   while (ymin > samp):
      samp = rd.normalvariate(mean,sig)
   return samp

# TEST IS PASSED !!
'''
import numpy as np
import random
y = np.zeros(200)
x = np.zeros(200)
for i in range(200):
   y[i] = random.normalvariate(15.,4.)
   x[i] = correlatedNormalDistribution(10.,15.,2.,4.,0.25,y[i])

print np.corrcoef(x,y)
print np.mean(x), np.std(x)
'''


