
import pyfits
fdir = '../../../Catalog/Output_File/'
fname = 'XCAT_RA_84.02_DEC_-45.00_XCAT_Aardvark_1.0.fit'
xcub = pyfits.open(fdir+fname)[1].data

xold = pyfits.open(fdir+'/old/'+fname)[1].data

fname = 'XCAT_RA_84.02_DEC_-45.00_XCAT_Aardvark_1.0_old.fit'
xnew = pyfits.open(fdir+fname)[1].data

print xcub['FLUX'][:4] * 1e44
print xnew['FLUX'][:4] * 1387.71 * 6.2415 * 10**8

a = xcub['FLUX'][:] 
b = xnew['FLUX'][:] * 1387.71 * 6.2415 * 10**8

print "1 - CUBE/ARYA: "
c = 1 - a/b
print c[:-10]
print xcub['T'][:-10]
print xcub['Z'][:-10]

#print xcub['FLUX'][:10] 
#print xnew['FLUX'][:10]
