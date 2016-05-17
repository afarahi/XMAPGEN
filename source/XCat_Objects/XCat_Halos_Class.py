# For reading the halos catalog fit file
# Return fals , none if it fails otherwise 
# true and halo infomations

from XCat_Utilities import read_data_string
from XCat_Objects import FDIR_HALOS
import numpy as np
import pyfits
import os

def Read_Halos_Catalog_fit_file(fname, fdir= './Catalog/Input_File/'):
   fdir = fdir + fname
   #print "Reading Halos Catalog :", fdir
   try:
      Halos_list = pyfits.open(fdir)
      #Halos_list.info()
      Halos_head = Halos_list[1].header
      Halos_data = Halos_list[1].data[::]
      Halos_data = Halos_data[ Halos_data.M500 > 8e12 ]
      #print "Reading Halos Catalog is done."
      Halos_list.close()
   except IOError:
      print "Error: can\'t find file or read data"
      raw_input("Press enter to continue. ")
      return (False, 0)
   return (fname[:-4], Halos_data,Halos_list[0].header)


# Halos Class 
class Halos_Class():

   def __init__(self):
      from XCat_Objects import Input_Parameters,\
                               Halo_Status_Class,\
                               Solvers_Class
      self.InputParam = Input_Parameters()
      self.Status     = Halo_Status_Class()
      self.Solvers    = Solvers_Class()
      #self.Solvers.ProperDistanceTabulate(self.InputParam,4.0)

      self.ID   = []; self.pd   = []; self.RA   = []; self.DEC   = []
      self.R500 = []; self.M500 = []; self.Z_red = []
      self.lgT  = []; self.lgLx = []; self.lgFx = []
      self.Rc   = []; self.beta = []
      self.XCatRA  = [];  self.XCatDEC = []
      self.RAref  = 0.0
      self.DECref = 0.0
      #self.OutputParam = Output_Parameters()

   def readPDistance(self,Z):
      import numpy as np
      z_tab, r_tab = np.loadtxt("./Output/tabulated_data/Proper_Distance.txt"\
                                , unpack=True)
      len_z_tab = len(z_tab); max_z_tab = max(z_tab)
      PDis = np.zeros(len(Z))
      for i in range(len(Z)):
         z_loc= int(len_z_tab*Z[i]/max_z_tab)
         PDis[i] = r_tab[z_loc] + (Z[i]-z_tab[z_loc])\
             *(r_tab[z_loc+1]-r_tab[z_loc])/(z_tab[z_loc+1]-r_tab[z_loc])
      return PDis

   def addHalosCatalog(self, fname=None):

      from numpy  import log, log10, zeros
      from XCat_Objects import R2D, D2R

      self.fname = fname

      fname,Halos,hdr = Read_Halos_Catalog_fit_file(fname=self.fname, fdir=FDIR_HALOS)


      if fname:

        self.RA = Halos['RA']   
        self.DEC= Halos['DEC']
        self.R500 = Halos['R500']
        self.M500 = Halos['M500']
        #! --- 1.55 is the mean extrated from arXiv:0705.0358 
        #self.M500.append(Halos_info[i][1]/1.5)
        self.Z_red = Halos['Z']
        self.pd = self.readPDistance(Halos['Z'])
           
        self.ID = Halos['HALOID'] 
        self.XCatRA = Halos['RA']
        self.XCatDEC = Halos['DEC']
             
        n = len(self.RA[:])
        self.lgT  = list(zeros(n))
        self.lgLx = list(zeros(n))
        self.lgFx = list(zeros(n))
        self.Rc   = list(zeros(n))
        self.beta = list(zeros(n))
        self.Z_min= min(self.Z_red)
        self.Z_max= max(self.Z_red)
        self.number_of_halos = n
        self.fname= fname
        self.Status.update(self)
        self.Status.XCatPreferedCoordinate = True

        del Halos, n

        print "Number of Halos = ", self.number_of_halos
        print "Halos class initialized successfully."



   #def creatAGN(self):
   #   from XCat_Objects import AGN_Class 
   #   self.AGNs = AGN_Class()
   #   self.Status.AGNsDataExist = self.AGNs.creatAGN(self)

   def SolveLxTxFlux(self):
      if (self.Status.HalosDataExist):
         self.Solvers.LxTxSolver(self)
         self.Status.LxTxSolved = True
      else:
         print "ERROR: Halos class is not initialized! "
         print "Please add a halos catalog! "
         raw_input("Press enter to continue ... ")


   def SaveHalosData(self):
     # Cheking for error
     if (self.Status.HalosDataExist == False):
        print "ERROR: Halos class is not initialized! "
        print "Please add a halos catalog! "
        raw_input("Press enter to continue ... ")
        return 0
     if (self.Status.LxTxSolved == False):
        print "ERROR: Lx, Tx, and flux was not created! "
        print "Please solve Lx, Tx, and flux for the halo catalog! "
        raw_input("Press enter to continue ... ")
        return 0

     fdir = './Catalog/Output_File/%s.fit'%(self.fname)
     print "Removing old halo catalog (if exist) :", fdir
     os.system('rm -r -f %s'%fdir)
     print "Saving new halo catalog :", fdir

     hdr =  pyfits.Header()      
     self._setHeader(hdr)
     PrimaryHDUList = np.arange(1)
     # HALOS INFO
     T = np.array(self.lgT)
     Lx= np.array(self.lgLx)
     F = self.InputParam.xray_band.fluxFac*np.power(10.,np.array(self.lgFx))
     FU= self.InputParam.xray_band.fluxFacStr
     col1 = pyfits.Column(name='HALOID',              
                          format='K', array=np.array(self.ID))
     col2 = pyfits.Column(name='RA',   unit='degree', 
                          format='E', array=np.array(self.RA))
     col3 = pyfits.Column(name='DEC',  unit='degree', 
                          format='E', array=np.array(self.DEC))
     col4 = pyfits.Column(name='R500', unit='Mpc/h',  
                          format='E', array=np.array(self.R500))
     col5 = pyfits.Column(name='M500', unit='Msun/h', 
                          format='E', array=np.array(self.M500))
     col6 = pyfits.Column(name='Z',                   
                          format='E', array=np.array(self.Z_red))
     col7 = pyfits.Column(name='T',    unit='keV',    
                          format='E', array=T)
     col8 = pyfits.Column(name='Lx',   unit='ergs/s', 
                          format='D', array=Lx)
     col9 = pyfits.Column(name='FLUX', unit=FU,       
                          format='E', array=F)
     col10= pyfits.Column(name='DIS',  unit='Mpc/h',  
                          format='E', array=np.array(self.pd))
     col11= pyfits.Column(name='XCAT_RA' ,unit='degree', 
                          format='E', array=np.array(self.XCatRA))
     col12= pyfits.Column(name='XCAT_DEC' ,unit='degree',
                          format='E', array=np.array(self.XCatDEC))
     col13= pyfits.Column(name='RC' ,unit='degree', 
                          format='E', array=np.array(self.Rc))
     col14= pyfits.Column(name='BETA',
                          format='E', array=np.array(self.beta))
     cols = pyfits.ColDefs([col1,  col2, col3, col4,  col5,  col6,
                            col7,  col8, col9, col10, col11, col12,
                            col13, col14])
 
     #tbhdu_halo = pyfits.new_table(cols)
     hdu = pyfits.PrimaryHDU(data=PrimaryHDUList,header=hdr)
     tbhdu_halo = pyfits.BinTableHDU.from_columns(cols)

     tbhdu_gals = pyfits.open(FDIR_HALOS+'%s.fit'%self.fname)[2]

     thdulist = pyfits.HDUList([hdu, tbhdu_halo, tbhdu_gals]) 
     thdulist.writeto(fdir)
     return 0

     if (self.Status.AGNsDataExist == False): 
        hdu = pyfits.PrimaryHDU(data=PrimaryHDUList,header=hdr)
        thdulist = pyfits.HDUList([hdu, tbhdu_halo]) 
        thdulist.writeto(fdir)
        return 0

     # AGNS INFO  
     col1 = pyfits.Column(name='AGN',              
                          format='K', array=np.array(self.AGNs.AGNID))
     col2 = pyfits.Column(name='HALOID',              
                          format='K', array=np.array(self.AGNs.HALOID))
     col3 = pyfits.Column(name='RA',   unit='degree', 
                          format='E', array=np.array(self.AGNs.RA))
     col4 = pyfits.Column(name='DEC',  unit='degree', 
                          format='E', array=np.array(self.AGNs.DEC))
     col5 = pyfits.Column(name='Z',                   
                          format='E', array=np.array(self.AGNs.Z))
     col6 = pyfits.Column(name='Lx',   unit='ergs/s', 
                          format='D', array=np.array(self.AGNs.LxAGN))
     col7 = pyfits.Column(name='FLUX', unit=FU,       
                          format='E', array=np.array(self.AGNs.fluxAGN))
     col8 = pyfits.Column(name='XCAT_RA' ,unit='degree', 
                          format='E', array=np.array(self.AGNs.XCatRA))
     col9 = pyfits.Column(name='XCAT_DEC' ,unit='degree',
                          format='E', array=np.array(self.AGNs.XCatDEC))
     cols = pyfits.ColDefs([col1, col2, col3, col4,  col5,
                            col6, col7, col8, col9])
     tbhdu_AGN = pyfits.TableHDU.from_columns(cols)
 
     hdu = pyfits.PrimaryHDU(data=PrimaryHDUList,header=hdr)
     thdulist = pyfits.HDUList([hdu, tbhdu_halo, tbhdu_AGN]) 
     thdulist.writeto(fdir)

   def _setHeader(self,hdr):
      if self.Status.XCatPreferedCoordinate:
         hdr.set('XCATPREF','TRUE','XCat prefered coordinate transformation')
         hdr.set('RA_REF',str(self.RAref),'XCat prefered coordinate transformation')
         hdr.set('DEC_REF',str(self.DECref),'XCat prefered coordinate transformation')
      else:
         hdr.set('XCATPREF','FALSE','XCat prefered coordinate transformation')
      if self.Status.LxTxSolved:
         hdr.set('XRAYINFO','TRUE','whether xray information exist.')
         hdr.set('XRAYBAND',self.InputParam.xray_band.Str,'unit:keV')
      else:
         hdr.set('XRAYINFO','FALSE','whether xray information exist.')
      hdr.set('XRAY_MAP','FALSE')
      hdr.set('SZ_MAP','FALSE')
      # COSMOLOGY 
      hdr.add_comment('------- COSMOLOGY --------')
      hdr.set('H0',str(self.InputParam.h_0))
      hdr.set('OMEGA_DE',str(self.InputParam.Omega_DE))
      hdr.set('OMEGA_M',str(self.InputParam.Omega_M))
      hdr.set('OMEGA_b',str(self.InputParam.Omega_b))
      hdr.set('OMEGA_R',str(self.InputParam.Omega_R))
      hdr.set('sigma_8',str(self.InputParam.sigma_8))
      hdr.set('w',str(self.InputParam.w))
      hdr.set('ns',str(self.InputParam.ns))
      hdr.add_comment('This table is created by XCAT version 0.0.3')
  




class Map_Halos_Class:

   def __init__(self):
      self.DIS  = []; self.DEC = []; self.F = []
      self.R500 = []; self.RA  = []; self.Z = []
      self.Rc   = []; self.beta = []
      self.RA_ref = 0.; self.DEC_ref = 0. 
      self.XCat_prefered = False; self.halosDataExist= False 
 
   def update(self,RACAT=0.0,DECCAT=0.0,Index=None,fname=None):
      from XCat_Utilities import read_data_string
      from XCat_index2coordinateReference import index2coordinateReference
      import pyfits

      if ( fname != None ): pass
      else:
         if ( Index != None ):
            RACAT,DECCAT = index2coordinateReference(index=Index)

         fname = read_data_string(tag_name='File_name',
                        file_name='./parameters/Input_Parameters.xml')
         fname = 'XCAT_RA_%0.2f_DEC_%0.2f_'%(RACAT,DECCAT) + fname

      Halos = pyfits.open('./Catalog/Output_File/'+fname)[1].data
      hdr = pyfits.open('./Catalog/Output_File/'+fname)[0].header

      try: 
        #self.XCat_prefered = bool(hdr['XCATPREF'])
        self.RA = Halos['XCAT_RA']
        self.DEC = Halos['XCAT_DEC']
        self.R500 = Halos['R500']
        self.Z = Halos['Z']
        #self.DIS = Halos['DIS']
        self.F = Halos['FLUX']
        self.Rc = Halos['RC']
        self.beta = Halos['BETA']
        self.RA_ref = float(hdr['RA_REF'])
        self.DEC_ref= float(hdr['DEC_REF'])
        self.fname = fname[:-4]
        self.halosDataExist= True
        self.XrayBandstr = hdr['XRAYBAND']
        print "Number of halos are : ",len(Halos[:]['Z'])
      except KeyError:
        print "The input table has some problem. "
        print "It is not consistent with XCat table structure."
        self.halosDataExist= False

