from XCat_Utilities     import *
from XCat_Properties    import *
from XCat_Objects       import REALIZATION_ID

def makingXRayCatalogs():
   Print_logo()
   from XCat_Objects import Halos_Class
   import sys; import time
   import random
   random.seed(int(sys.argv[1])+1000*REALIZATION_ID)

   ##################################################################
   #@@@@ STEP 1 [Initializing]
   #@@@@ 
   sTime = time.time()
   Halos = Halos_Class()
   #Halos.reInitilizeHalosClass(fname=sys.argv[2])
   eTime = time.time()-sTime
   print "Eplased time to initialize the code : ", eTime, " s"


   ##################################################################
   #@@@@ STEP 2 [Reading halos]
   #@@@@ 
   sTime = time.time()
   Halos.addHalosCatalog(fname=sys.argv[2])
   eTime = time.time()-sTime
   print "Eplased time to read halo catalogs : ", eTime, " s"

   ##################################################################
   #@@@@ STEP 3 [Solving for LX and Tx]
   #@@@@ 
   sTime = time.time()
   Halos.SolveLxTxFlux()
   eTime = time.time()-sTime
   print "Eplased time to solve X-ray properties : ", eTime, " s"


   ##################################################################
   #@@@@ STEP 4 [Saving clusters]
   #@@@@ 
   sTime = time.time()
   Halos.SaveHalosData()
   eTime = time.time()-sTime
   print "Eplased time to save outputs : ", eTime, " s"





def makingXRayRealization():
   from XCat_Objects import Halos_Class
   from XCat_Objects import Map_Halos_Class
   from XCat_Objects import Map_Class
   import sys; import time

   Print_logo()

   ##################################################################
   #@@@@ STEP 1 [Initializing and reading clusters catalogs]
   #@@@@ 
   sTime = time.time()
   Halos = Map_Halos_Class()
   Halos.update(fname=sys.argv[2])
   eTime = time.time()-sTime
   print "Eplased time to initialize and read halo catalog: ",eTime, " s"


   ##################################################################
   #@@@@ STEP 2 [Reading halos]
   #@@@@ 
   sTime = time.time()
   Map = Map_Class()
   Map.update()
   eTime = time.time()-sTime
   print "Eplased time to initialize the map object: ", eTime, " s"


   ##################################################################
   #@@@@ STEP 3 [Adding clusters to map]
   #@@@@ 
   sTime = time.time()
   Map.addHalos2Map(Halos)
   eTime = time.time()-sTime
   print "Eplased time to add clusters to map: ", eTime, " s"


   ##################################################################
   #@@@@ STEP 4 [Saving maps]
   #@@@@ 
   sTime = time.time()
   Map.saveMapPic()
   Map.saveMapFits()
   eTime = time.time()-sTime
   print "Eplased time to save map: ", eTime, " s"


