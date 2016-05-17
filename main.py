#! /usr/bin/env python
import sys
sys.path.insert(0,sys.path[0]+'/source')
sys.path.insert(0,sys.path[0]+'/XCat_Objects')
sys.path.insert(0,sys.path[1]+'/XCat_Calculator')
sys.path.insert(0,sys.path[2]+'/XCat_Utilities')
sys.path.insert(0,sys.path[3]+'/XCat_Solver')
sys.path.insert(0,sys.path[4]+'/XCat_Properties')
sys.path.insert(0,sys.path[5]+'/XCat_Models')
sys.path.insert(0,sys.path[6]+'/XCat_Models/AGN_Models')
sys.path.insert(0,sys.path[7]+'/XCat_Models/Surface_Brightness_Models')
sys.path.insert(0,sys.path[8]+'/XCat_Models/Event_Map_Models')

#print sys.path
FLUX_MODE = False

if FLUX_MODE:
   import matplotlib
   matplotlib.use('Agg')

from mainPipeline import makingXRayCatalogs,\
                         makingXRayRealization

print "(1) HALOS MODE"
print "(2) MAP MODE"
print "(3) TEST MODE"

ans = int(sys.argv[1]) #int(raw_input("Please enter the mode number : "))

if (ans == 1):
   makingXRayCatalogs()
if (ans == 2):
   makingXRayRealization()
#if (ans == 3):
#    test_tasks()
else:
    pass
