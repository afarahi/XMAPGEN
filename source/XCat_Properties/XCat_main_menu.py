import os
from XCat_Utilities import Read_Integer_Input

def Main_menu():
    os.system('clear')
    print " Main Menu : "
    print " (1)  Loading new input parameters. "
    print " (2)  Adding new halos catalog. " 
    print " (3)  Tranform halos into XCat prefered coordinate. " 
    print " (4)  Solve existing data for Lx, Tx, and flux. "
    print " (5)  Creating and adding new AGN Catalog. "
    print " (6)  Save existing data. "
    print " (7)  Plot Menue. "
    print " (0)  Exit. "
    print "-------------------------------------------------------"   
    ans = Read_Integer_Input(" Please enter a number : ")
    return ans

