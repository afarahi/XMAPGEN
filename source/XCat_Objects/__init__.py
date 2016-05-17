import warnings

try:
    ImportWarning
except NameError:
    class ImportWarning(Warning):
        pass


try:
    from XCat_Conversions_Constants import *
    from XCat_Global_Var import *
except ImportError:
    warnings.warn("Warning: Cannot import ...",
                  category=ImportWarning)

try:
    from XCat_Input_Parameters         import Input_Parameters
except ImportError:
    warnings.warn("Warning: Cannot import ...",
                  category=ImportWarning)

try:
    from XCat_Halos_Class              import Halos_Class, Map_Halos_Class
except ImportError:
    warnings.warn("Warning: Cannot import XCat_Halos_Data",
                  category=ImportWarning)

try: 
    from XCat_AGN_Class               import AGN_Class
except ImportError:
    warnings.warn("Warning: Cannot import XCat_AGN_Data",
                  category=ImportWarning)

try:
    from XCat_Map_Class                import Map_Class
except ImportError:
    warnings.warn("Warning: Cannot import XCat_Halos_Data",
                  category=ImportWarning)

try:
    from XCat_Solvers_Class         import Solvers_Class 
except ImportError:
    warnings.warn("Warning: Cannot import XCat_Solvers_Class",
                  category=ImportWarning)

try:
    from XCat_Xray_Band_Class       import Xray_Band_Class
except ImportError:
    warnings.warn("Warning: Cannot import XCat_Xray_Band_Option_Class",
                  category=ImportWarning)


try:
    from XCat_Scaling_Parameters import Temprature_scaling, Luminocity_scaling
except ImportError:
    warnings.warn("Warning: Cannot import ...",
                  category=ImportWarning)



try:
    from XCat_Plots_Class         import Plots_Class 
except ImportError:
    warnings.warn("Warning: Cannot import XCat_Solvers_Class",
                  category=ImportWarning)


try:
    from XCat_Halos_Data_type2         import ( Halo_Object_type2 )
except ImportError:
    warnings.warn("Warning: Cannot import XCat_Halos_Data_type2",
                  category=ImportWarning)

try:
    from XCat_Halos_General_Properties import Halo_Status_Class
except ImportError:
    warnings.warn("Warning: Cannot import XCat_Halos_Sample",
                  category=ImportWarning)

