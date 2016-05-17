import warnings

try:
    ImportWarning
except NameError:
    class ImportWarning(Warning):
        pass

try:
    from XCAT_Xray_Map_Class import Surface_Brightness_Class 
except ImportError:
    warnings.warn("Warning: Cannot import XCat_Xray_Map_Class",
                  category=ImportWarning)



