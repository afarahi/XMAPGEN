import warnings

try:
    ImportWarning
except NameError:
    class ImportWarning(Warning):
        pass

try:
    from XCAT_Event_Map_Class import Event_Class 
except ImportError:
    warnings.warn("Warning: Cannot import XCat_Event_Map_Class",
                  category=ImportWarning)



