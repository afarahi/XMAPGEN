import warnings

try:
    ImportWarning
except NameError:
    class ImportWarning(Warning):
        pass

from XCat_logo            import Print_logo
from XCat_main_menu       import Main_menu
