class Input_Parameters:

   def __init__(self):

      from XCat_Utilities  import read_data_string, read_data_float, read_data_int, read_data_bool 
      import sys 

      #Auto_mode
      self.Auto_mode = read_data_bool(tag_name = 'Auto_mode',file_name = 'parameters/Input_Parameters.xml')

      #Cosmology         
      self.h_0 = read_data_float(tag_name = 'Hubble_Parameter',file_name = 'parameters/Cosmological_Parameters.xml')
      self.Omega_DE = read_data_float(tag_name = 'Omega_DE',file_name = 'parameters/Cosmological_Parameters.xml')
      self.Omega_M = read_data_float(tag_name = 'Omega_M',file_name = 'parameters/Cosmological_Parameters.xml')
      self.Omega_b = read_data_float(tag_name = 'Omega_b',file_name = 'parameters/Cosmological_Parameters.xml')
      self.Omega_R = read_data_float(tag_name = 'Omega_R',file_name = 'parameters/Cosmological_Parameters.xml')
      self.Omega_k = read_data_float(tag_name = 'Omega_k',file_name = 'parameters/Cosmological_Parameters.xml')
      self.sigma_8 = read_data_float(tag_name = 'sigma_8',file_name = 'parameters/Cosmological_Parameters.xml')
      self.w = read_data_float(tag_name = 'w',file_name = 'parameters/Cosmological_Parameters.xml')
      self.ns = read_data_float(tag_name = 'ns',file_name = 'parameters/Cosmological_Parameters.xml')

      #Models
      self.Lxm_mod = read_data_string(tag_name = 'Lxm_model',file_name = 'parameters/Input_Parameters.xml')
      self.Txm_mod = read_data_string(tag_name = 'Txm_model',file_name = 'parameters/Input_Parameters.xml')
      xray_band_option = read_data_int(tag_name = 'XRay_band_Option',file_name = 'parameters/Input_Parameters.xml')
   
      from XCat_Objects import Xray_Band_Class
      from XCat_Objects import Temprature_scaling, Luminocity_scaling

      self.xray_band = Xray_Band_Class(option=xray_band_option)

      try:
         self.LxScaling = Luminocity_scaling(self.Lxm_mod)
      except KeyError:
         print "ERROR: Input Lxm Parameter method does not exist!"
         print "Please change XML file and try again."
         raw_input("Press enter to exit! ")
         sys.exit(2)
      try:
         self.TxScaling = Temprature_scaling(self.Txm_mod)
      except KeyError:
         print "ERROR: Input Txm Parameter method does not exist!"
         print "Please change XML file and try again."
         raw_input("Press enter to exit! ")
         sys.exit(2)


      self.SB_beta     = read_data_float(tag_name = 'SB_beta',file_name = 'parameters/Input_Parameters.xml')
  
      #Exctra Parameters
      self.rTL         = read_data_float(tag_name = 'rTL_correlation',file_name = 'parameters/Input_Parameters.xml')

      #Limits
      self.Flim        = read_data_float(tag_name = 'F_limit',file_name = 'parameters/Input_Parameters.xml')
      self.Mlim        = read_data_float(tag_name = 'Mass_limit',file_name = 'parameters/Input_Parameters.xml')

      
      print "COSMOLOGY and INPUTS are assigned SUCCESSFULLY."


