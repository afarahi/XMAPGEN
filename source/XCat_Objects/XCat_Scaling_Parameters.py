from XCat_Utilities  import read_data_float
import os.path

class Temprature_scaling:

   def __init__(self, label):

      fname = './parameters/Models/Txm/' + label + '_parameters.xml'

   

      if os.path.isfile(fname) == False:
         print fname, " does not exists it uses Tx scaling default parameters."
         fname = './parameters/Models/Txm/default_parameters.xml'

      #Parameters
      self.Norm    = read_data_float(tag_name = 'a',file_name = fname)
      self.M_slope = read_data_float(tag_name = 'M_slope',file_name = fname)
      self.E_slope = read_data_float(tag_name = 'E_slope',file_name = fname)

      self.M_p = read_data_float(tag_name = 'M_p',file_name = fname)
      self.z_p = read_data_float(tag_name = 'z_p',file_name = fname)

      self.sig  = read_data_float(tag_name = 'sig',file_name = fname)




class Luminocity_scaling:

   def __init__(self, label):

      fname = './parameters/Models/Lxm/' + label + '_parameters.xml'
   
      if os.path.isfile(fname) == False:
         print fname, " does not exists it uses Lx scaling default parameters."
         fname = './parameters/Models/Lxm/default_parameters.xml'

      #Parameters
      self.Norm    = read_data_float(tag_name = 'a',file_name = fname)
      self.M_slope = read_data_float(tag_name = 'M_slope',file_name = fname)
      self.E_slope = read_data_float(tag_name = 'E_slope',file_name = fname)

      self.M_p = read_data_float(tag_name = 'M_p',file_name = fname)
      self.z_p = read_data_float(tag_name = 'z_p',file_name = fname)

      self.sig  = read_data_float(tag_name = 'sig',file_name = fname)







