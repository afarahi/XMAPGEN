class Solvers_Class():

   def __init__(self):      
      pass

   def ProperDistanceTabulate(self,Input_Param,z_max):
      from XCat_Distance_Solver import Proper_Distance_Tabulate
      Proper_Distance_Tabulate(Input_Param,z_max)

   def LxTxSolver(self,Halos):
      from XCat_LxTx_Solver import LxTx_Solver
      solver = LxTx_Solver()
      solver.solve(Halos)
      #LxTx_Solver(Halos)


