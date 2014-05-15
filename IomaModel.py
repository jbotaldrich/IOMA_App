import numpy as np
import cobra.io.sbml as c_iosbml
import cobra.core as c_core

#class IomaModel:

#S = np.matrix([])
#C = np.Matrix([])    
class IOMAModel(object):        
    model = None
    def __init__(self, sbmlfile):
        temp = c_iosbml.create_cobra_model_from_sbml_file(sbmlfile)
        self.model = c_core.ArrayBasedModel(temp)
        
    def getS(self):
        return self.model.S
        
    def getb(self):
        return self.model.b 
     
    def getlb(self):
        return self.model.lower_bounds
        
    def getub(self):
        return self.model.upper_bounds
        
    def getGenes(self):
        return self.model.genes
        
    def getRxns(self):
        return self.model.reactions
        
    def getMetabolites(self):
        return self.model.metabolites

    def setUpper(previous, newer, self):
        self.model.upper_bounds[self.findIndices(self.model.upper_bounds, previous)] = newer
        
    def setLower(previous, newer, self):
        self.model.lower_bounds[self.findIndices(self.model.lower_bounds, previous)] = newer

    def findIndices(my_list, toFind):
        return[i for i, x in enumerate(my_list) if x == toFind]     

#    def _convertSBMLToMatrix(sbml):
        #do stuff
       # result = array(sbml)
    #    return result
        
    #def getS():
      #  return self.S
        
    #def getC():
  #      return C
        
   # def metabolite_num():
   #     return S.shape[0]
    
 #   def reaction_num():
  #      return self.S.shape[1]
    
    