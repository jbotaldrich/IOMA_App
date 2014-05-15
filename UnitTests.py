# -*- coding: utf-8 -*-
"""
Created on Sun Nov 03 14:36:39 2013

@author: aldr699
"""
import unittest
import IomaModel
from IomaModel import IOMAModel


class TestModelIOMAModel(unittest.TestCase):
    
    filename = "C:/Users/aldr699/Documents/MATLAB/Cobra/cobra/testing/testSBML/Ec_iJR904.xml"    
    
    def setUp(self):
        self.seq = range(10)
        
    def test_model(self):
        model = IomaModel.IOMAModel(self.filename)
        print(model.getS())
        print(model.getb())
        print(model.getGenes())
        print(model.getMetabolites())
        print(model.getRxns())
        
        
        
if __name__ == '__main__':
    unittest.main()