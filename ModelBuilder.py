# -*- coding: utf-8 -*-
import numpy as np
import scipy as sp
import math
import IomaModel


MAX_FLUX = 100.0
MIN_FLUX = -100.0

""" Main function for generating the model
model, is the input stoichiometric model
exp_mat, the experiment matrix
forward_met_rxn, fill this out
reverse_met_rxn, fill this out
rxn_idx, index of the reactions
ko_rxns, knocked out rxns
ex_flux, external flux 
ex_rxns, external reactions
GR, growth rate
"""


def BuildModel(model, exp_mat, forward_met_rxn,
               reverse_met_rxn, rxns_idx, ko_rxns,
               ex_flux, ex_rxns, GR):
   met_num = model.metabolite_num()
   rxn_num = model.reaction_num()
   ko_exp_num = exp_mat.shape[0]
   biomass_idx = findBiomass(model.GetC)
   
   reaction_with_v_max = findSumZeros(forward_met_rxn)
   v_max_num = len(reaction_with_v_max)
   
   ko_num = np.ones(v_max_num, 1) * ko_exp_num
   
   variables_number = 2 * v_max_num + ko_exp_num * (rxn_num + 2 * v_max_num)
   ko_vars = (rxn_num + 2* v_max_num)
   
   temp_size1 = 0
   temp_size2 = 0
   m = None
   
   for i in range(0, ko_exp_num):
       MBConMat = BuildMassBalanceConMatrix(i, met_num, temp_size1, v_max_num, 
                                            model, ko_exp_num, ko_vars)
                                            
       MMConMat = BuildMichaelisMentenConMatrix(i, exp_mat, temp_size1, v_max_num,
                                                model, forward_met_rxn, 
                                                reverse_met_rxn, rxn_num,
                                                rxns_idx)
       temp_size2 = temp_size2 + 1
       QPConMat = BuildVarianceConMatrix(i, exp_mat, temp_size1, temp_size2, 
                                         v_max_num, model, forward_Met_Rxn, 
                                         reverse_Met_Rxn, rxn_Num, rxns_idx, 
                                         ko_exp_num, ko_vars, ko_num)
                                         
       if m == None:
           m = np.vstack([MBconMat, MMConMat, QPConMat])
       else:
           m = np.vstack([m,MBConMat, MMConMat, QPConMat])
           
    count_y = rxn_num + 3*v_max_num #y variables
    count_v = 2*v_max_num #v variables
    count_e = rxn_num + 2*v_max_num #epsilon variables
    count_b = 2*v_max_num #biomass variables

    ub = np.zeros([variables_number,1])
    lb = np.zeros([variables_number,1])
    #v_max bounds
    
    ub[0:v_max_num] = MAX_FLUX
    lb[v_max_num+1:2*v_max_num+1] = MIN_FLUX
          
    
    
    #v bounds
    model.setLower(1000.0, MIN_FLUX)
    model.setUpper(1000.0, MAX_FLUX)  
          

    y_idx = []
    v_idx = []
    e_idx = []
    b_idx = []

    for i in range(0,ko_exp_num):
        y_idx = np.concatenate(y_idx, np.arange(count_y+1, count_y+v_max_num))
        e_idx = np.concatenate(e_idx, np.arange(count_e+1, count_e+v_max_num))
        v_idx = np.concatenate(v_idx, np.arange(count_v+1, count_v+rxn_num))
        b_idx = np.concatenate(b_idx, count_b+biomass_idx)
        

        ub[count_v+1:count_v+rxn_num] = model.getub()
        lb[count_v+1:count_v+rxn_num] = model.getlb()
        %exchange reaction bounds
        ub(ex_rxns'+count_v.*ones(length(ex_rxns),1)) = ex_flux(:,i);
        lb(ex_rxns'+count_v.*ones(length(ex_rxns),1)) = ex_flux(:,i);
        
        %set KO bounds
        ko = ko_rxns(i,:);
        ko = ko(ko~=0);
        ub(ko'+count_v.*ones(length(ko),1)) = 0;
        lb(ko'+count_v.*ones(length(ko),1)) = 0;
        
        %update counters
        count_y = count_y + rxn_num + 2*v_max_num;
        count_v = count_v + rxn_num + 2*v_max_num;
        count_e = count_e + rxn_num + 2*v_max_num;
        count_b = count_b + rxn_num + 2*v_max_num;
        

def findBiomass(CVector):
    try:
        index = np.where(CVector == 1)
    except ValueError:
        print "Incorrect biomass array length, should be 1"
        
    return index
    
def findSumZeros(inArray):
    return nonzero(inArray.sum(axis=1))
    

def BuildMassBalanceConMatrix(iteration, met_num, temp_size1,
         v_max_num, model, ko_exp_num, ko_vars):
    m_b_constraints = np.hstack(np.zeros([met_num, 2 * v_max_num]), 
                                np.zeros([met_num, temp_size1]),
                                model.getS(),
                                np.zeros([met_num, 2 * v_max_num]),
                                np.zeros([met_num, (ko_exp_num-iteration)*ko_vars]))
    return m_b_constraints
         
         
   
         
def BuildMichaelisMentenConMatrix(iteration, exp_mat, temp_size1, v_max_num, 
                                  model, met_frw_mat, met_bck_mat, rxn_num, 
                                  rxns_idx, ko_exp_num, ko_vars):
    exp_e = exp_mat[iteration,:]
    exp_ce = exp_mat[iterationm,:] * met_frw_mat[iteration,:]
    exp_de = exp_mat[iterationm,:] * met_bck_mat[iteration,:]
    s = np.argmax(rxns_idx)
    one_vec = np.ones([v_max_num, 1])
    t1 = arange(1.0, v_max_num, v_max_num - 1)
    m2f = np.column_stack([t1, t1, exp_ce])* (-1);
    m2b = np.column_stack([t1, t1, exp_de])
    m_zeros = np.zeros([v_max_num, temp_size1])
    m3 = np.column_stack([t1, rxns_idx, one_vec, np.zeros([v_max_num, rxn_num - s])])
    m4 = np.column_stack([t1, t1, exp_e]) * (-1.0)
    m5 = np.zeros([v_max_num, v_max_num])
    MM_constraints = np.hstack([m2f, m2b, m_zeros, m3, m4, m5, np.zeros([v_max_num, (ko_exp_num - iteration)*ko_vars])])
    
    return MM_constraints;
         
def BuildVarianceConMatrix(iteration, exp_mat,  temp_size1,  temp_size2,
                           v_max_num, model, met_frw_mat, met_bck_mat,  
                           rxn_num, rxns_idx,  ko_exp_num,  ko_vars, ko_num):

    e1_vec = (ko_num - np.ones([v_max_num,1])) / ko_num
    e2_vec = np.ones([v_max_num,1]) / ko_num

    y_vec = np.ones([v_max_num, 1]) / (math.sqrt(ko_exp_num))
#
    m_zeros = np.zeros([v_max_num, v_max_num])
#
    m_vmax = np.zeros([v_max_num, 2 * v_max_num])
#
    m_s = np.zeros([v_max_num, rxn_num])
#
    m_e1 = np.dot(np.identity(v_max_num),(e1_vec)) * (-1.0);
#
    m_y = np.dot(np.identity(v_max_num),y_vec);
#
    m_e2 = np.dot(np.identity(v_max_num),e2_vec);
#
    m_temp1 = np.hstack([m_s, m_e1, m_y]);
#
    m_temp2 = np.stack([m_s, m_e2, m_zeros]);
#
    QP_constraints = m_vmax;
    
    for i in range(1, temp_size2):
        QP_constraints = np.hstack([QP_constraints, m_temp2])
    
    QP_constraints = np.hstack([QP_constraints, m_temp1])
    j = ko_exp_num - iteration;
    while (j >= 1):
        QP_constraints = np.hstack([QP_constraints, m_temp2]);
        j = j - 1;

    return QP_constraints
}     
 
         