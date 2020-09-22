# -*- coding: utf-8 -*-
"""
Created on Thu Sep 26 11:06:23 2019
@author: Thomas Bernard (Python 3.6)
Description: misfit functions to access robustness of model results compare to observed data
"""

import numpy as np

def misfit_function_v1(result_Hr, Hr_max, Hr_t0, result_Hb, Hb_max, Hb_t0, result_wb, subsi_t23, result_etar, rhoc, rhom, rhoi, Lr, Lb, Lambda):
    
    x_range = np.linspace(0.,Lr, 1000)
    x_basin = np.linspace(Lr,(Lr+Lb) , 1000)
    a_range = x_range + Lr                      
    b_range = Lr-x_range                        
    b_basin = x_basin-Lr                        
    a_basin = 2*Lr+b_basin
    
    ## Calculation of range topographic misfit
    topo_r_misfit1 = abs(Hr_max-max(result_Hr))
    topo_r_misfit2 = abs(Hr_t0-result_Hr[6599])
    topo_r_misfit = topo_r_misfit1*(1/1) + topo_r_misfit2*(1/1)    
        
    topo_r_misfit1_norm = (np.sqrt(((Hr_max-max(result_Hr))/1)**2)/(Hr_max/1+np.sqrt(((Hr_max-max(result_Hr))/1)**2)))
    topo_r_misfit2_norm = (np.sqrt(((Hr_t0-result_Hr[6599])/1)**2)/(Hr_t0/1+np.sqrt(((Hr_t0-result_Hr[6599])/1)**2)))
    topo_r_misfit_norm = topo_r_misfit1_norm*(1/2) + topo_r_misfit2_norm*(1/2)
    
    topo_r_misfit_norm2 = (topo_r_misfit1/(Hr_max*0.1))*(1/2) + (topo_r_misfit2/(Hr_t0*0.1))*(1/2)
    
    ## Calculation of basin topographic misfit
    topo_b_misfit1 = abs(Hb_max-max(result_Hb))
    topo_b_misfit2 = abs(Hb_t0-result_Hb[6599])
    topo_b_misfit = topo_b_misfit1*(1/1) + topo_b_misfit2*(1/1)
    
    topo_b_misfit1_norm = (np.sqrt(((Hb_max-max(result_Hb))/1)**2)/(Hb_max/1+np.sqrt(((Hb_max-max(result_Hb)/1))**2)))
    topo_b_misfit2_norm = (np.sqrt(((Hb_t0-result_Hb[6599])/1)**2)/(Hb_t0/1+np.sqrt(((Hb_t0-result_Hb[6599])/1)**2)))
    topo_b_misfit_norm = topo_b_misfit1_norm*(1/2) + topo_b_misfit2_norm*(1/2)
    
    topo_b_misfit_norm2 = (topo_b_misfit1/(Hb_max*0.1))*(1/2) + (topo_b_misfit2/(Hb_t0*0.1))*(1/2)
    
    ## Calculation of basin subsidence misfit    
    #wxb_x0 = -(rhoc*result_etar[4300]/(2*(rhom-rhoi)))*(np.exp(-Lambda*a_basin[0])*np.cos(Lambda*a_basin[0])-np.exp(-Lambda*b_basin[0])*np.cos(Lambda*b_basin[0]));
    
    #subs_misfit1 = abs(subsi_t23-wxb_x0)
    subs_misfit2 = abs(subsi_t23-result_wb[4300])
    subs_misfit = subs_misfit2
    
    #subs_misfit1_norm = (np.sqrt((subsi_t23-wxb_x0)**2)/(subsi_t23+np.sqrt((subsi_t23-wxb_x0)**2)))
    subs_misfit2_norm = (np.sqrt(((subsi_t23-result_wb[4300])/1)**2)/(subsi_t23/1+np.sqrt(((subsi_t23-result_wb[4300])/1)**2)))
    subs_misfit_norm = subs_misfit2_norm
    
    subs_misfit_norm2 = (subs_misfit2/(subsi_t23*0.1))*(1/1)
  
    ## Calculation of global misfit
    global_misfit = topo_r_misfit*(1/1) + subs_misfit*(1/1) + topo_b_misfit*(1/1)
    global_misfit_norm = topo_r_misfit_norm*(1/3) + subs_misfit_norm*(1/3) + topo_b_misfit_norm*(1/3)
    global_misfit_norm2 = topo_r_misfit_norm2*(1/3) + subs_misfit_norm2*(1/3) + topo_b_misfit_norm2*(1/3)
    
    return topo_r_misfit_norm, topo_b_misfit_norm, subs_misfit_norm, global_misfit, global_misfit_norm,  global_misfit_norm2

def misfit_function_v2(result_Hr, Hr_max, Hr_t0, result_Hb, Hb_max, Hb_t0, result_wb, subsi_t23, result_etar, rhoc, rhom, rhoi, Lr, Lb, Lambda):
    
    topo_r_misfit1 = ((Hr_max - max(result_Hr))/(Hr_max*0.1))**2
    topo_r_misfit2 = ((Hr_t0 - result_Hr[-1])/(Hr_t0*0.1))**2
    
    topo_b_misfit1 = ((Hb_max - max(result_Hb))/(Hb_max*0.1))**2
    topo_b_misfit2 = ((Hb_t0 - result_Hb[-1])/(Hb_t0*0.1))**2
    
    depth_b_misfit = ((subsi_t23 - result_wb[4300])/(subsi_t23*0.1))**2
    
    global_misfit = topo_r_misfit1 + topo_r_misfit2 + topo_b_misfit1 + topo_b_misfit2 + depth_b_misfit
    
    tolerance = np.sqrt(2*5)
    
    return global_misfit, tolerance
    