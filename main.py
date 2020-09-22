"""
Created on Wed Mar 27 10:59:04 2019
@author: Thomas Bernard and Mark Naylor (Python 3.6)
Description:
Main file for post-orogenic development of a mountain range and its foreland basin based on
Tucker and Van der Beek Basin Research (2013) 24, 241–259, doi: 10.1111/j.1365-2117.2012.00559.x
and modified for the Northern Pyrenees case with inverse modelling Monte Carlo sampling algorithm 
"""
 
#import matplotlib;matplotlib.use("Agg")
from __future__ import division
#from scipy import stats
from datetime import datetime
import numpy as np
import sys
import functions_publication_figures
import functions_plotting_results
import function_TuckerCode
import functions_misfit
#from matplotlib import gridspec
#from matplotlib.colors import LogNorm
#import random as rd
import matplotlib.pyplot as plt
import matplotlib.colors
import matplotlib.ticker as ticker
import importlib
#import pandas as pd
#import seaborn as sns
#import matplotlib.patches as patches
#import matplotlib.lines as lines
#from matplotlib import gridspec
#from matplotlib.tri import Triangulation, TriAnalyzer, UniformTriRefiner
#from scipy.interpolate import griddata
#import scipy.ndimage
#from scipy.ndimage.filters import gaussian_filter
#from mpl_toolkits.mplot3d import axes3d

#%%=========================================================================%%#
#=================== DEFINE AND INITIATE FIXE PARAMETER ======================#

start_time = datetime.now()

## Set initial conditions:
Hr0 = 0.                    # range height
Hb0 = 0.                    # basin height
etar0 = 0.                  # range thickness
etab0 = 0.                  # basin thickness
tstart = 0.                 # start

## Define general variables:
Lt = 150000.                    # range width, m
alpha = 0.30                    # retro- vs pro-side ratio (*0.01)
alpha_min = 20                  # mimimum retro-side ratio
alpha_max = 40                  # maximum retro-side ratio 
Lr = Lt*alpha                   # retro-side range witdh
porosity = 0.4                  # porosity
dt = 10000.                     # Time step size, yr
tt = 66e6                       # total time , yr
nt = int( tt/dt )               # number of iterations

## Define flexure-related parameters: 
rhoc = 2700.                            # crust density, kg/m3
rhom = 3300.                            # mantle density, kg/m3
rhoi = 1620.                            # infill density, kg/m3
rhos = rhoc*(1-porosity)                # basin infill density (sediments)
Te = 21930.                             # North effective elastic thickness, m (Brunet et al., 1986)
Te_min = 15000.                         # minimum effective elastic thickness
Te_max = 40000.                         # maximum effective elastic thickness
E = 1.e11                               # Young's modulus, Pa
nu = 0.25                               # Poisson's ratio
g = 9.81                                # gravity constant

## Define key timing parameters:
T1 = 0                              # start of the model
T2 = 10e6
T3 = 25e6
T4 = 32e6                           # increase of convergence rate
T5 = 43e6                           # decrease convergence rate
T6_min = 45.5e6                     # minimum end of convergence
T6_max = 50.5e6                     # maximum end of convergence
T7 = 66e6                           # end of the model

T1_time_steps = int(round(T1/dt)); T2_time_steps = int(round(T2/dt)); T3_time_steps = int(round(T3/dt)); T4_time_steps = int(round(T4/dt)); T5_time_steps = int(round(T5/dt)); T7_time_steps = int(round(T7/dt));     # transform key timing for model calculation

## Define tectonic parameters:
Tc = 30000.                                                     # Convergence thickness, m
f = 0.0                                                         # fractional thickness of proximal basin (fraction of mean)
V_initial = 0                                                   # initial convergence velocity
V_max = 0.004                                                   # max convergence velocity
V1 = 0.0000; V2 = 0.0032; V3= 0.0024; V4 = 0.004; V5 = 0.0002   # realistic convergence velocity 

## Define transport coefficient parameters:
kappar_min = 100                                # minimum range transport coefficient
kappar_max = 5000                               # maximum range transport coefficient
kappab_min = 1000                               # minimum basin transport coefficient
kappab_max = 50000                              # maximum basin transport coefficient

## Define epei-orogeny parameters:
epei_crest = -0.004                       # uplift rate imposed at range crest (left side of "range box")
epeistart = int( (5.e6+3.*3.e5)/dt )
epeiend = int( (5.e6+40.*3.e5)/dt )

## State variables:
Hr = Hr0*np.ones(nt)            # This vector stores the range height at each time step ...
Hb = Hb0*np.ones(nt)            # ... same for basin height ...
etar = etar0*np.ones(nt)        # ... and range thickness ...
etab = etab0*np.ones(nt)        # ... and basin sediment thickness ...
wr = etar0*np.ones(nt)          # ... andflexure under the range ...
wb = etar0*np.ones(nt)          # ... and under the basin ...
Fe = np.zeros(nt)               # ... and erosional fluxes from the range to the basin ...
Fb = np.zeros(nt)               # ... and erosional fluxes from the basin to the base level ...
Fr = np.zeros(nt)               # ... and tectonic fluxes to the range ...
Hbl = np.zeros(nt)              # ... and base level thickness.

## Dataset for misfit calculations
Hr_max = 2000                       # mean topographic elevation of the range at 23 Ma
Hr_t0 = 1500                        # mean topographic elevation of the range at 0 Ma
Hb_t0 = 250                         # mean topographic elevation of the basin at 0 Ma
Hb_max = 350                        # mean maximum topographic elevation of the basin
subsi_t23 = 5000                    # maximum subsidence of the basin at 23 Ma
global_misfit_reference = 2.0       # maximum misfit for model acceptance

## Flexure equations:
#wxr = (rhoc*etar[0]/(2*(rhom-rhoi)))*(2-np.exp(-Lambda*b_range)*np.cos(Lambda*b_range)-np.exp(-Lambda*a_range)*np.cos(Lambda*a_range))
#wxb = -(rhoc*etar[0]/(2*(rhom-rhoi)))*(np.exp(-Lambda*a_basin)*np.cos(Lambda*a_basin)-np.exp(-Lambda*b_basin)*np.cos(Lambda*b_basin))

#%%=========================================================================%%#
#========================== INITIATE RESULTS STORING =========================#
                                        
result_Hr_list = np.empty((0, nt))                                            # Initiate an array for the range elevation accepted ...
result_Hb_list = np.empty((0, nt))                                            # ... and for the basin elevation accepted ...
result_etar_list = np.empty((0, nt))                                          # ... and for the range thickness accepted ...
result_etab_list = np.empty((0, nt))                                          # ... and for the basin thickness accepted ...
result_wr_list = np.empty((0, nt))                                            # ... and for range deflection accepted ...
result_wb_list = np.empty((0, nt))                                            # ... and for basin deflection accepted ...
result_Fe_list = np.empty((0, nt))                                            # ... and for the range sedimentary fluxes accepted ...
result_Fb_list = np.empty((0, nt))                                            # ... and for the basin sediemntary fluxes accepted ...
result_Fr_list = np.empty((0, nt))                                            # ... and for the tectonic flux accepted ...
result_time_list = np.empty((0, nt))                                          # ... and for the time accepted ...

result_global_misfit_list = np.empty((0))                                     # ... and for the global misfit accepted ...
result_global_misfit2_list = np.empty((0))
result_global_misfit3_list = np.empty((0))
result_kappar_list = np.empty((0))                                            # ... and for the ranhe transport coefficient accepted ...
result_kappab_list = np.empty((0))                                            # ... and for the basin transport coefficient accepted ...
result_Te_list = np.empty((0))                                                # ... and for the lithosphere elastic thickness accepted ...
result_alpha_list = np.empty((0))                                             # ... and for the alpha accepted ...
result_T4_list = np.empty((0))                                                # ... and for the time of convergence decrease accepted ...

#%%=========================================================================%%#
#====================== MC CALCULATION - FIX ITERATION =======================#
#
#i = 0;
#MC_iteration = 1000
#topo_r_misfit1, topo_r_misfit2, topo_r_misfit_list = (np.ones(MC_iteration) for k in range(3))     # initiate a list for the range topographic misfit ...
#topo_b_misfit1, topo_b_misfit2, topo_b_misfit_list = (np.ones(MC_iteration) for k in range(3))     # ... and for the basin topographic misfit ... 
#subs_misfit_list = np.ones(MC_iteration)                                                           # ... and for the subsidence misfit ...
#global_misfit_list = np.ones(MC_iteration) 
#
### Calculate arrays of random numbers (lithosphere elastic thickness, range and basin diffusivity coefficient, decrease time, alpha factor) for the MCMC computation
#Te_list = functions_MonteCarlo_sampling.lithosphere_elastic_thickness_MC_sampling(MC_iteration, Te_min, Te_max)
#kappar_list, kappab_list = functions_MonteCarlo_sampling.diffusivity_MC_sampling(MC_iteration, kappar_min, kappar_max, kappab_min, kappab_max)
#T4_time_steps_list = functions_MonteCarlo_sampling.decreasetime_MC_sampling(MC_iteration, dt, T4_min, T4_max)
##alpha_list = functions_MonteCarlo_sampling.alpha_MC_sampling(MC_iteration, alpha_min, alpha_max); alpha_list = alpha_list*0.01
#                         
#print("----------------------------------------------------------------------")
#print("magic start here")
#print("----------------------------------------------------------------------")
#
#while i < MC_iteration:
#
#  
#  ## define alpha factor parameter and update length of the wedge
#  #alpha = alpha_list[i]
#  #Lr = Lt*alpha
#  
#  ## define uplift deceleration parameter and update the derived tectonic parameters
#  T4_time_steps = T4_time_steps_list[i]                                            
#  
#  Vt = np.zeros(nt)
#  Vt[T1_time_steps:T2_time_steps] = V1
#  Vt[T2_time_steps:T3_time_steps] = V2
#  Vt[T3_time_steps:T6_time_steps] = V3
#  
#  V_T3 = V2
#  for v in range(T3_time_steps,T4_time_steps):
#    Vt[v] = V_T3 - (V2-V3)/(T4_time_steps-T3_time_steps)
#    V_T3 = Vt[v]
#       
#  mut = Vt*Tc/Lt
#  
#  ## define lithosphere elastic thickness parameter and update derived flexure parameters
#  Te = int(Te_list[i])
#  
#  
#  D = E*Te**3/(12*(1-nu**2))
#  Lambda = ((rhom-rhoi)*g/(4*D))**(1./4.)
#  tlaml = 2*Lambda*Lr
#  Psir = drat*np.exp(-tlaml)*(np.exp(tlaml)*(-1+2*tlaml)+np.cos(tlaml)-np.sin(tlaml))/(2*tlaml)
#  Lb = 3*np.pi/(4*Lambda)-Lr
#  Llm = Lambda*Lr
#  Psib = (rhoc/(2*(rhom-rhoi)))*(1/(-4*Llm+3*np.pi))*2*np.exp(-2*Llm-pi34)*(np.sqrt(2)*np.exp(Llm)*(-1+np.exp(2*Llm))*np.cos(Llm)+np.exp(pi34)*(np.exp(2*Llm)-np.cos(2*Llm)+np.sin(2*Llm)))  #average deflexion under the basin
# 
#  ## define response time parameters and update derived response time parameters
#  tr = (Lr**2)/kappar_list[i]; tb = (Lb**2)/kappab_list[i]                                
#  trt = np.zeros(nt); tbt = np.zeros(nt)
#  trt[T1_time_steps:T6_time_steps] = tr; 
#  tbt[T1_time_steps:T6_time_steps] = tb
#  
#  ## compute the model metric                                          
#  result_Hr, result_Hb, result_etar, result_etab, result_wr, result_wb, result_Fe, result_Fb, result_Fr, time = function_TuckerCode.Tucker_Code(dt, nt, tstart, trt, tbt, mut, Tc, rhos, rhoc, rhoi, Lr, Lb, f, Vt, epeirt, epeibt, Hr, Hb, etar, etab, wr, Psir, wb, Psib, Fe, Fb, Fr)
#  
#  ## compute range topographic misfit
#  topo_r_misfit1[i] = (np.sqrt(((Hr_max-max(result_Hr))/1)**2)/(Hr_max/1+np.sqrt(((Hr_max-max(result_Hr))/1)**2)))
#  topo_r_misfit2[i] = (np.sqrt(((Hr_t0-result_Hr[6599])/1)**2)/(Hr_t0/1+np.sqrt(((Hr_t0-result_Hr[6599])/1)**2)))  
#  
#  topo_r_misfit_list[i] = topo_r_misfit1[i]*(1/2) + topo_r_misfit2[i]*(1/2)
#  
#  ## compute basin topographic misfit
#  topo_b_misfit1[i] = (np.sqrt(((Hb_max-max(result_Hb))/1)**2)/(Hb_max/1+np.sqrt(((Hb_max-max(result_Hb)/1))**2)))
#  topo_b_misfit2[i] = (np.sqrt(((Hb_t0-result_Hb[6599])/1)**2)/(Hb_t0/1+np.sqrt(((Hb_t0-result_Hb[6599])/1)**2)))
#  
#  topo_b_misfit_list[i] = topo_b_misfit1[i]*(0/2) + topo_b_misfit2[i]*(2/2)
#  
#  ## compute subsidence misfit
#  #wxb_x0 = -(rhoc*result_etar[4300]/(2*(rhom-rhoi)))*(np.exp(-Lambda*a_basin[0])*np.cos(Lambda*a_basin[0])-np.exp(-Lambda*b_basin[0])*np.cos(Lambda*b_basin[0]));
#  #subs_misfit[i] = (np.sqrt((subsi_t23-wxb_x0)**2)/(subsi_t23+np.sqrt((subsi_t23-wxb_x0)**2)))
#  subs_misfit_list[i] = (np.sqrt(((subsi_t23-result_wb[4300])/1)**2)/(subsi_t23/1+np.sqrt(((subsi_t23-result_wb[4300])/1)**2)))
#  
#  #subs_misfit[i] = (np.sqrt(((subsi_t23-result_wb[4300])/200)**2))
#  
#  ## compute global misfit
#  global_misfit_list[i] = topo_r_misfit_list[i]*(1/3) + subs_misfit_list[i]*(1/3) + topo_b_misfit_list[i]*(1/3)
#  
#  ## store accepted result from a mimimum misfit
#  if global_misfit_list[i] <= global_misfit_reference:
#    ## store accepted metrics
#    result_Hr_list = np.append(result_Hr_list, result_Hr.reshape((1,nt)), axis=0)
#    result_Hb_list = np.append(result_Hb_list, result_Hb.reshape((1,nt)), axis=0)
#    result_etar_list = np.append(result_etar_list, result_etar.reshape((1,nt)), axis=0)
#    result_etab_list = np.append(result_etab_list, result_etab.reshape((1,nt)), axis=0)
#    result_wr_list = np.append(result_wr_list, result_wr.reshape((1,nt)), axis=0)
#    result_wb_list = np.append(result_wb_list, result_wb.reshape((1,nt)), axis=0)
#    result_time_list = np.append(result_time_list, time.reshape((1,nt)), axis=0)
#    result_Fe_list = np.append(result_Fe_list, result_Fe.reshape((1,nt)), axis=0)
#    result_Fb_list = np.append(result_Fb_list, result_Fb.reshape((1,nt)), axis=0)
#    result_Fr_list = np.append(result_Fr_list, result_Fr.reshape((1,nt)), axis=0)
#    
#    ## store accepted misfit
#    result_global_misfit_list = np.append(result_global_misfit_list, global_misfit_list[i])
#    
#    ## store accepted variable parameters
#    result_kappar_list = np.append(result_kappar_list, kappar_list[i])
#    result_kappab_list = np.append(result_kappab_list, kappab_list[i])
#    result_Te_list = np.append(result_Te_list, Te_list[i])
#    result_T4_list = np.append(result_T4_list, T4_time_steps_list[i])
#    #result_alpha_list = np.append(result_alpha_list, alpha_list[i])
#     
#  ## print the calcul progress (for fix number of iteration)
#  P = round(i/MC_iteration*100,3)                                            
#  if any([P==5, P==10, P==15, P==20, P==25, P==30, P==35, P==40, P==45, P==50, P==55, P==60, P==65, P==70, P==75, P==80, P==85, P==90, P==95, i==MC_iteration-1]):
#    print('percentage progress:', round(i/MC_iteration*100,3), '%, best misfit:', round(min(global_misfit_list),3), '...')
#    
#  ## update indexation
#  i += 1
#
#end_time = datetime.now()

#%%=========================================================================%%#
#================== MC CALCULATION - FIX ACCEPTED MODEL ======================#
importlib.reload(function_TuckerCode)
importlib.reload(functions_misfit)

i = 0
acm = 0; acm_update = 0
MC_result = 10

Te_list = np.empty((0))
T6_time_steps_list = np.empty((0))
kappar_list = np.empty((0))
kappab_list = np.empty((0))
global_misfit_list = np.empty((0))
global_misfit2_list = np.empty((0))
global_misfit3_list = np.empty((0))
topo_r_misfit_list = np.empty((0))
topo_b_misfit_list = np.empty((0))
subs_misfit_list = np.empty((0))

print("----------------------------------------------------------------------")
print("magic start here")
print("----------------------------------------------------------------------")

while acm < MC_result:
    
    ## define uplift deceleration parameter and update the derived tectonic parameters
    T6_time_steps = np.random.randint(low = int(round(T6_min/dt)), high = int(round(T6_max/dt)))
    T6_time_steps_list = np.append(T6_time_steps_list, T6_time_steps)
    
    Vt, mut = function_TuckerCode.Update_tectonic_parameter(nt, T1_time_steps, T2_time_steps, T3_time_steps, T4_time_steps, T5_time_steps, T6_time_steps, T7_time_steps, Tc, Lt, V1, V2, V3, V4, V5)
    
    ## define lithosphere elastic thickness parameter and update derived flexure parameters
    Te = int(np.random.randint(low = Te_min, high = Te_max))
    Te_list = np.append(Te_list, Te)
    
    Psir, Psib, Lb, Lambda, pi34 = function_TuckerCode.Update_flexural_parameter(Te, E, nu, g, Lr, rhom, rhoi, rhoc) 
    
    ## define transport coefficient parameters and update derived response time parameters
    kappar = np.random.randint(low = kappar_min, high = kappar_max); kappab = np.random.randint(low = kappab_min, high = kappab_max)
    kappar_list = np.append(kappar_list, kappar); kappab_list = np.append(kappab_list, kappab)
    
    tr, tb, trt, tbt = function_TuckerCode.Update_response_time_parameter(nt, Lr, Lb, kappar, kappab, T1_time_steps, T7_time_steps)
    
    ## update epeiorogeny parameters
    epeirt, epeibt = function_TuckerCode.Update_epeiorogenic_parameter(nt, pi34, Lr, Lambda, epei_crest, epeistart, epeiend)
    
    ## compute the model metric                                          
    result_Hr, result_Hb, result_etar, result_etab, result_wr, result_wb, result_Fe, result_Fb, result_Fr, time = function_TuckerCode.Tucker_Code(dt, nt, tstart, trt, tbt, mut, Tc, rhos, rhoc, rhoi, Lr, Lb, f, Vt, epeirt, epeibt, Hr, Hb, etar, etab, wr, Psir, wb, Psib, Fe, Fb, Fr)
  
    ## compute misfit for model acceptance
    topo_r_misfit, topo_b_misfit, subs_misfit, global_misfit3, global_misfit, global_misfit2  = functions_misfit.misfit_function_v1(result_Hr, Hr_max, Hr_t0, result_Hb, Hb_max, Hb_t0, result_wb, subsi_t23, result_etar, rhoc, rhom, rhoi, Lr, Lb, Lambda)    
    #global_misfit, tolerance = functions_misfit.misfit_function_v2(result_Hr, Hr_max, Hr_t0, result_Hb, Hb_max, Hb_t0, result_wb, subsi_t23, result_etar, rhoc, rhom, rhoi, Lr, Lb, Lambda)
    
    ## store misfit
    topo_r_misfit_list = np.append(topo_r_misfit_list, topo_r_misfit)
    topo_b_misfit_list = np.append(topo_b_misfit_list, topo_b_misfit)
    subs_misfit_list = np.append(subs_misfit_list, subs_misfit)
    global_misfit_list = np.append(global_misfit_list, global_misfit)
    
    global_misfit2_list = np.append(global_misfit2_list, global_misfit2)
    
    global_misfit3_list = np.append(global_misfit3_list, global_misfit3)
  
    ## store accepted result from a mimimum misfit
    if global_misfit2 <= global_misfit_reference:
        
        ## store accepted metrics
        result_Hr_list = np.append(result_Hr_list, result_Hr.reshape((1,nt)), axis=0)
        result_Hb_list = np.append(result_Hb_list, result_Hb.reshape((1,nt)), axis=0)
        result_etar_list = np.append(result_etar_list, result_etar.reshape((1,nt)), axis=0)
        result_etab_list = np.append(result_etab_list, result_etab.reshape((1,nt)), axis=0)
        result_wr_list = np.append(result_wr_list, result_wr.reshape((1,nt)), axis=0)
        result_wb_list = np.append(result_wb_list, result_wb.reshape((1,nt)), axis=0)
        result_time_list = np.append(result_time_list, time.reshape((1,nt)), axis=0)
        result_Fe_list = np.append(result_Fe_list, result_Fe.reshape((1,nt)), axis=0)
        result_Fb_list = np.append(result_Fb_list, result_Fb.reshape((1,nt)), axis=0)
        result_Fr_list = np.append(result_Fr_list, result_Fr.reshape((1,nt)), axis=0)
        
        ## store accepted misfit
        result_global_misfit_list = np.append(result_global_misfit_list, global_misfit)
        result_global_misfit2_list = np.append(result_global_misfit2_list, global_misfit2)
        result_global_misfit3_list = np.append(result_global_misfit3_list, global_misfit3)
        
        ## store accepted variable parameters
        result_kappar_list = np.append(result_kappar_list, kappar)
        result_kappab_list = np.append(result_kappab_list, kappab)
        result_Te_list = np.append(result_Te_list, Te)
        result_T4_list = np.append(result_T4_list, T4_time_steps)
        #result_alpha_list = np.append(result_alpha_list, alpha)
        
        ## store the number of accepted model
        acm += 1
        
    ## print the calcul progress (for fix number of accepted model)
    if acm != acm_update:
        P = round(acm/MC_result*100,3) 
        print('percentage progress:', P,'%; misfit for accepted model:', round(global_misfit2, 3),', iterations:', i)
        acm_update = acm
               
    ## update indexation
    i += 1
    
end_time = datetime.now()

#%%=========================================================================%%#
#============================= BEST MISFIT RESULT ============================#

ind_min_misfit_1 = np.argmin(topo_r_misfit_list)        # find location of the best range topographic misfit
ind_min_misfit_2 = np.argmin(topo_b_misfit_list)        # find location of the best basin topographic misfit
ind_min_misfit_3 = np.argmin(subs_misfit_list)          # find location of the best subsidence misfit
ind_min_misfit_4 = np.argmin(global_misfit_list)        # find location of the best global (range topogrphy + basin topography + subsidence) misfit

print('----------------------------------------------------------------------')
print('Duration: {}'.format(end_time - start_time))
print('Best range topographic misfit result is: ',min(topo_r_misfit_list))
print('with Te: ',Te_list[ind_min_misfit_1],'; range diffusivity: ',kappar_list[ind_min_misfit_1],'; basin diffusivity: ',kappab_list[ind_min_misfit_1])
print('Best basin topographic misfit result is: ',min(topo_b_misfit_list))
print('with Te: ',Te_list[ind_min_misfit_2],'; range diffusivity: ',kappar_list[ind_min_misfit_2],'; basin diffusivity: ',kappab_list[ind_min_misfit_2])
print('Best subsidence misfit result is: ',min(subs_misfit_list))
print('with Te: ',Te_list[ind_min_misfit_3],'; range diffusivity: ',kappar_list[ind_min_misfit_3],'; basin diffusivity: ',kappab_list[ind_min_misfit_3])
print('Best global misfit result is: ',min(global_misfit_list))
print('with Te: ',Te_list[ind_min_misfit_4],'; range diffusivity: ',kappar_list[ind_min_misfit_4],'; basin diffusivity: ',kappab_list[ind_min_misfit_4])
print('----------------------------------------------------------------------')
print('mean lithosphere elastic thickness accepted: ', np.mean(result_Te_list),' with 1 standard deviation: ',np.std(result_Te_list)*1)
print('mean range diffusivity accepted: ', np.mean(result_kappar_list),' with 1 standard deviation: ',np.std(result_kappar_list)*1)
print('mean basin diffusivity accepted: ', np.mean(result_kappab_list),' with 1 standard deviation: ',np.std(result_kappab_list)*1)

sys.exit("...")

#%%=========================================================================%%#
font = {'family':'arial','size':8}
matplotlib.rc('font',**font)            

#%%=========================================================================%%#
#================================= FIGURE 8 ==================================#
fig = plt.figure(figsize=(7.125, 5.25))
ax1 = fig.add_subplot(221)
ax2 = fig.add_subplot(222)
ax3 = fig.add_subplot(223)
ax4 = fig.add_subplot(224)

kappa_ratio = kappab_list/kappar_list
misfit = np.log10((global_misfit2_list)**(1/1))
order = np.argsort(misfit)[::-1]
ss = 0.75

sc1 = ax1.scatter(Te_list[order], kappar_list[order], c=misfit[order], vmin=min(misfit), vmax=max(misfit)*ss, cmap='jet_r', s=5)
fig.colorbar(sc1, label='Log(\u03C7)', ax=ax1)

sc2 = ax2.scatter(Te_list[order], kappab_list[order], c=misfit[order], vmin=min(misfit), vmax=max(misfit)*ss, cmap='jet_r', s=5)
fig.colorbar(sc2, label='Log(\u03C7)', ax=ax2)

sc3 = ax3.scatter(kappar_list[order], kappab_list[order], c=misfit[order], vmin=min(misfit), vmax=max(misfit)*ss, cmap='jet_r', s=5)
fig.colorbar(sc3, label='Log(\u03C7)', ax=ax3)

sc4 = ax4.scatter(Te_list[order], kappa_ratio[order], c=misfit[order], vmin=min(misfit), vmax=max(misfit)*ss, cmap='jet_r', s=5)
fig.colorbar(sc4, label='Log(\u03C7)', ax=ax4)

ax1.axis([Te_min, Te_max, kappar_min, kappar_max])
ax2.axis([Te_min, Te_max, kappab_min, kappab_max])
ax3.axis([kappar_min, kappar_max, kappab_min, kappab_max])
ax4.axis([Te_min, Te_max, 0, 200])

ax1.set_xlabel('Lithosphere elastic thickness (km)'); ax1.set_ylabel('Range transport coefficient (m$^{2}$.yr$^{-1}$)')
ax2.set_xlabel('Lithosphere elastic thickness (km)'); ax2.set_ylabel('Basin transport coefficient (m$^{2}$.yr$^{-1}$)')
ax3.set_xlabel('Range transport coefficient (m$^{2}$.yr$^{-1}$)'); ax3.set_ylabel('Basin transport coefficient (m$^{2}$.yr$^{-1}$)')
ax4.set_xlabel('Lithosphere elastic thickness (km)'); ax4.set_ylabel('Transport coefficient ratio')

ticks_x = ticker.FuncFormatter(lambda x, pos: '{0:g}'.format(x*1e-3))
ax1.xaxis.set_major_formatter(ticks_x)
ax2.xaxis.set_major_formatter(ticks_x)
ax4.xaxis.set_major_formatter(ticks_x)

fig.tight_layout()
#fig.savefig('R:\\Post-orogenic sediment flux to continental margins\\Publication_Report\\Paper_Post-orogen_Pyrenees+Aquitaine\\Review_v1\\Figure8.jpg', dpi=1200)
#fig.savefig('R:\\Post-orogenic sediment flux to continental margins\\Publication_Report\\Paper_Post-orogen_Pyrenees+Aquitaine\\Review_v1\\Figure8.pdf', dpi=1200)

#%%=========================================================================%%#
#======================== FIGURE PLOTTING ==================================%%#
importlib.reload(functions_plotting_results)

#functions_plotting_results.range_vs_basin_diffusivity_sampling_pattern(kappar_list, kappab_list)

#functions_plotting_results.range_vs_basin_responsetime_sampling_pattern(tr_list, tb_list)

#functions_plotting_results.convergence_decrease_sampling_pattern(T4_time_steps_list, T3_time_steps)

#functions_plotting_results.misfit_vs_parameters_scatter_version(global_misfit, 'Misfit', kappar_list, kappab_list, Te_list, kappar_min, kappar_max, kappab_min, kappab_max, Te_min, Te_max, MC_iteration)

#functions_plotting_results.misfit_vs_parameters_hexbin_version(global_misfit, 'Misfit', int(MC_iteration/1000), 0.1, kappar_list, kappab_list, Te_list, kappar_min, kappar_max, kappab_min, kappab_max, Te_min, Te_max, MC_iteration)

#functions_plotting_results.misfit_vs_parameters_interpolated_1_version(kappab_list, kappar_list, Te_list, global_misfit, Te_min, Te_max, kappar_min, kappar_max, kappab_min, kappab_max)

#functions_plotting_results.misfit_vs_parameters_interpolated_2_version(kappab_list, kappar_list, Te_list, global_misfit, Te_min, Te_max, kappar_min, kappar_max, kappab_min, kappab_max)

#functions_plotting_results.misfit_vs_parameters_interpolated_3_version(kappab_list, kappar_list, Te_list, global_misfit, Te_min, Te_max, kappar_min, kappar_max, kappab_min, kappab_max)

#functions_plotting_results.range_elevation_vs_time(result_Hr_list, time, tt)

#functions_plotting_results.basin_elevation_vs_time(result_Hb_list, time, tt)

#functions_plotting_results.basin_deflexion_vs_time(result_wb_list, time, tt)

#functions_plotting_results.density_range_elevation_vs_time(result_Hr_list, result_time_list, time, tt)

#functions_plotting_results.density_basin_elevation_vs_time(result_Hb_list, result_time_list, time, tt)

#functions_plotting_results.density_basin_deflexion_vs_time(result_wb_list, result_time_list, time, tt)

#functions_plotting_results.range_flux_elevation_vs_time(result_Fe_list, time, tt)

#functions_plotting_results.basin_flux_elevation_vs_time(result_Fb_list, time, tt)

#functions_plotting_results.density_range_flux_elevation_vs_time(result_Fe_list, result_time_list, time, tt)

#functions_plotting_results.density_basin_flux_elevation_vs_time(result_Fb_list, result_time_list, time, tt)

#functions_plotting_results.specific_elevation_range_basin(MC_iteration, result_Hr_list, result_Hb_list, T3_time_steps)

#functions_plotting_results.specific_timing_range_basin(MC_iteration, T3_time_steps, T6_time_steps, result_Hb_list, result_Fe_list, result_Fb_list)

#functions_plotting_results.density_relief_vs_time(result_Hr_list, result_Hb_list, result_time_list, time, tt)

#functions_plotting_results.volume_change_range_vs_time(result_Fr_list, result_Fe_list, rhoc, result_time_list, time, tt)

#functions_plotting_results.volume_change_basin_vs_time(result_Fe_list, result_Fb_list, rhos, result_time_list, time, tt)

#functions_plotting_results.cross_section_thickness_range_vs_time(result_Hr_list, result_wr_list, time, tt)

#functions_plotting_results.cross_section_thickness_basin_vs_time(result_Hb_list, result_wb_list, time, tt)

#functions_plotting_results.contribution_elevation_change_basin(result_Fe_list, result_Fb_list, result_wb_list, rhos, Lb, time, T3_time_steps)

#functions_plotting_results.model_cross_section_vs_time(nt, result_Hr_list, result_Hb_list, result_etar_list, result_etab_list, result_Te_list, result_T4_list, rhos, rhoc, rhom, rhoi, E, nu, g, Lr, T1_time_steps, T2_time_steps, T3_time_steps, T6_time_steps, V1, V2, V3)

#functions_plotting_results.model_cross_section(4300, result_Hr_list, result_Hb_list, result_etar_list, result_etab_list, result_Te_list, rhos, rhoc, rhom, rhoi, E, nu, g, Lr, Vt)

#functions_plotting_results.figure_misfit_vs_parameters_test(kappab_list, kappar_list, Te_list, global_misfit_list, Te_min, Te_max, kappar_min, kappar_max, kappab_min, kappab_max)

#%%=========================================================================%%#
#=================== FIGURE PLOTTING FOR PUBLICATION =======================%%#
importlib.reload(functions_publication_figures)

#functions_figure_geology_paper.figure_parameter_setting(kappar_list, kappab_list, tr_list, tb_list, T1_time_steps, T2_time_steps, T3_time_steps, T4_time_steps, T6_time_steps,V1, V2, V3, time, dt, nt, path='Z:\\Post-orogenic sediment flux to continental margins\\Publication_Report\\Paper_Post-orogen_Pyrenees+Aquitaine\\geology\\')

#functions_publication_figures.figure_topo_sed_fluxes_results(time, tt, result_time_list, result_Hr_list, result_Hb_list, result_Fe_list, result_Fb_list, path='R:\\Post-orogenic sediment flux to continental margins\\Publication_Report\\Paper_Post-orogen_Pyrenees+Aquitaine\\truc\\')

#functions_publication_figures.figure_ratio_positive_vs_negative_post_basin_elevation_contrib(result_Fe_list, result_Fb_list, result_wb_list, result_Hb_list, result_global_misfit_list, rhos, Lb, time, T5_time_steps, path='R:\\Post-orogenic sediment flux to continental margins\\Publication_Report\\Paper_Post-orogen_Pyrenees+Aquitaine\\truc\\')

#functions_publication_figures.misfit_vs_parameters_scatter_version(global_misfit, 'Misfit', kappar_list, kappab_list, Te_list, kappar_min, kappar_max, kappab_min, kappab_max, Te_min, Te_max, MC_iteration)

#functions_publication_figures.misfit_vs_parameters_hexbin_version(global_misfit, 'Misfit', int(MC_iteration/2000), 0.1, kappar_list, kappab_list, Te_list, kappar_min, kappar_max, kappab_min, kappab_max, Te_min, Te_max, MC_iteration, path='Z:\\Post-orogenic sediment flux to continental margins\\Publication_Report\\Paper_Post-orogen_Pyrenees+Aquitaine\\geology\\')

#functions_publication_figures.figure_misfit_vs_parameters_parameters_interpolated_1_version(kappab_list, kappar_list, Te_list, global_misfit_list, Te_min, Te_max, kappar_min, kappar_max, kappab_min, kappab_max, path='R:\\Post-orogenic sediment flux to continental margins\\Publication_Report\\Paper_Post-orogen_Pyrenees+Aquitaine\\truc\\')

#functions_publication_figures.figure_alpha_effect(nt, result_alpha_list, result_Hb_list, time, result_Te_list, result_global_misfit_list, result_kappab_list, result_kappar_list)

## SUPPLEMENTARY FIGURE 2
#functions_publication_figures.misfit_vs_parameters_3D_1(kappab_list, kappar_list, Te_list, np.log10(global_misfit2_list), Te_min, Te_max, kappar_min, kappar_max, kappab_min, kappab_max, path='R:\\Post-orogenic sediment flux to continental margins\\Publication_Report\\Paper_Post-orogen_Pyrenees+Aquitaine\\Review_v1\\')

## SUPPLEMENTARY FIGURE 3
#functions_publication_figures.misfit_vs_parameters_3D_2(kappab_list, kappar_list, Te_list, np.log10(global_misfit2_list), Te_min, Te_max, kappar_min, kappar_max, kappab_min, kappab_max, path='R:\\Post-orogenic sediment flux to continental margins\\Publication_Report\\Paper_Post-orogen_Pyrenees+Aquitaine\\Review_v1\\')

## SUPPLEMENTARY FIGURE 1
#functions_publication_figures.figure_misfit_vs_time_convergence_decrease(T6_time_steps_list, np.log10(global_misfit2_list), Te_list, kappar_list, kappab_list, path='R:\\Post-orogenic sediment flux to continental margins\\Publication_Report\\Paper_Post-orogen_Pyrenees+Aquitaine\\Review_v1\\')

#functions_publication_figures.figure_cross_section_range_results(result_Hr_list, result_wr_list, time, result_time_list, tt, path='R:\\Post-orogenic sediment flux to continental margins\\Publication_Report\\Paper_Post-orogen_Pyrenees+Aquitaine\\truc\\')

#functions_publication_figures.figure_cross_section_basin_results(result_Hb_list, result_wb_list, time, result_time_list, tt, path='R:\\Post-orogenic sediment flux to continental margins\\Publication_Report\\Paper_Post-orogen_Pyrenees+Aquitaine\\truc\\')

#functions_publication_figures.figure_range_volume_change_results(result_Hr_list, result_wr_list, Lr, result_Fr_list, result_Fe_list, result_time_list, time, tt, rhoc, path='R:\\Post-orogenic sediment flux to continental margins\\Publication_Report\\Paper_Post-orogen_Pyrenees+Aquitaine\\truc\\')

#functions_publication_figures.figure_basin_volume_change_results(result_Hb_list, result_wb_list, result_Fe_list, result_Fb_list, time, tt, result_time_list, Lb, rhos, path='R:\\Post-orogenic sediment flux to continental margins\\Publication_Report\\Paper_Post-orogen_Pyrenees+Aquitaine\\truc\\')