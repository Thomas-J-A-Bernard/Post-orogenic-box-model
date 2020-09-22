"""
Created on Tue Mar 19 09:29:44 2019
@author: Thomas Bernard (Python 3.6)
Description:
functions to plot results for publication
"""

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
from matplotlib.colors import LogNorm
#from matplotlib import gridspec
from matplotlib import gridspec
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator)
from scipy.interpolate import griddata
import scipy.ndimage
from scipy.ndimage.filters import gaussian_filter
from mpl_toolkits.mplot3d import axes3d
import matplotlib

#%%=========================================================================%%# 
def figure_parameter_setting(kappar_list, kappab_list, tr_list, tb_list, T1_time_steps, T2_time_steps, T3_time_steps, T4_time_steps, T6_time_steps,V1, V2, V3, time, dt, nt, path):
  
  Vt = np.zeros(nt)
  Vt[T1_time_steps:T2_time_steps] = V1
  Vt[T2_time_steps:T3_time_steps] = V2
  Vt[T3_time_steps:T6_time_steps] = V3
  T4_time_steps = int(round(45.5e6/dt))
  V_T3 = V2
  for v in range(T3_time_steps,T4_time_steps):
    Vt[v] = V_T3 - (V2-V3)/(T4_time_steps-T3_time_steps)
    V_T3 = Vt[v]
  Vt1 = Vt
  
  Vt = np.zeros(nt)
  Vt[T1_time_steps:T2_time_steps] = V1
  Vt[T2_time_steps:T3_time_steps] = V2
  Vt[T3_time_steps:T6_time_steps] = V3
  T4_time_steps = int(round(50.5e6/dt))
  V_T3 = V2
  for v in range(T3_time_steps,T4_time_steps):
    Vt[v] = V_T3 - (V2-V3)/(T4_time_steps-T3_time_steps)
    V_T3 = Vt[v]
  Vt2 = Vt
  
  fig = plt.figure(figsize = (7.125,2.375))
  ax1 = fig.add_subplot(132)
  ax2 = fig.add_subplot(133)
  ax3 = fig.add_subplot(131)
  
  ax1.scatter(kappar_list, kappab_list, s=15, marker='o', color='steelblue', linewidths=0.5, edgecolor='k',zorder=2)
  ax2.scatter(tr_list, tb_list, s=15, marker='o', color='darkred', linewidths=0.5, edgecolor='k',zorder=2)
  ax3.plot(Vt1, color='darkgreen', linestyle='-',zorder=2)
  ax3.plot(Vt2, color='darkgreen', linestyle='-',zorder=2)
  ax3.fill_between(time*100, Vt2, Vt1, facecolor = 'darkgreen', alpha='0.33')
  
  plt.text(0.9, 0.95, 'A', fontweight='bold', horizontalalignment='left', verticalalignment='top', transform=ax3.transAxes)
  plt.text(0.9, 0.95, 'B', fontweight='bold', horizontalalignment='left', verticalalignment='top', transform=ax1.transAxes)
  plt.text(0.9, 0.95, 'C', fontweight='bold', horizontalalignment='left', verticalalignment='top', transform=ax2.transAxes)
  
  ax1.set_xlim(0,5000); ax1.set_ylim(10000,50000)
  ax2.set_xlim(0,max(tr_list)); ax2.set_ylim(200000,max(tb_list))
  ax3.set_xlim(0,6599);
  
  ax1.yaxis.grid(which='major', linestyle='--', color='k', alpha=0.5, zorder=0); ax1.xaxis.grid(which='major', linestyle='--', color='k', alpha=0.5, zorder=1)
  ax2.yaxis.grid(which='major', linestyle='--', color='k', alpha=0.5, zorder=0); ax2.xaxis.grid(which='major', linestyle='--', color='k', alpha=0.5, zorder=1)
  ax3.yaxis.grid(which='major', linestyle='--', color='k', alpha=0.5, zorder=0); ax3.xaxis.grid(which='major', linestyle='--', color='k', alpha=0.5, zorder=1)
  
  ax1.set_xlabel('Range diffusivity\n (10$^{3}$ m$^{2}$.yr$^{-1}$)', fontsize=8, fontname='arial'); ax1.set_ylabel('Basin diffusivity\n (10$^{3}$ m$^{2}$.yr$^{-1}$)', fontsize=8, fontname='arial')
  ax2.set_xlabel('Range response time\n (10$^{5}$ yrs)', fontsize=8, fontname='arial'); ax2.set_ylabel('Basin response time\n (10$^{5}$ yrs)', fontsize=8, fontname='arial')
  ax3.set_xlabel('Time (Myrs)', fontsize=8, fontname='arial'); ax3.set_ylabel('Convergence velocity\n (mm.yr$^{-1}$)', fontsize=8, fontname='arial')
  
  ticks_x = ticker.FuncFormatter(lambda x, pos: '{0:g}'.format(x/1000))
  ax1.xaxis.set_major_formatter(ticks_x)
  ticks_y = ticker.FuncFormatter(lambda y, pos: '{0:g}'.format(y/1000))
  ax1.yaxis.set_major_formatter(ticks_y)
  
  ticks_x = ticker.FuncFormatter(lambda x, pos: '{0:g}'.format(x/100000))
  ax2.xaxis.set_major_formatter(ticks_x)
  ticks_y = ticker.FuncFormatter(lambda y, pos: '{0:g}'.format(y/100000))
  ax2.yaxis.set_major_formatter(ticks_y)
  
  ticks_x = ticker.FuncFormatter(lambda x, pos: '{0:g}'.format(x/100))
  ax3.xaxis.set_major_formatter(ticks_x)
  ticks_y = ticker.FuncFormatter(lambda y, pos: '{0:g}'.format(y/0.001))
  ax3.yaxis.set_major_formatter(ticks_y)
  
  ax1.xaxis.set_tick_params(labelsize=8); ax1.yaxis.set_tick_params(labelsize=8)
  ax2.xaxis.set_tick_params(labelsize=8); ax2.yaxis.set_tick_params(labelsize=8)
  ax3.xaxis.set_tick_params(labelsize=8); ax3.yaxis.set_tick_params(labelsize=8)
  
  fig.subplots_adjust(left=0.095, bottom=0.25, right=0.985, top=0.95, wspace=0.4)
  #fig.savefig(path+'FIGURE2_montecarlo_model_setting'+'.jpeg', dpi=1200)
  #fig.savefig(path+'FIGURE2_montecarlo_model_setting'+'.pdf', dpi=1200)
  
#%%=========================================================================%%#
  
def figure_topo_sed_fluxes_results(time, tt, result_time_list, result_Hr_list, result_Hb_list, result_Fe_list, result_Fb_list, path):
  
  meanElevation_Hr = result_Hr_list.mean(axis=0)
  percentile_4p55_Elevation_Hr = np.percentile(result_Hr_list, 4.55, axis=0)
  percentile_95p45_Elevation_Hr = np.percentile(result_Hr_list, 95.45, axis=0)

  meanElevation_Hb = result_Hb_list.mean(axis=0)
  percentile_4p55_Elevation_Hb = np.percentile(result_Hb_list, 4.55, axis=0)
  percentile_95p45_Elevation_Hb = np.percentile(result_Hb_list, 95.45, axis=0)

  meanFlux_Fe = result_Fe_list.mean(axis=0)
  percentile_4p55_Flux_Fe = np.percentile(result_Fe_list, 4.55, axis=0)
  percentile_95p45_Flux_Fe = np.percentile(result_Fe_list, 95.45, axis=0)

  meanFlux_Fb = result_Fb_list.mean(axis=0)
  percentile_4p55_Flux_Fb = np.percentile(result_Fb_list, 4.55, axis=0)
  percentile_95p45_Flux_Fb = np.percentile(result_Fb_list, 95.45, axis=0)
  
  fig = plt.figure(figsize = (7.125, 3.5625))
  ax1 = fig.add_subplot(221)
  ax2 = fig.add_subplot(222)
  ax3 = fig.add_subplot(223)
  ax4 = fig.add_subplot(224)
  
  ax1.hist2d(result_time_list.flatten(), result_Hr_list.flatten(), bins = (100, 100), cmap = plt.cm.jet, norm = LogNorm(), zorder = 2)
  ax1.plot(time, meanElevation_Hr, linewidth = 1, color = 'k', linestyle = '-', zorder = 3)
  ax1.plot(time, percentile_4p55_Elevation_Hr, linewidth = 1, color = 'k', linestyle = '--', zorder = 3)
  ax1.plot(time, percentile_95p45_Elevation_Hr, linewidth = 1, color = 'k', linestyle = '--', zorder = 3)
  
  ax2.hist2d(result_time_list.flatten(), result_Hb_list.flatten(), bins = (100, 100), cmap = plt.cm.jet, norm = LogNorm(), zorder = 2)
  ax2.plot(time, meanElevation_Hb, linewidth = 1, color = 'k', linestyle = '-', zorder = 3)
  ax2.plot(time, percentile_4p55_Elevation_Hb, linewidth = 1, color = 'k', linestyle = '--', zorder = 3)
  ax2.plot(time, percentile_95p45_Elevation_Hb, linewidth = 1, color = 'k', linestyle = '--', zorder = 3)
  
  ax3.hist2d(result_time_list.flatten(), result_Fe_list.flatten(), bins = (100, 100), cmap = plt.cm.jet, norm = LogNorm(), zorder = 2)
  ax3.plot(time, meanFlux_Fe, linewidth = 1, color = 'k', linestyle = '-', zorder = 3)
  ax3.plot(time, percentile_4p55_Flux_Fe, linewidth = 1, color = 'k', linestyle = '--', zorder = 3)
  ax3.plot(time, percentile_95p45_Flux_Fe, linewidth = 1, color = 'k', linestyle = '--', zorder = 3)
  
  ax4.hist2d(result_time_list.flatten(), result_Fb_list.flatten(), bins = (100, 100), cmap = plt.cm.jet, norm = LogNorm(), zorder = 2)
  ax4.plot(time, meanFlux_Fb, linewidth = 1, color = 'k', linestyle = '-', zorder = 3)
  ax4.plot(time, percentile_4p55_Flux_Fb, linewidth = 1, color = 'k', linestyle = '--', zorder = 3)
  ax4.plot(time, percentile_95p45_Flux_Fb, linewidth = 1, color = 'k', linestyle = '--', zorder = 3)
  
  plt.text(0.03, 0.95, 'A) Range elevation', fontweight='bold', fontsize=8, horizontalalignment='left', verticalalignment='top', transform=ax1.transAxes)
  plt.text(0.03, 0.95, 'B) Basin elevation', fontweight='bold', fontsize=8, horizontalalignment='left', verticalalignment='top', transform=ax2.transAxes)
  plt.text(0.03, 0.95, 'C) Sediment flux to basin', fontweight='bold', fontsize=8, horizontalalignment='left', verticalalignment='top', transform=ax3.transAxes)
  plt.text(0.03, 0.95, 'D) Sediment flux to sea', fontweight='bold', fontsize=8, horizontalalignment='left', verticalalignment='top', transform=ax4.transAxes)
  
  ax1.axis([0,tt*1e-6,result_Hr_list.min()*1.1,result_Hr_list.max()*1.1]); ax2.axis([0,tt*1e-6,result_Hb_list.min()*1.1,result_Hb_list.max()*1.1]); ax3.axis([0,tt*1e-6,0,result_Fe_list.max()*1.1]); ax4.axis([0,tt*1e-6,result_Fb_list.min()*1.1,result_Fb_list.max()*1.1])

  ax1.yaxis.grid(which = 'major', linestyle = '--', color = 'k', alpha = 0.5, zorder = 0); ax1.xaxis.grid(which = 'major', linestyle = '--', color = 'k', alpha = 0.5, zorder = 1)
  ax2.yaxis.grid(which = 'major', linestyle = '--', color = 'k', alpha = 0.5, zorder = 0); ax2.xaxis.grid(which = 'major', linestyle = '--', color = 'k', alpha = 0.5, zorder = 1)
  ax3.yaxis.grid(which = 'major', linestyle = '--', color = 'k', alpha = 0.5, zorder = 0); ax3.xaxis.grid(which = 'major', linestyle = '--', color = 'k', alpha = 0.5, zorder = 1)
  ax4.yaxis.grid(which = 'major', linestyle = '--', color = 'k', alpha = 0.5, zorder = 0); ax4.xaxis.grid(which = 'major', linestyle = '--', color = 'k', alpha = 0.5, zorder = 1)
  
  ax1.set_ylabel('Elevation (m)', fontsize=8, fontname='arial'); ax3.set_xlabel('Time (Ma)', fontsize=8, fontname='arial'); ax3.set_ylabel('Sediment Flux (kg.m$^{-1}$.yr$^{-1}$)', fontsize=8, fontname='arial'); ax4.set_xlabel('Time (Ma)', fontsize=8, fontname='arial')
 
  ax1.xaxis.set_tick_params(labelsize=8); ax1.yaxis.set_tick_params(labelsize=8)
  ax2.xaxis.set_tick_params(labelsize=8); ax2.yaxis.set_tick_params(labelsize=8)
  ax3.xaxis.set_tick_params(labelsize=8); ax3.yaxis.set_tick_params(labelsize=8)
  ax4.xaxis.set_tick_params(labelsize=8); ax4.yaxis.set_tick_params(labelsize=8)
  ax3.yaxis.offsetText.set_fontsize(8); ax4.yaxis.offsetText.set_fontsize(8)
  
  fig.subplots_adjust(left=0.095, bottom=0.125, right=0.985, top=0.985, wspace=0.15, hspace=0.25)
  #fig.savefig(path+'FIGURE7_topo+sed_fluxes_result'+'.jpeg', dpi=1200)
  #fig.savefig(path+'FIGURE7_topo+sed_fluxes_result-01'+'.pdf', dpi=1200)

#%%=========================================================================%%#
  
def figure_ratio_positive_vs_negative_post_basin_elevation_contrib(result_Fe_list, result_Fb_list, result_wb_list, result_Hb_list, result_global_misfit_list, rhos, Lb, time, T5_time_steps, path):
  
  Flux_Fe = result_Fe_list.mean(axis=0)
  Flux_Fb = result_Fb_list.mean(axis=0)
  Elevation_Hb = result_Hb_list.mean(axis=0)
  Deflexion_wb = result_wb_list.mean(axis=0)
  
  t_post = T5_time_steps; t_end = int(round(time[-1]*100, 0)); post = t_end - t_post
  Hb_4300 = Elevation_Hb[4300]
  
  post_Deflexion_wb_change = np.zeros(post)
  post_positive_elevation = Flux_Fe[t_post:t_end]/rhos/Lb
  post_negative_elevation = Flux_Fb[t_post:t_end]/rhos/Lb
  post_elevation_Hb = Elevation_Hb[t_post:t_end] - Hb_4300
  
  #max_index_post_Elevation_Hb = np.argmax(post_elevation_Hb)
  
  for i in range(0,post):
    post_Deflexion_wb_change[i] = Deflexion_wb[i+(t_post-1)]-Deflexion_wb[i+t_post]

  index_post_Deflexion_wb_change = np.where(post_Deflexion_wb_change >= 0); index_post_Deflexion_wb_change = index_post_Deflexion_wb_change[0][0]
  
  for i in range(0,post):
    if post_Deflexion_wb_change[i] < 0:
      post_negative_elevation[i] = post_negative_elevation[i] - post_Deflexion_wb_change[i]
    else:
      post_positive_elevation[i] = post_positive_elevation[i] + post_Deflexion_wb_change[i]
  
  post_basin_elevation = post_positive_elevation - post_negative_elevation
  post_basin_elevation = np.cumsum(post_basin_elevation)
  
  fig = plt.figure(figsize = (3.5, 4))  
  #gs = gridspec.GridSpec(2,2)
  #ax1 = fig.add_subplot(gs[0,:])
  #ax2 = fig.add_subplot(gs[1,0])
  #ax3 = fig.add_subplot(gs[1,1])
  ax1 = fig.add_subplot(311)
  ax2 = fig.add_subplot(312)
  ax3 = fig.add_subplot(313)
  
  ax1.plot(time[t_post:t_end], post_Deflexion_wb_change, color='darkgreen', linestyle='-')
  ax1.plot(time[t_post:t_end], -Flux_Fb[t_post:t_end]/rhos/Lb, color='steelblue', linestyle='-')
  ax1.plot(time[t_post:t_end], Flux_Fe[t_post:t_end]/rhos/Lb, color='darkred', linestyle='-')
  ax1.fill_between(time[t_post:t_end], post_Deflexion_wb_change, 0, facecolor='darkseagreen', alpha=0.5)
  ax1.fill_between(time[t_post:t_post+index_post_Deflexion_wb_change+1], Flux_Fe[t_post:t_post+index_post_Deflexion_wb_change+1]/rhos/Lb, 0, facecolor='indianred', alpha=0.5)
  ax1.fill_between(time[t_post+index_post_Deflexion_wb_change:t_end], Flux_Fe[t_post+index_post_Deflexion_wb_change:t_end]/rhos/Lb, post_Deflexion_wb_change[index_post_Deflexion_wb_change:post], facecolor='indianred', alpha=0.5)
  ax1.fill_between(time[t_post:t_post+index_post_Deflexion_wb_change+1], -Flux_Fb[t_post:t_post+index_post_Deflexion_wb_change+1]/rhos/Lb, post_Deflexion_wb_change[0:index_post_Deflexion_wb_change+1], facecolor='powderblue', alpha=0.5)
  ax1.fill_between(time[t_post+index_post_Deflexion_wb_change:t_end], -Flux_Fb[t_post+index_post_Deflexion_wb_change:t_end]/rhos/Lb, 0, facecolor='powderblue', alpha=0.5)
  
  ax2.plot(time[t_post:t_end], post_Deflexion_wb_change - Flux_Fb[t_post:t_end]/rhos/Lb + Flux_Fe[t_post:t_end]/rhos/Lb, color='k', linestyle='-')
  
  ax3.plot(time[t_post:t_end], post_elevation_Hb, color='k', linestyle='-', label = 'ratio')
  #ax3.plot([], [], color='powderblue', alpha=0.66, label='elevation gain')
  #ax3.plot([], [], color='indianred', alpha=0.66, label='elevation loss')
  #ax3.fill_between(time[t_post:t_post+max_index_post_Elevation_Hb], post_elevation_Hb[0:max_index_post_Elevation_Hb], 0, facecolor='powderblue', alpha=0.5)
  #ax3.fill_between(time[t_post+max_index_post_Elevation_Hb:t_end], post_elevation_Hb[max_index_post_Elevation_Hb:post], 0, facecolor='indianred', alpha=0.5)
  
  #plt.text(0.96, 0.95, 'A', fontweight='bold', horizontalalignment='left', verticalalignment='top', transform=ax1.transAxes)
  #plt.text(0.96, 0.95, 'B', fontweight='bold', horizontalalignment='left', verticalalignment='top', transform=ax2.transAxes)
  
  ax1.axis([time[t_post], time[t_end], np.max(Flux_Fb[t_post:t_end]/rhos/Lb)*(-1.1), np.max(Flux_Fe[t_post:t_end]/rhos/Lb)*1.1])
  ax1.yaxis.grid(which='major', linestyle='--', color='k', alpha=0.5, zorder=0); ax1.xaxis.grid(which='major', linestyle='--', color='k', alpha=0.5, zorder=0)
  ax1.set_ylabel('Average rate\n(m.yr$^{-1}$.10$^{-5}$)')
  
  ax2.set_xlim(time[t_post], time[t_end])
  ax2.yaxis.grid(which='major', linestyle='--', color='k', alpha=0.5, zorder=0); ax2.xaxis.grid(which='major', linestyle='--', color='k', alpha=0.5, zorder=0)
  ax2.set_ylabel('Average rate\n(m.yr$^{-1}$.10$^{-5}$)')
  
  ax3.axis([time[t_post], time[t_end], min(post_elevation_Hb)*1.1, max(post_elevation_Hb)*1.2])
  ax3.yaxis.grid(which='major', linestyle='--', color='k', alpha=0.5, zorder=0); ax3.xaxis.grid(which='major', linestyle='--', color='k', alpha=0.5, zorder=0)
  ax3.set_xlabel('Time (Ma)'); ax3.set_ylabel('Elevation\n(m)')  
  ax3.yaxis.set_major_locator(MultipleLocator(100))
  
  xlabel=['0','21','16','11','6','1']; ax1.set_xticklabels(xlabel); ax2.set_xticklabels(xlabel); ax3.set_xticklabels(xlabel)
  #xlabel=['0','16','6']; 
  
  fig.tight_layout()    
  #fig.savefig(path+'FIGURE8_post_rebound_contribution'+'.jpeg', dpi=1200)
  #fig.savefig(path+'FIGURE8_post_rebound_contribution'+'.pdf', dpi=1200)
  
#%%=========================================================================%%#
def misfit_vs_parameters_scatter_version(misfit, label_misfit, kappar_list, kappab_list, Te_list, kappar_min, kappar_max, kappab_min, kappab_max, Te_min, Te_max, MC_iteration):
  
  kappa_ratio_list = kappar_list/kappab_list
  kappa_ratio_list2 = (kappa_ratio_list).tolist()
  Te_list2 = Te_list.tolist()
  misfit2 = misfit.tolist()
  
  index = np.empty((0))
  for i in range (0, MC_iteration):
    if kappa_ratio_list2[i]>1:
      index = np.append(index, int(i))
  
  index = index.tolist()
  for i in range(0,np.size(index)):
    index[i] = int(index[i])
  
  for i in sorted(index, reverse=True):
    del kappa_ratio_list2[i]
    del Te_list2[i]
    del misfit2[i]
  
  kappa_ratio_list2 = np.array(kappa_ratio_list2)
  Te_list2 = np.array(Te_list2)
  misfit2 = np.array(misfit2)
  
  fig = plt.figure(figsize=(9,8))
  ax1 = fig.add_subplot(221)
  ax2 = fig.add_subplot(222)
  ax3 = fig.add_subplot(223)
  ax4 = fig.add_subplot(224)
  
  sc = ax1.scatter(kappar_list, kappab_list, c=misfit, cmap='rainbow_r', vmin=min(misfit), vmax=max(misfit), marker='o', linewidths=1, edgecolors='black')
  fig.colorbar(sc, label=label_misfit, aspect=30, ax=ax1)
  
  sc = ax2.scatter(Te_list, kappar_list, c=misfit, cmap='rainbow_r', vmin=min(misfit), vmax=max(misfit), marker='o', linewidths=1, edgecolors='black')
  fig.colorbar(sc, label=label_misfit, aspect=30, ax=ax2)
  
  sc = ax3.scatter(Te_list, kappab_list, c=misfit, cmap='rainbow_r', vmin=min(misfit), vmax=max(misfit), marker='o', linewidths=1, edgecolors='black')
  fig.colorbar(sc, label=label_misfit, aspect=30, ax=ax3)
  
  sc = ax4.scatter(Te_list2, kappa_ratio_list2, c=misfit2, cmap='rainbow_r', vmin=min(misfit), vmax=max(misfit), marker='o', linewidths=1, edgecolors='black')
  fig.colorbar(sc, label=label_misfit, aspect=30, ax=ax4)
  
  ax1.set_xlabel('Range diffusivity (m$^{2}$.yr$^{-1}$)'); ax1.set_ylabel('Basin diffusivity (m$^{2}$.yr$^{-1}$)')
  ax1.axis([kappar_min, kappar_max, kappab_min, kappab_max])

  ax2.set_xlabel('Lithosphere elastic thickness (m)'); ax2.set_ylabel('Range diffusivity (m$^{2}$.yr$^{-1}$)')
  ax2.axis([Te_min, Te_max, kappar_min, kappar_max])
  
  ax3.set_xlabel('Lithosphere elastic thickness (m)'); ax3.set_ylabel('Basin diffusivity (m$^{2}$.yr$^{-1}$)')
  ax3.axis([Te_min, Te_max, kappab_min, kappab_max])
  
  ax4.set_xlabel('Lithosphere elastic thickness (m)'); ax4.set_ylabel('Diffusivity ratio')
  ax4.axis([Te_min, Te_max, min(kappa_ratio_list2), max(kappa_ratio_list2)])
  
  fig.subplots_adjust(left=0.1, bottom=0.1, right=0.985, top=0.95, wspace=0.3, hspace=0.3)
  #fig.savefig('Z:\Post-orogenic sediment flux to continental margins\Modelling_Post_orogenic_system\Tucker_Code\Result_Northern_Pyrenees_InverseModelling_MC_approach\GBS_Hr=1690_wb=4434_t=23_Hr=1470_Hb=200_t=0_Hb=500_max_Te+tr+tb_var_i=100000.png', dpi=1200)

#%%=========================================================================%%#
def misfit_vs_parameters_hexbin_version(misfit, label_misfit, gridsize, linewidths, kappar_list, kappab_list, Te_list, kappar_min, kappar_max, kappab_min, kappab_max, Te_min, Te_max, MC_iteration, path):
  
  kappa_ratio_list = kappar_list/kappab_list
  kappa_ratio_list2 = (kappa_ratio_list).tolist()
  Te_list2 = Te_list.tolist()
  misfit2 = misfit.tolist()
  
  index = np.empty((0))
  for i in range (0, MC_iteration):
    if kappa_ratio_list2[i]>1:
      index = np.append(index, int(i))
  
  index = index.tolist()
  for i in range(0,np.size(index)):
    index[i] = int(index[i])
  
  for i in sorted(index, reverse=True):
    del kappa_ratio_list2[i]
    del Te_list2[i]
    del misfit2[i]
  
  kappa_ratio_list2 = np.array(kappa_ratio_list2)
  Te_list2 = np.array(Te_list2)
  misfit2 = np.array(misfit2)
  
  fig = plt.figure(figsize=(9,8))
  ax1 = fig.add_subplot(221)
  ax2 = fig.add_subplot(222)
  ax3 = fig.add_subplot(223)
  ax4 = fig.add_subplot(224)
  
  sc = ax1.hexbin(kappar_list, kappab_list, C=misfit, gridsize=gridsize, cmap='rainbow_r', vmin=min(misfit), vmax=max(misfit), linewidths=linewidths, edgecolors='black')
  fig.colorbar(sc, label=label_misfit, aspect=30, ax=ax1)
  
  sc = ax2.hexbin(Te_list, kappar_list, C=misfit, gridsize=gridsize, cmap='rainbow_r', vmin=min(misfit), vmax=max(misfit), linewidths=linewidths, edgecolors='black')
  fig.colorbar(sc, label=label_misfit, aspect=30, ax=ax2)
  
  sc = ax3.hexbin(Te_list, kappab_list, C=misfit, gridsize=gridsize, cmap='rainbow_r', vmin=min(misfit), vmax=max(misfit), linewidths=linewidths, edgecolors='black')
  fig.colorbar(sc, label=label_misfit, aspect=30, ax=ax3)
  
  sc = ax4.hexbin(Te_list2, kappa_ratio_list2, C=misfit2, gridsize=gridsize, cmap='rainbow_r', vmin=min(misfit), vmax=max(misfit), linewidths=linewidths, edgecolors='black')
  fig.colorbar(sc, label=label_misfit, aspect=30, ax=ax4)
  
  ax1.set_xlabel('Range diffusivity (m$^{2}$.yr$^{-1}$)'); ax1.set_ylabel('Basin diffusivity (m$^{2}$.yr$^{-1}$)')
  ax1.axis([kappar_min, kappar_max, kappab_min, kappab_max])
  
  ax2.set_xlabel('Lithosphere elastic thickness (m)'); ax2.set_ylabel('Range diffusivity (m$^{2}$.yr$^{-1}$)')
  ax2.axis([Te_min, Te_max, kappar_min, kappar_max])
  
  ax3.set_xlabel('Lithosphere elastic thickness (m)'); ax3.set_ylabel('Basin diffusivity (m$^{2}$.yr$^{-1}$)')
  ax3.axis([Te_min, Te_max, kappab_min, kappab_max])
  
  ax4.set_xlabel('Lithosphere elastic thickness (m)'); ax4.set_ylabel('Diffusivity ratio')
  ax4.axis([Te_min, Te_max, min(kappa_ratio_list2), max(kappa_ratio_list2)])
  
  plt.subplots_adjust(left=0.1, bottom=0.1, right=0.985, top=0.95, wspace=0.3, hspace=0.3)
  
  #fig.savefig(path+'SUP_FIG2_misfit_vs_variable_parameters'+'.jpeg', dpi=1200)
  #fig.savefig(path+'SUP_FIG2_misfit_vs_variable_parameters'+'.pdf', dpi=1200)
  
#%%=========================================================================%%#
def figure_misfit_vs_parameters_parameters_interpolated_1_version(kappab_list, kappar_list, Te_list, global_misfit, Te_min, Te_max, kappar_min, kappar_max, kappab_min, kappab_max, path):

  kappa_ratio_list = kappab_list/kappar_list
  
  x1 = Te_list; y1 = kappar_list; z1 = global_misfit
  grid_x1, grid_y1 = np.mgrid[Te_min:Te_max:100j, kappar_min:kappar_max:100j]
  grid_z1 = griddata((x1, y1), z1, (grid_x1, grid_y1), method='nearest')
  grid_z1 = scipy.ndimage.zoom(grid_z1, 3)
  
  x2 = Te_list; y2 = kappab_list; z2 = global_misfit
  grid_x2, grid_y2 = np.mgrid[Te_min:Te_max:100j, kappab_min:kappab_max:100j]
  grid_z2 = griddata((x2, y2), z2, (grid_x2, grid_y2), method='nearest')
  grid_z2 = scipy.ndimage.zoom(grid_z2, 3)
  
  x3 = kappar_list; y3 = kappab_list; z3 = global_misfit
  grid_x3, grid_y3 = np.mgrid[kappar_min:kappar_max:100j, kappab_min:kappab_max:100j]
  grid_z3 = griddata((x3, y3), z3, (grid_x3, grid_y3), method='nearest')
  grid_z3 = scipy.ndimage.zoom(grid_z3, 3)
  
  x4 = Te_list; y4 = kappa_ratio_list; z4 = global_misfit
  grid_x4, grid_y4 = np.mgrid[Te_min:Te_max:100j, 0:500:100j]
  grid_z4 = griddata((x4, y4), z4, (grid_x4, grid_y4), method='nearest')
  grid_z4 = scipy.ndimage.zoom(grid_z4, 3)
  
  fig = plt.figure(figsize=(7.125,5.5))
  ax1 = fig.add_subplot(221); ax2 = fig.add_subplot(222); ax3 = fig.add_subplot(223); ax4 = fig.add_subplot(224)
  
  ax1.imshow(grid_z1.T, extent=(0,100,0,100), origin='lower', cmap='gray_r')
  co1 = ax1.contour(grid_z1.T, levels=[0.15, 0.25, 0.35, 0.45, 0.55, 0.65], cmap='rainbow_r', extent=(0,100,0,100), linewidths=1.25)
  norm1= matplotlib.colors.Normalize(vmin=co1.cvalues.min(), vmax=co1.cvalues.max())
  sm1 = plt.cm.ScalarMappable(norm=norm1, cmap = co1.cmap)
  sm1.set_array([])
  fig.colorbar(sm1, ticks=co1.levels, label='Misfit', fraction=0.046, ax=ax1)
  
  ax2.imshow(grid_z2.T, extent=(0,100,0,100), origin='lower', cmap='gray_r')
  co2 = ax2.contour(grid_z2.T, levels=[0.15, 0.25, 0.35, 0.45, 0.55, 0.65], cmap='rainbow_r', extent=(0,100,0,100), linewidths=1.25)
  norm2= matplotlib.colors.Normalize(vmin=co2.cvalues.min(), vmax=co2.cvalues.max())
  sm2 = plt.cm.ScalarMappable(norm=norm2, cmap = co2.cmap)
  sm2.set_array([])
  fig.colorbar(sm2, ticks=co1.levels, label='Misfit', fraction=0.046, ax=ax2)
  
  ax3.imshow(grid_z3.T, extent=(0,100,0,100), origin='lower', cmap='gray_r')
  co3 = ax3.contour(grid_z3.T, levels=[0.15, 0.25, 0.35, 0.45, 0.55, 0.65], cmap='rainbow_r', extent=(0,100,0,100), linewidths=1.25)
  norm3= matplotlib.colors.Normalize(vmin=co3.cvalues.min(), vmax=co3.cvalues.max())
  sm3 = plt.cm.ScalarMappable(norm=norm3, cmap = co3.cmap)
  sm3.set_array([])
  fig.colorbar(sm3, ticks=co3.levels, label='Misfit', fraction=0.046, ax=ax3)
  
  ax4.imshow(grid_z4.T, extent=(0,100,0,100), origin='lower', cmap='gray_r')
  co4 = ax4.contour(grid_z4.T, levels=[0.15, 0.25, 0.35, 0.45, 0.55, 0.65], cmap='rainbow_r', extent=(0,100,0,100), linewidths=1.25)
  norm4= matplotlib.colors.Normalize(vmin=co4.cvalues.min(), vmax=co4.cvalues.max())
  sm4 = plt.cm.ScalarMappable(norm=norm4, cmap = co4.cmap)
  sm4.set_array([])
  fig.colorbar(sm4, ticks=co3.levels, label='Misfit', fraction=0.046, ax=ax4)
  
  ax1.set_xlabel('Lithosphere elastic thickness (km)'); ax1.set_ylabel('Range transport coefficient (m$^{2}$.yr$^{-1}$)')
  ax2.set_xlabel('Lithosphere elastic thickness (km)'); ax2.set_ylabel('Basin transport coefficient (m$^{2}$.yr$^{-1}$)')
  ax3.set_xlabel('Range transport coefficient (m$^{2}$.yr$^{-1}$)'); ax3.set_ylabel('Basin transport coefficient (m$^{2}$.yr$^{-1}$)')
  ax4.set_xlabel('Lithosphere elastic thickness (km)'); ax4.set_ylabel('transport coefficient ratio')
  
  Te_label = ['15','20','25','30','35','40']
  ax1.set_xticklabels(Te_label); ax2.set_xticklabels(Te_label); ax4.set_xticklabels(Te_label)
  kappar_label = ['100','1080','2060','3040','4020','5000']
  ax3.set_xticklabels(kappar_label); ax1.set_yticklabels(kappar_label)
  kappab_label = ['1000','10800','20600','30400','40200','50000']
  ax2.set_yticklabels(kappab_label); ax3.set_yticklabels(kappab_label)
  ratio_label = ['0','100','200','300','400','500']
  ax4.set_yticklabels(ratio_label)
  
  fig.subplots_adjust(left=0.1, bottom=0.085, right=0.925, top=0.975, wspace=0.5, hspace=0.25)
  #fig.savefig(path+'FIGURE6_misfit+contour_vs_parameters'+'.jpeg', dpi=1200)
  #fig.savefig(path+'FIGURE6_misfit+contour_vs_parameters-01'+'.pdf', dpi=1200)

#%%=========================================================================%%#
def figure_alpha_effect(nt, result_alpha_list, result_Hb_list, time, result_Te_list, result_global_misfit_list, result_kappab_list, result_kappar_list):

  B1 = np.empty((0, nt)); B2 = np.empty((0, nt)); B3 = np.empty((0, nt)); B4 = np.empty((0, nt)); B5 = np.empty((0, nt)); B6 = np.empty((0, nt)); B7 = np.empty((0, nt)); B8 = np.empty((0, nt)); B9 = np.empty((0, nt)); B10 = np.empty((0, nt))
  
  for i in range(0, np.shape(result_alpha_list)[0]):
    if result_alpha_list[i] >= 0.2 and result_alpha_list[i] < 0.25:
      B1 = np.append (B1, result_Hb_list[i].reshape((1,nt)), axis=0)
  for i in range(0, np.shape(result_alpha_list)[0]):
    if result_alpha_list[i] >= 0.25 and result_alpha_list[i] < 0.30:
      B2 = np.append (B2, result_Hb_list[i].reshape((1,nt)), axis=0)
  for i in range(0, np.shape(result_alpha_list)[0]):
    if result_alpha_list[i] >= 0.30 and result_alpha_list[i] < 0.35:
      B3 = np.append (B3, result_Hb_list[i].reshape((1,nt)), axis=0)
  for i in range(0, np.shape(result_alpha_list)[0]):
    if result_alpha_list[i] >= 0.35 and result_alpha_list[i] < 0.40:
      B4 = np.append (B4, result_Hb_list[i].reshape((1,nt)), axis=0)
  for i in range(0, np.shape(result_alpha_list)[0]):
    if result_alpha_list[i] >= 0.28 and result_alpha_list[i] < 0.30:
      B5 = np.append (B5, result_Hb_list[i].reshape((1,nt)), axis=0)
  for i in range(0, np.shape(result_alpha_list)[0]):
    if result_alpha_list[i] >= 0.30 and result_alpha_list[i] < 0.32:
      B6 = np.append (B6, result_Hb_list[i].reshape((1,nt)), axis=0)
  for i in range(0, np.shape(result_alpha_list)[0]):
    if result_alpha_list[i] >= 0.32 and result_alpha_list[i] < 0.34:
      B7 = np.append (B7, result_Hb_list[i].reshape((1,nt)), axis=0)
  for i in range(0, np.shape(result_alpha_list)[0]):
    if result_alpha_list[i] >= 0.34 and result_alpha_list[i] < 0.36:
      B8 = np.append (B8, result_Hb_list[i].reshape((1,nt)), axis=0)
  for i in range(0, np.shape(result_alpha_list)[0]):
    if result_alpha_list[i] >= 0.36 and result_alpha_list[i] < 0.38:
      B9 = np.append (B9, result_Hb_list[i].reshape((1,nt)), axis=0)
  for i in range(0, np.shape(result_alpha_list)[0]):
    if result_alpha_list[i] >= 0.38 and result_alpha_list[i] < 0.40:
      B10 = np.append (B10, result_Hb_list[i].reshape((1,nt)), axis=0)
  
  B1 = B1.mean(axis=0); B2 = B2.mean(axis=0); B3 = B3.mean(axis=0); B4 = B4.mean(axis=0); B5 = B5.mean(axis=0); B6 = B6.mean(axis=0); B7 = B7.mean(axis=0); B8 = B8.mean(axis=0); B9 = B9.mean(axis=0); B10 = B10.mean(axis=0)
  
  fig = plt.figure(figsize=(3.5,3))
  gs = gridspec.GridSpec(2,2)
  ax1 = fig.add_subplot(gs[0, :])
  ax2 = fig.add_subplot(gs[1, 0])
  ax3 = fig.add_subplot(gs[1, 1])
  
  ax1.plot(time, B1, label='0.20<=alpha<0.25'); ax1.plot(time, B2, label='0.25<=alpha<0.30'); ax1.plot(time, B3, label='0.30<=alpha<0.35'); ax1.plot(time, B4, label='0.35<=alpha<0.40');
  ax1.legend(facecolor='white', edgecolor='k', framealpha=1, borderpad=0.5)
  ax1.set_xlabel('Time (Ma)'); ax1.set_ylabel('Basin elevation (m)')
  ax1.yaxis.grid(which='major',linestyle='--',color='k',alpha=0.5,zorder=0); ax1.xaxis.grid(which='major',linestyle='--',color='k',alpha=0.5,zorder=0)
  ax1.axis([43, 66,50, 400])
  xlabel=['0','21','16','11','6','1']; ax1.set_xticklabels(xlabel)
  
  
  sc1 = ax2.scatter(result_alpha_list, result_Te_list, c=result_global_misfit_list, s=20, vmin=np.min(result_global_misfit_list), vmax=np.max(result_global_misfit_list), zorder=2)
  fig.colorbar(sc1, label='Misfit', ax=ax2)
  ax2.set_xlabel('alpha'); ax2.set_ylabel('lithosphere elastic thickness')
  ax2.yaxis.grid(which='major',linestyle='--',color='k',alpha=0.5,zorder=1); ax2.xaxis.grid(which='major',linestyle='--',color='k',alpha=0.5,zorder=1)
  ax2.axis([0.2, 0.4, 15000, 40000])
  
  sc2 = ax3.scatter(result_alpha_list, result_kappab_list/result_kappar_list, c=result_global_misfit_list, s=20, vmin=np.min(result_global_misfit_list), vmax=np.max(result_global_misfit_list), zorder=2)
  fig.colorbar(sc2, label='Misfit', ax=ax3)
  ax3.set_xlabel('alpha'); ax3.set_ylabel('transport coefficient ratio')
  ax3.yaxis.grid(which='major',linestyle='--',color='k',alpha=0.5,zorder=1); ax3.xaxis.grid(which='major',linestyle='--',color='k',alpha=0.5,zorder=1)
  ax3.axis([0.2, 0.4, 0, 100])
  
  fig.subplots_adjust(left=0.1, bottom=0.1, right=0.950, top=0.950, wspace=0.30, hspace=0.25)
  
  #fig.savefig('Z:\\Post-orogenic sediment flux to continental margins\\Publication_Report\\Paper_Post-orogen_Pyrenees+Aquitaine\\long_version\\FIGURE_alpha_effect.pdf', dpi=1200)
  #fig.savefig('Z:\\Post-orogenic sediment flux to continental margins\\Publication_Report\\Paper_Post-orogen_Pyrenees+Aquitaine\\long_version\\FIGURE_alpha_effect.jpg', dpi=1200)

#%%=========================================================================%%#
def misfit_vs_parameters_3D_1(kappab_list, kappar_list, Te_list, global_misfit_list, Te_min, Te_max, kappar_min, kappar_max, kappab_min, kappab_max, path):

    x1 = Te_list; y1 = kappar_list; z1 = kappab_list; s1 = global_misfit_list
    grid_x1, grid_y1, grid_z1 = np.mgrid[Te_min:Te_max:15j, kappar_min:kappar_max:15j, kappab_min:kappab_max:15j]
    
    M = griddata((x1, y1, z1), s1, (grid_x1, grid_y1, grid_z1), method='nearest')
      
    M2 = griddata((x1, y1, z1), s1, (grid_x1, grid_y1, grid_z1), method='nearest')
    M2[0:5,:,:] = np.nan; M2[6:15,:,:] = np.nan
      
    M3 = griddata((x1, y1, z1), s1, (grid_x1, grid_y1, grid_z1), method='nearest')
    M3[:,0:1,:] = np.nan; M3[:,2:15,:] = np.nan
      
    M4 = griddata((x1, y1, z1), s1, (grid_x1, grid_y1, grid_z1), method='nearest')
    M4[:,:,0:4] = np.nan; M4[:,:,5:15] = np.nan
      
    fig = plt.figure(figsize=(7.125,6.75))
    ax1 = fig.add_subplot(221, projection='3d')
    ax2 = fig.add_subplot(222, projection='3d')
    ax3 = fig.add_subplot(223, projection='3d')
    ax4 = fig.add_subplot(224, projection='3d')
      
    ax1.scatter(grid_x1, grid_y1, grid_z1, c=M.flatten(), vmin=np.min(M), vmax=np.max(M)*0.75, cmap='jet_r', alpha=0.5)
    ax1.set_xlabel('Lithosphere elastic thickness (km)'); ax1.set_ylabel('Range transport coefficient\n(m$^{2}$.yr$^{-1}$)', labelpad=10); ax1.set_zlabel('Basin transport coefficient\n(m$^{2}$.yr$^{-1}$)', labelpad=10)
      
    ax2.scatter(grid_x1, grid_y1, grid_z1, c=M2.flatten(), vmin=np.min(M), vmax=np.max(M)*0.75, cmap='jet_r', alpha=1)
    ax2.set_xlim(15000, 40000)
    ax2.set_xlabel('Lithosphere elastic thickness (km)'); ax2.set_ylabel('Range transport coefficient\n(m$^{2}$.yr$^{-1}$)', labelpad=10); ax2.set_zlabel('Basin transport coefficient\n(m$^{2}$.yr$^{-1}$)', labelpad=10)
      
    ax3.scatter(grid_x1, grid_y1, grid_z1, c=M3.flatten(), vmin=np.min(M), vmax=np.max(M)*0.75, cmap='jet_r', alpha=1)
    ax3.set_ylim(0, 5000)
    ax3.set_xlabel('Lithosphere elastic thickness (km)'); ax3.set_ylabel('Range transport coefficient\n(m$^{2}$.yr$^{-1}$)', labelpad=10); ax3.set_zlabel('Basin transport coefficient\n(m$^{2}$.yr$^{-1}$)', labelpad=10)
      
    S4 = ax4.scatter(grid_x1, grid_y1, grid_z1, c=M4.flatten(), vmin=np.min(M), vmax=np.max(M)*0.75, cmap='jet_r', alpha=1)
    ax4.set_zlim(0, 50000)
    ax4.set_xlabel('Lithosphere elastic thickness (km)'); ax4.set_ylabel('Range transport coefficient\n(m$^{2}$.yr$^{-1}$)', labelpad=10); ax4.set_zlabel('Basin transport coefficient\n(m$^{2}$.yr$^{-1}$)', labelpad=10)
      
    cax = fig.add_axes([0.1, 0.075, 0.8, 0.02])
    fig.colorbar(S4, label='Log(\u03C7)', orientation='horizontal', cax=cax)
    Te_label = ['15','20','25','30','35','40']
    ax1.set_xticklabels(Te_label); ax2.set_xticklabels(Te_label); ax3.set_xticklabels(Te_label); ax4.set_xticklabels(Te_label)
      
    fig.subplots_adjust(left=0.020, bottom=0.2, right=0.925, top=0.975, wspace=0.25, hspace=0.20)
    fig.savefig(path+'Sup_Figure2'+'.pdf', dpi=1200)
    fig.savefig(path+'Sup_Figure2'+'.jpg', dpi=1200)

#%%=========================================================================%%#
def misfit_vs_parameters_3D_2(kappab_list, kappar_list, Te_list, global_misfit_list, Te_min, Te_max, kappar_min, kappar_max, kappab_min, kappab_max, path):
  
  x1 = Te_list; y1 = kappar_list; z1 = kappab_list; s1 = global_misfit_list
  x2 = np.empty((0)); y2 = np.empty((0)); z2 = np.empty((0)); s2 = np.empty((0));
  x3 = np.empty((0)); y3 = np.empty((0)); z3 = np.empty((0)); s3 = np.empty((0));
  x4 = np.empty((0)); y4 = np.empty((0)); z4 = np.empty((0)); s4 = np.empty((0));
  x5 = np.empty((0)); y5 = np.empty((0)); z5 = np.empty((0)); s5 = np.empty((0));
  x6 = np.empty((0)); y6 = np.empty((0)); z6 = np.empty((0)); s6 = np.empty((0));
  x7 = np.empty((0)); y7 = np.empty((0)); z7 = np.empty((0)); s7 = np.empty((0));
  
  fig = plt.figure(figsize=(7.125,9))
  ax1 = fig.add_subplot(321, projection='3d')
  ax2 = fig.add_subplot(322, projection='3d')
  ax3 = fig.add_subplot(323, projection='3d')
  ax4 = fig.add_subplot(324, projection='3d')
  ax5 = fig.add_subplot(325, projection='3d')
  ax6 = fig.add_subplot(326, projection='3d')
  
  for i in range (0, np.shape(s1)[0]):
    if s1[i] <= 0.301:    
      x2 = np.append(x2, x1[i])
      y2 = np.append(y2, y1[i])
      z2 = np.append(z2, z1[i])
      s2 = np.append(s2, s1[i])
    if 0.301 < s1[i] <= 0.6:    
      x3 = np.append(x3, x1[i])
      y3 = np.append(y3, y1[i])
      z3 = np.append(z3, z1[i])
      s3 = np.append(s3, s1[i])
    if 0.6 < s1[i] <= 0.9:    
      x4 = np.append(x4, x1[i])
      y4 = np.append(y4, y1[i])
      z4 = np.append(z4, z1[i])
      s4 = np.append(s4, s1[i])
    if 0.9 < s1[i] <= 1.2:    
      x5 = np.append(x5, x1[i])
      y5 = np.append(y5, y1[i])
      z5 = np.append(z5, z1[i])
      s5 = np.append(s5, s1[i])
    if 1.2 < s1[i] <= 1.5:    
      x6 = np.append(x6, x1[i])
      y6 = np.append(y6, y1[i])
      z6 = np.append(z6, z1[i])
      s6 = np.append(s6, s1[i])
    if s1[i] > 1.5:    
      x7 = np.append(x7, x1[i])
      y7 = np.append(y7, y1[i])
      z7 = np.append(z7, z1[i])
      s7 = np.append(s7, s1[i])
             
  S = ax1.scatter(x2, y2, z2, c=s2, s=30, vmin=np.min(s1), vmax=np.max(s1), cmap='jet_r',alpha=0.66)
  ax1.set_xlim(15000, 40000); ax1.set_ylim(100, 5000); ax1.set_zlim(1000, 50000)
  ax1.set_xlabel('Lithosphere elastic thickness (km)'); ax1.set_ylabel('Range transport coefficient\n(m$^{2}$.yr$^{-1}$)', labelpad=10); ax1.set_zlabel('Basin transport coefficient\n(m$^{2}$.yr$^{-1}$)', labelpad=10)
  
  ax2.scatter(x3, y3, z3, c=s3, s=30, vmin=np.min(s1), vmax=np.max(s1), cmap='jet_r',alpha=0.66)
  ax2.set_xlim(15000, 40000); ax2.set_ylim(100, 5000); ax2.set_zlim(1000, 50000)
  ax2.set_xlabel('Lithosphere elastic thickness (km)'); ax2.set_ylabel('Range transport coefficient\n(m$^{2}$.yr$^{-1}$)', labelpad=10); ax2.set_zlabel('Basin transport coefficient\n(m$^{2}$.yr$^{-1}$)', labelpad=10)
  
  ax3.scatter(x4, y4, z4, c=s4, s=30, vmin=np.min(s1), vmax=np.max(s1), cmap='jet_r',alpha=0.66)
  ax3.set_xlim(15000, 40000); ax3.set_ylim(100, 5000); ax3.set_zlim(1000, 50000)
  ax3.set_xlabel('Lithosphere elastic thickness (km)'); ax3.set_ylabel('Range transport coefficient\n(m$^{2}$.yr$^{-1}$)', labelpad=10); ax3.set_zlabel('Basin transport coefficient\n(m$^{2}$.yr$^{-1}$)', labelpad=10)
  
  ax4.scatter(x5, y5, z5, c=s5, s=30, vmin=np.min(s1), vmax=np.max(s1), cmap='jet_r',alpha=0.66)
  ax4.set_xlim(15000, 40000); ax4.set_ylim(100, 5000); ax4.set_zlim(1000, 50000)
  ax4.set_xlabel('Lithosphere elastic thickness (km)'); ax4.set_ylabel('Range transport coefficient\n(m$^{2}$.yr$^{-1}$)', labelpad=10); ax4.set_zlabel('Basin transport coefficient\n(m$^{2}$.yr$^{-1}$)', labelpad=10)
  
  ax5.scatter(x6, y6, z6, c=s6, s=30, vmin=np.min(s1), vmax=np.max(s1), cmap='jet_r',alpha=0.66)
  ax5.set_xlim(15000, 40000); ax5.set_ylim(100, 5000); ax5.set_zlim(1000, 50000)
  ax5.set_xlabel('Lithosphere elastic thickness (km)'); ax5.set_ylabel('Range transport coefficient\n(m$^{2}$.yr$^{-1}$)', labelpad=10); ax5.set_zlabel('Basin transport coefficient\n(m$^{2}$.yr$^{-1}$)', labelpad=10)
  
  ax6.scatter(x7, y7, z7, c=s7, s=30, vmin=np.min(s1), vmax=np.max(s1), cmap='jet_r',alpha=0.66)
  ax6.set_xlim(15000, 40000); ax6.set_ylim(100, 5000); ax6.set_zlim(1000, 50000)
  ax6.set_xlabel('Lithosphere elastic thickness (km)'); ax6.set_ylabel('Range transport coefficient\n(m$^{2}$.yr$^{-1}$)', labelpad=10); ax6.set_zlabel('Basin transport coefficient\n(m$^{2}$.yr$^{-1}$)', labelpad=10)
  
  cax = fig.add_axes([0.1, 0.05, 0.8, 0.02])
  fig.colorbar(S, label='Log(\u03C7)', orientation='horizontal',cax=cax)
  Te_label = ['15','20','25','30','35','40']
  ax1.set_xticklabels(Te_label); ax2.set_xticklabels(Te_label); ax3.set_xticklabels(Te_label); ax4.set_xticklabels(Te_label); ax5.set_xticklabels(Te_label); ax6.set_xticklabels(Te_label)
  
  fig.subplots_adjust(left=0.020, bottom=0.15, right=0.925, top=0.975, wspace=0.25, hspace=0.20)
  fig.savefig(path+'Sup_Figure3'+'.pdf', dpi=1200)
  fig.savefig(path+'Sup_Figure3'+'.jpg', dpi=1200)

#%%=========================================================================%%#
def figure_misfit_vs_time_convergence_decrease(T4_time_steps_list, global_misfit, Te_list, kappar_list, kappab_list, path):
  
  fig = plt.figure(figsize=(7.125,2.75))
  ax1 = fig.add_subplot(121)
  ax2 = fig.add_subplot(122)
  
  order = np.argsort(global_misfit)[::-1]
  
  S1 = ax1.scatter(T4_time_steps_list, global_misfit, c=Te_list, s=10)
  cb1 = fig.colorbar(S1, label='Lithosphere elastic thickness (km)', ax=ax1)
  
  S2 = ax2.scatter(T4_time_steps_list, global_misfit, c=(np.log(kappab_list/kappar_list)), s=10)
  cb2 = fig.colorbar(S2, label='Log(transport coefficient ratio)', ax=ax2)
  
  cb1.set_ticks([15000, 20000, 25000, 30000, 35000])
  cb1.set_ticklabels(['15', '20', '25', '30', '35'])
  
  ax1.set_xlabel('Time of convergence velocity decrease (Myrs)'); ax2.set_xlabel('Time of convergence velocity decrease (Myrs)')
  ax1.set_ylabel('Log(\u03C7)'); ax2.set_ylabel('Log(\u03C7)')
  T4_label = ['2','3','4','5','6','7']
  ax1.set_xticklabels(T4_label); ax2.set_xticklabels(T4_label)
  
  fig.tight_layout()
  fig.savefig(path + 'Sup_Figure1.pdf', dpi=1200)
  fig.savefig(path + 'Sup_Figure1.jpg', dpi=1200)

#%%=========================================================================%%#  
def figure_cross_section_range_results(result_Hr_list, result_wr_list, time, result_time_list, tt, path):
    
    meanElevation_Hr = result_Hr_list.mean(axis=0)
    percentile_4p55_Elevation_Hr = np.percentile(result_Hr_list, 4.55, axis=0)
    percentile_95p45_Elevation_Hr = np.percentile(result_Hr_list, 95.45, axis=0)
    
    mean_Deflexion_wr = result_wr_list.mean(axis=0)
    percentile_4p55_Deflexion_wr = np.percentile(result_wr_list, 4.55, axis=0)
    percentile_95p45_Deflexion_wr = np.percentile(result_wr_list, 95.45, axis=0)
    
    fig = plt.figure(figsize = (4.82, 3))
    ax1 = fig.add_subplot(111)
 
    ax1.hist2d(result_time_list.flatten(), result_Hr_list.flatten(), bins = (100, 100), cmap = plt.cm.jet, norm = LogNorm(), zorder = 2)
    ax1.plot(time, meanElevation_Hr, linewidth = 2, color = 'k', linestyle = '-', zorder = 3)
    ax1.plot(time, percentile_4p55_Elevation_Hr, linewidth = 2, color = 'k', linestyle = '--', zorder = 3)
    ax1.plot(time, percentile_95p45_Elevation_Hr, linewidth = 2, color = 'k', linestyle = '--', zorder = 3)    
    
    ax1.hist2d(result_time_list.flatten(), -result_wr_list.flatten(), bins = (100, 100), cmap = plt.cm.jet, norm = LogNorm(), zorder = 2)
    ax1.plot(time, -mean_Deflexion_wr, linewidth = 2, color = 'k', linestyle = '-', zorder = 3)
    ax1.plot(time, -percentile_4p55_Deflexion_wr, linewidth = 2, color = 'k', linestyle = '--', zorder = 3)
    ax1.plot(time, -percentile_95p45_Deflexion_wr, linewidth = 2, color = 'k', linestyle = '--', zorder = 3)
    
    ax1.yaxis.grid(which='major', linestyle='--', color='k', alpha=0.5, zorder=0); ax1.xaxis.grid(which='major', linestyle='--', color='k', alpha=0.5, zorder=0)
    ax1.axis([0, tt*1e-6, -percentile_95p45_Deflexion_wr.max()*1.15, percentile_95p45_Elevation_Hr.max()*1.5])
    ax1.set_xlabel('Time (Ma)'); ax1.set_ylabel('Range elevation and depth (m)')
    xlabel=['66','56','46','36','26','16','6','0']; ax1.set_xticklabels(xlabel);
    
    fig.tight_layout()
    #fig.savefig(path+'SUP_FIGURE4_range_cross-section_evolution.pdf', dpi=1200)
    #fig.savefig(path+'SUP_FIGURE4_range_cross-section_evolution-01.jpeg', dpi=1200)

#%%=========================================================================%%#  
def figure_cross_section_basin_results(result_Hb_list, result_wb_list, time, result_time_list, tt, path):
    
    meanElevation_Hb = result_Hb_list.mean(axis=0)
    percentile_4p55_Elevation_Hb = np.percentile(result_Hb_list, 4.55, axis=0)
    percentile_95p45_Elevation_Hb = np.percentile(result_Hb_list, 95.45, axis=0)
    
    mean_Deflexion_wb = result_wb_list.mean(axis=0)
    percentile_4p55_Deflexion_wb = np.percentile(result_wb_list, 4.55, axis=0)
    percentile_95p45_Deflexion_wb = np.percentile(result_wb_list, 95.45, axis=0)
    
    fig = plt.figure(figsize = (4.82, 3))
    ax1 = fig.add_subplot(111)
 
    ax1.hist2d(result_time_list.flatten(), result_Hb_list.flatten(), bins = (100, 100), cmap = plt.cm.jet, norm = LogNorm(), zorder = 2)
    ax1.plot(time, meanElevation_Hb, linewidth = 2, color = 'k', linestyle = '-', zorder = 3)
    ax1.plot(time, percentile_4p55_Elevation_Hb, linewidth = 2, color = 'k', linestyle = '--', zorder = 3)
    ax1.plot(time, percentile_95p45_Elevation_Hb, linewidth = 2, color = 'k', linestyle = '--', zorder = 3)    
    
    ax1.hist2d(result_time_list.flatten(), -result_wb_list.flatten(), bins = (100, 100), cmap = plt.cm.jet, norm = LogNorm(), zorder = 2)
    ax1.plot(time, -mean_Deflexion_wb, linewidth = 2, color = 'k', linestyle = '-', zorder = 3)
    ax1.plot(time, -percentile_4p55_Deflexion_wb, linewidth = 2, color = 'k', linestyle = '--', zorder = 3)
    ax1.plot(time, -percentile_95p45_Deflexion_wb, linewidth = 2, color = 'k', linestyle = '--', zorder = 3)
    
    ax1.yaxis.grid(which='major', linestyle='--', color='k', alpha=0.5, zorder=0); ax1.xaxis.grid(which='major', linestyle='--', color='k', alpha=0.5, zorder=0)
    ax1.axis([0, tt*1e-6, -percentile_95p45_Deflexion_wb.max()*1.15, percentile_95p45_Elevation_Hb.max()*1.5])
    ax1.set_xlabel('Time (Ma)'); ax1.set_ylabel('Range elevation and depth (m)')
    xlabel=['66','56','46','36','26','16','6','0']; ax1.set_xticklabels(xlabel)
    
    fig.tight_layout()
    #fig.savefig(path+'SUP_FIGURE5_basin_cross-section_evolution'+'.pdf', dpi=1200)
    #fig.savefig(path+'SUP_FIGURE5_basin_cross-section_evolution-01.jpeg', dpi=1200)

#%%=========================================================================%%#
def figure_range_volume_change_results(result_Hr_list, result_wr_list, Lr, result_Fr_list, result_Fe_list, result_time_list, time, tt, rhoc, path):
    
    massRange = (result_Hr_list + result_wr_list) * Lr * rhoc
    mean_massRange = massRange.mean(axis=0)
    percentile_4p55_massRange = np.percentile(massRange, 4.55, axis=0)
    percentile_95p45_massRange = np.percentile(massRange, 95.45, axis=0)
    
    massRangeChange = (result_Fr_list - result_Fe_list) * 1e-5
    mean_massRangeChange = massRangeChange.mean(axis=0)
    percentile_4p55_massRangeChange = np.percentile(massRangeChange, 4.55, axis=0) 
    percentile_95p45_massRangeChange = np.percentile(massRangeChange, 95.45, axis=0)
    
    fig = plt.figure(figsize = (4.82, 6))
    ax1 = fig.add_subplot(211)
    ax2 = fig.add_subplot(212)
    
    ax1.hist2d(result_time_list.flatten(), massRange.flatten(), bins = (100, 100), cmap = plt.cm.jet, norm = LogNorm(), zorder = 2)
    ax1.plot(time, mean_massRange, linewidth = 2, color = 'k', linestyle = '-', zorder = 3)
    ax1.plot(time, percentile_4p55_massRange, linewidth = 2, color = 'k', linestyle = '--', zorder = 3)
    ax1.plot(time, percentile_95p45_massRange, linewidth = 2, color = 'k', linestyle = '--', zorder = 3)
    
    ax2.hist2d(result_time_list.flatten(), massRangeChange.flatten(), bins = (100, 100), cmap = plt.cm.jet, norm = LogNorm(), zorder = 2)
    ax2.plot(time, mean_massRangeChange, linewidth = 2, color = 'k', linestyle = '-', zorder = 3)
    ax2.plot(time, percentile_4p55_massRangeChange, linewidth = 2, color = 'k', linestyle = '--', zorder = 3)
    ax2.plot(time, percentile_95p45_massRangeChange, linewidth = 2, color = 'k', linestyle = '--', zorder = 3)
    
    ax1.axis([0, tt*1e-6, 0, percentile_95p45_massRange.max()*1.05])
    ax1.set_xlabel('Time (Ma)'); ax1.set_ylabel('Range mass (kg)')
    ax1.yaxis.grid(which='major',linestyle='--',color='k',alpha=0.5,zorder=0); ax1.xaxis.grid(which='major',linestyle='--',color='k',alpha=0.5,zorder=0)
    xlabel=['66','56','46','36','26','16','6','0']; ax1.set_xticklabels(xlabel)

    ax2.axis([0, tt*1e-6, percentile_4p55_massRangeChange.min()*1.05, percentile_95p45_massRangeChange.max()*1.05])
    ax2.set_xlabel('Time (Ma)'); ax2.set_ylabel('Range mass change (kg.yr$^{-1}$)')
    ax2.yaxis.grid(which='major',linestyle='--',color='k',alpha=0.5,zorder=0); ax2.xaxis.grid(which='major',linestyle='--',color='k',alpha=0.5,zorder=0)
    xlabel=['66','56','46','36','26','16','6','0']; ax2.set_xticklabels(xlabel)  
    
    fig.tight_layout()
    #fig.savefig(path+'SUP_FIGURE6_range_mass_change_evolution'+'.pdf', dpi=1200)
    #fig.savefig(path+'SUP_FIGURE6_range_mass_change_evolution-01'+'.jpeg', dpi=1200)

#%%=========================================================================%%#
def figure_basin_volume_change_results(result_Hb_list, result_wb_list, result_Fe_list, result_Fb_list, time, tt, result_time_list, Lb, rhos, path):
    
    massBasin = (result_Hb_list + result_wb_list) * Lb * rhos
    mean_massBasin = massBasin.mean(axis=0)
    percentile_4p55_massBasin = np.percentile(massBasin, 4.55, axis=0)
    percentile_95p45_massBasin = np.percentile(massBasin, 95.45, axis=0)
    
    massBasinChange = (result_Fe_list - result_Fb_list) * 1e-5
    mean_massBasinChange = massBasinChange.mean(axis=0)
    percentile_4p55_massBasinChange = np.percentile(massBasinChange, 4.55, axis=0) 
    percentile_95p45_massBasinChange = np.percentile(massBasinChange, 95.45, axis=0)
    
    fig = plt.figure(figsize = (4.82, 6))
    ax1 = fig.add_subplot(211)
    ax2 = fig.add_subplot(212)
    
    ax1.hist2d(result_time_list.flatten(), massBasin.flatten(), bins = (100, 100), cmap = plt.cm.jet, norm = LogNorm(), zorder = 2)
    ax1.plot(time, mean_massBasin, linewidth = 2, color = 'k', linestyle = '-', zorder = 3)
    ax1.plot(time, percentile_4p55_massBasin, linewidth = 2, color = 'k', linestyle = '--', zorder = 3)
    ax1.plot(time, percentile_95p45_massBasin, linewidth = 2, color = 'k', linestyle = '--', zorder = 3)
    
    ax2.hist2d(result_time_list.flatten(), massBasinChange.flatten(), bins = (100, 100), cmap = plt.cm.jet, norm = LogNorm(), zorder = 2)
    ax2.plot(time, mean_massBasinChange, linewidth = 2, color = 'k', linestyle = '-', zorder = 3)
    ax2.plot(time, percentile_4p55_massBasinChange, linewidth = 2, color = 'k', linestyle = '--', zorder = 3)
    ax2.plot(time, percentile_95p45_massBasinChange, linewidth = 2, color = 'k', linestyle = '--', zorder = 3)
    
    ax1.axis([0, tt*1e-6, 0, percentile_95p45_massBasin.max()*1.05])
    ax1.set_xlabel('Time (Ma)'); ax1.set_ylabel('Basin mass (kg)')
    ax1.yaxis.grid(which='major',linestyle='--',color='k',alpha=0.5,zorder=0); ax1.xaxis.grid(which='major',linestyle='--',color='k',alpha=0.5,zorder=0)
    xlabel=['66','56','46','36','26','16','6','0']; ax1.set_xticklabels(xlabel)

    ax2.axis([0, tt*1e-6, percentile_4p55_massBasinChange.min()*1.05, percentile_95p45_massBasinChange.max()*1.05])
    ax2.set_xlabel('Time (Ma)'); ax2.set_ylabel('Basin mass change (kg.yr$^{-1}$)')
    ax2.yaxis.grid(which='major',linestyle='--',color='k',alpha=0.5,zorder=0); ax2.xaxis.grid(which='major',linestyle='--',color='k',alpha=0.5,zorder=0)
    xlabel=['66','56','46','36','26','16','6','0']; ax2.set_xticklabels(xlabel)
    
    fig.tight_layout()
    #fig.savefig(path+'SUP_FIGURE7_basin_mass_change_evolution'+'.pdf', dpi=1200)
    #fig.savefig(path+'SUP_FIGURE7_basin_mass_change_evolution-01'+'.jpeg', dpi=1200)
