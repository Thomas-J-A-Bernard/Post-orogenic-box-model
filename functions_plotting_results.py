"""
Created on Wed Mar 20 09:45:08 2019
@author: Thomas Bernard (Python 3.6)
Description:
functions to plot results from the Tucker code with Monte Carlo sampling approch
"""

import matplotlib.pyplot as plt
import numpy as np
#import seaborn as sns
#import pandas as pd
from matplotlib.colors import LogNorm
import matplotlib.patches as patches
import matplotlib.lines as lines
#from scipy.interpolate import griddata
#import scipy.ndimage
#from scipy.ndimage.filters import gaussian_filter
from mpl_toolkits.mplot3d import axes3d
import matplotlib

#%%=========================================================================%%#
def range_vs_basin_diffusivity_sampling_pattern(kappar_list, kappab_list):
  
  df0 = pd.DataFrame({"Range diffusivity (m$^{2}$.yr$^{-1}$)": kappar_list, "Basin diffusivity (m$^{2}$.yr$^{-1}$)": kappab_list})
  sns.jointplot("Range diffusivity (m$^{2}$.yr$^{-1}$)", "Basin diffusivity (m$^{2}$.yr$^{-1}$)", data = df0, marginal_kws = dict(bins = 15, rug = False), annot_kws = dict(stat = "r"), kind = "reg", color = "steelblue", height = 7.5, xlim = (100,5000), ylim=(10000,50000))
  plt.subplots_adjust(left=0.15, bottom=0.1, right=0.95, top=0.95, wspace=0.2)

  #plt.savefig("Z:\Post-orogenic sediment flux to continental margins\Modelling_Post_orogenic_system\Tucker_Code\Result_NorthernPyrenees_MC_Approach\Figure_1_range_vs_basin_diffusivity.png",dpi=1200, transparent = True)
  #plt.savefig("Z:\Post-orogenic sediment flux to continental margins\Modelling_Post_orogenic_system\Tucker_Code\Result_NorthernPyrenees_MC_Approach\Figure_1_range_vs_basin_diffusivity.pdf",dpi=1200, transparent = False)

#%%=========================================================================%%#
def range_vs_basin_responsetime_sampling_pattern(tr_list, tb_list):
  
  df1 = pd.DataFrame({"Range response time (yr)": tr_list, "Basin response time (yr)": tb_list})
  sns.jointplot("Range response time (yr)", "Basin response time (yr)", data = df1, marginal_kws = dict(bins = 15, rug = False), annot_kws = dict(stat = "r"), kind = "reg", color = "steelblue", height = 7.5, xlim = (0,max(tr_list)), ylim=(0,max(tb_list)))
  plt.subplots_adjust(left=0.15, bottom=0.1, right=0.95, top=0.95, wspace=0.2)

  #plt.savefig("Z:\Post-orogenic sediment flux to continental margins\Modelling_Post_orogenic_system\Tucker_Code\Result_NorthernPyrenees_MC_Approach\Figure_2_range_vs_basin_response_time.png",dpi=1200, transparent = True)
  #plt.savefig("Z:\Post-orogenic sediment flux to continental margins\Modelling_Post_orogenic_system\Tucker_Code\Result_NorthernPyrenees_MC_Approach\Figure_2_range_vs_basin_response_time.pdf",dpi=1200, transparent = False)

#%%=========================================================================%%#
def convergence_decrease_sampling_pattern(T4_time_steps_list, T3_time_steps):
  
  df2 = pd.DataFrame({ "Convergence decrease time (Myrs)": (T4_time_steps_list - T3_time_steps)/100})
  sns.jointplot("Convergence decrease time (Myrs)", "Convergence decrease time (Myrs)", data = df2, marginal_kws = dict(bins = 15, rug = False), annot_kws = dict(stat = "r"), kind = "scatter", color = "steelblue", height = 7.5, xlim = (2,8), ylim=(2,8))
  plt.subplots_adjust(left=0.15, bottom=0.1, right=0.95, top=0.95, wspace=0.2)

  #plt.savefig("Z:\Post-orogenic sediment flux to continental margins\Modelling_Post_orogenic_system\Tucker_Code\Result_NorthernPyrenees_MC_Approach\Figure_3_convergence_velocity_decrease.png",dpi=1200, transparent = True)
  #plt.savefig("Z:\Post-orogenic sediment flux to continental margins\Modelling_Post_orogenic_system\Tucker_Code\Result_NorthernPyrenees_MC_Approach\Figure_3_convergence_velocity_decrease.pdf",dpi=1200, transparent = False)

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
def misfit_vs_parameters_hexbin_version(misfit, label_misfit, gridsize, linewidths, kappar_list, kappab_list, Te_list, kappar_min, kappar_max, kappab_min, kappab_max, Te_min, Te_max, MC_iteration):
  
  kappa_ratio_list = kappab_list/kappar_list
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
  
  sc = ax4.hexbin(Te_list, kappa_ratio_list, C=misfit, gridsize=gridsize, cmap='rainbow_r', vmin=min(misfit), vmax=max(misfit), linewidths=linewidths, edgecolors='black')
  fig.colorbar(sc, label=label_misfit, aspect=30, ax=ax4)
  
  ax1.set_xlabel('Range diffusivity (m$^{2}$.yr$^{-1}$)'); ax1.set_ylabel('Basin diffusivity (m$^{2}$.yr$^{-1}$)')
  ax1.axis([kappar_min, kappar_max, kappab_min, kappab_max])
  
  ax2.set_xlabel('Lithosphere elastic thickness (m)'); ax2.set_ylabel('Range diffusivity (m$^{2}$.yr$^{-1}$)')
  ax2.axis([Te_min, Te_max, kappar_min, kappar_max])
  
  ax3.set_xlabel('Lithosphere elastic thickness (m)'); ax3.set_ylabel('Basin diffusivity (m$^{2}$.yr$^{-1}$)')
  ax3.axis([Te_min, Te_max, kappab_min, kappab_max])
  
  ax4.set_xlabel('Lithosphere elastic thickness (m)'); ax4.set_ylabel('Diffusivity ratio')
  ax4.axis([Te_min, Te_max, min(kappa_ratio_list), max(kappa_ratio_list)])
  
  fig.subplots_adjust(left=0.1, bottom=0.1, right=0.985, top=0.95, wspace=0.3, hspace=0.3)

#%%=========================================================================%%#
def misfit_vs_parameters_interpolated_1_version(kappab_list, kappar_list, Te_list, global_misfit, Te_min, Te_max, kappar_min, kappar_max, kappab_min, kappab_max):
 
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
  
  fig = plt.figure(figsize=(9,7))
  ax1 = fig.add_subplot(221); ax2 = fig.add_subplot(222); ax3 = fig.add_subplot(223); ax4 = fig.add_subplot(224)
  
  im1 = ax1.imshow(grid_z1.T, extent=(0,100,0,100), origin='lower', cmap='rainbow_r')
  ax1.contour(grid_z1.T, levels=[0.15, 0.25, 0.35, 0.45, 0.6], extent=(0,100,0,100), linewidths=1.0, colors='k')
  fig.colorbar(im1, label='Misfit', fraction=0.046, ax=ax1)
  
  im2 = ax2.imshow(grid_z2.T, extent=(0,100,0,100), origin='lower', cmap='rainbow_r')
  ax2.contour(grid_z2.T, levels=[0.15, 0.25, 0.35, 0.45, 0.6], extent=(0,100,0,100), linewidths=1.0, colors='k')
  fig.colorbar(im2, label='Misfit', fraction=0.046, ax=ax2)
  
  im3 = ax3.imshow(grid_z3.T, extent=(0,100,0,100), origin='lower', cmap='rainbow_r')
  ax3.contour(grid_z3.T, levels=[0.15, 0.25, 0.35, 0.45, 0.6], extent=(0,100,0,100), linewidths=1.0, colors='k')
  fig.colorbar(im3, label='Misfit', fraction=0.046, ax=ax3)
  
  im4 = ax4.imshow(grid_z4.T, extent=(0,100,0,100), origin='lower', cmap='rainbow_r')
  ax4.contour(grid_z4.T, levels=[0.15, 0.25, 0.35, 0.45, 0.6], extent=(0,100,0,100), linewidths=1.0, colors='k')
  fig.colorbar(im4, label='Misfit', fraction=0.046, ax=ax4)
  
  ax1.set_xlabel('Lithosphere elastic thickness (km)'); ax1.set_ylabel('Range diffusivity (m$^{2}$.yr$^{-1}$)')
  ax2.set_xlabel('Lithosphere elastic thickness (km)'); ax2.set_ylabel('Basin diffusivity (m$^{2}$.yr$^{-1}$)')
  ax3.set_xlabel('Range diffusivity (m$^{2}$.yr$^{-1}$)'); ax3.set_ylabel('Basin diffusivity (m$^{2}$.yr$^{-1}$)')
  ax4.set_xlabel('Lithosphere elastic thickness (km)'); ax4.set_ylabel('Diffusivity ratio')
  
  Te_label = ['15','20','25','30','35','40']
  ax1.set_xticklabels(Te_label); ax2.set_xticklabels(Te_label); ax4.set_xticklabels(Te_label)
  kappar_label = ['100','1080','2060','3040','4020','5000']
  ax3.set_xticklabels(kappar_label); ax1.set_yticklabels(kappar_label)
  kappab_label = ['1000','10800','20600','30400','40200','50000']
  ax2.set_yticklabels(kappab_label); ax3.set_yticklabels(kappab_label)
  ratio_label = ['0','100','200','300','400','500']
  ax4.set_yticklabels(ratio_label)
  
  fig.subplots_adjust(left=0.1, bottom=0.085, right=0.925, top=0.975, wspace=0.45, hspace=0.25)

#%%=========================================================================%%#
def misfit_vs_parameters_interpolated_2_version(kappab_list, kappar_list, Te_list, global_misfit, Te_min, Te_max, kappar_min, kappar_max, kappab_min, kappab_max):
  
  kappa_ratio_list = kappab_list/kappar_list
  
  x1 = Te_list; y1 = kappar_list; z1 = global_misfit
  grid_x1, grid_y1 = np.mgrid[Te_min:Te_max:100j, kappar_min:kappar_max:100j]
  grid_z1 = griddata((x1, y1), z1, (grid_x1, grid_y1), method='nearest')
  grid_z1 = gaussian_filter(grid_z1, 4)
  
  x2 = Te_list; y2 = kappab_list; z2 = global_misfit
  grid_x2, grid_y2 = np.mgrid[Te_min:Te_max:100j, kappab_min:kappab_max:100j]
  grid_z2 = griddata((x2, y2), z2, (grid_x2, grid_y2), method='nearest')
  grid_z2 = gaussian_filter(grid_z2, 4)
  
  x3 = kappar_list; y3 = kappab_list; z3 = global_misfit
  grid_x3, grid_y3 = np.mgrid[kappar_min:kappar_max:100j, kappab_min:kappab_max:100j]
  grid_z3 = griddata((x3, y3), z3, (grid_x3, grid_y3), method='nearest')
  grid_z3 = gaussian_filter(grid_z3, 4)
  
  x4 = Te_list; y4 = kappa_ratio_list; z4 = global_misfit
  grid_x4, grid_y4 = np.mgrid[Te_min:Te_max:100j, 0:500:100j]
  grid_z4 = griddata((x4, y4), z4, (grid_x4, grid_y4), method='nearest')
  grid_z4 = gaussian_filter(grid_z4, 4)
  
  fig = plt.figure(figsize=(9,7))
  ax1 = fig.add_subplot(221); ax2 = fig.add_subplot(222); ax3 = fig.add_subplot(223); ax4 = fig.add_subplot(224)
  
  im1 = ax1.imshow(grid_z1.T, extent=(0,100,0,100), origin='lower', cmap='rainbow_r')
  ax1.contour(grid_z1.T, levels=[0.25, 0.275, 0.3, 0.325, 0.36, 0.39, 0.475, 0.55], extent=(0,100,0,100), linewidths=1.0, colors='k')
  fig.colorbar(im1, label='Misfit', fraction=0.046, ax=ax1)
  
  im2 = ax2.imshow(grid_z2.T, extent=(0,100,0,100), origin='lower', cmap='rainbow_r')
  ax2.contour(grid_z2.T, levels=[0.25, 0.275, 0.3, 0.325, 0.36, 0.39, 0.475, 0.55], extent=(0,100,0,100), linewidths=1.0, colors='k')
  fig.colorbar(im2, label='Misfit', fraction=0.046, ax=ax2)
  
  im3 = ax3.imshow(grid_z3.T, extent=(0,100,0,100), origin='lower', cmap='rainbow_r')
  ax3.contour(grid_z3.T, levels=[0.25, 0.275, 0.3, 0.325, 0.36, 0.39, 0.475, 0.55], extent=(0,100,0,100), linewidths=1.0, colors='k')
  fig.colorbar(im3, label='Misfit', fraction=0.046, ax=ax3)
  
  im4 = ax4.imshow(grid_z4.T, extent=(0,100,0,100), origin='lower', cmap='rainbow_r')
  ax4.contour(grid_z4.T, levels=[0.25, 0.275, 0.3, 0.325, 0.36, 0.39, 0.475, 0.55], extent=(0,100,0,100), linewidths=1.0, colors='k')
  fig.colorbar(im4, label='Misfit', fraction=0.046, ax=ax4)
  
  ax1.set_xlabel('Lithosphere elastic thickness (km)'); ax1.set_ylabel('Range diffusivity (m$^{2}$.yr$^{-1}$)')
  ax2.set_xlabel('Lithosphere elastic thickness (km)'); ax2.set_ylabel('Basin diffusivity (m$^{2}$.yr$^{-1}$)')
  ax3.set_xlabel('Range diffusivity (m$^{2}$.yr$^{-1}$)'); ax3.set_ylabel('Basin diffusivity (m$^{2}$.yr$^{-1}$)')
  ax4.set_xlabel('Lithosphere elastic thickness (km)'); ax4.set_ylabel('Diffusivity ratio')
  
  Te_label = ['15','20','25','30','35','40']
  ax1.set_xticklabels(Te_label); ax2.set_xticklabels(Te_label); ax4.set_xticklabels(Te_label)
  kappar_label = ['100','1080','2060','3040','4020','5000']
  ax3.set_xticklabels(kappar_label); ax1.set_yticklabels(kappar_label)
  kappab_label = ['1000','10800','20600','30400','40200','50000']
  ax2.set_yticklabels(kappab_label); ax3.set_yticklabels(kappab_label)
  ratio_label = ['0','100','200','300','400','500']
  ax4.set_yticklabels(ratio_label)
  
  fig.subplots_adjust(left=0.1, bottom=0.085, right=0.925, top=0.975, wspace=0.45, hspace=0.25)

#%%=========================================================================%%#
def misfit_vs_parameters_interpolated_3_version(kappab_list, kappar_list, Te_list, global_misfit, Te_min, Te_max, kappar_min, kappar_max, kappab_min, kappab_max):

  kappa_ratio_list = kappab_list/kappar_list
  
  x1 = Te_list; y1 = kappar_list; z1 = global_misfit
  grid_x1, grid_y1 = np.mgrid[Te_min:Te_max:100j, kappar_min:kappar_max:100j]
  grid_z1 = griddata((x1, y1), z1, (grid_x1, grid_y1), method='nearest')
  grid_z1 = gaussian_filter(grid_z1, 4)
  
  x2 = Te_list; y2 = kappab_list; z2 = global_misfit
  grid_x2, grid_y2 = np.mgrid[Te_min:Te_max:100j, kappab_min:kappab_max:100j]
  grid_z2 = griddata((x2, y2), z2, (grid_x2, grid_y2), method='nearest')
  grid_z2 = gaussian_filter(grid_z2, 4)
  
  x3 = kappar_list; y3 = kappab_list; z3 = global_misfit
  grid_x3, grid_y3 = np.mgrid[kappar_min:kappar_max:100j, kappab_min:kappab_max:100j]
  grid_z3 = griddata((x3, y3), z3, (grid_x3, grid_y3), method='nearest')
  grid_z3 = gaussian_filter(grid_z3, 4)
  
  x4 = Te_list; y4 = kappa_ratio_list; z4 = global_misfit
  grid_x4, grid_y4 = np.mgrid[Te_min:Te_max:100j, 0:500:100j]
  grid_z4 = griddata((x4, y4), z4, (grid_x4, grid_y4), method='nearest')
  grid_z4 = gaussian_filter(grid_z4, 4)
  
  fig = plt.figure(figsize=(9,7))
  ax1 = fig.add_subplot(221); ax2 = fig.add_subplot(222); ax3 = fig.add_subplot(223); ax4 = fig.add_subplot(224)
  
  ax1.imshow(grid_z1.T, extent=(0,100,0,100), origin='lower', cmap='gray_r')
  co1 = ax1.contour(grid_z1.T, levels=[0.25, 0.275, 0.3, 0.325, 0.36, 0.39, 0.475, 0.55], cmap='rainbow_r', extent=(0,100,0,100), linewidths=1.25)
  norm1= matplotlib.colors.Normalize(vmin=co1.cvalues.min(), vmax=co1.cvalues.max())
  sm1 = plt.cm.ScalarMappable(norm=norm1, cmap = co1.cmap)
  sm1.set_array([])
  fig.colorbar(sm1, ticks=co1.levels, label='Misfit', fraction=0.046, ax=ax1)
  
  ax2.imshow(grid_z2.T, extent=(0,100,0,100), origin='lower', cmap='gray_r')
  co2 = ax2.contour(grid_z2.T, levels=[0.25, 0.275, 0.3, 0.325, 0.36, 0.39, 0.475, 0.55], cmap='rainbow_r', extent=(0,100,0,100), linewidths=1.25)
  norm2= matplotlib.colors.Normalize(vmin=co2.cvalues.min(), vmax=co2.cvalues.max())
  sm2 = plt.cm.ScalarMappable(norm=norm2, cmap = co2.cmap)
  sm2.set_array([])
  fig.colorbar(sm2, ticks=co1.levels, label='Misfit', fraction=0.046, ax=ax2)
  
  ax3.imshow(grid_z3.T, extent=(0,100,0,100), origin='lower', cmap='gray_r')
  co3 = ax3.contour(grid_z3.T, levels=[0.25, 0.275, 0.3, 0.325, 0.36, 0.39, 0.475, 0.55], cmap='rainbow_r', extent=(0,100,0,100), linewidths=1.25)
  norm3= matplotlib.colors.Normalize(vmin=co3.cvalues.min(), vmax=co3.cvalues.max())
  sm3 = plt.cm.ScalarMappable(norm=norm3, cmap = co3.cmap)
  sm3.set_array([])
  fig.colorbar(sm3, ticks=co3.levels, label='Misfit', fraction=0.046, ax=ax3)
  
  ax4.imshow(grid_z4.T, extent=(0,100,0,100), origin='lower', cmap='gray_r')
  co4 = ax4.contour(grid_z4.T, levels=[0.25, 0.275, 0.3, 0.325, 0.36, 0.39, 0.475, 0.55], cmap='rainbow_r', extent=(0,100,0,100), linewidths=1.25)
  norm4= matplotlib.colors.Normalize(vmin=co4.cvalues.min(), vmax=co4.cvalues.max())
  sm4 = plt.cm.ScalarMappable(norm=norm4, cmap = co4.cmap)
  sm4.set_array([])
  fig.colorbar(sm4, ticks=co3.levels, label='Misfit', fraction=0.046, ax=ax4)
  
  ax1.set_xlabel('Lithosphere elastic thickness (km)'); ax1.set_ylabel('Range diffusivity (m$^{2}$.yr$^{-1}$)')
  ax2.set_xlabel('Lithosphere elastic thickness (km)'); ax2.set_ylabel('Basin diffusivity (m$^{2}$.yr$^{-1}$)')
  ax3.set_xlabel('Range diffusivity (m$^{2}$.yr$^{-1}$)'); ax3.set_ylabel('Basin diffusivity (m$^{2}$.yr$^{-1}$)')
  ax4.set_xlabel('Lithosphere elastic thickness (km)'); ax4.set_ylabel('Diffusivity ratio')
  
  Te_label = ['15','20','25','30','35','40']
  ax1.set_xticklabels(Te_label); ax2.set_xticklabels(Te_label); ax4.set_xticklabels(Te_label)
  kappar_label = ['100','1080','2060','3040','4020','5000']
  ax3.set_xticklabels(kappar_label); ax1.set_yticklabels(kappar_label)
  kappab_label = ['1000','10800','20600','30400','40200','50000']
  ax2.set_yticklabels(kappab_label); ax3.set_yticklabels(kappab_label)
  ratio_label = ['0','100','200','300','400','500']
  ax4.set_yticklabels(ratio_label)
  
  fig.subplots_adjust(left=0.1, bottom=0.085, right=0.925, top=0.975, wspace=0.45, hspace=0.25)
  
#%%=========================================================================%%#
def range_elevation_vs_time(result_Hr_list, time, tt):
  
  meanElevation_Hr = result_Hr_list.mean(axis=0)
  percentile_4p55_Elevation_Hr = np.percentile(result_Hr_list, 4.55, axis=0)
  percentile_95p45_Elevation_Hr = np.percentile(result_Hr_list, 95.45, axis=0)
  
  fig = plt.figure(figsize = (8,4))
  ax1 = fig.add_subplot(111)

  ax1.plot(time, result_Hr_list.T, color = 'steelblue', linestyle = '-.')
  ax1.plot(time, meanElevation_Hr, linewidth = 2, color = 'k', linestyle = '-')
  ax1.plot(time, percentile_4p55_Elevation_Hr, linewidth = 2, color = 'k', linestyle = '--')
  ax1.plot(time, percentile_95p45_Elevation_Hr, linewidth = 2, color = 'k', linestyle = '--')

  ax1.axis([0, tt*1e-6, result_Hr_list.min()*1.05, result_Hr_list.max()*1.05])
  ax1.set_xlabel('Time (Myr)'); ax1.set_ylabel('Range elevation (m)')
  ax1.yaxis.grid(which='major',linestyle='--',color='k',alpha=0.5,zorder=0); ax1.xaxis.grid(which='major',linestyle='--',color='k',alpha=0.5,zorder=0)
  fig.subplots_adjust(left=0.1, bottom=0.15, right=0.95, top=0.95, wspace=0.2)

  #fig.savefig("Z:\Post-orogenic sediment flux to continental margins\Modelling_Post_orogenic_system\Tucker_Code\Result_NorthernPyrenees_MC_Approach\Figure_4_range_elevation_evolution.png",dpi=1200, transparent = True)
  #fig.savefig("Z:\Post-orogenic sediment flux to continental margins\Modelling_Post_orogenic_system\Tucker_Code\Result_NorthernPyrenees_MC_Approach\Figure_4_range_elevation_evolution.pdf",dpi=1200, transparent = False)

#%%=========================================================================%%#
def basin_elevation_vs_time(result_Hb_list, time, tt):
  
  meanElevation_Hb = result_Hb_list.mean(axis=0)
  percentile_4p55_Elevation_Hb = np.percentile(result_Hb_list, 4.55, axis=0)
  percentile_95p45_Elevation_Hb = np.percentile(result_Hb_list, 95.45, axis=0)
  
  fig = plt.figure(figsize = (8,4))
  ax1 = fig.add_subplot(111)

  ax1.plot(time, result_Hb_list.T, color = 'steelblue', linestyle = '-.')
  ax1.plot(time, meanElevation_Hb, linewidth = 2, color = 'k', linestyle = '-')
  ax1.plot(time, percentile_4p55_Elevation_Hb, linewidth = 2, color = 'k', linestyle = '--')
  ax1.plot(time, percentile_95p45_Elevation_Hb, linewidth = 2, color = 'k', linestyle = '--')

  ax1.axis([0, tt*1e-6, result_Hb_list.min()*1.05, result_Hb_list.max()*1.05])
  ax1.set_xlabel('Time (Myr)'); ax1.set_ylabel('Basin elevation (m)')
  ax1.yaxis.grid(which = 'major', linestyle = '--', color = 'k', alpha = 0.5, zorder=0); ax1.xaxis.grid(which = 'major', linestyle = '--', color = 'k', alpha = 0.5, zorder = 1)
  fig.subplots_adjust(left = 0.1, bottom = 0.15, right = 0.95, top = 0.95, wspace = 0.2)

  #fig.savefig("Z:\Post-orogenic sediment flux to continental margins\Modelling_Post_orogenic_system\Tucker_Code\Result_NorthernPyrenees_MC_Approach\Figure_5_basin_elevation_evolution.png",dpi=1200, transparent = True)
  #fig.savefig("Z:\Post-orogenic sediment flux to continental margins\Modelling_Post_orogenic_system\Tucker_Code\Result_NorthernPyrenees_MC_Approach\Figure_5_basin_elevation_evolution.pdf",dpi=1200, transparent = False)

#%%=========================================================================%%#
def basin_deflexion_vs_time(result_wb_list, time, tt):
  
  meanDeflexion_wb = result_wb_list.mean(axis=0)
  percentile_4p55_Deflexion_wb = np.percentile(result_wb_list, 4.55, axis=0)
  percentile_95p45_Deflexion_wb = np.percentile(result_wb_list, 95.45, axis=0)
  
  fig = plt.figure(figsize = (8,4))
  ax1 = fig.add_subplot(111)

  ax1.plot(time, result_wb_list.T, color = 'steelblue', linestyle = '-.')
  ax1.plot(time, meanDeflexion_wb, linewidth = 2, color = 'k', linestyle = '-')
  ax1.plot(time, percentile_4p55_Deflexion_wb, linewidth = 2, color = 'k', linestyle = '--')
  ax1.plot(time, percentile_95p45_Deflexion_wb, linewidth = 2, color = 'k', linestyle = '--')

  ax1.axis([0, tt*1e-6, result_wb_list.min()*1.05, result_wb_list.max()*1.05])
  ax1.set_xlabel('Time (Myr)'); ax1.set_ylabel('Basin deflexion (m)')
  ax1.yaxis.grid(which = 'major', linestyle = '--', color = 'k', alpha = 0.5, zorder=0); ax1.xaxis.grid(which = 'major', linestyle = '--', color = 'k', alpha = 0.5, zorder = 1)
  plt.subplots_adjust(left = 0.1, bottom = 0.15, right = 0.95, top = 0.95, wspace = 0.2)

  #fig.savefig("Z:\Post-orogenic sediment flux to continental margins\Modelling_Post_orogenic_system\Tucker_Code\Result_NorthernPyrenees_MC_Approach\Figure_5_basin_elevation_evolution.png",dpi=1200, transparent = True)
  #fig.savefig("Z:\Post-orogenic sediment flux to continental margins\Modelling_Post_orogenic_system\Tucker_Code\Result_NorthernPyrenees_MC_Approach\Figure_5_basin_elevation_evolution.pdf",dpi=1200, transparent = False)

#%%=========================================================================%%# 
def density_range_elevation_vs_time(result_Hr_list, result_time_list, time, tt):
  
  meanElevation_Hr = result_Hr_list.mean(axis=0)
  percentile_4p55_Elevation_Hr = np.percentile(result_Hr_list, 4.55, axis=0)
  percentile_95p45_Elevation_Hr = np.percentile(result_Hr_list, 95.45, axis=0)
  
  fig = plt.figure(figsize = (8,4))
  ax1 = fig.add_subplot(111)
  
  ax1.hist2d(result_time_list.flatten(), result_Hr_list.flatten(), bins = (100, 100), cmap = plt.cm.jet, norm = LogNorm(), zorder = 2)
  ax1.plot(time, meanElevation_Hr, linewidth = 2, color = 'k', linestyle = '-', zorder = 3)
  ax1.plot(time, percentile_4p55_Elevation_Hr, linewidth = 2, color = 'k', linestyle = '--', zorder = 3)
  ax1.plot(time, percentile_95p45_Elevation_Hr, linewidth = 2, color = 'k', linestyle = '--', zorder = 3)
  
  ax1.axis([0, tt*1e-6, result_Hr_list.min()*1.05, result_Hr_list.max()*1.05])
  ax1.set_xlabel('Time (Myr)'); ax1.set_ylabel('Range elevation (m)')
  ax1.yaxis.grid(which = 'major', linestyle = '--', color = 'k', alpha = 0.5, zorder = 0); ax1.xaxis.grid(which = 'major', linestyle = '--', color = 'k', alpha = 0.5, zorder = 1)
  fig.subplots_adjust(left=0.1, bottom=0.15, right=0.95, top=0.95, wspace=0.2)
  
  #fig.savefig("Z:\Post-orogenic sediment flux to continental margins\Modelling_Post_orogenic_system\Tucker_Code\Result_NorthernPyrenees_MC_Approach\Figure_6_range_elevation_density.png", dpi=1200, transparent = True)
  #fig.savefig("Z:\Post-orogenic sediment flux to continental margins\Modelling_Post_orogenic_system\Tucker_Code\Result_NorthernPyrenees_MC_Approach\Figure_6_range_elevation_density.pdf", dpi=1200, transparant = False)

#%%=========================================================================%%#
def density_basin_elevation_vs_time(result_Hb_list, result_time_list, time, tt):
  
  meanElevation_Hb = result_Hb_list.mean(axis=0)
  percentile_4p55_Elevation_Hb = np.percentile(result_Hb_list, 4.55, axis=0)
  percentile_95p45_Elevation_Hb = np.percentile(result_Hb_list, 95.45, axis=0)
  
  fig = plt.figure(6, figsize = (8,4))
  ax1 = fig.add_subplot(111)
  
  ax1.hist2d(result_time_list.flatten(), result_Hb_list.flatten(), bins = (100, 100), cmap = plt.cm.jet, norm = LogNorm(), zorder = 2)
  ax1.plot(time, meanElevation_Hb, linewidth = 2, color = 'k', linestyle = '-', zorder = 3)
  ax1.plot(time, percentile_4p55_Elevation_Hb, linewidth = 2, color = 'k', linestyle = '--', zorder = 3)
  ax1.plot(time, percentile_95p45_Elevation_Hb, linewidth = 2, color = 'k', linestyle = '--', zorder = 3)
  
  ax1.axis([0, tt*1e-6, result_Hb_list.min()*1.05, result_Hb_list.max()*1.05])
  ax1.set_xlabel('Time (Myr)'); ax1.set_ylabel('Basin elevation (m)')
  ax1.yaxis.grid(which = 'major', linestyle = '--', color = 'k', alpha = 0.5, zorder = 0); ax1.xaxis.grid(which = 'major', linestyle = '--', color = 'k', alpha = 0.5, zorder = 1)
  fig.subplots_adjust(left=0.1, bottom=0.15, right=0.95, top=0.95, wspace=0.2)
  
  #fig.savefig("Z:\Post-orogenic sediment flux to continental margins\Modelling_Post_orogenic_system\Tucker_Code\Result_NorthernPyrenees_MC_Approach\Figure_7_basin_elevation_density.png",dpi=1200, transparent = True)
  #fig.savefig("Z:\Post-orogenic sediment flux to continental margins\Modelling_Post_orogenic_system\Tucker_Code\Result_NorthernPyrenees_MC_Approach\Figure_7_basin_elevation_density.pdf",dpi=1200, transparent = False)

#%%=========================================================================%%#
def density_basin_deflexion_vs_time(result_wb_list, result_time_list, time, tt):
  
  meanDeflexion_wb = result_wb_list.mean(axis=0)
  percentile_4p55_Deflexion_wb = np.percentile(result_wb_list, 4.55, axis=0)
  percentile_95p45_Deflexion_wb = np.percentile(result_wb_list, 95.45, axis=0)
  
  fig = plt.figure(figsize = (8,4))
  ax1 = fig.add_subplot(111)

  ax1.hist2d(result_time_list.flatten(), result_wb_list.flatten(), bins = (100, 100), cmap = plt.cm.jet, norm = LogNorm(), zorder = 2)
  ax1.plot(time, meanDeflexion_wb, linewidth = 2, color = 'k', linestyle = '-', zorder = 3)
  ax1.plot(time, percentile_4p55_Deflexion_wb, linewidth = 2, color = 'k', linestyle = '--', zorder = 3)
  ax1.plot(time, percentile_95p45_Deflexion_wb, linewidth = 2, color = 'k', linestyle = '--', zorder = 3)

  ax1.axis([0, tt*1e-6, result_wb_list.min()*1.05, result_wb_list.max()*1.05])
  ax1.set_xlabel('Time (Myr)'); ax1.set_ylabel('Basin deflexion (m)')
  ax1.yaxis.grid(which = 'major', linestyle = '--', color = 'k', alpha = 0.5, zorder=0); ax1.xaxis.grid(which = 'major', linestyle = '--', color = 'k', alpha = 0.5, zorder = 1)
  fig.subplots_adjust(left = 0.1, bottom = 0.15, right = 0.95, top = 0.95, wspace = 0.2)

  #fig.savefig("Z:\Post-orogenic sediment flux to continental margins\Modelling_Post_orogenic_system\Tucker_Code\Result_NorthernPyrenees_MC_Approach\Figure_5_basin_elevation_evolution.png",dpi=1200, transparent = True)
  #fig.savefig("Z:\Post-orogenic sediment flux to continental margins\Modelling_Post_orogenic_system\Tucker_Code\Result_NorthernPyrenees_MC_Approach\Figure_5_basin_elevation_evolution.pdf",dpi=1200, transparent = False)

#%%=========================================================================%%#
def range_flux_vs_time(result_Fe_list, time, tt):
  
  meanFlux_Fe = result_Fe_list.mean(axis=0)
  percentile_4p55_Flux_Fe = np.percentile(result_Fe_list, 4.55, axis=0)
  percentile_95p45_Flux_Fe = np.percentile(result_Fe_list, 95.45, axis=0)
  
  fig = plt.figure(figsize = (8,4))
  ax1 = fig.add_subplot(111)
  
  ax1.plot(time, result_Fe_list.T, color = 'steelblue', linestyle = '-.')
  ax1.plot(time, meanFlux_Fe, linewidth = 2, color = 'k', linestyle = '-')
  ax1.plot(time, percentile_4p55_Flux_Fe, linewidth = 2, color = 'k', linestyle = '--')
  ax1.plot(time, percentile_95p45_Flux_Fe, linewidth = 2, color = 'k', linestyle = '--')
  
  ax1.axis([0, tt*1e-6, result_Fe_list.min()*1.05, result_Fe_list.max()*1.05])
  ax1.set_xlabel('Time (Myr)'); ax1.set_ylabel('Range Sedimentary Fluxes (kg.m$^{-1}$.yr$^{-1}$)')
  ax1.yaxis.grid(which='major',linestyle='--',color='k',alpha=0.5,zorder=0); ax1.xaxis.grid(which='major',linestyle='--',color='k',alpha=0.5,zorder=0)
  fig.subplots_adjust(left=0.125, bottom=0.15, right=0.95, top=0.95, wspace=0.2)
  
  #fig.savefig("Z:\Post-orogenic sediment flux to continental margins\Modelling_Post_orogenic_system\Tucker_Code\Result_NorthernPyrenees_MC_Approach\Figure_8_range_flux_evolution.png",dpi=1200, transparent = True)
  #fig.savefig("Z:\Post-orogenic sediment flux to continental margins\Modelling_Post_orogenic_system\Tucker_Code\Result_NorthernPyrenees_MC_Approach\Figure_8_range_flux_evolution.pdf",dpi=1200, transparent = False)

#%%=========================================================================%%#
def basin_flux_vs_time(result_Fb_list, time, tt):
  
  meanFlux_Fb = result_Fb_list.mean(axis=0)
  percentile_4p55_Flux_Fb = np.percentile(result_Fb_list, 4.55, axis=0)
  percentile_95p45_Flux_Fb = np.percentile(result_Fb_list, 95.45, axis=0)  
  
  fig = plt.figure(figsize = (8,4))
  ax1 = fig.add_subplot(111)
  
  ax1.plot(time, result_Fb_list.T, color = 'steelblue', linestyle = '-.')
  ax1.plot(time, meanFlux_Fb, linewidth = 2, color = 'k', linestyle = '-')
  ax1.plot(time, percentile_4p55_Flux_Fb, linewidth = 2, color = 'k', linestyle = '--')
  ax1.plot(time, percentile_95p45_Flux_Fb, linewidth = 2, color = 'k', linestyle = '--')
  
  ax1.axis([0, tt*1e-6, result_Fb_list.min()*1.05, result_Fb_list.max()*1.05])
  ax1.set_xlabel('Time (Myr)'); ax1.set_ylabel('Basin Sedimentary Fluxes (kg.m$^{-1}$.yr$^{-1}$)')
  ax1.yaxis.grid(which='major',linestyle='--',color='k',alpha=0.5,zorder=0); ax1.xaxis.grid(which='major',linestyle='--',color='k',alpha=0.5,zorder=0)
  fig.subplots_adjust(left=0.125, bottom=0.15, right=0.95, top=0.95, wspace=0.2)
  
  #fig.savefig("Z:\Post-orogenic sediment flux to continental margins\Modelling_Post_orogenic_system\Tucker_Code\Result_NorthernPyrenees_MC_Approach\Figure_9_basin_flux_evolution.png",dpi=1200, transparent = True)
  #fig.savefig("Z:\Post-orogenic sediment flux to continental margins\Modelling_Post_orogenic_system\Tucker_Code\Result_NorthernPyrenees_MC_Approach\Figure_9_basin_flux_evolution.pdf",dpi=1200, transparent = False)

#%%=========================================================================%%#
def density_range_flux_elevation_vs_time(result_Fe_list, result_time_list, time, tt):
  
  meanFlux_Fe = result_Fe_list.mean(axis=0)
  percentile_4p55_Flux_Fe = np.percentile(result_Fe_list, 4.55, axis=0)
  percentile_95p45_Flux_Fe = np.percentile(result_Fe_list, 95.45, axis=0)
  
  fig = plt.figure(figsize = (8,4))
  ax1 = fig.add_subplot(111)
  
  ax1.hist2d(result_time_list.flatten(), result_Fe_list.flatten(), bins = (100, 100), cmap = plt.cm.jet, norm = LogNorm(), zorder = 2)
  ax1.plot(time, meanFlux_Fe, linewidth = 2, color = 'k', linestyle = '-', zorder = 3)
  ax1.plot(time, percentile_4p55_Flux_Fe, linewidth = 2, color = 'k', linestyle = '--', zorder = 3)
  ax1.plot(time, percentile_95p45_Flux_Fe, linewidth = 2, color = 'k', linestyle = '--', zorder = 3)
  
  ax1.axis([0, tt*1e-6, result_Fe_list.min()*1.05, result_Fe_list.max()*1.05])
  ax1.set_xlabel('Time (Myr)'); ax1.set_ylabel('Range Sedimentary Fluxes (kg.m$^{-1}$.yr$^{-1}$)')
  ax1.yaxis.grid(which = 'major', linestyle = '--', color = 'k', alpha = 0.5, zorder = 0); ax1.xaxis.grid(which = 'major', linestyle = '--', color = 'k', alpha = 0.5, zorder = 1)
  fig.subplots_adjust(left=0.125, bottom=0.15, right=0.95, top=0.95, wspace=0.2)
  
  #fig.savefig("Z:\Post-orogenic sediment flux to continental margins\Modelling_Post_orogenic_system\Tucker_Code\Result_NorthernPyrenees_MC_Approach\Figure_10_range_flux_density.png",dpi=1200, transparent = True)
  #fig.savefig("Z:\Post-orogenic sediment flux to continental margins\Modelling_Post_orogenic_system\Tucker_Code\Result_NorthernPyrenees_MC_Approach\Figure_10_range_flux_density.pdf",dpi=1200, transparent = False)

#%%=========================================================================%%#
def density_basin_flux_elevation_vs_time(result_Fb_list, result_time_list, time, tt):
  
  meanFlux_Fb = result_Fb_list.mean(axis=0)
  percentile_4p55_Flux_Fb = np.percentile(result_Fb_list, 4.55, axis=0)
  percentile_95p45_Flux_Fb = np.percentile(result_Fb_list, 95.45, axis=0)
  
  fig = plt.figure(10, figsize = (8,4))
  ax1 = fig.add_subplot(111)
  
  ax1.hist2d(result_time_list.flatten(), result_Fb_list.flatten(), bins = (100, 100), cmap = plt.cm.jet, norm = LogNorm(), zorder = 2)
  ax1.plot(time, meanFlux_Fb, linewidth = 2, color = 'k', linestyle = '-', zorder = 3)
  ax1.plot(time, percentile_4p55_Flux_Fb, linewidth = 2, color = 'k', linestyle = '--', zorder = 3)
  ax1.plot(time, percentile_95p45_Flux_Fb, linewidth = 2, color = 'k', linestyle = '--', zorder = 3)
  
  ax1.axis([0, tt*1e-6, result_Fb_list.min()*1.05, result_Fb_list.max()*1.05])
  ax1.set_xlabel('Time (Myr)'); ax1.set_ylabel('Basin Sedimentary Fluxes (kg.m$^{-1}$.yr$^{-1}$)')
  ax1.yaxis.grid(which = 'major', linestyle = '--', color = 'k', alpha = 0.5, zorder = 0); ax1.xaxis.grid(which = 'major', linestyle = '--', color = 'k', alpha = 0.5, zorder = 1)
  fig.subplots_adjust(left=0.125, bottom=0.15, right=0.95, top=0.95, wspace=0.2)
  
  #fig.savefig("Z:\Post-orogenic sediment flux to continental margins\Modelling_Post_orogenic_system\Tucker_Code\Result_NorthernPyrenees_MC_Approach\Figure_11_basin_flux_density.png",dpi=1200, transparent = True)
  #fig.savefig("Z:\Post-orogenic sediment flux to continental margins\Modelling_Post_orogenic_system\Tucker_Code\Result_NorthernPyrenees_MC_Approach\Figure_11_basin_flux_density.pdf",dpi=1200, transparent = False)

#%%=========================================================================%%#
def specific_elevation_range_basin(MC_iteration, result_Hr_list, result_Hb_list, T3_time_steps):
  
  result_max_Hr_list = np.empty(MC_iteration)
  result_max_Hb_list = np.empty(MC_iteration)
  result_max_Hb_rebound_list = np.empty(MC_iteration)
  
  for i in range(0, MC_iteration):
    result_max_Hr_list[i] = np.max(result_Hr_list[i])
    result_max_Hb_list[i] = np.max(result_Hb_list[i])
    result_max_Hb_rebound_list[i] = np.max(result_Hb_list[i]) - result_Hb_list[i,T3_time_steps]
  
  max_Hr_mn = np.mean(result_max_Hr_list); max_Hr_max = np.max(result_max_Hr_list); max_Hr_min = np.min(result_max_Hr_list); max_Hr_2SDinf = np.percentile(result_max_Hr_list, 4.55); max_Hr_2SDsup = np.percentile(result_max_Hr_list, 95.45)
  max_Hb_mn = np.mean(result_max_Hb_list); max_Hb_max = np.max(result_max_Hb_list); max_Hb_min = np.min(result_max_Hb_list); max_Hb_2SDinf = np.percentile(result_max_Hb_list, 4.55); max_Hb_2SDsup = np.percentile(result_max_Hb_list, 95.45)
  max_Hb_rebound_mn = np.mean(result_max_Hb_rebound_list); max_Hb_rebound_max = np.max(result_max_Hb_rebound_list); max_Hb_rebound_min = np.min(result_max_Hb_rebound_list); max_Hb_rebound_2SDinf = np.percentile(result_max_Hb_rebound_list, 4.55); max_Hb_rebound_2SDsup = np.percentile(result_max_Hb_rebound_list, 95.45)
  
  fig = plt.figure(figsize = (5,3.5))
  ax1 = fig.add_subplot(131)
  ax2 = fig.add_subplot(132)
  ax3 = fig.add_subplot(133)
  
  ax1.errorbar(x = 1, y = max_Hr_mn, yerr = [[max_Hr_mn - max_Hr_2SDinf],[max_Hr_2SDsup - max_Hr_mn]], fmt = 'ok', markersize = 6.5, markeredgecolor = 'k', markerfacecolor = 'darkred', elinewidth = 2, capsize= 4, capthick = 2, zorder = 3)
  ax1.errorbar(x = 1, y = max_Hr_mn, yerr = [[max_Hr_mn - max_Hr_min],[max_Hr_max - max_Hr_mn]], fmt = 'ok', ecolor='darkred', markersize = 7, markeredgecolor = 'k', markerfacecolor = 'darkred', elinewidth = 1.5, capsize= 5, capthick = 1.5, zorder=2)
  
  ax2.errorbar(x = 1, y = max_Hb_mn, yerr = [[max_Hb_mn - max_Hb_2SDinf],[max_Hb_2SDsup - max_Hb_mn]], fmt = 'ok', markersize = 6.5, markeredgecolor = 'k', markerfacecolor = 'darkred', elinewidth = 2, capsize= 4, capthick = 2, zorder = 3)
  ax2.errorbar(x = 1, y = max_Hb_mn, yerr = [[max_Hb_mn - max_Hb_min],[max_Hb_max - max_Hb_mn]], fmt = 'ok', ecolor='darkred', markersize = 7, markeredgecolor = 'k', markerfacecolor = 'darkred', elinewidth = 1.5, capsize= 5, capthick = 1.5, zorder=2)
  
  ax3.errorbar(x = 1, y = max_Hb_rebound_mn, yerr = [[max_Hb_rebound_mn - max_Hb_rebound_2SDinf],[max_Hb_rebound_2SDsup - max_Hb_rebound_mn]], fmt = 'ok', markersize = 6.5, markeredgecolor = 'k', markerfacecolor = 'darkred', elinewidth = 2, capsize= 4, capthick = 2, zorder = 3)
  ax3.errorbar(x = 1, y = max_Hb_rebound_mn, yerr = [[max_Hb_rebound_mn - max_Hb_rebound_min],[max_Hb_rebound_max - max_Hb_rebound_mn]], fmt = 'ok', ecolor='darkred', markersize = 7, markeredgecolor = 'k', markerfacecolor = 'darkred', elinewidth = 1.5, capsize= 5, capthick = 1.5, zorder=2)
  
  ax1.yaxis.grid(which = 'major', linestyle = '--', color = 'k', alpha = 0.5, zorder = 0); ax1.text(0.15, 0.95, "A", fontsize = 12, horizontalalignment='center', verticalalignment='center', transform=ax1.transAxes); ax1.set_ylabel('Elevation (m)'); ax1.axis([0.95,1.05,max_Hr_min*0.9,max_Hr_max*1.1])
  ax2.yaxis.grid(which = 'major', linestyle = '--', color = 'k', alpha = 0.5, zorder = 0); ax2.text(0.15, 0.95, "B", fontsize = 12, horizontalalignment='center', verticalalignment='center', transform=ax2.transAxes); ax2.axis([0.95,1.05,max_Hb_min*0.9,max_Hb_max*1.1])
  ax3.yaxis.grid(which = 'major', linestyle = '--', color = 'k', alpha = 0.5, zorder = 0); ax3.text(0.15, 0.95, "C", fontsize = 12, horizontalalignment='center', verticalalignment='center', transform=ax3.transAxes); ax3.axis([0.95,1.05,max_Hb_rebound_min*0.9,max_Hb_rebound_max*1.1])
  fig.subplots_adjust(left=0.15, bottom=0.1, right=0.95, top=0.95, wspace=0.5)
  
  #fig.savefig("Z:\Post-orogenic sediment flux to continental margins\Modelling_Post_orogenic_system\Tucker_Code\Result_NorthernPyrenees_MC_Approach\Figure_12_spatial_estimation.png",dpi=1200, transparent = True)
  #fig.savefig("Z:\Post-orogenic sediment flux to continental margins\Modelling_Post_orogenic_system\Tucker_Code\Result_NorthernPyrenees_MC_Approach\Figure_12_spatial_estimation.pdf",dpi=1200, transparent = False)

#%%=========================================================================%%#
def specific_timing_range_basin(MC_iteration, T3_time_steps, T6_time_steps, result_Hb_list, result_Fe_list, result_Fb_list):
  
  time_aggradation = np.zeros(MC_iteration)
  for i in range(0, MC_iteration):
    for j in range(T3_time_steps,T6_time_steps):
      if result_Fb_list[i,j] > result_Fe_list[i,j]:
        time_aggradation[i] = j - T3_time_steps
        break
  time_aggradation = time_aggradation/100
  
  ## calcul the time of elevation rebound in the basin 
  time_rebound = np.zeros(MC_iteration)
  for i in range(0, MC_iteration):
     for j in range(T3_time_steps,T6_time_steps):
       if result_Hb_list[i,j] < result_Hb_list[i,j-1]:
         time_rebound[i] = j - T3_time_steps
         break
  time_rebound = time_rebound/100
  
  ## calcul the time of more important sedimentary fuxes to the base level after decrease of convergence
  time_sedfluxbl = np.zeros(MC_iteration)
  for i in range(0, MC_iteration):
    for j in range(T3_time_steps+1, T6_time_steps):
      if j == T6_time_steps-1:
        time_sedfluxbl[i] = T6_time_steps - T3_time_steps
        break
      if result_Fb_list[i,j] < result_Fb_list[i,T3_time_steps]:
        time_sedfluxbl[i] = j - T3_time_steps
        break
  time_sedfluxbl = time_sedfluxbl/100
  
  time_aggradation_mn = np.mean(time_aggradation); time_aggradation_max = np.max(time_aggradation); time_aggradation_min = np.min(time_aggradation); time_aggradation_2SDinf = np.percentile(time_aggradation, 4.55); time_aggradation_2SDsup = np.percentile(time_aggradation, 95.45)
  time_rebound_mn = np.mean(time_rebound); time_rebound_max = np.max(time_rebound); time_rebound_min = np.min(time_rebound); time_rebound_2SDinf = np.percentile(time_rebound, 4.55); time_rebound_2SDsup = np.percentile(time_rebound, 95.45)
  time_sedfluxbl_mn = np.mean(time_sedfluxbl); time_sedfluxbl_max = np.max(time_sedfluxbl); time_sedfluxbl_min = np.min(time_sedfluxbl); time_sedfluxbl_2SDinf = np.percentile(time_sedfluxbl, 4.55); time_sedfluxbl_2SDsup = np.percentile(time_sedfluxbl, 95.45)
  
  fig = plt.figure(figsize = (5,4))
  ax1 = fig.add_subplot(311)
  ax2 = fig.add_subplot(312)
  ax3 = fig.add_subplot(313)
  
  ax1.errorbar(x = time_aggradation_mn, y = 1, xerr = [[time_aggradation_mn - time_aggradation_2SDinf],[time_aggradation_2SDsup - time_aggradation_mn]], fmt = 'ok', markersize = 6.5, markeredgecolor = 'k', markerfacecolor = 'darkred', elinewidth = 2, capsize= 4, capthick = 2, zorder = 3)
  ax1.errorbar(x = time_aggradation_mn, y = 1, xerr = [[time_aggradation_mn - time_aggradation_min],[time_aggradation_max - time_aggradation_mn]], fmt = 'ok', ecolor='darkred', markersize = 7, markeredgecolor = 'k', markerfacecolor = 'darkred', elinewidth = 1.5, capsize= 5, capthick = 1.5, zorder=2)
  
  ax2.errorbar(x = time_rebound_mn, y = 1, xerr = [[time_rebound_mn - time_rebound_2SDinf],[time_rebound_2SDsup - time_rebound_mn]], fmt = 'ok', markersize = 6.5, markeredgecolor = 'k', markerfacecolor = 'darkred', elinewidth = 2, capsize= 4, capthick = 2, zorder = 3)
  ax2.errorbar(x = time_rebound_mn, y = 1, xerr = [[time_rebound_mn - time_rebound_min],[time_rebound_max - time_rebound_mn]], fmt = 'ok', ecolor='darkred', markersize = 7, markeredgecolor = 'k', markerfacecolor = 'darkred', elinewidth = 1.5, capsize= 5, capthick = 1.5, zorder=2)
  
  ax3.errorbar(x = time_sedfluxbl_mn, y = 1, xerr = [[time_sedfluxbl_mn - time_sedfluxbl_2SDinf],[time_sedfluxbl_2SDsup - time_sedfluxbl_mn]], fmt = 'ok', markersize = 6.5, markeredgecolor = 'k', markerfacecolor = 'darkred', elinewidth = 2, capsize= 4, capthick = 2, zorder = 3)
  ax3.errorbar(x = time_sedfluxbl_mn, y = 1, xerr = [[time_sedfluxbl_mn - time_sedfluxbl_min],[time_sedfluxbl_max - time_sedfluxbl_mn]], fmt = 'ok', ecolor='darkred', markersize = 7, markeredgecolor = 'k', markerfacecolor = 'darkred', elinewidth = 1.5, capsize= 5, capthick = 1.5, zorder=2)
  
  ax1.xaxis.grid(which = 'major', linestyle = '--', color = 'k', alpha = 0.5, zorder = 0); ax1.text(0.04, 0.85, "D", fontsize = 12, horizontalalignment='center', verticalalignment='center', transform=ax1.transAxes); #ax1.axis([1.0,2.6,0.95,1.05])
  ax2.xaxis.grid(which = 'major', linestyle = '--', color = 'k', alpha = 0.5, zorder = 0); ax2.text(0.04, 0.85, "E", fontsize = 12, horizontalalignment='center', verticalalignment='center', transform=ax2.transAxes); #ax2.axis([4.0,4.9,0.95,1.05])
  ax3.xaxis.grid(which = 'major', linestyle = '--', color = 'k', alpha = 0.5, zorder = 0);  ax3.text(0.04, 0.85, "F", fontsize = 12, horizontalalignment='center', verticalalignment='center', transform=ax3.transAxes); ax3.set_xlabel('Time (Myr)') #ax3.axis([8,24,0.95,1.05]);
  
  fig.subplots_adjust(left=0.10, bottom=0.16, right=0.95, top=0.95, hspace=0.35)
  
  #fig.savefig("Z:\Post-orogenic sediment flux to continental margins\Modelling_Post_orogenic_system\Tucker_Code\Result_NorthernPyrenees_MC_Approach\Figure_13_timing_estimation.png",dpi=1200, transparent = True)
  #fig.savefig("Z:\Post-orogenic sediment flux to continental margins\Modelling_Post_orogenic_system\Tucker_Code\Result_NorthernPyrenees_MC_Approach\Figure_13_timing_estimation.pdf",dpi=1200, transparent = False)

#%%=========================================================================%%#
def density_relief_vs_time(result_Hr_list, result_Hb_list, result_time_list, time, tt):

  result_relief_rb_list = result_Hr_list - result_Hb_list
  meanrelief_rb = result_relief_rb_list.mean(axis=0)
  percentile_4p55_relief_rb = np.percentile(result_relief_rb_list, 4.55, axis=0)
  percentile_95p45_relief_rb = np.percentile(result_relief_rb_list, 95.45, axis=0)
  
  fig = plt.figure(figsize = (8,4))
  ax1 = fig.add_subplot(111)
  
  ax1.hist2d(result_time_list.flatten(), result_relief_rb_list.flatten(), bins = (100, 100), cmap = plt.cm.jet, norm = LogNorm(), zorder = 2)
  ax1.plot(time, meanrelief_rb, linewidth = 2, color = 'k', linestyle = '-', zorder = 3)
  ax1.plot(time, percentile_4p55_relief_rb, linewidth = 2, color = 'k', linestyle = '--', zorder = 3)
  ax1.plot(time, percentile_95p45_relief_rb, linewidth = 2, color = 'k', linestyle = '--', zorder = 3)
  
  ax1.axis([0,tt*1e-6,0,result_relief_rb_list.max()*1.05])
  ax1.set_xlabel('Time (Myr)'); ax1.set_ylabel('Relief between range and basin (m)')
  ax1.yaxis.grid(which = 'major', linestyle = '--', color = 'k', alpha = 0.5, zorder = 0); ax1.xaxis.grid(which = 'major', linestyle = '--', color = 'k', alpha = 0.5, zorder = 1)
  fig.subplots_adjust(left=0.1, bottom=0.15, right=0.95, top=0.95, wspace=0.2)
  
  #fig.savefig("Z:\Post-orogenic sediment flux to continental margins\Modelling_Post_orogenic_system\Tucker_Code\Result_NorthernPyrenees_MC_Approach\Figure_14_relief_density.png",dpi=1200, transparent = True)
  #fig.savefig("Z:\Post-orogenic sediment flux to continental margins\Modelling_Post_orogenic_system\Tucker_Code\Result_NorthernPyrenees_MC_Approach\Figure_14_relief_density.pdf",dpi=1200, transparent = False)

#%%=========================================================================%%#
def volume_change_range_vs_time(result_Fr_list, result_Fe_list, rhoc, result_time_list, time, tt):

  volumechangerange = (result_Fr_list - result_Fe_list)/rhoc
  mean_volumechangerange = volumechangerange.mean(axis=0)
  percentile_4p55_volumechangerange = np.percentile(volumechangerange, 4.55, axis=0)
  percentile_95p45_volumechangerange = np.percentile(volumechangerange, 95.45, axis=0)
  
  fig = plt.figure(figsize =(8,4))
  ax1 = fig.add_subplot(111)
  
  ax1.hist2d(result_time_list.flatten(), volumechangerange.flatten(), bins = (100, 100), cmap = plt.cm.jet, norm = LogNorm(), zorder = 2)
  ax1.plot(time, mean_volumechangerange, linewidth = 2, color = 'k', linestyle = '-', zorder = 3)
  ax1.plot(time, percentile_4p55_volumechangerange, linewidth = 2, color = 'k', linestyle = '--', zorder = 3)
  ax1.plot(time, percentile_95p45_volumechangerange, linewidth = 2, color = 'k', linestyle = '--', zorder = 3)
  
  ax1.axis([0,tt*1e-6,volumechangerange.min()*1.05,volumechangerange.max()*1.05])
  ax1.set_xlabel('Time (Myr)'); ax1.set_ylabel('Range volume change (m$^{2}$.yr$^{-1}$)')
  ax1.yaxis.grid(which='major',linestyle='--',color='k',alpha=0.5,zorder=0); ax1.xaxis.grid(which='major',linestyle='--',color='k',alpha=0.5,zorder=0)
  fig.subplots_adjust(left=0.15, bottom=0.15, right=0.95, top=0.95, wspace=0.2)
  
  #fig.savefig("Z:\Post-orogenic sediment flux to continental margins\Modelling_Post_orogenic_system\Tucker_Code\Result_NorthernPyrenees_MC_Approach\Figure_15_range_mass_loss_density.png",dpi=1200, transparent = True)
  #fig.savefig("Z:\Post-orogenic sediment flux to continental margins\Modelling_Post_orogenic_system\Tucker_Code\Result_NorthernPyrenees_MC_Approach\Figure_15_range_mass_loss_density.pdf",dpi=1200, transparent = False)

#%%=========================================================================%%#
def volume_change_basin_vs_time(result_Fe_list, result_Fb_list, rhos, result_time_list, time, tt):

  volumechangebasin = (result_Fe_list - result_Fb_list)/rhos
  mean_volumechangebasin = volumechangebasin.mean(axis=0)
  percentile_4p55_volumechangebasin = np.percentile(volumechangebasin, 4.55, axis=0)
  percentile_95p45_volumechangebasin = np.percentile(volumechangebasin, 95.45, axis=0)
  
  fig = plt.figure(figsize =(8,4))
  ax1 = fig.add_subplot(111)
  
  ax1.hist2d(result_time_list.flatten(), volumechangebasin.flatten(), bins = (100, 100), cmap = plt.cm.jet, norm = LogNorm(), zorder = 2)
  ax1.plot(time, mean_volumechangebasin, linewidth = 2, color = 'k', linestyle = '-', zorder = 3)
  ax1.plot(time, percentile_4p55_volumechangebasin, linewidth = 2, color = 'k', linestyle = '--', zorder = 3)
  ax1.plot(time, percentile_95p45_volumechangebasin, linewidth = 2, color = 'k', linestyle = '--', zorder = 3)
  
  ax1.axis([0,tt*1e-6,volumechangebasin.min()*1.05,volumechangebasin.max()*1.05])
  ax1.set_xlabel('Time (Myr)'); ax1.set_ylabel('Basin volume change (m$^{2}$.yr$^{-1}$)')
  ax1.yaxis.grid(which='major',linestyle='--',color='k',alpha=0.5,zorder=0); ax1.xaxis.grid(which='major',linestyle='--',color='k',alpha=0.5,zorder=0)
  fig.subplots_adjust(left=0.15, bottom=0.15, right=0.95, top=0.95, wspace=0.2)
  
  #fig.savefig("Z:\Post-orogenic sediment flux to continental margins\Modelling_Post_orogenic_system\Tucker_Code\Result_NorthernPyrenees_MC_Approach\Figure_16_basin_mass_loss_density.png",dpi=1200, transparent = True)
  #fig.savefig("Z:\Post-orogenic sediment flux to continental margins\Modelling_Post_orogenic_system\Tucker_Code\Result_NorthernPyrenees_MC_Approach\Figure_16_basin_mass_loss_density.pdf",dpi=1200, transparent = False)

#%%=========================================================================%%#
def cross_section_thickness_range_vs_time(result_Hr_list, result_wr_list, time, tt):
  
  meanElevation_Hr = result_Hr_list.mean(axis=0)
  percentile_4p55_Elevation_Hr = np.percentile(result_Hr_list, 4.55, axis=0)
  percentile_95p45_Elevation_Hr = np.percentile(result_Hr_list, 95.45, axis=0)
  
  mean_Deflexion_wr = result_wr_list.mean(axis=0)
  percentile_4p55_Deflexion_wr = np.percentile(result_wr_list, 4.55, axis=0)
  percentile_95p45_Deflexion_wr = np.percentile(result_wr_list, 95.45, axis=0)
  
  fig = plt.figure(figsize =(7.125,4))
  ax1 = fig.add_subplot(111)
  
  ax1.plot(time, meanElevation_Hr , color='k')
  ax1.plot(time, -mean_Deflexion_wr , color='k')
  ax1.plot(time, -percentile_95p45_Deflexion_wr , color='k', linestyle = '--', zorder=2)
  ax1.plot(time, -percentile_4p55_Deflexion_wr , color='k', linestyle = '--', zorder=2)
  ax1.plot(time, percentile_4p55_Elevation_Hr , color='k', linestyle = '--', zorder=2)
  ax1.plot(time, percentile_95p45_Elevation_Hr , color='k', linestyle = '--', zorder=2)
  
  #ax1.fill_between(time, -percentile_4p55_Deflexion_wr, percentile_4p55_Elevation_Hr, facecolor='powderblue',alpha=0.66, zorder=1)
  ax1.fill_between(time, -percentile_4p55_Deflexion_wr, -percentile_95p45_Deflexion_wr, facecolor='steelblue',alpha=0.33, zorder=1)
  ax1.fill_between(time, percentile_4p55_Elevation_Hr, percentile_95p45_Elevation_Hr, facecolor='steelblue',alpha=0.33, zorder=1)
  
  ax1.yaxis.grid(which='major', linestyle='--', color='k', alpha=0.5, zorder=0); ax1.xaxis.grid(which='major', linestyle='--', color='k', alpha=0.5, zorder=0)
  ax1.axis([0, tt*1e-6, -percentile_95p45_Deflexion_wr.max()*1.05, percentile_95p45_Elevation_Hr.max()*1.5])
  ax1.set_xlabel('Time (Ma)'); ax1.set_ylabel('Range cross-section thickness (m)')
  xlabel=['66','56','46','36','26','26','16','6']; ax1.set_xticklabels(xlabel);
  fig.tight_layout()
  
  #fig.savefig('Z:\\Post-orogenic sediment flux to continental margins\\Publication_Report\\Paper_Post-orogen_Pyrenees+Aquitaine\\long_version\\SUP_FIGURE4_range_cross-section_through_time.pdf', dpi=1200)
  #fig.savefig('Z:\\Post-orogenic sediment flux to continental margins\\Publication_Report\\Paper_Post-orogen_Pyrenees+Aquitaine\\long_version\\SUP_FIGURE4_range_cross-section_through_time.jpeg', dpi=1200)

#%%=========================================================================%%# 
def cross_section_thickness_basin_vs_time(result_Hb_list, result_wb_list, time, tt):
  
  meanElevation_Hb = result_Hb_list.mean(axis=0)
  percentile_4p55_Elevation_Hb = np.percentile(result_Hb_list, 4.55, axis=0)
  percentile_95p45_Elevation_Hb = np.percentile(result_Hb_list, 95.45, axis=0)
  
  mean_Deflexion_wb = result_wb_list.mean(axis=0)
  percentile_4p55_Deflexion_wb = np.percentile(result_wb_list, 4.55, axis=0)
  percentile_95p45_Deflexion_wb = np.percentile(result_wb_list, 95.45, axis=0)
    
  fig = plt.figure(figsize =(7.125,4))
  ax1 = fig.add_subplot(111)
  
  ax1.plot(time, meanElevation_Hb , color='k')
  ax1.plot(time, -mean_Deflexion_wb , color='k')
  ax1.plot(time, -percentile_95p45_Deflexion_wb , color='k', linestyle = '--', zorder=2)
  ax1.plot(time, -percentile_4p55_Deflexion_wb , color='k', linestyle = '--', zorder=2)
  ax1.plot(time, percentile_4p55_Elevation_Hb , color='k', linestyle = '--', zorder=2)
  ax1.plot(time, percentile_95p45_Elevation_Hb , color='k', linestyle = '--', zorder=2)
  
  #ax1.fill_between(time, -percentile_4p55_Deflexion_wb, percentile_4p55_Elevation_Hb, facecolor='powderblue',alpha=0.66, zorder=1)
  ax1.fill_between(time, -percentile_4p55_Deflexion_wb, -percentile_95p45_Deflexion_wb, facecolor='steelblue',alpha=0.33, zorder=1)
  ax1.fill_between(time, percentile_4p55_Elevation_Hb, percentile_95p45_Elevation_Hb, facecolor='steelblue',alpha=0.33, zorder=1)
  
  ax1.yaxis.grid(which='major', linestyle='--', color='k', alpha=0.5, zorder=0); ax1.xaxis.grid(which='major', linestyle='--', color='k', alpha=0.5, zorder=0)
  ax1.axis([0, tt*1e-6, -percentile_95p45_Deflexion_wb.max()*1.05, percentile_95p45_Elevation_Hb.max()*1.5])
  ax1.set_xlabel('Time (Ma)'); ax1.set_ylabel('Basin cross-section thickness (m)')
  xlabel=['66','56','46','36','26','26','16','6']; ax1.set_xticklabels(xlabel);
  fig.tight_layout()
  
  #fig.savefig('Z:\\Post-orogenic sediment flux to continental margins\\Publication_Report\\Paper_Post-orogen_Pyrenees+Aquitaine\\long_version\\SUP_FIGURE5_basin_cross-section_through_time.pdf', dpi=1200)
  #fig.savefig('Z:\\Post-orogenic sediment flux to continental margins\\Publication_Report\\Paper_Post-orogen_Pyrenees+Aquitaine\\long_version\\SUP_FIGURE5_basin_cross-section_through_time.jpeg', dpi=1200)

#%%=========================================================================%%#  
def contribution_elevation_change_basin(result_Fe_list, result_Fb_list, result_wb_list, rhos, Lb, time, T3_time_steps):

  meanFlux_Fe = result_Fe_list.mean(axis=0)
  meanFlux_Fb = result_Fb_list.mean(axis=0)
  mean_Deflexion_wb = result_wb_list.mean(axis=0)
  
  t_post = T3_time_steps
  t_end = int(round(time[-1]*100, 0))
  
  post = t_end - t_post
  post_mean_Deflexion_wb_change = np.zeros(post)
  post_mean_positive_elevation = meanFlux_Fe[t_post:t_end]/rhos/Lb
  post_mean_negative_elevation = meanFlux_Fb[t_post:t_end]/rhos/Lb
  
  for i in range(0,post):
    post_mean_Deflexion_wb_change[i] = mean_Deflexion_wb[i+(t_post-1)]-mean_Deflexion_wb[i+t_post]
  
  for i in range(0,post):
    if post_mean_Deflexion_wb_change[i] < 0:
      post_mean_negative_elevation[i] = post_mean_negative_elevation[i] - post_mean_Deflexion_wb_change[i]
    else:
      post_mean_positive_elevation[i] = post_mean_positive_elevation[i] + post_mean_Deflexion_wb_change[i]
      
  fig = plt.figure(figsize =(8,4))
  ax1 = fig.add_subplot(111)
  
  ax1.plot(time[t_post:t_end], post_mean_Deflexion_wb_change, color='darkseagreen', linestyle='-.')
  ax1.plot(time[t_post:t_end], -meanFlux_Fb[t_post:t_end]/rhos/Lb, color='steelblue', linestyle='-.')
  ax1.plot(time[t_post:t_end], meanFlux_Fe[t_post:t_end]/rhos/Lb, color='darkred', linestyle='-.')
  ax1.plot(time[t_post:t_end], -post_mean_negative_elevation, color='black', linestyle='--')
  ax1.plot(time[t_post:t_end], +post_mean_positive_elevation, color='black', linestyle='--')
  ax1.plot(time[t_post:t_end], post_mean_positive_elevation - post_mean_negative_elevation, color='black',)
  ax1.fill_between(time[t_post:t_end], (post_mean_positive_elevation - post_mean_negative_elevation), 0, facecolor='powderblue',alpha=0.66)
  
  ax1.axis([time[t_post], time[t_end], -15, 15])
  ax1.yaxis.grid(which='major', linestyle='--', color='k', alpha=0.5, zorder=0); ax1.xaxis.grid(which='major', linestyle='--', color='k', alpha=0.5, zorder=0)
  ax1.set_xlabel('Time (Myr)'); ax1.set_ylabel('Elevation contribution (m)')
  fig.subplots_adjust(left=0.125, bottom=0.15, right=0.975, top=0.95)
  
  #fig.savefig("Z:\Post-orogenic sediment flux to continental margins\Modelling_Post_orogenic_system\Tucker_Code\Result_NorthernPyrenees_MC_Approach\Figure_19_basin_contribution.png",dpi=1200, transparent = True)
  #fig.savefig("Z:\Post-orogenic sediment flux to continental margins\Modelling_Post_orogenic_system\Tucker_Code\Result_NorthernPyrenees_MC_Approach\Figure_19_basin_contribution.pdf",dpi=1200, transparent = False)

#%%=========================================================================%%#
def model_cross_section_vs_time(nt, result_Hr_list, result_Hb_list, result_etar_list, result_etab_list, result_Te_list, result_T4_list, rhos, rhoc, rhom, rhoi, E, nu, g, Lr, T1_time_steps, T2_time_steps, T3_time_steps, T6_time_steps, V1, V2, V3):
 
  meanElevation_Hr = result_Hr_list.mean(axis=0)
  percentile_4p55_Elevation_Hr = np.percentile(result_Hr_list, 4.55, axis=0)
  percentile_95p45_Elevation_Hr = np.percentile(result_Hr_list, 95.45, axis=0)
  
  meanElevation_Hb = result_Hb_list.mean(axis=0)
  percentile_4p55_Elevation_Hb = np.percentile(result_Hb_list, 4.55, axis=0)
  percentile_95p45_Elevation_Hb = np.percentile(result_Hb_list, 95.45, axis=0)
    
  mean_thickness_etar = result_etar_list.mean(axis=0)
  percentile_4p55_thickness_etar = np.percentile(result_etar_list, 4.55, axis=0)
  percentile_95p45_thickness_etar = np.percentile(result_etar_list, 95.45, axis=0)
  
  mean_thickness_etab = result_etab_list.mean(axis=0)
  #percentile_4p55_thickness_etab = np.percentile(result_etab_list, 4.55, axis=0)
  #percentile_95p45_thickness_etab = np.percentile(result_etab_list, 95.45, axis=0)
  
  Te = np.mean(result_Te_list)
  D = E*Te**3/(12*(1-nu**2))
  Lambda = ((rhom-rhoi)*g/(4*D))**(1./4.)
  #tlaml = 2*Lambda*Lr
  Lb = 3*np.pi/(4*Lambda)-Lr
  #Llm = Lambda*Lr
  
  x_range = np.linspace(0.,Lr, 1000)
  x_basin = np.linspace(Lr,(Lr+Lb) , 1000)
  a_range = x_range + Lr                      
  b_range = Lr-x_range                        
  b_basin = x_basin-Lr                        
  a_basin = 2*Lr+b_basin
  
  T4_time_steps = int(round(np.mean(result_T4_list)))
  Vt = np.zeros(nt)
  Vt[T1_time_steps:T2_time_steps] = V1
  Vt[T2_time_steps:T3_time_steps] = V2
  Vt[T3_time_steps:T6_time_steps] = V3
  
  V_T3 = V2
  for v in range(T3_time_steps,T4_time_steps):
    Vt[v] = V_T3 - (V2-V3)/(T4_time_steps-T3_time_steps)
    V_T3 = Vt[v]
  
  time3 = 0
  path_to_save1 = 'Z:\\Post-orogenic sediment flux to continental margins\\Modelling_Post_orogenic_system\\Tucker_Code\\Result_NorthernPyrenees_InverseModelling+AlphaEffect\\range+basin_cross-section_time\\'
  figure_name = 'cross-section_t'
  #general_view = [-10, 170, -35000, 5000]
  topo_view = [-5, 170, -500, 4500]
  n = 10
  
  for i in range(0,int(nt/n)):
  
    time1 = i*n
    time2 = i/n
    
    mean_wxb = -(rhoc*mean_thickness_etar[time1]/(2*(rhom-rhoi)))*(np.exp(-Lambda*a_basin)*np.cos(Lambda*a_basin)-np.exp(-Lambda*b_basin)*np.cos(Lambda*b_basin))
    percentile_4p55_wxb = -(rhoc*percentile_4p55_thickness_etar[time1]/(2*(rhom-rhoi)))*(np.exp(-Lambda*a_basin)*np.cos(Lambda*a_basin)-np.exp(-Lambda*b_basin)*np.cos(Lambda*b_basin))
    percentile_95p45_wxb = -(rhoc*percentile_95p45_thickness_etar[time1]/(2*(rhom-rhoi)))*(np.exp(-Lambda*a_basin)*np.cos(Lambda*a_basin)-np.exp(-Lambda*b_basin)*np.cos(Lambda*b_basin))
    
    mean_wxr = (rhoc*mean_thickness_etar[time1]/(2*(rhom-rhoi)))*(2-np.exp(-Lambda*b_range)*np.cos(Lambda*b_range)-np.exp(-Lambda*a_range)*np.cos(Lambda*a_range))
    percentile_4p55_wxr = (rhoc*percentile_4p55_thickness_etar[time1]/(2*(rhom-rhoi)))*(2-np.exp(-Lambda*b_range)*np.cos(Lambda*b_range)-np.exp(-Lambda*a_range)*np.cos(Lambda*a_range))
    percentile_95p45_wxr = (rhoc*percentile_95p45_thickness_etar[time1]/(2*(rhom-rhoi)))*(2-np.exp(-Lambda*b_range)*np.cos(Lambda*b_range)-np.exp(-Lambda*a_range)*np.cos(Lambda*a_range))
    
    mean_Hxb = -(meanElevation_Hb[time1]*2)*(b_basin/max(b_basin))+(meanElevation_Hb[time1]*2)
    percentile_4p55_Hxb = -(percentile_4p55_Elevation_Hb[time1]*2)*(b_basin/max(b_basin))+(percentile_4p55_Elevation_Hb[time1]*2)
    percentile_95p45_Hxb = -(percentile_95p45_Elevation_Hb[time1]*2)*(b_basin/max(b_basin))+(percentile_95p45_Elevation_Hb[time1]*2)
    
    mean_Hxr = max(mean_Hxb)+(((meanElevation_Hr[time1]-max(mean_Hxb))*2)*(b_range/max(b_range)))
    percentile_4p55_Hxr = max(percentile_4p55_Hxb)+(((percentile_4p55_Elevation_Hr[time1]-max(percentile_4p55_Hxb))*2)*(b_range/max(b_range)))
    percentile_95p45_Hxr = max(percentile_95p45_Hxb)+(((percentile_95p45_Elevation_Hr[time1]-max(percentile_95p45_Hxb))*2)*(b_range/max(b_range)))
    
    fig = plt.figure(figsize =(8,4))
    ax1 = fig.add_subplot(111)
    
    ax1.plot(x_range/1000, -mean_wxr, color='steelblue'); ax1.plot(x_basin/1000, -mean_wxb, color='darkorange')
    ax1.plot(x_range/1000, mean_Hxr, color='steelblue'); ax1.plot(x_basin/1000, mean_Hxb, color='darkorange')
    
    ax1.plot(x_range/1000, -percentile_4p55_wxr, color='steelblue', linestyle='--'); ax1.plot(x_basin/1000, -percentile_4p55_wxb, color='darkorange', linestyle='--')
    ax1.plot(x_range/1000, percentile_4p55_Hxr, color='steelblue', linestyle='--'); ax1.plot(x_basin/1000, percentile_4p55_Hxb, color='darkorange', linestyle='--')
   
    ax1.plot(x_range/1000, -percentile_95p45_wxr, color='steelblue', linestyle='--'); ax1.plot(x_basin/1000, -percentile_95p45_wxb, color='darkorange', linestyle='--')
    ax1.plot(x_range/1000, percentile_95p45_Hxr, color='steelblue', linestyle='--'); ax1.plot(x_basin/1000, percentile_95p45_Hxb, color='darkorange', linestyle='--')
    
    ax1.fill_between(x_range/1000, -percentile_4p55_wxr, -percentile_95p45_wxr, facecolor = 'steelblue', alpha='0.33')
    ax1.fill_between(x_basin/1000, -percentile_4p55_wxb, -percentile_95p45_wxb, facecolor = 'darkorange', alpha='0.33')
    ax1.fill_between(x_range/1000, percentile_4p55_Hxr, percentile_95p45_Hxr, facecolor = 'steelblue', alpha='0.33')
    ax1.fill_between(x_basin/1000, percentile_4p55_Hxb, percentile_95p45_Hxb, facecolor = 'darkorange', alpha='0.33')
    
    lineBL = lines.Line2D((0,(Lr+Lb)/1000), (0,0), linestyle='-.', color='k')
    patchR = patches.Rectangle((0, meanElevation_Hr[time1]), Lr/1000, -mean_thickness_etar[time1], fill=False, edgecolor='k', linewidth=1.4)
    patchB = patches.Rectangle((Lr/1000, meanElevation_Hb[time1]), Lb/1000, -mean_thickness_etab[time1], fill=False, edgecolor='k', linewidth=1.4)
    ax1.add_patch(patchR); ax1.add_patch(patchB); ax1.add_artist(lineBL)
    
    plt.text(0.95, 0.9, 'Time = %s Ma'%round((time2-66),1), horizontalalignment='right', verticalalignment='top', transform=ax1.transAxes)
    plt.text(0.95, 0.85, 'Convergence = %s mm/yrs'%round((Vt[time1]*1000),1), horizontalalignment='right', verticalalignment='top', transform=ax1.transAxes)
    if time1 < 4300:
      plt.text(0.95, 0.8, 'Syn-orogenesis', horizontalalignment='right', verticalalignment='top', transform=ax1.transAxes)
    if time1 >= 4300:
      plt.text(0.95, 0.8, 'Post-orogenesis', horizontalalignment='right', verticalalignment='top', transform=ax1.transAxes)
    
    ax1.axis(topo_view)
    ax1.set_xlabel('Distance (km)'); ax1.set_ylabel('Elevation (m)')
    plt.ioff()
    
    if time3 < 10:
      plt.savefig(path_to_save1 + figure_name + '000'+ str(time3) +'.png', dpi=300)
      plt.close(fig)
    if time3 >= 10 and time3 < 100:
      plt.savefig(path_to_save1 + figure_name + '00'+ str(time3) +'.png', dpi=300)
      plt.close(fig)
    if time3 >= 100:
      plt.savefig(path_to_save1 + figure_name + '0'+ str(time3) +'.png', dpi=300)
      plt.close(fig)
    
    time3 += 1

#%%=========================================================================%%#
def model_cross_section(time1, result_Hr_list, result_Hb_list, result_etar_list, result_etab_list, result_Te_list, rhos, rhoc, rhom, rhoi, E, nu, g, Lr, Vt):
    
  general_view = [-5, 170, -500, 2500]
  
  meanElevation_Hr = result_Hr_list.mean(axis=0)
  percentile_4p55_Elevation_Hr = np.percentile(result_Hr_list, 4.55, axis=0)
  percentile_95p45_Elevation_Hr = np.percentile(result_Hr_list, 95.45, axis=0)
  
  meanElevation_Hb = result_Hb_list.mean(axis=0)
  percentile_4p55_Elevation_Hb = np.percentile(result_Hb_list, 4.55, axis=0)
  percentile_95p45_Elevation_Hb = np.percentile(result_Hb_list, 95.45, axis=0)
    
  mean_thickness_etar = result_etar_list.mean(axis=0)
  percentile_4p55_thickness_etar = np.percentile(result_etar_list, 4.55, axis=0)
  percentile_95p45_thickness_etar = np.percentile(result_etar_list, 95.45, axis=0)
  
  mean_thickness_etab = result_etab_list.mean(axis=0)
  #percentile_4p55_thickness_etab = np.percentile(result_etab_list, 4.55, axis=0)
  #percentile_95p45_thickness_etab = np.percentile(result_etab_list, 95.45, axis=0)  
  
  Te = np.mean(result_Te_list)
  D = E*Te**3/(12*(1-nu**2))
  Lambda = ((rhom-rhoi)*g/(4*D))**(1./4.)
  #tlaml = 2*Lambda*Lr
  Lb = 3*np.pi/(4*Lambda)-Lr
  #Llm = Lambda*Lr
  
  x_range = np.linspace(0.,Lr, 1000)
  x_basin = np.linspace(Lr,(Lr+Lb) , 1000)
  a_range = x_range + Lr                      
  b_range = Lr-x_range                        
  b_basin = x_basin-Lr                        
  a_basin = 2*Lr+b_basin
  
  mean_wxb = -(rhoc*mean_thickness_etar[time1]/(2*(rhom-rhoi)))*(np.exp(-Lambda*a_basin)*np.cos(Lambda*a_basin)-np.exp(-Lambda*b_basin)*np.cos(Lambda*b_basin))
  percentile_4p55_wxb = -(rhoc*percentile_4p55_thickness_etar[time1]/(2*(rhom-rhoi)))*(np.exp(-Lambda*a_basin)*np.cos(Lambda*a_basin)-np.exp(-Lambda*b_basin)*np.cos(Lambda*b_basin))
  percentile_95p45_wxb = -(rhoc*percentile_95p45_thickness_etar[time1]/(2*(rhom-rhoi)))*(np.exp(-Lambda*a_basin)*np.cos(Lambda*a_basin)-np.exp(-Lambda*b_basin)*np.cos(Lambda*b_basin))
  
  mean_wxr = (rhoc*mean_thickness_etar[time1]/(2*(rhom-rhoi)))*(2-np.exp(-Lambda*b_range)*np.cos(Lambda*b_range)-np.exp(-Lambda*a_range)*np.cos(Lambda*a_range))
  percentile_4p55_wxr = (rhoc*percentile_4p55_thickness_etar[time1]/(2*(rhom-rhoi)))*(2-np.exp(-Lambda*b_range)*np.cos(Lambda*b_range)-np.exp(-Lambda*a_range)*np.cos(Lambda*a_range))
  percentile_95p45_wxr = (rhoc*percentile_95p45_thickness_etar[time1]/(2*(rhom-rhoi)))*(2-np.exp(-Lambda*b_range)*np.cos(Lambda*b_range)-np.exp(-Lambda*a_range)*np.cos(Lambda*a_range))
  
  mean_Hxb = -(meanElevation_Hb[time1]*2)*(b_basin/max(b_basin))+(meanElevation_Hb[time1]*2)
  percentile_4p55_Hxb = -(percentile_4p55_Elevation_Hb[time1]*2)*(b_basin/max(b_basin))+(percentile_4p55_Elevation_Hb[time1]*2)
  percentile_95p45_Hxb = -(percentile_95p45_Elevation_Hb[time1]*2)*(b_basin/max(b_basin))+(percentile_95p45_Elevation_Hb[time1]*2)
  
  mean_Hxr = max(mean_Hxb)+(((meanElevation_Hr[time1]-max(mean_Hxb))*2)*(b_range/max(b_range)))
  percentile_4p55_Hxr = max(percentile_4p55_Hxb)+(((percentile_4p55_Elevation_Hr[time1]-max(percentile_4p55_Hxb))*2)*(b_range/max(b_range)))
  percentile_95p45_Hxr = max(percentile_95p45_Hxb)+(((percentile_95p45_Elevation_Hr[time1]-max(percentile_95p45_Hxb))*2)*(b_range/max(b_range)))
  
  fig = plt.figure(figsize =(8,4))
  ax1 = fig.add_subplot(111)
  
  ax1.plot(x_range/1000, -mean_wxr, color='steelblue'); ax1.plot(x_basin/1000, -mean_wxb, color='darkorange')
  ax1.plot(x_range/1000, mean_Hxr, color='steelblue'); ax1.plot(x_basin/1000, mean_Hxb, color='darkorange')
  
  ax1.plot(x_range/1000, -percentile_4p55_wxr, color='steelblue', linestyle='--'); ax1.plot(x_basin/1000, -percentile_4p55_wxb, color='darkorange', linestyle='--')
  ax1.plot(x_range/1000, percentile_4p55_Hxr, color='steelblue', linestyle='--'); ax1.plot(x_basin/1000, percentile_4p55_Hxb, color='darkorange', linestyle='--')
 
  ax1.plot(x_range/1000, -percentile_95p45_wxr, color='steelblue', linestyle='--'); ax1.plot(x_basin/1000, -percentile_95p45_wxb, color='darkorange', linestyle='--')
  ax1.plot(x_range/1000, percentile_95p45_Hxr, color='steelblue', linestyle='--'); ax1.plot(x_basin/1000, percentile_95p45_Hxb, color='darkorange', linestyle='--')
  
  ax1.fill_between(x_range/1000, -percentile_4p55_wxr, -percentile_95p45_wxr, facecolor = 'steelblue', alpha='0.33')
  ax1.fill_between(x_basin/1000, -percentile_4p55_wxb, -percentile_95p45_wxb, facecolor = 'darkorange', alpha='0.33')
  ax1.fill_between(x_range/1000, percentile_4p55_Hxr, percentile_95p45_Hxr, facecolor = 'steelblue', alpha='0.33')
  ax1.fill_between(x_basin/1000, percentile_4p55_Hxb, percentile_95p45_Hxb, facecolor = 'darkorange', alpha='0.33')
  
  lineBL = lines.Line2D((0,(Lr+Lb)/1000), (0,0), linestyle='-.', color='k')
  patchR = patches.Rectangle((0, meanElevation_Hr[time1]), Lr/1000, -mean_thickness_etar[time1], fill=False, edgecolor='k', linewidth=1.4)
  patchB = patches.Rectangle((Lr/1000, meanElevation_Hb[time1]), Lb/1000, -mean_thickness_etab[time1], fill=False, edgecolor='k', linewidth=1.4)
  ax1.add_patch(patchR); ax1.add_patch(patchB); ax1.add_artist(lineBL)
  
  ax1.axis(general_view)

#%%=========================================================================%%#
  
def figure_misfit_vs_parameters_test(kappab_list, kappar_list, Te_list, global_misfit, Te_min, Te_max, kappar_min, kappar_max, kappab_min, kappab_max):

  #kappa_ratio_list = kappab_list/kappar_list
  
  x1 = Te_list; y1 = kappar_list; z1 = global_misfit
  grid_x1, grid_y1 = np.mgrid[Te_min:Te_max:100j, kappar_min:kappar_max:100j]
  grid_z1 = griddata((x1, y1), z1, (grid_x1, grid_y1), method='nearest')
  grid_z1 = scipy.ndimage.zoom(grid_z1, 3)
  
  x2 = Te_list; y2 = kappar_list; z2 = global_misfit
  grid_x2, grid_y2 = np.mgrid[Te_min:Te_max:100j, kappar_min:kappar_max:100j]
  grid_z2 = griddata((x2, y2), z2, (grid_x2, grid_y2), method='nearest')
  grid_z2 = scipy.ndimage.spline_filter(grid_z2, 5)
  
  x3 = Te_list; y3 = kappar_list; z3 = global_misfit
  grid_x3, grid_y3 = np.mgrid[Te_min:Te_max:100j, kappar_min:kappar_max:100j]
  grid_z3 = griddata((x3, y3), z3, (grid_x3, grid_y3), method='nearest')
  grid_z3 = scipy.ndimage.zoom(grid_z3, 3)
  
  x4 = Te_list; y4 = kappar_list; z4 = global_misfit
  grid_x4, grid_y4 = np.mgrid[Te_min:Te_max:100j, kappar_min:kappar_max:100j]
  grid_z4 = griddata((x4, y4), z4, (grid_x4, grid_y4), method='nearest')
  grid_z4 = scipy.ndimage.gaussian_filter(grid_z4, 1)
  
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
  ax2.set_xlabel('Lithosphere elastic thickness (km)'); ax2.set_ylabel('Range transport coefficient (m$^{2}$.yr$^{-1}$)')
  ax3.set_xlabel('Lithosphere elastic thickness (km)'); ax3.set_ylabel('Range transport coefficient (m$^{2}$.yr$^{-1}$)')
  ax4.set_xlabel('Lithosphere elastic thickness (km)'); ax4.set_ylabel('Range transport coefficient (m$^{2}$.yr$^{-1}$)')
  
  Te_label = ['15','20','25','30','35','40']
  ax1.set_xticklabels(Te_label); ax2.set_xticklabels(Te_label); ax3.set_xticklabels(Te_label); ax4.set_xticklabels(Te_label)
  kappar_label = ['100','1080','2060','3040','4020','5000']
  ax1.set_yticklabels(kappar_label); ax2.set_yticklabels(kappar_label); ax3.set_yticklabels(kappar_label); ax4.set_yticklabels(kappar_label)
  
  fig.subplots_adjust(left=0.1, bottom=0.085, right=0.925, top=0.975, wspace=0.5, hspace=0.25)
  #fig.savefig('Z:\\Post-orogenic sediment flux to continental margins\\Publication_Report\\Paper_Post-orogen_Pyrenees+Aquitaine\\long_version\\FIGURE3_misfit+contour_vs_parameters.pdf', dpi=1200)
