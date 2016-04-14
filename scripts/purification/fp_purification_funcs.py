''' 
Module with auxilairy functions for the analysis of fingerprint data, THT
'''
import numpy as np
from analysis.lib.tools import toolbox
from analysis.lib.m2.ssro import mbi
import copy as cp

import hyperfine_params as module_hyperfine_params; reload(module_hyperfine_params) 
hf = module_hyperfine_params.hyperfine_params


def load_mult_dat(timestamp,   
      ssro_calib_folder =''):
  ''' 
  function to load and combine multiple msmts. 
  '''
  cum_pts = 0
  kk = 0
  while True:
    folder = toolbox.data_from_time(timestamp)
    a = mbi.MBIAnalysis(folder)
    a.get_sweep_pts()
    try:
      a.get_readout_results(name='measurement' + str(kk))
    except:
      print 'found no more data, stopping the loop after {} datasets'.format(kk)
      break
    a.get_electron_ROC(ssro_calib_folder=ssro_calib_folder)
    cum_pts += a.pts

    if kk == 0:
      cum_p0 = a.p0
      cum_u_p0 = a.u_p0
    else:
      cum_p0 = np.concatenate((cum_p0, a.p0))
      cum_u_p0 = np.concatenate((cum_u_p0, a.u_p0))
    kk+=1

  print 'Number of datasets is',kk

  step = a.sweep_pts[-1]-a.sweep_pts[0]
  pts_per_msmt = len(a.sweep_pts)
  x_axis_start = a.sweep_pts[-1]-step*(kk)
  cum_sweep_pts = np.linspace(x_axis_start,a.sweep_pts[-1],cum_pts)

  a.pts   = cum_pts
  a.sweep_pts = cum_sweep_pts
  a.p0    = cum_p0
  a.u_p0  = cum_u_p0


  return a, folder

# def get_hyperfine_params(ms = 'plus', carbon_spins = 'all'):
  
#   if carbon_spins == 'all':

#     HF_perp = []
#     HF_par  = []

#     for kk in range(len(hf)):
#       carbon_string     = 'C' + str(kk+1) 
#       HF_perp.append(hf[carbon_string]['perp'])
#       if ms == 'plus':
#         HF_par.append(hf[carbon_string]['par'])
#       elif ms == 'min':
#         HF_par.append(-1*hf[carbon_string]['par'])
  
#   return HF_perp, HF_par

def get_hyperfine_params(ms = 'plus', carbon_spins = 'all',NV = None):
    
    HF_perp = []
    HF_par 	= []
    print 'NV (get hyperfine_params): ' + str(NV)
    if NV =='Pippin_SIL1':
      hf = module_hyperfine_params.hyperfine_params_pippin_SIL1_msm1
    elif NV =='Pippin_SIL3':
      hf = module_hyperfine_params.hyperfine_params_pippin_SIL3_msm1
    elif NV == 'hans':
      hf = module_hyperfine_params.hyperfine_params_hans_SIL1_msm1
    else: 
      hf = module_hyperfine_params.hyperfine_params
    print 'hf = ' + str(hf)
    if carbon_spins == 'all':
      for kk in range(len(hf)):
        carbon_string     = 'C' + str(kk+1) 
        HF_perp.append(hf[carbon_string]['perp'])
        if ms == 'plus':
          HF_par.append(hf[carbon_string]['par'])
        elif ms == 'min':
          HF_par.append(-1*hf[carbon_string]['par'])

    else:
      for carbon in carbon_spins:
        carbon_string     = 'C' + str(carbon) 
        HF_perp.append(hf[carbon_string]['perp'])
        if ms == 'plus':
          HF_par.append(hf[carbon_string]['par'])
        elif ms == 'min':
          HF_par.append(-1*hf[carbon_string]['par'])


    return HF_perp, HF_par
