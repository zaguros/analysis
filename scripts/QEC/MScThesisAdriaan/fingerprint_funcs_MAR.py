'''
Module with auxilairy functions for the analysis of fingerprint data
Copy of Fingerprint funcs edited to be used for MSc Thesis Adriaan and remain unchanged
'''
import numpy as np
from analysis.lib.tools import toolbox
from analysis.lib.m2.ssro import mbi

import hyperfine_params_MAR as module_hyperfine_params; reload(module_hyperfine_params)
hf = module_hyperfine_params.hyperfine_params


def load_mult_dat(timestamp,
			number_of_msmts,
			x_axis_start 		= 2.0,
			x_axis_step 		= 1.0,
			x_axis_pts_per_msmnt= 101,
			ssro_calib_folder	=''):
   '''
   function to load and combine multiple msmts.
   '''
   cum_pts = 0
   for kk in range(number_of_msmts):
       folder = toolbox.data_from_time(timestamp)
       a = mbi.MBIAnalysis(folder)
       a.get_sweep_pts()
       a.get_readout_results(name='measurement' + str(kk))
       a.get_electron_ROC(ssro_calib_folder=ssro_calib_folder)
       cum_pts += a.pts

       if kk == 0:
           cum_sweep_pts = np.linspace(x_axis_start, x_axis_start+x_axis_step, x_axis_pts_per_msmnt)#a.sweep_pts
           cum_p0 = a.p0
           cum_u_p0 = a.u_p0
       else:
           cum_sweep_pts = np.concatenate((cum_sweep_pts, np.linspace(x_axis_start+kk*x_axis_step, x_axis_start+(kk+1)*x_axis_step, x_axis_pts_per_msmnt)))
           cum_p0 = np.concatenate((cum_p0, a.p0))
           cum_u_p0 = np.concatenate((cum_u_p0, a.u_p0))

   a.pts   = cum_pts
   a.sweep_pts = cum_sweep_pts
   a.p0    = cum_p0
   a.u_p0  = cum_u_p0

   return a, folder

def load_mult_dat_2(timestamps, number_of_msmts, ssro_calib_folders=['']):
   cum_pts = 0
   for tt in range(len(timestamps)):
    for kk in range(number_of_msmts[tt]):
      folder = toolbox.data_from_time(timestamps[tt])
      a = mbi.MBIAnalysis(folder)
      a.get_sweep_pts()
      a.get_readout_results(name='measurement' + str(kk))
      a.get_electron_ROC(ssro_calib_folder=ssro_calib_folders[tt])
      cum_pts += a.pts

      if kk == 0 and tt==0:
        cum_sweep_pts = np.linspace(2.0, 2.5, 51)#a.sweep_pts
        cum_p0 = a.p0
        cum_u_p0 = a.u_p0
      elif tt==0:
        cum_sweep_pts = np.concatenate((cum_sweep_pts, np.linspace(2.0+kk*(50)*10e-3, 2.5+kk*(50)*10e-3, 51)))
        cum_p0 = np.concatenate((cum_p0, a.p0))
        cum_u_p0 = np.concatenate((cum_u_p0, a.u_p0))
      elif tt==1:
        cum_sweep_pts = np.concatenate((cum_sweep_pts, np.linspace(2.0+(kk+90)*(50)*10e-3, 2.5+(kk+90)*(50)*10e-3, 51)))
        cum_p0 = np.concatenate((cum_p0, a.p0))
        cum_u_p0 = np.concatenate((cum_u_p0, a.u_p0))

   a.pts   = cum_pts
   a.sweep_pts = cum_sweep_pts
   a.p0    = cum_p0
   a.u_p0  = cum_u_p0

   return a, folder

def get_hyperfine_params(ms = 'plus', carbon_spins = 'all'):

	if carbon_spins == 'all':

		HF_perp = []
		HF_par 	= []

		for kk in range(len(hf)):
			carbon_string 		= 'C' + str(kk+1)
			HF_perp.append(hf[carbon_string]['perp'])
			if ms == 'plus':
				HF_par.append(hf[carbon_string]['par'])
			elif ms == 'min':
				HF_par.append(-1*hf[carbon_string]['par'])

	return HF_perp, HF_par



