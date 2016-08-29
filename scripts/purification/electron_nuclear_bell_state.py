"""
Evaluates the density matrix of a combined electron nuclear spin state.
Based upon 12 measurements for different initial measurements of the electron

TODO: proper error propagation onto elements of the dnesity matrix.
"""


import numpy as np
import os
from analysis.lib.tools import toolbox as tb
from analysis.lib.m2.ssro import mbi; reload(mbi)
from matplotlib import pyplot as plt

import h5py



####################
#                  # 
# Helper functions #
#                  #
####################

def generate_pauli_matrices():
    return [np.matrix([[1,0],[0,1]],dtype=complex),np.matrix([[0,1],[1,0]],dtype=complex),np.matrix([[0,-1j],[1j,0]],dtype=complex),np.matrix([[1,0],[0,-1]],dtype=complex)]

def get_RO_results(folder,ssro_calib_folder):
    a = mbi.MBIAnalysis(folder)
    a.get_sweep_pts()
    a.get_readout_results(name='adwindata')
    a.get_electron_ROC(ssro_calib_folder)
    y_a= ((a.p0.reshape(-1)[:])-0.5)*2
    y_err_a = 2*a.u_p0.reshape(-1)[:] 

    return y_a,y_err_a

def get_correlations(**kw):

    ### pull data
    ssro_calib_timestamp = kw.pop('ssro_calib_timestamp',None)

    if ssro_calib_timestamp == None: 
        ssro_calib_folder = tb.latest_data('SSROCalibration')
    else:
        ssro_dstmp, ssro_tstmp = tb.verify_timestamp(ssro_calib_timestamp)
        ssro_calib_folder = tb.datadir + '\\'+ssro_dstmp+'\\'+ssro_tstmp+'_AdwinSSRO_SSROCalibration_Pippin_SIL2'
        print ssro_calib_folder

    tomo_pulses_p = ['none','x','y'] ### used when looking for folders
    tomo_pulses_m = ['X','mx','my'] ### used when looking for folders
    f_list_pos_p,f_list_neg_p,f_list_pos_m,f_list_neg_m = [],[],[],[]
    exp_values,exp_vals_u = np.array([]),np.array([]) ### this will be a list with all combined expectation values such as XX or YZ
    sweep_pts = [] ### this will be a list of the bases associated expectation values

    ### carbon 1-qubit correlations
    c_x,c_y,c_z = 0.,0.,0.
    c_x_u,c_y_u,c_z_u = 0.,0.,0.

    ### this dictionary is used to translate the electron RO pulse into a measurement basis
    tomo_pulse_translation_dict = {'none': 'Z','X':'Z','x':'X','mx':'X','y':'Y','my':'Y'}

    for p,m in zip(tomo_pulses_p,tomo_pulses_m):
        f_list_pos_p.append(tb.latest_data('el_13C_dm_'+p+'_positive' ,return_timestamp = False,**kw))
        f_list_neg_p.append(tb.latest_data('el_13C_dm_'+p+'_negative' ,return_timestamp = False,**kw))
        f_list_pos_m.append(tb.latest_data('el_13C_dm_'+m+'_positive' ,return_timestamp = False,**kw))
        f_list_neg_m.append(tb.latest_data('el_13C_dm_'+m+'_negative' ,return_timestamp = False,**kw))

        #### now also calculate the measured contrast
        y_a,y_err_a = get_RO_results(f_list_pos_p[-1],ssro_calib_folder)
        y_b,y_err_b= get_RO_results(f_list_neg_p[-1],ssro_calib_folder)


        y_c,y_err_c = get_RO_results(f_list_pos_m[-1],ssro_calib_folder)
        y_d,y_err_d= get_RO_results(f_list_neg_m[-1],ssro_calib_folder)

        y_pos_electron   = (y_a - y_b)/2.
        y_pos_electron_u = 1./2*(y_err_a**2 + y_err_b**2)**0.5
        y_neg_electron   = (y_c - y_d)/2.
        y_neg_electron_u =  1./2*(y_err_c**2 + y_err_d**2)**0.5

        #### now do the final combination of results regarding the electron RO basis to get the combined expectation value
        y   = (y_pos_electron - y_neg_electron)/2.
        y_u = 1./2*(y_pos_electron_u**2 + y_neg_electron_u**2)**0.5

        ### Combine data
        exp_values = np.append(exp_values,y)
        exp_vals_u = np.append(exp_vals_u,y_u)

        #### we also need to calculate the single qubit expectation values.
        #### easy for the electron spin. One simply adds the electron outcomes and divides by 2.
        #### for the nuclear spin one has to add up several RO bases and look at the correlations there:
        exp_val_sum = (y_pos_electron + y_neg_electron)/2.
        exp_val_sum_u = 1./2*(y_pos_electron_u**2 + y_neg_electron_u**2)


        ### update the running average of the nuclear spin expectation values:
        c_x,c_y,c_z = c_x + exp_val_sum[0]/3.,c_y + exp_val_sum[1]/3., c_z + exp_val_sum[2]/3.
        c_x_u,c_y_u,c_z_u = c_x_u + exp_val_sum_u[0]/9.,c_y_u + exp_val_sum_u[1]/9., c_z_u + exp_val_sum_u[2]/9.

        ### add the electron values
        y_electron = np.sum(exp_val_sum)/3.
        y_electron_u = np.sum(exp_val_sum_u/9.)**0.5
        exp_values = np.append(exp_values,y_electron)
        exp_vals_u = np.append(exp_vals_u,y_electron_u)




        ### write up the combined analysis bases:
        a = mbi.MBIAnalysis(f_list_pos_p[-1])
        a.get_sweep_pts()
        a.get_readout_results(name='adwindata')
        e_base = np.array([tomo_pulse_translation_dict[p]]*4)
        c_base = a.adgrp.attrs['Tomography Bases'].flatten()
        comb_base = np.core.defchararray.add(e_base, np.append(c_base,'I')) ### add identity to 13C bases as this is evaluated on the fly
        sweep_pts = np.append(sweep_pts,comb_base)


    ### the nuclear spin single qubit correlations have been esimated in a running average (transposed correlations)
    ### and are now added to the results (square root still needs to be taken for the uncertainties)
    sweep_pts = np.append(sweep_pts,['IX','IY','IZ'])
    exp_values = np.append(exp_values,np.array([c_x,c_y,c_z]))
    exp_vals_u = np.append(exp_vals_u,np.array([c_x_u,c_y_u,c_z_u]))



    return f_list_pos_p[0],sweep_pts,exp_values,exp_vals_u

def add_single_qubit_correlations(sweep_pts,exp_values,exp_vals_u):
    """
    gets two qubit correlations as input, traces over one of the two and adds individual correlations.
    returns a dictionary where the keys are te measurement basis and the outcomes are the associated values
    """
    vals_dict,vals_dict_u = {},{}
    ### load up dictionaries
    # for t,exp,exp_u in sweep_pts,exp_values,exp_vals_u:
    #     vals_dict.update(t:exp)
    #     vals_dict_u.update(t:exp_u)


def carbon_ROC(exp,exp_u,folder):
    """
    setup dependent carbon read-correction coefficients are implemented here
    """
    if 'Pippin' in folder:
        ROC_coeff = 
        ROC_coeff_u = 
    else:
        ROC_coeff = 
        ROC_coeff_u = 

    return exp/ROC_coeff,exp_u

def electron_carbon_density_matrix(**kw):
    """
    This function currently neglects error propagation
    """

    folder,sweep_pts,exp_values,exp_vals_u = get_correlations(**kw)
    

    paulis = generate_pauli_matrices()
    ### initialize the dm via the identity correlations
    dm = np.kron(paulis[0],paulis[0])/4.

    ### basis definition
    t_dict = {'I' : 0, 'X':1, 'Y':2, 'Z':3}

    ### put dm together
    for t,exp,exp_u in zip(sweep_pts,exp_values,exp_vals_u)
            if t+t2 == 'II':
                continue
                
            sigma_kron = np.kron(paulis[t_dict[t[0]]],paulis[t_dict[t[1]]])


            ### carbon ROC necessary?
            if t[1] =='I':
                dm +=exp*sigma_kron/4.

            else:
                exp,exp_u = carbon_ROC(exp,exp_u,folder)

                dm +=exp*sigma_kron/4.

