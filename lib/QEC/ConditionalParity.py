from analysis.lib.m2 import m2
from analysis.lib.m2.ssro import ssro
from analysis.lib.math import error
from analysis.lib.tools import toolbox
from analysis.lib.m2.ssro import mbi
import numpy as np

#NOTE: Function could be moved to analysis.lib.m2.ssro folder when complete and checked for errors.

class ConditionalParityAnalysis(mbi.MBIAnalysis):
    '''
    Class used to extract and analyze C13 experiments with conditional outcomes.
    Child-class of MBI analysis (as all C13 experiments are based on the MBI architecture).
    '''

    def get_readout_results(self, name='',post_select = False,post_select_QEC = False):
        '''
        Get the readout results.
        self.ssro_results contains the readout results (sum of the photons for
            a given sweep point and readout in a sequence)
        self.normalized_ssro contains the normalized result (i.e., probability
            for getting a photon)
        
        if post_select      is true, select on the true or false condition (single parity measurement)
        if post_select_QEC  is true, select on the syndrome outcome (two parity measurements)
        else                use the original get_readout_results from the MBI class
 
        '''

        if post_select == True:
            
            self.post_select        = True
            self.result_corrected   = False

            adwingrp    = self.adwingrp(name)
            self.adgrp  = adwingrp

            self.pts      = adwingrp.attrs['sweep_length']
            self.reps     = adwingrp.attrs['reps_per_ROsequence']
            self.readouts = adwingrp.attrs['nr_of_ROsequences']

            ### Step 0 extract data from hdf5 file
            self.parity_result = adwingrp['parity_RO_results'].value 
            self.ssro_results  = adwingrp['ssro_results'].value       

            ### Step 1 Multiply results with post selection parameter
            ssro_results_0 = self.parity_result     * self.ssro_results
            ssro_results_1 = (1-self.parity_result) * self.ssro_results

            ### Step 2 reshape
            self.parity_result  = self.parity_result.reshape((-1,self.pts,self.readouts)).sum(axis=0)
            self.ssro_results_0 = ssro_results_0.reshape((-1,self.pts,self.readouts)).sum(axis=0)
            self.ssro_results_1 = ssro_results_1.reshape((-1,self.pts,self.readouts)).sum(axis=0)

            ### Step 3 normalization and uncertainty
            self.normalized_ssro_0 = self.ssro_results_0/(self.parity_result).astype('float')
            self.u_normalized_ssro_0 = (self.normalized_ssro_0*(1-self.normalized_ssro_0)/(self.parity_result))**0.5
            
            self.normalized_ssro_1 = self.ssro_results_1/(self.reps-self.parity_result).astype('float')
            self.u_normalized_ssro_1 = (self.normalized_ssro_1*(1-self.normalized_ssro_1)/(self.reps-self.parity_result))**0.5

            # print 'Probabilities ms=0 and ms=-1'
            # print np.average(self.parity_result/self.reps.astype('float'))
            # print np.average((self.reps-self.parity_result)/self.reps.astype('float'))

        elif post_select_QEC == True:
            self.post_select = True
            self.result_corrected = False

            adwingrp = self.adwingrp(name)
            self.adgrp = adwingrp

            self.pts        = adwingrp.attrs['sweep_length']
            self.reps       = adwingrp.attrs['reps_per_ROsequence']
            self.readouts   = adwingrp.attrs['nr_of_ROsequences']

            ### Step 0 extract data from hdf5 file
            self.parity_result = adwingrp['parity_RO_results'].value
            self.ssro_results  = adwingrp['ssro_results'].value 

            parity_a_result = self.parity_result[0::2]  ### The two parity outcomes are stored sequentially in an array 
            parity_b_result = self.parity_result[1::2] 
            

            ### Step 1 Multiply results with post selection parameter
            ssro_results_00 = parity_a_result*parity_b_result*self.ssro_results
            ssro_results_01 = parity_a_result*(1-parity_b_result)*self.ssro_results
            ssro_results_10 = (1-parity_a_result)*parity_b_result*self.ssro_results
            ssro_results_11 = (1-parity_a_result)*(1-parity_b_result)*self.ssro_results
           
            # print'((1-parity_a_result)*parity_b_result)'
            # print ((1-parity_a_result)*parity_b_result).reshape((-1,self.pts,self.readouts)).sum(axis=0)


            ### Step 2 reshape
            parity_result_00 = (parity_a_result*parity_b_result).reshape((-1,self.pts,self.readouts)).sum(axis=0)
            parity_result_01 = (parity_a_result*(1-parity_b_result)).reshape((-1,self.pts,self.readouts)).sum(axis=0)
            parity_result_10 = ((1-parity_a_result)*parity_b_result).reshape((-1,self.pts,self.readouts)).sum(axis=0)
            parity_result_11 = ((1-parity_a_result)*(1-parity_b_result)).reshape((-1,self.pts,self.readouts)).sum(axis=0)



            self.ssro_results_00 = ssro_results_00.reshape((-1,self.pts,self.readouts)).sum(axis=0)
            self.ssro_results_01 = ssro_results_01.reshape((-1,self.pts,self.readouts)).sum(axis=0)
            self.ssro_results_10 = ssro_results_10.reshape((-1,self.pts,self.readouts)).sum(axis=0)
            self.ssro_results_11 = ssro_results_11.reshape((-1,self.pts,self.readouts)).sum(axis=0)


            ### Step 3 normalization and uncertainty
            self.normalized_ssro_00 = self.ssro_results_00/(parity_result_00).astype('float')
            self.u_normalized_ssro_00 = (self.normalized_ssro_00*(1-self.normalized_ssro_00)/(parity_result_00))**0.5
            
            self.normalized_ssro_01 = self.ssro_results_01/(parity_result_01).astype('float')
            self.u_normalized_ssro_01 = (self.normalized_ssro_01*(1-self.normalized_ssro_01)/(parity_result_01))**0.5
            
            self.normalized_ssro_10 = self.ssro_results_10/(parity_result_10 ).astype('float')
            self.u_normalized_ssro_10 = (self.normalized_ssro_10*(1-self.normalized_ssro_10)/(parity_result_10))**0.5
            
            self.normalized_ssro_11 = self.ssro_results_11/(parity_result_11).astype('float')
            self.u_normalized_ssro_11 = (self.normalized_ssro_11*(1-self.normalized_ssro_11)/(parity_result_11))**0.5
            # print 'test'
            ### 'Probabilities 00, 01, 10, 11'
            # print self.reps.astype('float')
            # print float(len(self.ssro_results))
            self.p00 = (parity_result_00/float(len(self.ssro_results)/self.pts))
            # print  self.p00
            self.p01 = ((parity_result_01)/float(len(self.ssro_results)/self.pts))
            # print parity_result_01
            # print  self.p01
            self.p10 = (parity_result_10/float(len(self.ssro_results)/self.pts))
            # print parity_result_10
            # print  self.p10
            self.p11 = ((parity_result_11)/float(len(self.ssro_results)/self.pts))
            # print  self.p11

        else:
            mbi.MBIAnalysis.get_readout_results(self,name) #NOTE: super cannot be used as this is an "old style class"
            self.post_select = False

    def get_electron_ROC(self, ssro_calib_folder='',post_select_QEC = False):
        '''
        Performs Readout Correction, needs to be updated to correct for post selected results and apply error-bars correctly.
        '''
        
        if post_select_QEC == True and self.post_select == True:
            if ssro_calib_folder == '':
                ssro_calib_folder = toolbox.latest_data('SSRO')

            self.p0_00 = np.zeros(self.normalized_ssro_00.shape)
            self.u_p0_00 = np.zeros(self.normalized_ssro_00.shape)
            self.p0_01 = np.zeros(self.normalized_ssro_01.shape)
            self.u_p0_01 = np.zeros(self.normalized_ssro_01.shape)
            self.p0_10 = np.zeros(self.normalized_ssro_10.shape)
            self.u_p0_10 = np.zeros(self.normalized_ssro_10.shape)
            self.p0_11 = np.zeros(self.normalized_ssro_11.shape)
            self.u_p0_11 = np.zeros(self.normalized_ssro_11.shape)
            
            ro_durations = self.g.attrs['E_RO_durations']

            roc = error.SingleQubitROC()

            for i in range(len(self.normalized_ssro_00[0])):
                roc.F0, roc.u_F0, roc.F1, roc.u_F1 = \
                    ssro.get_SSRO_calibration(ssro_calib_folder,
                            ro_durations[i])
                p0_00, u_p0_00 = roc.num_eval(self.normalized_ssro_00[:,i],
                        self.u_normalized_ssro_00[:,i])
                p0_01, u_p0_01 = roc.num_eval(self.normalized_ssro_01[:,i],
                        self.u_normalized_ssro_01[:,i])
                p0_10, u_p0_10 = roc.num_eval(self.normalized_ssro_10[:,i],
                        self.u_normalized_ssro_10[:,i])
                p0_11, u_p0_11 = roc.num_eval(self.normalized_ssro_11[:,i],
                        self.u_normalized_ssro_11[:,i])

                self.p0_00[:,i] = p0_00
                self.u_p0_00[:,i] = u_p0_00
                self.p0_01[:,i] = p0_01
                self.u_p0_01[:,i] = u_p0_01
                self.p0_10[:,i] = p0_10
                self.u_p0_10[:,i] = u_p0_10
                self.p0_11[:,i] = p0_11
                self.u_p0_11[:,i] = u_p0_11


        elif self.post_select == True:
            if ssro_calib_folder == '':
                ssro_calib_folder = toolbox.latest_data('SSRO')

            self.p0_0 = np.zeros(self.normalized_ssro_0.shape)
            self.u_p0_0 = np.zeros(self.normalized_ssro_0.shape)
            self.p0_1 = np.zeros(self.normalized_ssro_1.shape)
            self.u_p0_1 = np.zeros(self.normalized_ssro_1.shape)
            ro_durations = self.g.attrs['E_RO_durations']

            roc = error.SingleQubitROC()

            for i in range(len(self.normalized_ssro_0[0])):
                roc.F0, roc.u_F0, roc.F1, roc.u_F1 = \
                    ssro.get_SSRO_calibration(ssro_calib_folder,
                            ro_durations[i])
                p0_0, u_p0_0 = roc.num_eval(self.normalized_ssro_0[:,i],
                        self.u_normalized_ssro_0[:,i])
                p0_1, u_p0_1 = roc.num_eval(self.normalized_ssro_1[:,i],
                        self.u_normalized_ssro_1[:,i])

                self.p0_0[:,i] = p0_0
                self.u_p0_0[:,i] = u_p0_0
                self.p0_1[:,i] = p0_1
                self.u_p0_1[:,i] = u_p0_1

            self.result_corrected = True
        


            self.result_corrected = True
        else:
            mbi.MBIAnalysis.get_electron_ROC(self,ssro_calib_folder)  #NOTE: super cannot be used as this is an "old style class"

    def convert_fidelity_to_contrast(self,p0,u_p0):
        '''
        takes data and corresponding uncertainty and converts it to contrast
        '''
        c0= ((p0.reshape(-1)[:])-0.5)*2
        u_c0 = 2*u_p0.reshape(-1)[:]
        # y_a= ((self.p0.reshape(-1)[:])-0.5)*2
        # y_err_a = 2*a.u_p0.reshape(-1)[:]

        return c0, u_c0


