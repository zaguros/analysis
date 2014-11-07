from analysis.lib.m2 import m2
from analysis.lib.m2.ssro import ssro
from analysis.lib.math import error
from analysis.lib.tools import toolbox
from analysis.lib.m2.ssro import mbi
import numpy as np

#NOTE: Function could be moved to analysis.lib.m2.ssro folder when complete and checked for errors.

class ConditionalParityAnalysis(mbi.MBIAnalysis):
    '''
    Class used for extracting data from conditional C13 experiments.
    Child-class of MBI analysis (as conditional C13 experiments are based on the MBI architecture).
    '''

    def get_readout_results(self, name='',post_select = False):
        '''
        Get the readout results.
        self.ssro_results contains the readout results (sum of the photons for
            a given sweep point and readout in a sequence)
        self.normalized_ssro contains the normalized result (i.e., probability
            for getting a photon)
        if Post selection is false, use the get_readout results from the MBI class
        if post selection is true, select on the true or false condition
        '''

        if post_select == True:
            self.post_select = True
            self.result_corrected = False

            adwingrp = self.adwingrp(name)
            self.adgrp = adwingrp

            self.pts = adwingrp.attrs['sweep_length']
            self.reps = adwingrp.attrs['reps_per_ROsequence']
            self.readouts= adwingrp.attrs['nr_of_ROsequences']


            #Step 0 extract data from hdf5 file
            self.parity_result = adwingrp['parity_RO_results'].value #Creates a list of 0 and 1 's for when the parity measurement was success
            self.ssro_results = adwingrp['ssro_results'].value #Extracts all the SSRO data

            # #Test cases, condition on self, condition on random successes, only zeros
            # self.parity_result = adwingrp['ssro_results'].value #Creates a list of 0 and 1 's for when the parity measurement was success
            # self.parity_result = np.random.randint(2, size= len(self.ssro_results))
            # self.parity_result = np.random.randint(1, size= len(self.ssro_results))


            #Step 1 Multiply results with post selection parameter
            ssro_results_0 = self.parity_result*self.ssro_results
            ssro_results_1 = (1-self.parity_result)*self.ssro_results

            # Step 2 reshape
            self.parity_result = self.parity_result.reshape((-1,self.pts,self.readouts)).sum(axis=0)
            self.ssro_results_0 = ssro_results_0.reshape((-1,self.pts,self.readouts)).sum(axis=0)
            self.ssro_results_1 = ssro_results_1.reshape((-1,self.pts,self.readouts)).sum(axis=0)

            #Step 3 normalization and uncertainty, different per column
            self.normalized_ssro_0 = self.ssro_results_0/(self.parity_result ).astype('float')
            self.u_normalized_ssro_0 = (self.normalized_ssro_0*(1-self.normalized_ssro_0)/(self.parity_result))**0.5
            self.normalized_ssro_1 = self.ssro_results_1/(self.reps-self.parity_result).astype('float')
            self.u_normalized_ssro_1 = (self.normalized_ssro_1*(1-self.normalized_ssro_1)/(self.reps-self.parity_result))**0.5


            print 'Probabilities ms=0 and ms=-1'
            print np.average(self.parity_result/self.reps.astype('float'))
            print np.average((self.reps-self.parity_result)/self.reps.astype('float'))

        else:
            mbi.MBIAnalysis.get_readout_results(self,name) #NOTE: super cannot be used as this is an "old style class"
            self.post_select = False

    def get_electron_ROC(self, ssro_calib_folder=''):
        '''
        Performs Readout Correction, needs to be updated to correct for post selected results and apply error-bars correctly.
        '''
        if self.post_select == True:
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


