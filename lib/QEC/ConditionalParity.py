from analysis.lib.m2 import m2
from analysis.lib.m2.ssro import ssro
from analysis.lib.tools import toolbox
from analysis.lib.m2.ssro import mbi
from analysis.lib.math import error
import numpy as np

#NOTE: Function could be moved to analysis.lib.m2.ssro folder when complete and checked for errors.

class ConditionalParityAnalysis(mbi.MBIAnalysis):
    '''
    Class used to extract and analyze C13 experiments with conditional outcomes.
    Child-class of MBI analysis (as all C13 experiments are based on the MBI architecture).
    '''

    def get_readout_results(self, name='',post_select = False,post_select_QEC = False,
                                post_select_multiple_rounds = False, post_select_GHZ = False,orientation_correct=False):
        '''
        Get the readout results.
        self.ssro_results contains the readout results (sum of the photons for
            a given sweep point and readout in a sequence)
        self.normalized_ssro contains the normalized result (i.e., probability
            for getting a photon)
        
        if post_select                  is true, select on the true or false condition (single parity measurement)
        if post_select_QEC              is true, select on the syndrome outcome (two parity measurements)
        if post_select_multiple_rounds  is true, select on the syndrome outcome (4 parity measurements)
        if post_select_GHZ              is true, select on the first 3 measurement rounds
        else                            use the original get_readout_results from the MBI class
 
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

            if orientation_correct:
                orientation_a = adwingrp.attrs['Parity_A_RO_orientation']
                orientation_b = adwingrp.attrs['Tomo_RO_orientation']
                self.orientations = (orientation_a,orientation_b)

                #take into account the orientations of the first measurement for the postselection
                if orientation_a == 'negative':
                    self.parity_result = (1-self.parity_result)

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

            #Probabilities ms=0 and ms=-1
            self.p0 = np.average(self.parity_result/self.reps.astype('float'))
            self.p1 = np.average((self.reps-self.parity_result)/self.reps.astype('float'))

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

        elif post_select_multiple_rounds == True: ### Work in progress, THT
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

            parity_a_result = self.parity_result[0::4]  ### The four parity outcomes are stored sequentially in an array 
            parity_b_result = self.parity_result[1::4] 
            parity_c_result = self.parity_result[2::4]  
            parity_d_result = self.parity_result[3::4] 

            # print 'here they come'
            # print parity_a_result
            # print parity_b_result
            # print parity_c_result
            # print parity_d_result

            ### Step 1 Multiply results with post selection parameter
            ssro_results_0000 = parity_a_result     * parity_b_result   * parity_c_result * parity_d_result * self.ssro_results
            ssro_results_0100 = parity_a_result     *(1-parity_b_result)* parity_c_result * parity_d_result * self.ssro_results
            ssro_results_1000 = (1-parity_a_result) *parity_b_result    * parity_c_result * parity_d_result * self.ssro_results
            ssro_results_1100 = (1-parity_a_result) *(1-parity_b_result)* parity_c_result * parity_d_result * self.ssro_results

            ssro_results_0010 = parity_a_result     * parity_b_result   * (1-parity_c_result) * parity_d_result * self.ssro_results
            ssro_results_0110 = parity_a_result     *(1-parity_b_result)* (1-parity_c_result) * parity_d_result * self.ssro_results
            ssro_results_1010 = (1-parity_a_result) *parity_b_result    * (1-parity_c_result) * parity_d_result * self.ssro_results
            ssro_results_1110 = (1-parity_a_result) *(1-parity_b_result)* (1-parity_c_result) * parity_d_result * self.ssro_results

            ssro_results_0001 = parity_a_result     * parity_b_result   * parity_c_result * (1-parity_d_result) * self.ssro_results
            ssro_results_0101 = parity_a_result     *(1-parity_b_result)* parity_c_result * (1-parity_d_result) * self.ssro_results
            ssro_results_1001 = (1-parity_a_result) *parity_b_result    * parity_c_result * (1-parity_d_result) * self.ssro_results
            ssro_results_1101 = (1-parity_a_result) *(1-parity_b_result)* parity_c_result * (1-parity_d_result) * self.ssro_results

            ssro_results_0011 = parity_a_result     * parity_b_result   * (1-parity_c_result) * (1-parity_d_result) * self.ssro_results
            ssro_results_0111 = parity_a_result     *(1-parity_b_result)* (1-parity_c_result) * (1-parity_d_result) * self.ssro_results
            ssro_results_1011 = (1-parity_a_result) *parity_b_result    * (1-parity_c_result) * (1-parity_d_result) * self.ssro_results
            ssro_results_1111 = (1-parity_a_result) *(1-parity_b_result)* (1-parity_c_result) * (1-parity_d_result) * self.ssro_results

            # print'((1-parity_a_result)*parity_b_result)'
            # print ((1-parity_a_result)*parity_b_result).reshape((-1,self.pts,self.readouts)).sum(axis=0)

            ### Step 2 reshape (integrate the total of occurences for each results)
            parity_result_0000 = (parity_a_result*parity_b_result         * parity_c_result * parity_d_result  ).reshape((-1,self.pts,self.readouts)).sum(axis=0)
            parity_result_0100 = (parity_a_result*(1-parity_b_result)     * parity_c_result * parity_d_result  ).reshape((-1,self.pts,self.readouts)).sum(axis=0)
            parity_result_1000 = ((1-parity_a_result)*parity_b_result     * parity_c_result * parity_d_result  ).reshape((-1,self.pts,self.readouts)).sum(axis=0)
            parity_result_1100 = ((1-parity_a_result)*(1-parity_b_result) * parity_c_result * parity_d_result  ).reshape((-1,self.pts,self.readouts)).sum(axis=0)

            parity_result_0010 = (parity_a_result*parity_b_result         * (1-parity_c_result) * parity_d_result  ).reshape((-1,self.pts,self.readouts)).sum(axis=0)
            parity_result_0110 = (parity_a_result*(1-parity_b_result)     * (1-parity_c_result) * parity_d_result  ).reshape((-1,self.pts,self.readouts)).sum(axis=0)
            parity_result_1010 = ((1-parity_a_result)*parity_b_result     * (1-parity_c_result) * parity_d_result  ).reshape((-1,self.pts,self.readouts)).sum(axis=0)
            parity_result_1110 = ((1-parity_a_result)*(1-parity_b_result) * (1-parity_c_result) * parity_d_result  ).reshape((-1,self.pts,self.readouts)).sum(axis=0)

            parity_result_0001 = (parity_a_result*parity_b_result         * parity_c_result * (1-parity_d_result)  ).reshape((-1,self.pts,self.readouts)).sum(axis=0)
            parity_result_0101 = (parity_a_result*(1-parity_b_result)     * parity_c_result * (1-parity_d_result)  ).reshape((-1,self.pts,self.readouts)).sum(axis=0)
            parity_result_1001 = ((1-parity_a_result)*parity_b_result     * parity_c_result * (1-parity_d_result)  ).reshape((-1,self.pts,self.readouts)).sum(axis=0)
            parity_result_1101 = ((1-parity_a_result)*(1-parity_b_result) * parity_c_result * (1-parity_d_result)  ).reshape((-1,self.pts,self.readouts)).sum(axis=0)

            parity_result_0011 = (parity_a_result*parity_b_result         * (1-parity_c_result) * (1-parity_d_result)  ).reshape((-1,self.pts,self.readouts)).sum(axis=0)
            parity_result_0111 = (parity_a_result*(1-parity_b_result)     * (1-parity_c_result) * (1-parity_d_result)  ).reshape((-1,self.pts,self.readouts)).sum(axis=0)
            parity_result_1011 = ((1-parity_a_result)*parity_b_result     * (1-parity_c_result) * (1-parity_d_result)  ).reshape((-1,self.pts,self.readouts)).sum(axis=0)
            parity_result_1111 = ((1-parity_a_result)*(1-parity_b_result) * (1-parity_c_result) * (1-parity_d_result)  ).reshape((-1,self.pts,self.readouts)).sum(axis=0)


            self.ssro_results_0000 = ssro_results_0000.reshape((-1,self.pts,self.readouts)).sum(axis=0)
            self.ssro_results_0100 = ssro_results_0100.reshape((-1,self.pts,self.readouts)).sum(axis=0)
            self.ssro_results_1000 = ssro_results_1000.reshape((-1,self.pts,self.readouts)).sum(axis=0)
            self.ssro_results_1100 = ssro_results_1100.reshape((-1,self.pts,self.readouts)).sum(axis=0)

            self.ssro_results_0010 = ssro_results_0010.reshape((-1,self.pts,self.readouts)).sum(axis=0)
            self.ssro_results_0110 = ssro_results_0110.reshape((-1,self.pts,self.readouts)).sum(axis=0)
            self.ssro_results_1010 = ssro_results_1010.reshape((-1,self.pts,self.readouts)).sum(axis=0)
            self.ssro_results_1110 = ssro_results_1110.reshape((-1,self.pts,self.readouts)).sum(axis=0)

            self.ssro_results_0001 = ssro_results_0001.reshape((-1,self.pts,self.readouts)).sum(axis=0)
            self.ssro_results_0101 = ssro_results_0101.reshape((-1,self.pts,self.readouts)).sum(axis=0)
            self.ssro_results_1001 = ssro_results_1001.reshape((-1,self.pts,self.readouts)).sum(axis=0)
            self.ssro_results_1101 = ssro_results_1101.reshape((-1,self.pts,self.readouts)).sum(axis=0)

            self.ssro_results_0011 = ssro_results_0011.reshape((-1,self.pts,self.readouts)).sum(axis=0)
            self.ssro_results_0111 = ssro_results_0111.reshape((-1,self.pts,self.readouts)).sum(axis=0)
            self.ssro_results_1011 = ssro_results_1011.reshape((-1,self.pts,self.readouts)).sum(axis=0)
            self.ssro_results_1111 = ssro_results_1111.reshape((-1,self.pts,self.readouts)).sum(axis=0)

            ### Step 3 normalization and uncertainty
            self.normalized_ssro_0000     = self.ssro_results_0000/(parity_result_0000).astype('float')
            self.u_normalized_ssro_0000   = (self.normalized_ssro_0000*(1-self.normalized_ssro_0000)/(parity_result_0000))**0.5
            
            self.normalized_ssro_0100     = self.ssro_results_0100/(parity_result_0100).astype('float')
            self.u_normalized_ssro_0100   = (self.normalized_ssro_0100*(1-self.normalized_ssro_0100)/(parity_result_0100))**0.5
            
            self.normalized_ssro_1000     = self.ssro_results_1000/(parity_result_1000 ).astype('float')
            self.u_normalized_ssro_1000   = (self.normalized_ssro_1000*(1-self.normalized_ssro_1000)/(parity_result_1000))**0.5
            
            self.normalized_ssro_1100     = self.ssro_results_1100/(parity_result_1100).astype('float')
            self.u_normalized_ssro_1100   = (self.normalized_ssro_1100*(1-self.normalized_ssro_1100)/(parity_result_1100))**0.5
            

            self.normalized_ssro_0010     = self.ssro_results_0010/(parity_result_0010).astype('float')
            self.u_normalized_ssro_0010   = (self.normalized_ssro_0010*(1-self.normalized_ssro_0010)/(parity_result_0010))**0.5
            
            self.normalized_ssro_0110     = self.ssro_results_0110/(parity_result_0110).astype('float')
            self.u_normalized_ssro_0110   = (self.normalized_ssro_0110*(1-self.normalized_ssro_0110)/(parity_result_0110))**0.5
            
            self.normalized_ssro_1010     = self.ssro_results_1010/(parity_result_1010 ).astype('float')
            self.u_normalized_ssro_1010   = (self.normalized_ssro_1010*(1-self.normalized_ssro_1010)/(parity_result_1010))**0.5
            
            self.normalized_ssro_1110     = self.ssro_results_1110/(parity_result_1110).astype('float')
            self.u_normalized_ssro_1110   = (self.normalized_ssro_1110*(1-self.normalized_ssro_1110)/(parity_result_1110))**0.5
            

            self.normalized_ssro_0001     = self.ssro_results_0001/(parity_result_0001).astype('float')
            self.u_normalized_ssro_0001   = (self.normalized_ssro_0001*(1-self.normalized_ssro_0001)/(parity_result_0001))**0.5
            
            self.normalized_ssro_0101     = self.ssro_results_0101/(parity_result_0101).astype('float')
            self.u_normalized_ssro_0101   = (self.normalized_ssro_0101*(1-self.normalized_ssro_0101)/(parity_result_0101))**0.5
            
            self.normalized_ssro_1001     = self.ssro_results_1001/(parity_result_1001 ).astype('float')
            self.u_normalized_ssro_1001   = (self.normalized_ssro_1001*(1-self.normalized_ssro_1001)/(parity_result_1001))**0.5
            
            self.normalized_ssro_1101     = self.ssro_results_1101/(parity_result_1101).astype('float')
            self.u_normalized_ssro_1101   = (self.normalized_ssro_1101*(1-self.normalized_ssro_1101)/(parity_result_1101))**0.5
            

            self.normalized_ssro_0011     = self.ssro_results_0011/(parity_result_0011).astype('float')
            self.u_normalized_ssro_0011   = (self.normalized_ssro_0011*(1-self.normalized_ssro_0011)/(parity_result_0011))**0.5
            
            self.normalized_ssro_0111     = self.ssro_results_0111/(parity_result_0111).astype('float')
            self.u_normalized_ssro_0111   = (self.normalized_ssro_0111*(1-self.normalized_ssro_0111)/(parity_result_0111))**0.5
            
            self.normalized_ssro_1011     = self.ssro_results_1011/(parity_result_1011 ).astype('float')
            self.u_normalized_ssro_1011   = (self.normalized_ssro_1011*(1-self.normalized_ssro_1011)/(parity_result_1011))**0.5
            
            self.normalized_ssro_1111     = self.ssro_results_1111/(parity_result_1111).astype('float')
            self.u_normalized_ssro_1111   = (self.normalized_ssro_1111*(1-self.normalized_ssro_1111)/(parity_result_1111))**0.5
            
            ### 'Probabilities for different outcomes'
            
            self.p0000 = (parity_result_0000/float(len(self.ssro_results)/self.pts))
            self.p0100 = ((parity_result_0100)/float(len(self.ssro_results)/self.pts))
            self.p1000 = (parity_result_1000/float(len(self.ssro_results)/self.pts))
            self.p1100 = ((parity_result_1100)/float(len(self.ssro_results)/self.pts))

            self.p0010 = ( parity_result_0010/float(len(self.ssro_results)/self.pts))
            self.p0110 = ((parity_result_0110)/float(len(self.ssro_results)/self.pts))
            self.p1010 = ( parity_result_1010/float(len(self.ssro_results)/self.pts))
            self.p1110 = ((parity_result_1110)/float(len(self.ssro_results)/self.pts))

            self.p0001 = ( parity_result_0001/float(len(self.ssro_results)/self.pts)) 
            self.p0101 = ((parity_result_0101)/float(len(self.ssro_results)/self.pts))
            self.p1001 = ( parity_result_1001/float(len(self.ssro_results)/self.pts))
            self.p1101 = ((parity_result_1101)/float(len(self.ssro_results)/self.pts))

            self.p0011 = ( parity_result_0011/float(len(self.ssro_results)/self.pts)) 
            self.p0111 = ((parity_result_0111)/float(len(self.ssro_results)/self.pts))
            self.p1011 = ( parity_result_1011/float(len(self.ssro_results)/self.pts))
            self.p1111 = ((parity_result_1111)/float(len(self.ssro_results)/self.pts))
                       
        elif post_select_GHZ == True:
            self.post_select = True
            self.result_corrected = False

            adwingrp = self.adwingrp(name)
            self.adgrp = adwingrp

            self.pts        = adwingrp.attrs['sweep_length']
            self.reps       = adwingrp.attrs['reps_per_ROsequence']
            self.readouts   = adwingrp.attrs['nr_of_ROsequences']

            orientation_a = adwingrp.attrs['Parity_xyy_RO_orientation']
            orientation_b = adwingrp.attrs['Parity_yxy_RO_orientation']
            orientation_c = adwingrp.attrs['Parity_yyx_RO_orientation']
            orientation_d = adwingrp.attrs['Tomo_RO_orientation']
            self.orientations = (orientation_a,orientation_b,orientation_c,orientation_d)

            ### Step 0 extract data from hdf5 file
            self.parity_result = adwingrp['parity_RO_results'].value
            self.ssro_results = adwingrp['ssro_results'].value 

            #######take into account the readout orientations???
            parity_a_result = self.parity_result[0::3]  ### The three parity outcomes are stored sequentially in an array 
            parity_b_result = self.parity_result[1::3] 
            parity_c_result = self.parity_result[2::3]  
            parity_d_result = self.ssro_results 

            #take into account the orientations of the first three measurements for the postselection
            if orientation_a == 'negative':
                parity_a_result = (1-parity_a_result)
            if orientation_b == 'negative':
                parity_b_result = (1-parity_b_result)
            if orientation_c == 'negative':
                parity_c_result = (1-parity_c_result)
                        # if orientation_d == 'negative':
            #     parity_d_result = (1-parity_d_result)

            ### number of times the first 3 parity mmt outcomes occurred
            parity_result_000 = (parity_a_result*     parity_b_result*     parity_c_result     ).reshape((-1,self.pts,self.readouts)).sum(axis=0)
            parity_result_001 = (parity_a_result*     parity_b_result*     (1-parity_c_result) ).reshape((-1,self.pts,self.readouts)).sum(axis=0)
            parity_result_010 = (parity_a_result*     (1-parity_b_result)* parity_c_result     ).reshape((-1,self.pts,self.readouts)).sum(axis=0)
            parity_result_011 = (parity_a_result*     (1-parity_b_result)* (1-parity_c_result) ).reshape((-1,self.pts,self.readouts)).sum(axis=0)
            parity_result_100 = ((1-parity_a_result)* parity_b_result*     parity_c_result     ).reshape((-1,self.pts,self.readouts)).sum(axis=0)
            parity_result_101 = ((1-parity_a_result)* parity_b_result*     (1-parity_c_result) ).reshape((-1,self.pts,self.readouts)).sum(axis=0)
            parity_result_110 = ((1-parity_a_result)* (1-parity_b_result)* parity_c_result     ).reshape((-1,self.pts,self.readouts)).sum(axis=0)
            parity_result_111 = ((1-parity_a_result)* (1-parity_b_result)* (1-parity_c_result) ).reshape((-1,self.pts,self.readouts)).sum(axis=0)

            ssro_result_000 = (parity_a_result     *parity_b_result     *parity_c_result     *parity_d_result).reshape((-1,self.pts,self.readouts)).sum(axis=0)
            ssro_result_001 = (parity_a_result     *parity_b_result     *(1-parity_c_result) *parity_d_result).reshape((-1,self.pts,self.readouts)).sum(axis=0)
            ssro_result_010 = (parity_a_result     *(1-parity_b_result) *parity_c_result     *parity_d_result).reshape((-1,self.pts,self.readouts)).sum(axis=0)
            ssro_result_011 = (parity_a_result     *(1-parity_b_result) *(1-parity_c_result) *parity_d_result).reshape((-1,self.pts,self.readouts)).sum(axis=0)
            ssro_result_100 = ((1-parity_a_result) *parity_b_result     *parity_c_result     *parity_d_result).reshape((-1,self.pts,self.readouts)).sum(axis=0)
            ssro_result_101 = ((1-parity_a_result) *parity_b_result     *(1-parity_c_result) *parity_d_result).reshape((-1,self.pts,self.readouts)).sum(axis=0)
            ssro_result_110 = ((1-parity_a_result) *(1-parity_b_result) *parity_c_result     *parity_d_result).reshape((-1,self.pts,self.readouts)).sum(axis=0)
            ssro_result_111 = ((1-parity_a_result) *(1-parity_b_result) *(1-parity_c_result) *parity_d_result).reshape((-1,self.pts,self.readouts)).sum(axis=0)

            if parity_result_000 != 0:
                self.normalized_ssro_000 = ssro_result_000/(parity_result_000).astype('float')
                self.u_normalized_ssro_000 = ((self.normalized_ssro_000*(1-self.normalized_ssro_000))/(parity_result_000))**(0.5)
            else:
                self.normalized_ssro_000 = np.array([[0.]])
                self.u_normalized_ssro_000 = np.array([[1.]])
            if parity_result_001 != 0:
                self.normalized_ssro_001 = ssro_result_001/(parity_result_001).astype('float')
                self.u_normalized_ssro_001 = ((self.normalized_ssro_001*(1-self.normalized_ssro_001))/(parity_result_001))**(0.5)
            else:
                self.normalized_ssro_001 = np.array([[0.]])
                self.u_normalized_ssro_001 = np.array([[1.]])
            if parity_result_010 != 0:
                self.normalized_ssro_010 = ssro_result_010/(parity_result_010).astype('float')
                self.u_normalized_ssro_010 = ((self.normalized_ssro_010*(1-self.normalized_ssro_010))/(parity_result_010))**(0.5)
            else:
                self.normalized_ssro_010 = np.array([[0.]])
                self.u_normalized_ssro_010 = np.array([[1.]])
            if parity_result_011 != 0:
                self.normalized_ssro_011 = ssro_result_011/(parity_result_011).astype('float')
                self.u_normalized_ssro_011 = ((self.normalized_ssro_011*(1-self.normalized_ssro_011))/(parity_result_011))**(0.5)
            else:
                self.normalized_ssro_011 = np.array([[0.]])
                self.u_normalized_ssro_011 = np.array([[1.]])
            if parity_result_100 != 0:
                self.normalized_ssro_100 = ssro_result_100/(parity_result_100).astype('float')
                self.u_normalized_ssro_100 = ((self.normalized_ssro_100*(1-self.normalized_ssro_100))/(parity_result_100))**(0.5)
            else:
                self.normalized_ssro_100 = np.array([[0.]])
                self.u_normalized_ssro_100 = np.array([[1.]])
            if parity_result_101 != 0:
                self.normalized_ssro_101 = ssro_result_101/(parity_result_101).astype('float')
                self.u_normalized_ssro_101 = ((self.normalized_ssro_101*(1-self.normalized_ssro_101))/(parity_result_101))**(0.5)
            else:
                self.normalized_ssro_101 = np.array([[0.]])
                self.u_normalized_ssro_101 = np.array([[1.]])
            if parity_result_110 != 0:
                self.normalized_ssro_110 = ssro_result_110/(parity_result_110).astype('float')
                self.u_normalized_ssro_110 = ((self.normalized_ssro_110*(1-self.normalized_ssro_110))/(parity_result_110))**(0.5)
            else:
                self.normalized_ssro_110 = np.array([[0.]])
                self.u_normalized_ssro_110 = np.array([[1.]])
            if parity_result_111 != 0:
                self.normalized_ssro_111 = ssro_result_111/(parity_result_111).astype('float')
                self.u_normalized_ssro_111 = ((self.normalized_ssro_111*(1-self.normalized_ssro_111))/(parity_result_111))**(0.5)
            else:
                self.normalized_ssro_111 = np.array([[0.]])
                self.u_normalized_ssro_111 = np.array([[1.]])
 
            self.p000 = (parity_result_000/float(len(self.ssro_results)/self.pts))
            self.p001 = (parity_result_001/float(len(self.ssro_results)/self.pts))
            self.p010 = (parity_result_010/float(len(self.ssro_results)/self.pts))
            self.p011 = (parity_result_011/float(len(self.ssro_results)/self.pts))
            self.p100 = (parity_result_100/float(len(self.ssro_results)/self.pts))
            self.p101 = (parity_result_101/float(len(self.ssro_results)/self.pts))
            self.p110 = (parity_result_110/float(len(self.ssro_results)/self.pts))
            self.p111 = (parity_result_111/float(len(self.ssro_results)/self.pts))

        else:
            mbi.MBIAnalysis.get_readout_results(self,name) #NOTE: super cannot be used as this is an "old style class"
            self.post_select = False

    def get_electron_ROC(self, ssro_calib_folder='',post_select_QEC = False, post_select_multiple_rounds = False, post_select_GHZ = False):
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

        elif post_select_multiple_rounds == True and self.post_select == True:
              
            if ssro_calib_folder == '':
                ssro_calib_folder = toolbox.latest_data('SSRO')

            self.p0_0000   = np.zeros(self.normalized_ssro_0000.shape)
            self.u_p0_0000 =  np.zeros(self.normalized_ssro_0000.shape)
            self.p0_0100   = np.zeros(self.normalized_ssro_0100.shape)
            self.u_p0_0100 =  np.zeros(self.normalized_ssro_0100.shape)
            self.p0_1000   = np.zeros(self.normalized_ssro_1000.shape)
            self.u_p0_1000 =  np.zeros(self.normalized_ssro_1000.shape)
            self.p0_1100   = np.zeros(self.normalized_ssro_1100.shape)
            self.u_p0_1100 =  np.zeros(self.normalized_ssro_1100.shape)

            self.p0_0010   = np.zeros(self.normalized_ssro_0010.shape)
            self.u_p0_0010 =  np.zeros(self.normalized_ssro_0010.shape)
            self.p0_0110   = np.zeros(self.normalized_ssro_0110.shape)
            self.u_p0_0110 =  np.zeros(self.normalized_ssro_0110.shape)
            self.p0_1010   = np.zeros(self.normalized_ssro_1010.shape)
            self.u_p0_1010 =  np.zeros(self.normalized_ssro_1010.shape)
            self.p0_1110   = np.zeros(self.normalized_ssro_1110.shape)
            self.u_p0_1110 =  np.zeros(self.normalized_ssro_1110.shape)

            self.p0_0001   = np.zeros(self.normalized_ssro_0001.shape)
            self.u_p0_0001 =  np.zeros(self.normalized_ssro_0001.shape)
            self.p0_0101   = np.zeros(self.normalized_ssro_0101.shape)
            self.u_p0_0101 =  np.zeros(self.normalized_ssro_0101.shape)
            self.p0_1001   = np.zeros(self.normalized_ssro_1001.shape)
            self.u_p0_1001 =  np.zeros(self.normalized_ssro_1001.shape)
            self.p0_1101   = np.zeros(self.normalized_ssro_1101.shape)
            self.u_p0_1101 =  np.zeros(self.normalized_ssro_1101.shape)

            self.p0_0011   = np.zeros(self.normalized_ssro_0011.shape)
            self.u_p0_0011 =  np.zeros(self.normalized_ssro_0011.shape)
            self.p0_0111   = np.zeros(self.normalized_ssro_0111.shape)
            self.u_p0_0111 =  np.zeros(self.normalized_ssro_0111.shape)
            self.p0_1011   = np.zeros(self.normalized_ssro_1011.shape)
            self.u_p0_1011 =  np.zeros(self.normalized_ssro_1011.shape)
            self.p0_1111   = np.zeros(self.normalized_ssro_1111.shape)
            self.u_p0_1111 =  np.zeros(self.normalized_ssro_1111.shape)
            
            ro_durations = self.g.attrs['E_RO_durations']

            roc = error.SingleQubitROC()

            for i in range(len(self.normalized_ssro_0000[0])):
                roc.F0, roc.u_F0, roc.F1, roc.u_F1 = \
                    ssro.get_SSRO_calibration(ssro_calib_folder,
                            ro_durations[i])
                p0_0000, u_p0_0000 = roc.num_eval(self.normalized_ssro_0000[:,i],
                        self.u_normalized_ssro_0000[:,i])
                p0_0100, u_p0_0100 = roc.num_eval(self.normalized_ssro_0100[:,i],
                        self.u_normalized_ssro_0100[:,i])
                p0_1000, u_p0_1000 = roc.num_eval(self.normalized_ssro_1000[:,i],
                        self.u_normalized_ssro_1000[:,i])
                p0_1100, u_p0_1100 = roc.num_eval(self.normalized_ssro_1100[:,i],
                        self.u_normalized_ssro_1100[:,i])

                p0_0010, u_p0_0010 = roc.num_eval(self.normalized_ssro_0010[:,i],
                        self.u_normalized_ssro_0010[:,i])
                p0_0110, u_p0_0110 = roc.num_eval(self.normalized_ssro_0110[:,i],
                        self.u_normalized_ssro_0110[:,i])
                p0_1010, u_p0_1010 = roc.num_eval(self.normalized_ssro_1010[:,i],
                        self.u_normalized_ssro_1010[:,i])
                p0_1110, u_p0_1110 = roc.num_eval(self.normalized_ssro_1110[:,i],
                        self.u_normalized_ssro_1110[:,i])

                p0_0001, u_p0_0001 = roc.num_eval(self.normalized_ssro_0001[:,i],
                        self.u_normalized_ssro_0001[:,i])
                p0_0101, u_p0_0101 = roc.num_eval(self.normalized_ssro_0101[:,i],
                        self.u_normalized_ssro_0101[:,i])
                p0_1001, u_p0_1001 = roc.num_eval(self.normalized_ssro_1001[:,i],
                        self.u_normalized_ssro_1001[:,i])
                p0_1101, u_p0_1101 = roc.num_eval(self.normalized_ssro_1101[:,i],
                        self.u_normalized_ssro_1101[:,i])

                p0_0011, u_p0_0011 = roc.num_eval(self.normalized_ssro_0011[:,i],
                        self.u_normalized_ssro_0011[:,i])
                p0_0111, u_p0_0111 = roc.num_eval(self.normalized_ssro_0111[:,i],
                        self.u_normalized_ssro_0111[:,i])
                p0_1011, u_p0_1011 = roc.num_eval(self.normalized_ssro_1011[:,i],
                        self.u_normalized_ssro_1011[:,i])
                p0_1111, u_p0_1111 = roc.num_eval(self.normalized_ssro_1111[:,i],
                        self.u_normalized_ssro_1111[:,i])

                self.p0_0000[:,i]   = p0_0000
                self.u_p0_0000[:,i] = u_p0_0000
                self.p0_0100[:,i]   = p0_0100
                self.u_p0_0100[:,i] = u_p0_0100
                self.p0_1000[:,i]   = p0_1000
                self.u_p0_1000[:,i] = u_p0_1000
                self.p0_1100[:,i]   = p0_1100
                self.u_p0_1100[:,i] = u_p0_1100

                self.p0_0010[:,i]   = p0_0010
                self.u_p0_0010[:,i] = u_p0_0010
                self.p0_0110[:,i]   = p0_0110
                self.u_p0_0110[:,i] = u_p0_0110
                self.p0_1010[:,i]   = p0_1010
                self.u_p0_1010[:,i] = u_p0_1010
                self.p0_1110[:,i]   = p0_1110
                self.u_p0_1110[:,i] = u_p0_1110

                self.p0_0001[:,i]   = p0_0001
                self.u_p0_0001[:,i] = u_p0_0001
                self.p0_0101[:,i]   = p0_0101
                self.u_p0_0101[:,i] = u_p0_0101
                self.p0_1001[:,i]   = p0_1001
                self.u_p0_1001[:,i] = u_p0_1001
                self.p0_1101[:,i]   = p0_1101
                self.u_p0_1101[:,i] = u_p0_1101

                self.p0_0011[:,i]   = p0_0011
                self.u_p0_0011[:,i] = u_p0_0011
                self.p0_0111[:,i]   = p0_0111
                self.u_p0_0111[:,i] = u_p0_0111
                self.p0_1011[:,i]   = p0_1011
                self.u_p0_1011[:,i] = u_p0_1011
                self.p0_1111[:,i]   = p0_1111
                self.u_p0_1111[:,i] = u_p0_1111

        elif post_select_GHZ == True and self.post_select==True:
            if ssro_calib_folder == '':
                ssro_calib_folder = toolbox.latest_data('SSRO')
            
            self.p0_000 = np.zeros(self.normalized_ssro_000.shape)
            self.u_p0_000 = np.zeros(self.normalized_ssro_000.shape)
            self.p0_001 = np.zeros(self.normalized_ssro_001.shape)
            self.u_p0_001 = np.zeros(self.normalized_ssro_001.shape)
            self.p0_010 = np.zeros(self.normalized_ssro_010.shape)
            self.u_p0_010 = np.zeros(self.normalized_ssro_010.shape)
            self.p0_011 = np.zeros(self.normalized_ssro_011.shape)
            self.u_p0_011 = np.zeros(self.normalized_ssro_011.shape)
            self.p0_100 = np.zeros(self.normalized_ssro_100.shape)
            self.u_p0_100 = np.zeros(self.normalized_ssro_100.shape)
            self.p0_101 = np.zeros(self.normalized_ssro_101.shape)
            self.u_p0_101 = np.zeros(self.normalized_ssro_101.shape)
            self.p0_110 = np.zeros(self.normalized_ssro_110.shape)
            self.u_p0_110 = np.zeros(self.normalized_ssro_110.shape)
            self.p0_111 = np.zeros(self.normalized_ssro_111.shape)
            self.u_p0_111 = np.zeros(self.normalized_ssro_111.shape)
            
            ro_durations = self.g.attrs['E_RO_durations']

            roc = error.SingleQubitROC()

            for i in range(len(self.normalized_ssro_000[0])):
                roc.F0, roc.u_F0, roc.F1, roc.u_F1 = \
                    ssro.get_SSRO_calibration(ssro_calib_folder,
                            ro_durations[i])
                p0_000, u_p0_000 = roc.num_eval(self.normalized_ssro_000[:,i],
                        self.u_normalized_ssro_000[:,i])
                p0_001, u_p0_001 = roc.num_eval(self.normalized_ssro_001[:,i],
                        self.u_normalized_ssro_001[:,i])
                p0_010, u_p0_010 = roc.num_eval(self.normalized_ssro_010[:,i],
                        self.u_normalized_ssro_010[:,i])
                p0_011, u_p0_011 = roc.num_eval(self.normalized_ssro_011[:,i],
                        self.u_normalized_ssro_011[:,i])
                p0_100, u_p0_100 = roc.num_eval(self.normalized_ssro_100[:,i],
                        self.u_normalized_ssro_100[:,i])
                p0_101, u_p0_101 = roc.num_eval(self.normalized_ssro_101[:,i],
                        self.u_normalized_ssro_101[:,i])
                p0_110, u_p0_110 = roc.num_eval(self.normalized_ssro_110[:,i],
                        self.u_normalized_ssro_110[:,i])
                p0_111, u_p0_111 = roc.num_eval(self.normalized_ssro_111[:,i],
                        self.u_normalized_ssro_111[:,i])

                self.p0_000[:,i]   = p0_000
                self.u_p0_000[:,i]   = u_p0_000
                self.p0_001[:,i]   = p0_001
                self.u_p0_001[:,i]   = u_p0_001
                self.p0_010[:,i]   = p0_010
                self.u_p0_010[:,i]   = u_p0_010
                self.p0_011[:,i]   = p0_011
                self.u_p0_011[:,i]   = u_p0_011
                self.p0_100[:,i]   = p0_100
                self.u_p0_100[:,i]   = u_p0_100
                self.p0_101[:,i]   = p0_101
                self.u_p0_101[:,i]   = u_p0_101
                self.p0_110[:,i]   = p0_110
                self.u_p0_110[:,i]   = u_p0_110
                self.p0_111[:,i]   = p0_111
                self.u_p0_111[:,i]   = u_p0_111

            self.result_corrected = True

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


