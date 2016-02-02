import numpy as np
import os
import datetime
import h5py
import analysis.lib.tools.toolbox as tb
from analysis.lib.lde import ro_c_err
from analysis.lib.bell import bell_events as be
from analysis.lib.math import error

   #corr defs
rnd_corr=[[0,0],[0,1],[1,0],[1,1]] #'RND [LT3,LT4] rnd 00,01,10,11'
ro_corr =[[1,1],[1,0],[0,1],[0,0]] #'RO [LT3,LT4] ms  00, 01, 10, 11'

def correlator_error(N00,N01,N10,N11):
    n=float(N00+N01+N10+N11)
    p_est=(N00+N11)/n
    return 2*np.sqrt(p_est*(1-p_est)/n) #Wald method 68% confidence interval (z=1) for binomial parameter estimation. Factor 2 comes from E=2p-1
    #this one was bs. #return 2*np.sqrt((N00**3*(N01+N10)+N00**2*(2*N01*N10+3*(N01+N10)*N11)+N11*((N01+N10)**3+2*N01*N10*N11+(N01+N10)*N11**2)+N00*((N01+N10)**3+2*(N01**2+4*N01*N10+N10**2)*N11+3*(N01+N10)*N11**2))/(N00+N01+N10+N11)**5)

def RO_correct(RO_norm, F0A, F1A, F0B, F1B, reverse=False):
    RO_norm =np.asmatrix(RO_norm)
    U = np.matrix([[F0A*F0B,         F0A*(1-F1B),     F0B*(1-F1A),     (1-F1A)*(1-F1B)],
                   [F0A*(1-F0B),     F0A*F1B,         (1-F1A)*(1-F0B), F1B*(1-F1A)],
                   [F0B*(1-F0A),     (1-F0A)*(1-F1B), F1A*F0B,         F1A*(1-F1B)],
                   [(1-F0A)*(1-F0B), F1B*(1-F0A),     F1A*(1-F0B),     F1A*F1B]])
    if reverse:
        Uinv = U
    else:
        Uinv=U.I
    return np.matrix.dot(Uinv,RO_norm.T)

def RO_correct_err(RO_norm, F0A, F1A, F0B, F1B):
    (N00, N01, N10, N11) = RO_norm
    return np.array(ro_c_err.get_readout_correction_errors(F0a=F0A,F1a=F1A, F0b=F0B, F1b=F1B, dF0a = 0.005,dF1a = 0.005,dF0b = 0.005,dF1b = 0.005,N00=N00, N01=N01,N10=N10, N11=N11))

def RO_correct_single_qubit(p0, u_p0, F0, F1, u_F0 , u_F1):
    roc = error.SingleQubitROC()
    roc.F0, roc.u_F0, roc.F1, roc.u_F1 = (F0,u_F0,F1,u_F1)
    return roc.num_eval(p0,u_p0)


def merge_dicts(*dict_args):
    '''
    Given any number of dicts, shallow copy and merge into a new dict,
    precedence goes to key value pairs in latter dicts.
    '''
    result = {}
    for dictionary in dict_args:
        result.update(dictionary)
    return result

def print_correlators(corr_mats, VERBOSE=True, psis=['psi_min', 'psi_plus']):
    # print the results:
    for psi in psis:
        corr_mat=np.zeros((4,4))
        for k in corr_mats:
            if psi in k:
                corr_mat+=corr_mats[k]
        noof_ev_fltr = int(np.sum(corr_mat))
                
        Es=np.zeros(4)
        dEs=np.zeros(4)

        for i,rnd in enumerate(rnd_corr):
            Es[i] = (corr_mat[i,0] - corr_mat[i,1] - corr_mat[i,2] + corr_mat[i,3])/float(np.sum(corr_mat[i,:]))
            dEs[i] = correlator_error(corr_mat[i,0],corr_mat[i,1],corr_mat[i,2],corr_mat[i,3])
            
        if psi == 'psi_min':
            expected_Es= (-0.42,0.42,0.67,0.67)
            CHSH  = -Es[0] + Es[1] + Es[2] + Es[3]
        elif psi == 'psi_plus':
            expected_Es= (0.42,-0.42,0.67,0.67)
            CHSH  = Es[0] - Es[1] + Es[2] + Es[3]
        dCHSH = np.sqrt(dEs[0]**2 + dEs[1]**2 + dEs[2]**2 + dEs[3]**2)
        if VERBOSE:
            print '-'*40
            print 'FILTERED EVENTS {}: Number of events {}'.format(psi,noof_ev_fltr)
            print 'RO ms   00, 01, 10, 11'
            print 'RND00', corr_mat[0], '  +pi/2, +3pi/4'
            print 'RND01', corr_mat[1], '  +pi/2, -3pi/4'
            print 'RND10', corr_mat[2], '  0,     +3pi/4'
            print 'RND11', corr_mat[3], '  0,     -3pi/4\n'

            print ' E (RND00  RND01  RND10  RND11 )'
            print '   ({:+.3f}, {:+.3f}, {:+.3f}, {:+.3f}) expected'.format(*expected_Es)
            print '   ({:+.3f}, {:+.3f}, {:+.3f}, {:+.3f}) measured'.format(*Es)
            print '+/-( {:.3f},  {:.3f},  {:.3f},  {:.3f} )'.format(*dEs)

            print 'CHSH : {:.3f} +- {:.3f}'.format(CHSH, dCHSH)
    return CHSH, dCHSH, Es, dEs

def C_val(a,b,x,y,psi):#expects a,b,psi in {0,1} and x,y in {+1,-1}
        return ((-1)**(a*(b+psi))*x*y+1)/2
    
def get_p_val(K,N):
    from scipy.stats import binom
    tau = 5.4e-6*2
    p_lhv = 3./4+3*(tau+tau**2) #1-(0.5-eps)**2*(1-12*tau*(1+tau))#eps and tau correspond to the bias and the predictability.
    return 1-binom.cdf(K-1, N, p_lhv)

def calculate_p_lhv(corr_mats, VERBOSE=True):
    K = 0
    N = 0
    Kxx = 0
    Nxx = 0
    Kzz = 0
    Nzz = 0
    for i,rnd in enumerate(rnd_corr):
        for j,ro in enumerate(ro_corr):
            for psi in corr_mats:
                C=C_val(rnd[0]+1, rnd[1], 2*ro[0]-1, 2*ro[1]-1, 'psi_min' in psi)
                k=(C==1)*corr_mats[psi][i,j]
                n=corr_mats[psi][i,j]
                K+=k
                N+=n
                if rnd[0] == 0: #LT3 did a pi/2 pulse
                    Kxx+=k
                    Nxx+=n
                elif rnd[0] == 1: #LT3 did no pi/2 pulse
                    Kzz+=k
                    Nzz+=n      
    p_val = get_p_val(K,N)
    if VERBOSE:
        print 'All: {}/{} = {:.2f}'.format(K, N, K/N)
        print 'Probability of LHV model: {:.4f}%'.format(p_val*100)
        print 'XX: {}/{} = {:.2f}'.format(Kxx, Nxx, Kxx/Nxx)
        print 'ZZ: {}/{} = {:.2f}'.format(Kzz, Nzz, Kzz/Nzz)
    return K,N,Kxx,Nxx,Kzz,Nzz,p_val

def get_sp_corrs(db,dlt,db_fps, analysis_params, lt3, VERBOSE=False):
    st_start_ch0  = analysis_params['st_start_ch0']
    st_len        = analysis_params['st_len'] #50 ns
    st_len_w2     = st_len
    p_sep         = analysis_params['pulse_sep'] #600 ns
    st_start_ch1  = analysis_params['st_start_ch1']
            
    st_fltr_w1 =  (((st_start_ch0 <= db[:,be._cl_st_w1]) & (db[:,be._cl_st_w1] < (st_start_ch0 + st_len)) & (db[:,be._cl_ch_w1] == 0)) \
                 | ((st_start_ch1 <= db[:,be._cl_st_w1]) & (db[:,be._cl_st_w1] < (st_start_ch1 + st_len)) & (db[:,be._cl_ch_w1] == 1)) )  
    st_fltr_w2 =  (((st_start_ch0 + p_sep <= db[:,be._cl_st_w2]) & (db[:,be._cl_st_w2] < (st_start_ch0 + p_sep + st_len_w2)) & (db[:,be._cl_ch_w2] == 0)) \
                 | ((st_start_ch1 + p_sep <= db[:,be._cl_st_w2]) & (db[:,be._cl_st_w2] < (st_start_ch1 + p_sep + st_len_w2)) & (db[:,be._cl_ch_w2] == 1)) )   

    #no_invalid_mrkr_fltr = (dlt[:,be._cl_inv_mrkr]==0)

    sp_names=['w1','w2']
    sp_fltrs = [st_fltr_w1,st_fltr_w2]
    valid_event_fltr_SP = (db[:,be._cl_type] == 2) | (db[:,be._cl_type] == 3)
    st_fltr = st_fltr_w1 | st_fltr_w2
    dt_fltr = True
    rnd_fltr = (dlt[:,be._cl_noof_rnd_0] + dlt[:,be._cl_noof_rnd_1] == 1 )
    corr_mats={}

    ### reinserted because stuff broke: norbert
    F0 =analysis_params['F0A'] if lt3 else analysis_params['F0B']
    F1 =analysis_params['F1A'] if lt3 else analysis_params['F1B']
    print 'I changed stuff in get_sp_corrs!!! Norbert'

    for psi_name,psi_fltr in zip(sp_names,sp_fltrs):
        fltr = valid_event_fltr_SP & rnd_fltr & psi_fltr #& no_invalid_mrkr_fltr
        db_fltr = db[fltr]
        dlt_fltr = dlt[fltr]
        noof_ev_fltr = np.sum(fltr)
        if VERBOSE:
            print psi_name, 'noof events: ', noof_ev_fltr
        if noof_ev_fltr > 0:
            p0 = float(np.sum(dlt_fltr[:,be._cl_noof_ph_ro]>0))/noof_ev_fltr
            u_p0 = np.sqrt(p0*(1-p0)/noof_ev_fltr)
        else: p0, u_p0=0,0

        p0_corr, u_p0_corr = RO_correct_single_qubit(p0,u_p0, F0, F1,0.005,0.005) #XXXXXXXXXXXXXXXXXX
        corr_mats[psi_name] = [p0,u_p0,p0_corr,u_p0_corr,noof_ev_fltr]
    return corr_mats



def get_corr_mats(db,d3,d4, db_fps, analysis_params, bad_time_ranges, ret_fltr_psi_name='psi_min', VERBOSE=True):
    #Windows & other filtes:
    st_start_ch0  = analysis_params['st_start_ch0']
    st_len        = analysis_params['st_len'] #50 ns
    st_len_w2_00  = analysis_params['st_len_w2_00']
    st_len_w2_11  = analysis_params['st_len_w2_11']
    p_sep         = analysis_params['pulse_sep'] #600 ns
    st_start_ch1  = analysis_params['st_start_ch1']
    ro_start      = analysis_params['ro_start']
    ro_length     = analysis_params['ro_length']
    max_inv_distance = analysis_params['invalid_marker_max_sn_diff']
        
    #bad times due to lights on at EWI
    event_times = []
    bad_time_fltr = np.ones(len(db), dtype=np.bool)

    for i in range(len(db)):
        event_time = tb.get_datetime_from_folder(os.path.split(db_fps[i])[0]) \
                    + datetime.timedelta(seconds = db[i,be._cl_tt_w1]/1e12)
        event_times.append(event_time)
        for bad_time_range in bad_time_ranges:
            bad_time_fltr[i] = (bad_time_fltr[i]) and ((event_time <=  bad_time_range[0]) or (event_time >  bad_time_range[1]))
    if VERBOSE:
        print 'Events in bad time ranges: {}/{}'.format(len(db)-np.sum(bad_time_fltr), len(db))

    #invalid data marker filter & BS invalid event filter
    no_invalid_mrkr_fltr =  (((d3[:,be._cl_last_inv_mrkr]==0) | (d3[:,be._cl_last_inv_mrkr]>=max_inv_distance)) \
                          & ((d4[:,be._cl_last_inv_mrkr]==0) | (d4[:,be._cl_last_inv_mrkr]>=max_inv_distance)))
    valid_event_fltr = db[:,be._cl_type] == 1
    rnd_fltr = (d3[:,be._cl_noof_rnd_0] + d3[:,be._cl_noof_rnd_1] == 1 ) \
                 & (d4[:,be._cl_noof_rnd_0] + d4[:,be._cl_noof_rnd_1] == 1 ) 
    psb_tail_fltr = (d3[:,be._cl_first_tail_st] == 0) & (d4[:,be._cl_first_tail_st] == 0)

    #BS event channels
    psi_plus_00_fltr = (db[:,be._cl_ch_w1] == 0) & (db[:,be._cl_ch_w2] == 0)
    psi_min_01_fltr  = (db[:,be._cl_ch_w1] == 0) & (db[:,be._cl_ch_w2] == 1)
    psi_min_10_fltr  = (db[:,be._cl_ch_w1] == 1) & (db[:,be._cl_ch_w2] == 0)
    psi_plus_11_fltr = (db[:,be._cl_ch_w1] == 1) & (db[:,be._cl_ch_w2] == 1)
    psi_filters = [psi_plus_00_fltr,psi_min_01_fltr,psi_min_10_fltr,psi_plus_11_fltr]
    psi_names = ['psi_plus_00','psi_min_01','psi_min_10','psi_plus_11']
    ret_fltr = False
    corr_mats={}
    for psi_name,psi_fltr in zip(psi_names,psi_filters):
        if 'psi_min' in psi_name:
            st_len_w2 = st_len
        elif 'psi_plus_00' in psi_name:
            st_len_w2 = st_len_w2_00
        elif 'psi_plus_11' in psi_name:
             st_len_w2 = st_len_w2_11
                
        st_fltr_w1 =  (((st_start_ch0 <= db[:,be._cl_st_w1]) & (db[:,be._cl_st_w1] < (st_start_ch0 + st_len)) & (db[:,be._cl_ch_w1] == 0)) \
                     | ((st_start_ch1 <= db[:,be._cl_st_w1]) & (db[:,be._cl_st_w1] < (st_start_ch1 + st_len)) & (db[:,be._cl_ch_w1] == 1)) )  
        st_fltr_w2 =  (((st_start_ch0 + p_sep <= db[:,be._cl_st_w2]) & (db[:,be._cl_st_w2] < (st_start_ch0 + p_sep + st_len_w2)) & (db[:,be._cl_ch_w2] == 0)) \
                     | ((st_start_ch1 + p_sep <= db[:,be._cl_st_w2]) & (db[:,be._cl_st_w2] < (st_start_ch1 + p_sep + st_len_w2)) & (db[:,be._cl_ch_w2] == 1)) )   
        st_fltr = st_fltr_w1 & st_fltr_w2

        fltr = psi_fltr & valid_event_fltr
        
        noof_ev_fltr = np.sum(fltr) 
        fltr = fltr & st_fltr
        if VERBOSE:
            print '-'*40
            print 'Sync time filter {} : {}/{}'.format(psi_name,np.sum(fltr),noof_ev_fltr)
        
        noof_ev_fltr = np.sum(fltr)
        fltr = fltr & no_invalid_mrkr_fltr  
        if VERBOSE:
            print 'no_invalid_mrkr_fltr {} : {}/{}'.format(psi_name,np.sum(fltr),noof_ev_fltr)

        noof_ev_fltr = np.sum(fltr)
        fltr = fltr & rnd_fltr  
        if VERBOSE:
            print np.sum(np.logical_not(rnd_fltr))
            print 'rnd_fltr {} : {}/{}'.format(psi_name,np.sum(fltr),noof_ev_fltr)

        noof_ev_fltr = np.sum(fltr)
        fltr = fltr & psb_tail_fltr  
        if VERBOSE:
            print 'psb_tail_fltr {} : {}/{}'.format(psi_name,np.sum(fltr),noof_ev_fltr)

        noof_ev_fltr = np.sum(fltr)
        fltr = fltr & bad_time_fltr  
        if VERBOSE:
            print 'bad_time_fltr {} : {}/{}'.format(psi_name,np.sum(fltr),noof_ev_fltr)

        noof_ev_fltr = np.sum(fltr)
        if VERBOSE:
            lt_valid_fltr = (d3[:,be._cl_sn_ma] > 0) & (d4[:,be._cl_sn_ma] > 0)
            noof_lt_inv = noof_ev_fltr - np.sum(fltr & lt_valid_fltr)
            print 'LT invalid events {} : {}/{}'.format(psi_name,noof_lt_inv,noof_ev_fltr)

        
        db_fltr = db[fltr]
        d3_fltr = d3[fltr]
        d4_fltr = d4[fltr]
        noof_ev_fltr = np.sum(fltr)   
        
        corr_mat=np.zeros((4,4))

        for i,rnd in enumerate(rnd_corr):
            for j,ro in enumerate(ro_corr):
                corr_mat[i,j] =np.sum( 
                                       (d3_fltr[:,be._cl_noof_rnd_0]      == rnd[0]) \
                                     & (d4_fltr[:,be._cl_noof_rnd_1]      == rnd[1]) \
                                     & (((d3_fltr[:,be._cl_first_ph_st] > ro_start) & (d3_fltr[:,be._cl_first_ph_st] < ro_start+ro_length)) == ro[0] )\
                                     & (((d4_fltr[:,be._cl_first_ph_st] > ro_start) & (d4_fltr[:,be._cl_first_ph_st] < ro_start+ro_length)) == ro[1])
                                     )
        corr_mats[psi_name] = corr_mat
        if ret_fltr_psi_name in psi_name:
            ret_fltr = ret_fltr | fltr
    
    return corr_mats, ret_fltr

def print_ZZ_fidelity(corr_mats, analysis_params):
    for psi in ['psi_min', 'psi_plus']:
        corr_mat=np.zeros((4,4))
        for k in corr_mats:
            if psi in k:
                corr_mat+=corr_mats[k]
        noof_ev_fltr = int(np.sum(corr_mat))

        Ng_ZZ = corr_mat[0,1] + corr_mat[0,2] +corr_mat[1,0] + corr_mat[1,3] + corr_mat[2,0] + corr_mat[2,3] + corr_mat[3,1] + corr_mat[3,2] #for ZZ with RND
        corr_mat_RO_corr=np.zeros((4,4))
        for i in range(4) : 
            RO_row = corr_mat[i]
            RO_row_corr= RO_correct(RO_row/float(noof_ev_fltr),
                                    F0A =analysis_params['F0A'], F1A =analysis_params['F1A'],
                                    F0B =analysis_params['F0B'], F1B =analysis_params['F1B'] )
            for j in range(4):
                corr_mat_RO_corr[i,j] = RO_row_corr[j]
        Ng_ZZ_RO_corr = corr_mat_RO_corr[0,1] + corr_mat_RO_corr[0,2] +corr_mat_RO_corr[1,0] + corr_mat_RO_corr[1,3] \
                        + corr_mat_RO_corr[2,0] + corr_mat_RO_corr[2,3] + corr_mat_RO_corr[3,1] + corr_mat_RO_corr[3,2]
        RO_corr = np.zeros(4)
        RO_corr[0] = corr_mat_RO_corr[0,0] + corr_mat_RO_corr[1,1] +corr_mat_RO_corr[2,1] + corr_mat_RO_corr[3,0]
        RO_corr[1] = corr_mat_RO_corr[0,1] + corr_mat_RO_corr[1,0] +corr_mat_RO_corr[2,0] + corr_mat_RO_corr[3,1]
        RO_corr[2] = corr_mat_RO_corr[0,2] + corr_mat_RO_corr[1,3] +corr_mat_RO_corr[2,3] + corr_mat_RO_corr[3,2]
        RO_corr[3] = corr_mat_RO_corr[0,3] + corr_mat_RO_corr[1,2] +corr_mat_RO_corr[2,1] + corr_mat_RO_corr[3,3]

        Fid = Ng_ZZ/noof_ev_fltr
        dFid = np.sqrt(Fid*(1-Fid)/noof_ev_fltr)
        Fid_corr= Ng_ZZ_RO_corr
        
        print '-'*40
        print 'FILTERED EVENTS {}: Number of events {}'.format(psi,noof_ev_fltr)
        print 'Fidelity: {:.2f}% +- {:.1f}%'.format(Fid*100,dFid*100)
        print 'Corrected fidelity: {:.2f}%'.format(Fid_corr*100)




def print_XX_fidelity(corr_mats, analysis_params):
    for psi in ['psi_min', 'psi_plus']:
        corr_mat=np.zeros((4,4))
        for k in corr_mats:
            if psi in k:
                corr_mat+=corr_mats[k]
        noof_ev_fltr = int(np.sum(corr_mat))

        
        corr_mat_RO_corr=np.zeros((4,4))
        for i in range(4) : 
            RO_row = corr_mat[i]
            RO_row_corr= RO_correct(RO_row/float(noof_ev_fltr),
                                    F0A =analysis_params['F0A'], F1A =analysis_params['F1A'],
                                    F0B =analysis_params['F0B'], F1B =analysis_params['F1B'] )
            for j in range(4):
                corr_mat_RO_corr[i,j] = RO_row_corr[j]
        
        RO_corr = np.zeros(4)
        RO_corr[0] = corr_mat_RO_corr[0,0] + corr_mat_RO_corr[1,1] +corr_mat_RO_corr[2,1] + corr_mat_RO_corr[3,0]
        RO_corr[1] = corr_mat_RO_corr[0,1] + corr_mat_RO_corr[1,0] +corr_mat_RO_corr[2,0] + corr_mat_RO_corr[3,1]
        RO_corr[2] = corr_mat_RO_corr[0,2] + corr_mat_RO_corr[1,3] +corr_mat_RO_corr[2,3] + corr_mat_RO_corr[3,2]
        RO_corr[3] = corr_mat_RO_corr[0,3] + corr_mat_RO_corr[1,2] +corr_mat_RO_corr[2,1] + corr_mat_RO_corr[3,3]

        if psi == 'psi_min' :
            Ng_XX = corr_mat[0,0] + corr_mat[0,3] +corr_mat[1,1] + corr_mat[1,2] + corr_mat[2,1] + corr_mat[2,2] + corr_mat[3,0] + corr_mat[3,3]
            Ng_XX_RO_corr = corr_mat_RO_corr[0,0] + corr_mat_RO_corr[0,3] +corr_mat_RO_corr[1,1] + corr_mat_RO_corr[1,2] \
                        + corr_mat_RO_corr[2,1] + corr_mat_RO_corr[2,2] + corr_mat_RO_corr[3,0] + corr_mat_RO_corr[3,3]
        else : 
            Ng_XX = corr_mat[0,1] + corr_mat[0,2] +corr_mat[1,0] + corr_mat[1,3] + corr_mat[2,0] + corr_mat[2,3] + corr_mat[3,1] + corr_mat[3,2]
            Ng_XX_RO_corr = corr_mat_RO_corr[0,1] + corr_mat_RO_corr[0,2] +corr_mat_RO_corr[1,0] + corr_mat_RO_corr[1,3] \
                        + corr_mat_RO_corr[2,0] + corr_mat_RO_corr[2,3] + corr_mat_RO_corr[3,1] + corr_mat_RO_corr[3,2]

        Fid = Ng_XX/noof_ev_fltr
        dFid = np.sqrt(Fid*(1-Fid)/noof_ev_fltr)
        Fid_corr= Ng_XX_RO_corr
        
        print '-'*40
        print 'FILTERED EVENTS {}: Number of events {}'.format(psi,noof_ev_fltr)
        print 'Fidelity: {:.2f}% +- {:.1f}%'.format(Fid*100,dFid*100)
        print 'Corrected fidelity: {:.2f}%'.format(Fid_corr*100)



def plot_title(fp):
    return os.path.splitext(os.path.split(fp)[1])[0]

def save_fp(folder,analysis_fp):
    return os.path.join(folder,plot_title(analysis_fp))

def save_figure(name, ax,output_folder,analysis_fp):
    ax.figure.savefig(save_fp(output_folder,analysis_fp)+'_'+name+'.jpg')

def sweep_analysis_param(analysis_fp,settings,bad_time_ranges,sweep_params,xmin,xmax,pts, relative):

    f = h5py.File(analysis_fp,'r')
    db_fps=f['analysis']['total_ent_events_fps'].value
    db = f['analysis']['total_ent_events'].value
    d3 = f['analysis']['total_lt3_ssro'].value
    d4 = f['analysis']['total_lt4_ssro'].value
    f.close()
    X=np.linspace(xmin,xmax,pts)
    S=np.zeros(pts)
    dS=np.zeros(pts)
    Ns=np.zeros(pts)
    Ks=np.zeros(pts)
    P=np.zeros(pts)
    for i,x in enumerate(X):
        ap=settings.analysis_params.copy()
        #ap['st_start_ch0']=x
        #ap['st_start_ch1']=x-900
        for sweep_param in sweep_params:
            ap[sweep_param] = ap[sweep_param] + x if relative else x
        corr_mats, _tmp = get_corr_mats(db,d3,d4,db_fps, ap, bad_time_ranges, VERBOSE=False)
        K,N,Kxx,Nxx,Kzz,Nzz,p_val = calculate_p_lhv(corr_mats, VERBOSE=False)
        CHSH, dCHSH, Es, dEs = print_correlators(corr_mats, VERBOSE=False, psis=['psi_min'])
        S[i]  = CHSH
        dS[i] = dCHSH
        Ns[i] = N
        Ks[i] = K
        P[i]  = p_val

    return X, S, dS, P, Ns, Ks