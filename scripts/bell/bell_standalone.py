#this analyses all the raw bell data.
import os
import numpy as np
import h5py
#lets collect the data

def get_all_measurement_filepaths(folder, suffix='hdf5', pattern=''):
    filepaths = []
    suffixlen = len(suffix)
    
    for root,dirs,files in os.walk(folder):
        for f in files:
            if len(f) > suffixlen and f[-suffixlen:] == suffix and pattern in f:
                filepaths.append(os.path.join(root, f))
                
    filepaths = sorted(filepaths)
    
    return filepaths

folder_charlie = r'D:\measuring\data\2015-03-18-lhfbt2\BS'
measurement_pattern = 'second_ever'
folder_alice        = os.path.join(os.path.split(folder_charlie)[0],'LT3')
folder_bob          = os.path.join(os.path.split(folder_charlie)[0],'LT4')

filepaths_charlie   = get_all_measurement_filepaths(folder_charlie, pattern=measurement_pattern)[:]
filepaths_alice     = get_all_measurement_filepaths(folder_alice,   pattern=measurement_pattern)[:]
filepaths_bob       = get_all_measurement_filepaths(folder_bob,     pattern=measurement_pattern)[:]

start_coarse        = 5430000#5415000   
window_coarse       = 250000
start_P0            = 5443600#5426600 #ps after sync-pulse
start_P1            = 5443600+1000#5426600 - 900 #ps after sync-pulse
window_length       = 50000  #ps after start
window_separation   = 350000#250000 #ps between early and late windows

readout_start_P     = 10620 #ns after sync-pulse
readout_length_P    = 3700  #ns after sync-pulse
psb_start_P_alice   = 7480  #ns after sync-pulse
psb_start_P_bob     = 5350  #ns after sync-pulse
psb_window_length_P = 200   #ns after sync-pulse

total_event_types_charlie = -1*np.ones((0,1), dtype=np.int) # Charlie event type, 
total_random_input_alice  = -1*np.ones((0,1), dtype=np.int) # Alice random input,
total_random_input_bob    = -1*np.ones((0,1), dtype=np.int) # Bob random input,
total_readout_alice       = -1*np.ones((0,1), dtype=np.int) # Alice readout result,
total_readout_bob         = -1*np.ones((0,1), dtype=np.int) # Bob readout result,
total_invalid_alice       = -1*np.ones((0,1), dtype=np.int) # Alice not ready,
total_invalid_bob         = -1*np.ones((0,1), dtype=np.int) # Bob not ready. 

ii=0
for fp_c,fp_a,fp_b in zip(filepaths_charlie, filepaths_alice, filepaths_bob):
    if fp_c[-22:] != fp_a[-22:] or fp_c[-22:] != fp_b[-22:]:
        raise Exception('Combining wrong data files!: {}, {}, {}'.format(fp_c,fp_a,fp_b))
    print ii,
    
    ################# Get the heralded entanglement events from Charlie: #################
    f_c = h5py.File(fp_c,'r')
    channel     = f_c['/PQ_channel-1'].value       # defines event channel together with special_bit
    special     = f_c['/PQ_special-1'].value     # defines event channel together with channel
    sync_number = f_c['/PQ_sync_number-1'].value # number of sync-pulses received since start of run
    sync_time   = f_c['/PQ_sync_time-1'].value   # time since last sync-pulse in picoseconds
    total_time  = f_c['/PQ_time-1'].value          # time since start of run in picoseconds
    f_c.close()

    E_event_filter  = (special == 1) & ((channel & 1) == 1) # entanglement markers are defined by a high special bit, and a channel 1
    number_of_E_events   = np.sum(E_event_filter)
    E_sync_numbers  = sync_number[E_event_filter]

    #clicks on photon event inputs have special_bit zero
    filter_early_P0     = (special == 0) & (channel == 0) & (start_P0 <= sync_time)                       & (sync_time < (start_P0 + window_length))
    filter_early_P1     = (special == 0) & (channel == 1) & (start_P1 <= sync_time)                       & (sync_time < (start_P1 + window_length))
    filter_late_P0      = (special == 0) & (channel == 0) & ((start_P0 + window_separation) <= sync_time) & (sync_time < (start_P0 + window_separation + window_length))
    filter_late_P1      = (special == 0) & (channel == 1) & ((start_P1 + window_separation) <= sync_time) & (sync_time < (start_P1 + window_separation + window_length))
    
    filter_early = filter_early_P0 | filter_early_P1
    filter_late  = filter_late_P0  | filter_late_P1

    filter_early_coarse = (special == 0) & (start_coarse <= sync_time)                        & (sync_time < (start_coarse + window_coarse))
    filter_late_coarse  = (special == 0) & ((start_coarse + window_separation) <= sync_time)  & (sync_time < (start_coarse + window_separation + window_coarse))
    
    #lets mak a list of event types for every E-event: -1: not a valid event, 0: 00 event (psi-plus) 1: 01 event (psi_min) 2: 10 event (psi-min) 3: 11 event (psi_plus)
    event_types_charlie = -1*np.ones((number_of_E_events,1), dtype=np.int) 
    for i,E_sync_number in enumerate(E_sync_numbers):
        E_sync_number_filter = (sync_number == E_sync_number)
        if np.sum(filter_early_coarse & E_sync_number_filter) == 1 and np.sum(filter_late_coarse & E_sync_number_filter) == 1: #exactly one early and one late click
            if np.sum(filter_early & E_sync_number_filter) == 1 and np.sum(filter_late & E_sync_number_filter) == 1: 
                event_types_charlie[i] = 2*np.sum(filter_late_P1 & E_sync_number_filter) + np.sum(filter_early_P1 & E_sync_number_filter)

    ################# Get the corresponding correlation results from Alice and Bob #################
    random_input_alice = -1*np.ones((number_of_E_events,1), dtype=np.int)
    random_input_bob   = -1*np.ones((number_of_E_events,1), dtype=np.int)
    readout_alice      = -1*np.ones((number_of_E_events,1), dtype=np.int)
    readout_bob        = -1*np.ones((number_of_E_events,1), dtype=np.int)
    invalid_alice      = -1*np.ones((number_of_E_events,1), dtype=np.int)
    invalid_bob        = -1*np.ones((number_of_E_events,1), dtype=np.int)
    for location, fp_ab,  random_input, readout, invalid in zip(['Alice','Bob'], [fp_a,fp_b], [random_input_alice,random_input_bob], [readout_alice,readout_bob], [invalid_alice,invalid_bob]):
        f_ab = h5py.File(fp_ab,'r')
        channel     = f_ab['/PQ_channel-1'].value       # defines event channel together with special_bit
        special     = f_ab['/PQ_special-1'].value      # defines event channel together with channel
        sync_number = f_ab['/PQ_sync_number-1'].value  # number of sync-pulses received since start of run
        sync_time   = f_ab['/PQ_sync_time-1'].value    # time since last sync-pulse in picoseconds
        total_time  = f_ab['/PQ_time-1'].value          # time since start of run in picoseconds
        f_ab.close()

        E_event_filter_ab = (special == 1) & ((channel & 4) == 4) # entanglement markers at Alice & Bob are defined by a high special bit, and channel 4
        E_sync_numbers_ab = sync_number[E_event_filter_ab]

        I_event_filter_ab = (special == 1) & ((channel & 8) == 8)# entanglement markers at Alice & Bob are defined by a high special bit, and channel 8
        I_sync_numbers_ab = sync_number[I_event_filter_ab]

        psb_start_P = psb_start_P_alice if location == 'Alice' else  psb_start_P_bob
        filter_early_P    = (special == 0) & (channel == 0)  & (sync_time > psb_start_P)  & (sync_time < (psb_start_P  + psb_window_length_P))
        filter_late_P     = (special == 0) & (channel == 0)  & (sync_time > psb_start_P + window_separation)  & (sync_time < (psb_start_P  + psb_window_length_P + window_separation))
        filter_psb_P      = filter_early_P | filter_late_P
        filter_readout_P  = (special == 0) & (channel == 0)  & (sync_time > readout_start_P)  & (sync_time < (readout_start_P  + readout_length_P))
        filter_random_0   = (special == 1) & (channel == 1)
        filter_random_1   = (special == 1) & (channel == 2)
        filter_random     = filter_random_0 | filter_random_1

        for i,sn in enumerate(E_sync_numbers_ab):
            tigger_charlie = np.float(E_sync_numbers[i])/251.
            trigger_ab = np.float(sn)/250. if location == 'Alice' else  np.float(sn)/251. #Alice records one sync-pulse less per start
            if np.abs(tigger_charlie - trigger_ab) > 1.: #check that this E_event indeed corresponds to the same excitation round.
                #if not 'day1_run7' in fp_ab:
                print 'WARNING: sn diff error in {}: sn_ab-sn_c = {:.2f}'.format(fp_ab, trigger_ab - tigger_charlie)
                #print 's',
            entanglement_sync_number_filter = (sync_number == sn)
            E_event_sync_time  = sync_time[ entanglement_sync_number_filter & E_event_filter_ab]
            E_event_total_time = total_time[entanglement_sync_number_filter & E_event_filter_ab]
            if E_event_sync_time < 10000:
                entanglement_sync_number_filter = (sync_number == sn-1) #the entanglement signal arrives just after the next sync-pulse, 
                                                                        #unless it was the last in the 250 sync-pulses per CR check
            if np.int64(E_event_total_time) - np.int64(total_time[filter_random & entanglement_sync_number_filter][0]) > 7000 :
                raise Exception('Unexpected time difference between arrival of entanglement event and R event of previous round')

            #### the random input ####
            if np.sum(filter_random & entanglement_sync_number_filter ) == 1:
                random_input[i] = np.sum(filter_random_1 & entanglement_sync_number_filter)
                
            #### the readout result ####
            readout[i] = (np.sum(filter_readout_P & entanglement_sync_number_filter) == 0) #readout ms=1 if no photons, ms=0 otherwise

            #### the invalid marker ####
            sync_number_difference_E_vs_I_event = np.array(sn - I_sync_numbers_ab.astype(np.int64), dtype=np.int64)
            #if an I event occured less than 251 syncs in the past, we mark this run as invalid
            invalid_event_filter = (sync_number_difference_E_vs_I_event > 0) & (sync_number_difference_E_vs_I_event <= 251)
            #or if a PSB photon was detected during the entanglement generation
            invalid[i] = np.sum(invalid_event_filter) + np.sum(filter_psb_P & entanglement_sync_number_filter)


    total_event_types_charlie = np.vstack((total_event_types_charlie,   event_types_charlie))
    total_random_input_alice  = np.vstack((total_random_input_alice,    random_input_alice))
    total_random_input_bob    = np.vstack((total_random_input_bob,      random_input_bob))
    total_readout_alice       = np.vstack((total_readout_alice,         readout_alice))
    total_readout_bob         = np.vstack((total_readout_bob,           readout_bob))
    total_invalid_alice       = np.vstack((total_invalid_alice,         invalid_alice))
    total_invalid_bob         = np.vstack((total_invalid_bob,           invalid_bob))
    
    ii+=1




############## Calculate S for psi_min ################
psi_min_filter = (total_event_types_charlie == 1) | (total_event_types_charlie == 2)
not_invalid = (total_invalid_alice <= 0) & (total_invalid_bob <= 0)
fltr = psi_min_filter & not_invalid
errors_in_dataset = np.sum(total_random_input_alice[fltr]==-1)\
                  + np.sum(total_random_input_bob[fltr]==-1)\
                  + np.sum(total_readout_alice[fltr]==-1)\
                  + np.sum(total_readout_bob[fltr]==-1)
if errors_in_dataset > 0:
    raise Exception('errors_in_dataset: {:d}'.format(errors_in_dataset))


corr_mat=np.zeros((4,4))
rnd_corr=[[0,0],[0,1],[1,0],[1,1]] #'RND [LT3,LT4] rnd 00,01,10,11'
ro_corr =[[1,1],[1,0],[0,1],[0,0]] #'RO [LT3,LT4] ms  00, 01, 10, 11'

Es=np.zeros(4)
dEs=np.zeros(4)
K=0
N=0
for i,rnd in enumerate(rnd_corr):
    for j,ro in enumerate(ro_corr):
        corr_mat[i,j] =np.sum( (total_random_input_alice[fltr] == rnd[0]) \
                             & (total_random_input_bob[fltr]   == rnd[1]) \
                             & (total_readout_alice[fltr]      == ro[0] ) \
                             & (total_readout_bob[fltr]        == ro[1]))
        if rnd[0] ==0 and rnd[1]==0:
            C_val=(ro[0]+ro[1])%2
        else:
            C_val=(ro[0]+ro[1]+1)%2
        K+=(C_val==1)*corr_mat[i,j]
        N+= corr_mat[i,j]
    Es[i] = (corr_mat[i,0] - corr_mat[i,1] - corr_mat[i,2] + corr_mat[i,3])/float(np.sum(corr_mat[i,:]))
    n=np.float(np.sum(corr_mat[i,:]))
    p_est=(corr_mat[i,0]+corr_mat[i,3])/n
    dEs[i] = 2*np.sqrt(p_est*(1-p_est)/n)

expected_Es= (-0.48,0.48,0.6,0.6)
CHSH  = -Es[0] + Es[1] + Es[2] + Es[3]
dCHSH = np.sqrt(dEs[0]**2 + dEs[1]**2 + dEs[2]**2 + dEs[3]**2)

print ''        
print '-'*40
print 'FILTERED EVENTS {}: Number of events {}'.format('psi_min',N)
print 'RO ms   00, 01, 10, 11'
print 'RND00', corr_mat[0], '  +pi/2, +3pi/4'
print 'RND01', corr_mat[1], '  +pi/2, -3pi/4'
print 'RND10', corr_mat[2], '  0,     +3pi/4'
print 'RND11', corr_mat[3], '  0,     -3pi/4\n'

print ' E (RND00  RND01  RND10  RND11 )'
print '   ({:+.2f}, {:+.2f}, {:+.2f}, {:+.2f}) expected'.format(*expected_Es)
print '   ({:+.2f}, {:+.2f}, {:+.2f}, {:+.2f}) measured'.format(*Es)
print '+/-( {:.2f},  {:.2f},  {:.2f},  {:.2f} )'.format(*dEs)

print 'CHSH : {:.2f} +- {:.2f}'.format(CHSH, dCHSH)

from scipy.stats import binom
tau = 1e-25
eps = 1e-5
p_lhv = 1-(0.5-eps)**2*(1-12*tau*(1+tau))#eps and tau correspond to the bias and the predictability.
p_val = 1- binom.cdf(K-1, N, p_lhv)
print 'All: {}/{} = {:.2f}'.format(K, N, K/N)
print 'Probability of LHV model: {:.1f}%'.format(p_val*100)