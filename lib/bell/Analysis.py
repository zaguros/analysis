import numpy as np
import h5py
import Settings, files, Filter

def analyze_raw_data(Total_entanglement_events):
    # Makes boolean filters which determine if an event is psiminus or psiplus (see Filter.for_loop_sync_num for definition of both)
    is_psiminus = Total_entanglement_events[:,12] == 1
    is_psiplus = Total_entanglement_events[:,12] == 0
    
    if Settings.VERBOSE:
        print 'Number of psiminus events:', sum(is_psiminus)
        print 'Number of psiplus event:', sum(is_psiplus)

    # Makes boolean filters which determine if the SSRO correspond with the up (photons are emitted) or down (no photons are emitted) state
    is_up_LT1 = Total_entanglement_events[:,6] > 0
    is_down_LT1 = Total_entanglement_events[:,6] == 0
    is_up_LT3 = Total_entanglement_events[:,7] > 0
    is_down_LT3 = Total_entanglement_events[:,7] == 0

    # Makes boolean filters for up, up; up, down; down, up; down, down events
    is_upLT1_upLT3 = is_up_LT1 & is_up_LT3
    is_upLT1_downLT3 = is_up_LT1 & is_down_LT3
    is_downLT1_upLT3 = is_down_LT1 & is_up_LT3
    is_downLT1_downLT3 = is_down_LT1 & is_down_LT3

    # Makes boolean filters for psiplus state
    is_psiplus_up_up = is_upLT1_upLT3 & is_psiplus
    is_psiplus_up_down = is_upLT1_downLT3 & is_psiplus
    is_psiplus_down_up = is_downLT1_upLT3 & is_psiplus
    is_psiplus_down_down = is_downLT1_downLT3 & is_psiplus

    # Makes boolean filters for psiminus state
    is_psiminus_up_up = is_upLT1_upLT3 & is_psiminus
    is_psiminus_up_down = is_upLT1_downLT3 & is_psiminus
    is_psiminus_down_up = is_downLT1_upLT3 & is_psiminus
    is_psiminus_down_down = is_downLT1_downLT3 & is_psiminus

    # Determines the amount of up, up; up, down; down, up and down, down events for the psiplus state
    psiplus_up_up = sum(is_psiplus_up_up)
    psiplus_up_down = sum(is_psiplus_up_down)
    psiplus_down_up = sum(is_psiplus_down_up)
    psiplus_down_down = sum(is_psiplus_down_down)

    psiplus_bars = np.array([psiplus_up_up, psiplus_up_down, psiplus_down_up, psiplus_down_down])
    psiplus_bars_norm = psiplus_bars/float(sum(psiplus_bars))


    # Determines the amount of up, up; up, down; down, up and down, down events for the psiminus state
    psiminus_up_up = sum(is_psiminus_up_up)
    psiminus_up_down = sum(is_psiminus_up_down)
    psiminus_down_up = sum(is_psiminus_down_up)
    psiminus_down_down = sum(is_psiminus_down_down)

    psiminus_bars = np.array([psiminus_up_up, psiminus_up_down, psiminus_down_up, psiminus_down_down])
    psiminus_bars_norm = psiminus_bars/float(sum(psiminus_bars))
    
    return psiplus_bars, psiplus_bars_norm, psiminus_bars, psiminus_bars_norm