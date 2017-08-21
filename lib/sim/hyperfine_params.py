#########################################################
######### No1 SIL 18#####################################
#########################################################

# in Rad/Hz


### fitted to the ms = -1 data

hyperfine_params = {}
hyperfine_params['espin_trans'] = '-1'
hyperfine_params['C1']  = {'par' : 33.0e3, 'perp':35.0e3}
hyperfine_params['C2']  = {'par' : 26.5e3, 'perp':30.0e3}
# hyperfine_params['C3']  = {'par' : -365.0e3, 'perp':80.0e3}
# hyperfine_params['C4']  = {'par' : 8.1e3, 'perp':21.0e3}
# hyperfine_params['C5']  = {'par' : 24.7e3, 'perp':26.0e3}
# hyperfine_params['C6']  = {'par' : -48.7e3, 'perp':12.0e3}
# hyperfine_params['C7']  = {'par' : 14.5e3, 'perp':11.0e3}
# hyperfine_params['C8']  = {'par' : -20.5e3, 'perp':21.0e3}
# hyperfine_params['C8']  = {'par' : 7.32e3, 'perp':5.0e3}

### Adding some random test spins: spins 9 and 10 are not "real" spins obtained from data
hyperfine_params['C9']  = {'par' : 0.32e3, 'perp':0.2e3}
hyperfine_params['C10']  = {'par' : 3.32e3, 'perp':1.2e3}


#########################################################
######### HANS SIL 1#####################################
#########################################################

### fitted to the ms = +1 data
# hyperfine_params_hans_SIL1_msp1 = {}
# hyperfine_params_hans_SIL1_msp1['C1']  = {'par' : 30.0e3, 'perp':80.0e3}
# hyperfine_params_hans_SIL1_msp1['C2']  = {'par' : 27.0e3, 'perp':28.5e3}
# hyperfine_params_hans_SIL1_msp1['C3']  = {'par' :-51.0e3, 'perp':105.0e3}
# hyperfine_params_hans_SIL1_msp1['C4']  = {'par' : 45.1e3,   'perp':20e3}   # set paralel term from 45 to 45.1 (makes gates work a little better) (140728 - JULIA)
# hyperfine_params_hans_SIL1_msp1['C5']  = {'par' : 17e3,   'perp':10e3}
# hyperfine_params_hans_SIL1_msp1['C6']  = {'par' :-15e3,   'perp':12e3}
# hyperfine_params_hans_SIL1_msp1['C7']  = {'par' :-23e3,   'perp':12e3}
# hyperfine_params_hans_SIL1_msp1['C8']  = {'par' : 10e3,   'perp':8.0e3}
# hyperfine_params_hans_SIL1_msp1['C9']  = {'par' : 8.0e3,  'perp':12e3}
# hyperfine_params_hans_SIL1_msp1['C10'] = {'par' :-9.3e3,  'perp':13e3}
# hyperfine_params_hans_SIL1_msp1['C11'] = {'par' :-10.0e3, 'perp':5e3}
# hyperfine_params_hans_SIL1_msp1['C12'] = {'par' :-30.0e3, 'perp':35e3}
# hyperfine_params_hans_SIL1_msp1['C13'] = {'par' :-32.0e3, 'perp':20e3}

### fitted to the ms = -1 data
hyperfine_params_hans_SIL1_msm1 = {}
hyperfine_params_hans_SIL1_msm1['C1']  = {'par' : 30.7e3, 'perp':80.0e3}
hyperfine_params_hans_SIL1_msm1['C2']  = {'par' : 29.0e3, 'perp':28.5e3}
hyperfine_params_hans_SIL1_msm1['C3']  = {'par' :-51e3, 'perp':105.0e3}  # 'par' :-53.5e3 fits better for the ms = -1 msmts. But there is always the field uncertainty
hyperfine_params_hans_SIL1_msm1['C4']  = {'par' : 45.1e3,   'perp':25e3}   # set paralel term from 45 to 45.1 (makes gates work a little better) (140728 - JULIA)
hyperfine_params_hans_SIL1_msm1['C5']  = {'par' : 17e3,   'perp':10e3}
hyperfine_params_hans_SIL1_msm1['C6']  = {'par' :-15e3,   'perp':14e3}
hyperfine_params_hans_SIL1_msm1['C7']  = {'par' :-23e3,   'perp':16e3}
hyperfine_params_hans_SIL1_msm1['C8']  = {'par' : 10e3,   'perp':8.0e3}
hyperfine_params_hans_SIL1_msm1['C9']  = {'par' : 8.0e3,  'perp':12e3}
hyperfine_params_hans_SIL1_msm1['C10'] = {'par' :-9.3e3,  'perp':13e3}
hyperfine_params_hans_SIL1_msm1['C11'] = {'par' :-10.0e3, 'perp':5e3}
hyperfine_params_hans_SIL1_msm1['C12'] = {'par' :-30.0e3, 'perp':35e3}
hyperfine_params_hans_SIL1_msm1['C13'] = {'par' :-32.0e3, 'perp':20e3}


# ### fitted to the ms = -1 data
# hyperfine_params_hans_SIL1_msm1 = {}
# hyperfine_params_hans_SIL1_msm1['C1']  = {'par' : 30.7e3, 'perp':80.0e3}
# hyperfine_params_hans_SIL1_msm1['C2']  = {'par' : 29.0e3, 'perp':28.5e3}
# # hyperfine_params_hans_SIL1_msm1['C3']  = {'par' :-53.5e3, 'perp':105.0e3}  # 'par' :-53.5e3 fits better for the ms = -1 msmts. But there is always the field uncertainty
# hyperfine_params_hans_SIL1_msm1['C3']  = {'par' : 45.1e3,   'perp':25e3}   # set paralel term from 45 to 45.1 (makes gates work a little better) (140728 - JULIA)
# hyperfine_params_hans_SIL1_msm1['C4']  = {'par' : 17e3,   'perp':10e3}
# hyperfine_params_hans_SIL1_msm1['C5']  = {'par' :-15e3,   'perp':14e3}
# hyperfine_params_hans_SIL1_msm1['C6']  = {'par' :-23e3,   'perp':16e3}
# hyperfine_params_hans_SIL1_msm1['C7']  = {'par' : 10e3,   'perp':8.0e3}
# hyperfine_params_hans_SIL1_msm1['C8']  = {'par' : 8.0e3,  'perp':12e3}
# hyperfine_params_hans_SIL1_msm1['C9'] = {'par' :-9.3e3,  'perp':13e3}
# hyperfine_params_hans_SIL1_msm1['C10'] = {'par' :-10.0e3, 'perp':5e3}
# hyperfine_params_hans_SIL1_msm1['C11'] = {'par' :-30.0e3, 'perp':35e3}
# hyperfine_params_hans_SIL1_msm1['C12'] = {'par' :-32.0e3, 'perp':20e3}

# ###TEST FOR BUGFINDGING in theory
# hyperfine_params_hans_SIL1_msm1['C1']  = {'par' : 29.3e3, 'perp':80.0e3}
# hyperfine_params_hans_SIL1_msm1['C2']  = {'par' : 27.0e3, 'perp':28.5e3}
# hyperfine_params_hans_SIL1_msm1['C3']  = {'par' :-51.0e3, 'perp':105.0e3}
# hyperfine_params_hans_SIL1_msm1['C4']  = {'par' : 44.2e3,   'perp':23.5e3}   # set paralel term from 45 to 45.1 (makes gates work a little better) (140728 - JULIA)
# hyperfine_params_hans_SIL1_msm1['C5']  = {'par' : 17e3,   'perp':10e3}
# hyperfine_params_hans_SIL1_msm1['C6']  = {'par' :-15e3,   'perp':12e3}
# hyperfine_params_hans_SIL1_msm1['C7']  = {'par' :-23e3,   'perp':12e3}
