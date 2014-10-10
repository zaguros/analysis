''' Theoretically calibrated values for gate parameters. These differ slightly
from the experimental parameters due to uncertainties in the field and hyperfine params'''
import numpy as np
### for B_field = 304.74
gp = {} #gate params

gp['C1_freq']       =   345.124e3
gp['C1_freq_0']     =   325.787e3
gp['C1_freq_1']     =   364.570e3
gp['C1_freq_dec']   =   345.124e3
gp['C1_Ren_extra_phase_correction_list'] = np.array([0]*3 + [-132] + [0]*6)
gp['C1_Ren_tau']    =   [9.420e-6]
gp['C1_Ren_N']      =   [16]

gp['C2_freq']       =   339.955e3
gp['C2_Ren_tau']    =   [6.62e-6]
gp['C2_Ren_N']      =   [26]

gp['C3_freq']       =   302.521e3
gp['C3_freq_0']     =   325.775e3
gp['C3_freq_1']     =   293.888e3
gp['C3_freq_dec']   =   302.521e3
gp['C3_Ren_extra_phase_correction_list'] = np.array([0]*10)
gp['C3_Ren_tau']    =   [18.542e-6]
gp['C3_Ren_N']      =   [16]

gp['C4_freq']       =   348.574e3
gp['C4_freq_0']     =   325.787e3
gp['C4_freq_1']     =   370.115e3
gp['C4_freq_dec']   =   348.574e3
gp['C4_Ren_extra_phase_correction_list'] = np.array([0] +[-90] + [0]*8)
gp['C4_Ren_tau']    =   [6.456e-6   ]
gp['C4_Ren_N']      =   [40         ]
