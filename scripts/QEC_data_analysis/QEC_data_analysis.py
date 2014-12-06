'''script for analysis of all QEC data'''


import sys
import numpy as np

sys.path.append(r'D:/measuring')
sys.path.append(r'D:/measuring/analysis')

from analysis.scripts.QEC import QEC_analysis as QEC
from analysis.scripts.QEC import three_qubit_QEC_analysis as three_QEC
reload(QEC)
reload(three_QEC)

''' save the datasets '''

three_QEC.save_QEC_dataset(older_than = '20141120_103734', no_error = '00')

# three_QEC.save_QEC_dataset(older_than = '20141121_141906', no_error = '11')

# three_QEC.save_QEC_dataset(older_than = '20141122_192627', no_error = '01')

# three_QEC.save_QEC_dataset(older_than = '20141124_130000', no_error = '10_!')


''' plot single figures '''

# for state in ['Z','mZ','Y','mY','X','mX']:
#     for RO in range(6):
#         three_QEC.QEC_plot_single_state_RO(date = '20141120', no_error = '00',state = state,RO = RO, load_set = True)

# for state in ['Z','mZ','Y','mY','X','mX']:
#     for RO in range(6):
#         three_QEC.QEC_plot_single_state_RO(date = '20141121', no_error = '11',state = state,RO = RO, load_set = True)

# for state in ['Z','mZ','Y','mY','X','mX']:
#     for RO in range(6):
#         three_QEC.QEC_plot_single_state_RO(date = '20141122', no_error = '01',state = state,RO = RO, load_set = True)

# for state in ['Z','mZ','Y','mY','X','mX']:
#     for RO in range(6):
#         three_QEC.QEC_plot_single_state_RO(date = '20141124', no_error = '10_!',state = state,RO = RO, load_set = True)

''' plot process fidelity figures '''

three_QEC.QEC_plot_process_fids(date = '20141120',no_error = '00')

three_QEC.QEC_plot_process_fids(date = '20141121',no_error = '11')

three_QEC.QEC_plot_process_fids(date = '20141122',no_error = '01')

three_QEC.QEC_plot_process_fids(date = '20141124',no_error = '10_!')


''' plot 3 qubit and toffoli fidelity figures '''

# three_QEC.plot_QEC_sum_fidelities(date = '20141120',state = 'Z',no_error = '00')

# three_QEC.plot_QEC_sum_fidelities(date = '20141121',state = 'Z',no_error = '11')

# three_QEC.plot_QEC_sum_fidelities(date = '20141122',state = 'Z',no_error = '01')

# three_QEC.plot_QEC_sum_fidelities(date = '20141124',state = 'Z',no_error = '10_!')

''' test plots '''



