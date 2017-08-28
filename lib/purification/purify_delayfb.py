'''
Analyzes all non PQ data of the purification project.
All data are analyzed according to the MBI analysis class with varying input names
2016 NK
'''

import numpy as np
import os
from analysis.lib.tools import toolbox, plot;
from analysis.lib.m2.ssro import mbi, sequence, pqsequence
from matplotlib import pyplot as plt
# import matplotlib as mpl
from analysis.lib.fitting import fit, common
import copy as cp

reload(fit)
reload(mbi)
reload(common)
reload(toolbox)

CR_after_check = True  # global variable that let's us post select whether or not the NV was ionized

class PurificationDelayFBAnalysis(mbi.MBIAnalysis):
    max_nuclei = 6

    def get_phase_errors(self, name='', return_aggregate_results=True, carbon_idx=0, all_carbons_saved=True, **kw):
        agrp = self.adwingrp(name)
        reps = agrp.attrs['reps_per_ROsequence']
        pts = agrp.attrs['sweep_length']
        sweep_pts = agrp.attrs['sweep_pts']
        sweep_name = agrp.attrs['sweep_name']

        if all_carbons_saved:
            feedback_delay_cycles = agrp['feedback_delay_cycles'].value.reshape((self.max_nuclei, pts, reps), order='F')[carbon_idx]
            input_phases = agrp['compensated_phase'].value.reshape((self.max_nuclei,pts, reps), order='F')[carbon_idx]
        else:
            feedback_delay_cycles = agrp['feedback_delay_cycles'].value.reshape((pts, reps), order='F')
            input_phases = agrp['compensated_phase'].value.reshape((pts, reps), order='F')


        delay_feedback_N = agrp.attrs['delay_feedback_N']
        delay_clock_cycle_time = agrp.attrs['delay_clock_cycle_time']
        delay_feedback_target_phase = agrp.attrs['delay_feedback_target_phase']
        delay_time_offset = agrp.attrs['delay_time_offset']

        nuclear_frequency = agrp.attrs['nuclear_frequencies'][carbon_idx]
        # print "Carbon idx: ", carbon_idx

        feedback_delay_times = 2.0 * delay_feedback_N * (feedback_delay_cycles * delay_clock_cycle_time + delay_time_offset)

        feedback_phases = feedback_delay_times * nuclear_frequency * 360.0

        final_phases = feedback_phases + input_phases

        final_phase_errors = delay_feedback_target_phase - final_phases

        if return_aggregate_results:
            x = sweep_pts
            y = np.mean(final_phase_errors, axis=1)
            y_u = np.std(final_phase_errors, axis=1)

            return x, y, y_u
        else:
            return final_phase_errors

class PurificationDelayFBPQAnalysis(PurificationDelayFBAnalysis, pqsequence.PQSequenceAnalysis):
    def __init__(self, folder, **kw):
        PurificationDelayFBAnalysis.__init__(self, folder, **kw)
        self.pqf = self.f

    def select_dataset(self, name):
        # the pq_device attribute gets prepended to all group selections, perfect for us
        self.pq_device = "pq_data/" + name
        self.selected_dataset = name

    def extract_pulse_data(self, pulse_pts = None):
        sweep_pts = self.pts
        if pulse_pts is None:
            agrp = self.adwingrp(self.selected_dataset)
            pulse_pts = agrp.attrs['delay_feedback_N'] * 6

        self.pulse_pts = pulse_pts

        self.pulse_sweep_idxs = self.sweep_idxs.reshape((pulse_pts, sweep_pts, -1), order='F')

        self.pulse_sync_times = self.pqf[self.pq_device + '/PQ_sync_time-1'].value.reshape((pulse_pts, sweep_pts, -1), order='F')
        self.pulse_channels = self.pqf[self.pq_device + '/PQ_channel-1'].value.reshape((pulse_pts, sweep_pts, -1), order='F')
        self.pulse_sync_number = self.pqf[self.pq_device + '/PQ_sync_number-1'].value.reshape((pulse_pts, sweep_pts, -1),
                                                                                   order='F')


def get_tstamp_from_folder(folder):
    measurementstring = os.path.split(folder)[1]
    timestamp = os.path.split(os.path.split(folder)[0])[1] \
                + '/' + measurementstring[:6]
    return timestamp


def create_plot(f, **kw):
    ylabel = kw.pop('ylabel', None)
    xlabel = kw.pop('xlabel', None)
    title = kw.pop('title', None)

    fig = plt.figure()
    ax = plt.subplot(111)

    if xlabel != None:
        plt.xlabel(xlabel)
    if ylabel != None:
        plt.ylabel(ylabel)
    else:
        plt.ylabel('Contrast')
    if title != None:
        plt.title(title + ' ' + get_tstamp_from_folder(f))
    else:
        plt.title(get_tstamp_from_folder(f))

    return fig, ax


def save_and_close_plot(f, show=True):
    plt.savefig(os.path.join(f, 'Results.pdf'), format='pdf')
    plt.savefig(os.path.join(f, 'Results.png'), format='png')
    if show:
        plt.show()
    plt.close('all')


def plot_data(x, y, **kw):
    label = kw.pop('label', None)
    y_u = kw.pop('y_u', None)
    if not y_u is None:
        plt.errorbar(x, y, y_u, fmt='o', label=label, **kw)
    else:
        plt.plot(x, y)


def quadratic_addition(X, Y, X_u, Y_u):
    '''
    takes RO results and averages them according to theiry geometric average.
    inputs are numpy arrays and are assumed to have the same shape
    '''

    res = np.sqrt(Y ** 2 + X ** 2)
    res_u = np.sqrt((X * X_u) ** 2 + (Y * Y_u) ** 2) / np.sqrt((X ** 2 + Y ** 2))

    return res, res_u


def get_pos_neg_data(a, adwindata_str='', ro_array=['positive', 'negative'], **kw):
    '''
    Input: a : a data object of the class MBIAnalysis
    averages positive and negative data
    returns the sweep points, measured contrast and uncertainty
    '''

    ### get SSRO
    ssro_calib_folder = kw.pop('ssro_calib_folder', None)
    if ssro_calib_folder is None:
        use_preceding_ssro_calib = kw.pop('use_preceding_ssro_calib', False)
        if use_preceding_ssro_calib:
            kw['older_than'] = a.timestamp.replace('/', '')
            ssro_calib_folder = toolbox.latest_data('SSROCalib', **kw)
        else:
            ssro_calib_timestamp = kw.pop('ssro_calib_timestamp', None)

            if ssro_calib_timestamp == None:
                ssro_calib_folder = toolbox.latest_data('SSROCalib', **kw)  # , older_than = older_than)
            else:
                ssro_dstmp, ssro_tstmp = toolbox.verify_timestamp(ssro_calib_timestamp)
                ssro_calib_folder = toolbox.data_from_time(ssro_calib_timestamp)

    # if adwindata_str == '':
    #   return

    ##acquire pos_neg data
    for i, ro in enumerate(ro_array):
        a.get_sweep_pts()
        a.get_readout_results(name=adwindata_str + ro, CR_after_check=CR_after_check, **kw)
        a.get_electron_ROC(ssro_calib_folder, **kw)

        x_labels = a.sweep_pts
        if i == 0:
            res = ((a.p0.reshape(-1)) - 0.5) * 2
            res_u = 2 * a.u_p0.reshape(-1)

        else:
            y = ((a.p0.reshape(-1)) - 0.5) * 2  # Contrast
            y_u = 2 * a.u_p0.reshape(-1)  # contrast
            res = [y0 / 2 - y[ii] / 2 for ii, y0 in enumerate(res)]
            res_u = [np.sqrt(y0 ** 2 + y_u[ii] ** 2) / 2 for ii, y0 in enumerate(res_u)]

    return np.array(x_labels), np.array(res), np.array(res_u)


def average_repump_time(contains='', do_fit=False, **kw):
    '''
    gets data from a folder whose name contains the contains variable.
    Does or does not fit the data with a gaussian function
    '''

    ### kw for fitting

    fit_offset = kw.pop('fit_offset', 0)
    fit_amplitude = kw.pop('fit_amplitude', 1)
    fit_x0 = kw.pop('fit_x0', None)
    fit_sigma = kw.pop('fit_sigma', 0.2)
    fixed = kw.pop('fixed', [])
    show_guess = kw.pop('show_guess', False)

    ### folder choice
    if contains == '':
        contains = 'Sweep_Repump_time'
    elif len(contains) == 2:
        contains = 'Sweep_Repump_time' + contains
    elif len(contains) == 1:
        contains = 'Sweep_Repump_time_' + contains

    # older_than = kw.get('older_than',None) automatically handled by kws
    ### acquire data
    f = toolbox.latest_data(contains, **kw)
    a = mbi.MBIAnalysis(f)

    if '_Z' in f:
        x, y, y_u = get_pos_neg_data(a, adwindata_str='Z_', **kw)
        ylabel = 'Z'
    else:
        x, y1, y1_u = get_pos_neg_data(a, adwindata_str='X_', **kw)
        x2, y2, y2_u = get_pos_neg_data(a, adwindata_str='Y_', **kw)
        y, y_u = quadratic_addition(y1, y2, y1_u, y2_u)
        # y=y1
        # y_u = y1_u
        ylabel = 'Bloch vector length'

    ### create a plot
    xlabel = a.g.attrs['sweep_name']
    x = a.g.attrs['sweep_pts']  # could potentially be commented out?
    title = kw.pop("title", "avg repump time")
    fig, ax = create_plot(f, xlabel=xlabel, ylabel=ylabel, title=title)

    if fit_x0 is None:
        max_idx = np.argmax(y)
        fit_x0 = x[max_idx]

    ## plot data
    plot_data(x, y, y_u=y_u)

    ### fitting if you feel like it
    if do_fit:

        p0, fitfunc, fitfunc_str = common.fit_gauss(fit_offset, fit_amplitude, fit_x0, fit_sigma)

        if show_guess:
            # print decay
            ax.plot(np.linspace(x[0], x[-1], 201), fitfunc(np.linspace(x[0], x[-1], 201)), ':', lw=2)

        fit_result = fit.fit1d(x, y, None, p0=p0, fitfunc=fitfunc, do_print=True, fixed=fixed, ret=True)
        plot.plot_fit1d(fit_result, np.linspace(x[0], x[-1], 100), ax=ax, plot_data=False)

    ## save and close plot. We are done.
    save_and_close_plot(f)

    if kw.get('ret_data_fit'):
        return x, y, y_u, fit_result


def average_repump_time(contains='', do_fit=False, **kw):
    '''
    gets data from a folder whose name contains the contains variable.
    Does or does not fit the data with a gaussian function
    '''

    ### kw for fitting

    fit_offset = kw.pop('fit_offset', 0)
    fit_amplitude = kw.pop('fit_amplitude', 1)
    fit_x0 = kw.pop('fit_x0', None)
    fit_sigma = kw.pop('fit_sigma', 0.2)
    fixed = kw.pop('fixed', [])
    show_guess = kw.pop('show_guess', False)

    ### folder choice
    if contains == '':
        contains = 'Sweep_Repump_time'
    elif len(contains) == 2:
        contains = 'Sweep_Repump_time' + contains
    elif len(contains) == 1:
        contains = 'Sweep_Repump_time_' + contains

    # older_than = kw.get('older_than',None) automatically handled by kws
    ### acquire data
    f = toolbox.latest_data(contains, **kw)
    a = mbi.MBIAnalysis(f)

    if '_Z' in f:
        x, y, y_u = get_pos_neg_data(a, adwindata_str='Z_', **kw)
        ylabel = 'Z'
    else:
        x, y1, y1_u = get_pos_neg_data(a, adwindata_str='X_', **kw)
        x2, y2, y2_u = get_pos_neg_data(a, adwindata_str='Y_', **kw)
        y, y_u = quadratic_addition(y1, y2, y1_u, y2_u)
        # y=y1
        # y_u = y1_u
        ylabel = 'Bloch vector length'

    ### create a plot
    xlabel = a.g.attrs['sweep_name']
    x = a.g.attrs['sweep_pts']  # could potentially be commented out?
    title = kw.pop("title", "avg repump time")
    fig, ax = create_plot(f, xlabel=xlabel, ylabel=ylabel, title=title)

    if fit_x0 is None:
        max_idx = np.argmax(y)
        fit_x0 = x[max_idx]

    ## plot data
    plot_data(x, y, y_u=y_u)

    ### fitting if you feel like it
    if do_fit:

        p0, fitfunc, fitfunc_str = common.fit_gauss(fit_offset, fit_amplitude, fit_x0, fit_sigma)

        if show_guess:
            # print decay
            ax.plot(np.linspace(x[0], x[-1], 201), fitfunc(np.linspace(x[0], x[-1], 201)), ':', lw=2)

        fit_result = fit.fit1d(x, y, None, p0=p0, fitfunc=fitfunc, do_print=True, fixed=fixed, ret=True)
        plot.plot_fit1d(fit_result, np.linspace(x[0], x[-1], 100), ax=ax, plot_data=False)

        try:
            fit_result['carbon_id'] = a.g.attrs['carbons'][0]
        except:
            pass

    ## save and close plot. We are done.
    save_and_close_plot(f)

    if kw.get('ret_data_fit', False):
        return x, y, y_u, fit_result

    if kw.get('ret', False):
        return fit_result


def number_of_repetitions(contains='', do_fit=False, **kw):
    '''
    gets data from a folder whose name contains the contains variable.
    Does or does not fit the data with a gaussian function
    '''

    ### kw for fitting

    g_a = kw.pop('fit_a', 0)
    g_A = kw.pop('fit_A', 1)
    g_x0 = kw.pop('fit_x0', 0)
    g_T = kw.pop('fit_T', 500)
    g_n = kw.pop('fit_n', 1)
    g_f = kw.pop('fit_f', 0.0000)
    g_phi = kw.pop('fit_phi', 0)
    fixed = kw.pop('fixed', [])
    show_guess = kw.pop('show_guess', False)
    x_only = kw.pop('x_only', False)

    ### folder choice
    if contains == '':
        contains = '_sweep_number_of_reps'
    elif len(contains) == 2:
        contains = '_sweep_number_of_reps' + contains
    elif len(contains) == 1:
        contains = '_sweep_number_of_reps_' + contains

    # older_than = kw.get('older_than',None) automatically handled by kws
    ### acquire data
    f = toolbox.latest_data(contains, **kw)
    a = mbi.MBIAnalysis(f)

    if '_Z' in f and x_only == False:
        x, y, y_u = get_pos_neg_data(a, adwindata_str='Z_', **kw)
        ylabel = 'Z'
    else:
        x, y1, y1_u = get_pos_neg_data(a, adwindata_str='X_', **kw)
        if x_only:
            y = y1;
            y_u = y1_u;
            ylabel = 'X'
        else:
            x2, y2, y2_u = get_pos_neg_data(a, adwindata_str='Y_', **kw)
            y, y_u = quadratic_addition(y1, y2, y1_u, y2_u)
            ylabel = 'Bloch vector length'

    ### create a plot
    xlabel = a.g.attrs['sweep_name']
    x = a.g.attrs['sweep_pts']  # could potentially be commented out?
    fig, ax = create_plot(f, xlabel=xlabel, ylabel=ylabel, title='Number of repetitions')

    ## plot data
    plot_data(x, y, y_u=y_u)

    ### fitting if you feel like it
    if do_fit:

        p0, fitfunc, fitfunc_str = common.fit_exp_cos(g_a, g_A, g_x0, g_T, g_n, g_f, g_phi)

        if show_guess:
            # print decay
            ax.plot(np.linspace(x[0], x[-1], 201), fitfunc(np.linspace(x[0], x[-1], 201)), ':', lw=2)

        fit_result = fit.fit1d(x, y, None, p0=p0, fitfunc=fitfunc, do_print=True, fixed=fixed, ret=True)

        if isinstance(fit_result, int):
            print "Fit failed!"
        else:
            plot.plot_fit1d(fit_result, np.linspace(x[0], x[-1], 100), ax=ax, plot_data=False)

    ## save and close plot. We are done.
    save_and_close_plot(f)

    if kw.get('ret_data', False):
        return x, y, y_u

    if kw.get('ret', False):
        return fit_result

    if kw.get('ret_data_fit', False):
        return x, y, y_u, fit_result


def el_to_c_swap(contains='', input_el=['Z'], do_fit=False, **kw):
    '''
    gets data from a folder whose name contains the contains variable.
    Does or does not fit the data with a gaussian function
    '''

    ### folder choice
    if contains == '':
        contains = 'Swap_el_to_C'

    # older_than = kw.get('older_than',None) automatically handled by kws
    ### acquire data
    f = toolbox.latest_data(contains, **kw)
    a = mbi.MBIAnalysis(f)
    print 'this is the timestamp ', get_tstamp_from_folder(f)

    # data = np.empty([3,len(input_el)],dtype=str)
    data = []
    for i in range(len(input_el)):
        data.append([0, 0, 0])
    for ii, el in enumerate(input_el):
        # data.append([0,0,0])
        data_strings = []
        ro_str = 'el_state_' + el + '_'

        ro_array = ['positive', 'negative']
        x, y, y_u = get_pos_neg_data(a, adwindata_str=ro_str, ro_array=ro_array, **kw)
        y = np.round(y, decimals=2)
        y_u = np.round(y_u, decimals=2)
        # print el,x,y,y_u ### for debugging

        ### put output string together
        for jj, res, res_u in zip(range(3), y, y_u):
            data[ii][jj] = cp.deepcopy(str(res) + " +/- " + str(res_u))

    row_format = "{:>18}" * (len(x) + 1)
    headline_format = "{:>12}" + "{:>18}" * len(x)
    print headline_format.format("", *x)
    for el, row in zip(input_el, data):
        print "--------------------------------------------------------------------------------------------------"
        print row_format.format(el + ' |', *row)


def el_to_c_swap_success(contains='', input_el=['Z'], do_fit=False, **kw):
    '''
    gets data from a folder whose name contains the contains variable.
    Does or does not fit the data with a gaussian function
    '''

    ### folder choice
    if contains == '':
        contains = 'SwapSuccess_el_to_C'

    # older_than = kw.get('older_than',None) automatically handled by kws
    ### acquire data
    f = toolbox.latest_data(contains, **kw)
    a = mbi.MBIAnalysis(f)
    print 'this is the timestamp ', get_tstamp_from_folder(f)

    # data = np.empty([3,len(input_el)],dtype=str)
    data = []
    for i in range(len(input_el)):
        data.append([0, 0, 0])
    for ii, el in enumerate(input_el):
        # data.append([0,0,0])
        data_strings = []
        ro_str = 'el_state_' + el

        a.get_readout_results(name=ro_str, CR_after_check=CR_after_check)
        y = a.normalized_ssro[0]
        y_u = a.u_normalized_ssro[0]
        y = np.round(y, decimals=2)
        y_u = np.round(y_u, decimals=2)

        ### put output string together
        for jj, res, res_u in zip(range(3), y, y_u):
            data[ii][jj] = cp.deepcopy(str(res) + " +/- " + str(res_u))

    row_format = "{:>18}" * (2)
    for el, row in zip(input_el, data):
        print "--------------------------------------------------------------------------------------------------"
        print row_format.format(el + ' |', *row)


def phase_offset_after_LDE(contains='', **kw):
    '''
    gets data from a folder whose name contains the contains variable.
    Does or does not fit the data with a gaussian function
    '''

    ### folder choice
    if contains == '':
        contains = 'phase_offset_after_LDE'

    # older_than = kw.get('older_than',None) automatically handled by kws
    ### acquire data
    f = toolbox.latest_data(contains, **kw)
    a = mbi.MBIAnalysis(f)
    print 'this is the timestamp ', get_tstamp_from_folder(f)

    ro_array = ['positive', 'negative']
    x, y, y_u = get_pos_neg_data(a, adwindata_str='', ro_array=ro_array, **kw)
    y = np.round(y, decimals=2)
    y_u = np.round(y_u, decimals=2)
    # print el,x,y,y_u ### for debugging

    data = [0, 0, 0]
    ### put output string together
    for jj, res, res_u in zip(range(3), y, y_u):
        data[jj] = cp.deepcopy(str(res) + " +/- " + str(res_u))

    row_format = "{:>18}" * (len(x))

    headline_format = "{:>18}" * len(x)
    print headline_format.format(*x)
    print "-------------------------------------------------------------------------"
    print row_format.format(*data)


def calibrate_LDE_phase(contains='', do_fit=False, **kw):
    '''
    gets data from a folder whose name contains the contains variable.
    Does or does not fit the data with a gaussian function
    '''

    ### folder choice
    if contains == '':
        contains = 'LDE_phase_calibration'

    # tomography
    tomo = kw.pop('tomo_basis', 'X')
    # return fit
    ret = kw.get('ret', False)
    # for fitting
    freq = kw.pop('freq', None)
    decay = kw.pop('decay', 50)
    phi0 = kw.pop('phi0', 0)
    offset = kw.pop('offset', 0)
    A0 = kw.pop('A0', None)
    show_plot = kw.pop('show_plot', True)

    fixed = kw.pop('fixed', [1])
    show_guess = kw.pop('show_guess', False)

    # older_than = kw.get('older_than',None) automatically handled by kws
    ### acquire data
    if 'folder' in kw:
        f = kw.pop('folder')
    else:
        f = toolbox.latest_data(contains, **kw)
    a = mbi.MBIAnalysis(f)

    if freq is None:
        try:
            freq = a.g.attrs['phase_detuning'] / 360.0
        except:
            freq = 1./12. # voll auf die zwoelf.

    ro_array = ['positive', 'negative']
    # print ro_array
    if tomo == '':
        adwindata_str = tomo
    else:
        adwindata_str = tomo + '_'
    x, y, y_u = get_pos_neg_data(a, adwindata_str=adwindata_str, ro_array=ro_array, **kw)
    ylabel = tomo

    ### create a plot
    xlabel = a.g.attrs['sweep_name']
    x = a.g.attrs['sweep_pts']  # could potentially be commented out?
    fig, ax = create_plot(f, xlabel=xlabel, ylabel=ylabel, title='Acquired phase')

    ## plot data
    plot_data(x, y, y_u=y_u)

    ### fitting if you feel like it / still needs implementation
    if do_fit:
        if A0 is None:
            A0 = max(y)

        p0, fitfunc, fitfunc_str = common.fit_decaying_cos(freq, offset, A0, phi0, decay)

        if show_guess:
            # print decay
            ax.plot(np.linspace(x[0], x[-1], 201), fitfunc(np.linspace(x[0], x[-1], 201)), ':', lw=2)

        fit_result = fit.fit1d(x, y, None, p0=p0, fitfunc=fitfunc, do_print=True, ret=True, VERBOSE=True, fixed=fixed)

        plot.plot_fit1d(fit_result, np.linspace(x[0], x[-1], 1001), ax=ax, color='r', plot_data=False, add_txt=True,
                        lw=2)

        p_dict = fit_result['params_dict']
        e_dict = fit_result['error_dict']

        if 'carbons' in a.g.attrs:
            fit_result['carbon_id'] = a.g.attrs['carbons'][0]
        elif 'dps_carbons' in a.g.attrs:
            fit_result['carbon_id'] = a.g.attrs['dps_carbons'][0]
        elif 'carbon' in a.g.attrs:
            fit_result['carbon_id'] = a.g.attrs['carbon']
        else:
            fit_result['carbon_id'] = None

        if p_dict['A'] < 0:
            p_dict['phi'] = p_dict['phi'] + 180
            p_dict['A'] = p_dict['A'] * (-1)

        try:
            detuning = a.g.attrs['phase_detuning']
            fit_result['detuning'] = detuning
            fit_result['acq_phase_per_rep'] = 360 * (p_dict['f']) - detuning
            fit_result['u_acq_phase_per_rep'] = 360 * (e_dict['f'])
            print 'This is the phase detuning', detuning
            print 'Acquired phase per repetition (compensating for phase_detuning=) {:3.3f} +/- {:3.3f}'.format(
                round(fit_result['acq_phase_per_rep'], 3), round(fit_result['u_acq_phase_per_rep'], 3))
            print 'phase offset ', round(p_dict['phi'], 3)
        except:
            print 'no phase detuning found'
            ## save and close plot. We are done.

    save_and_close_plot(f, show=show_plot)

    if kw.get('ret', False):
        return fit_result

    if kw.get('ret_data', False):
        return x, y, y_u

    if kw.get('ret_fit_data', False):
        return fit_result, x, y, y_u

    if kw.get('ret_data_fit', False):
        return x, y, y_u, fit_result

def analyse_sequence_phase(contains='phase_fb_delayline', do_fit=False, **kw):
    '''
    gets data from a folder whose name contains the contains variable.
    Does or does not fit the data with a gaussian function
    '''

    # tomography
    tomo = kw.pop('tomo_basis', 'X')
    # return fit
    ret = kw.get('ret', False)
    # for fitting
    freq = kw.pop('freq', 1 / 12.)  # voll auf die zwoelf.
    decay = kw.pop('decay', 50)
    phi0 = kw.pop('phi0', 0)
    offset = kw.pop('offset', 0)
    A0 = kw.pop('A0', None)

    fixed = kw.pop('fixed', [1])
    show_guess = kw.pop('show_guess', False)

    # older_than = kw.get('older_than',None) automatically handled by kws
    ### acquire data
    f = toolbox.latest_data(contains, **kw)
    a = PurificationDelayFBAnalysis(f)

    ro_array = ['positive', 'negative']
    # print ro_array
    if tomo == '':
        adwindata_str = tomo
    else:
        adwindata_str = tomo + '_'
    x, y, y_u = get_pos_neg_data(a, adwindata_str=adwindata_str, ro_array=ro_array, **kw)
    ylabel = tomo

    ### create a plot
    xlabel = a.g.attrs['sweep_name']
    x = a.g.attrs['sweep_pts']  # could potentially be commented out?
    fig, ax = create_plot(f, xlabel=xlabel, ylabel=ylabel, title='Acquired phase')

    ## plot data
    plot_data(x, y, y_u=y_u)

    ### fitting if you feel like it / still needs implementation
    if do_fit:
        if A0 is None:
            A0 = max(y)

        _, phase_errors, _ = a.get_phase_errors(name=adwindata_str+ro_array[0], **kw)
        print phase_errors

        p0, fitfunc, fitfunc_str = common.fit_decaying_cos_with_phase_errors(freq, offset, A0, phi0, decay, phase_errors)

        if show_guess:
            # print decay
            ax.plot(np.linspace(x[0], x[-1], 201), fitfunc(np.linspace(x[0], x[-1], 201), phi_err=0.0), ':', lw=2)

        fit_result = fit.fit1d(x, y, None, p0=p0, fitfunc=fitfunc, do_print=True, ret=True, VERBOSE=True, fixed=fixed)

        plot.plot_fit1d(fit_result, np.linspace(x[0], x[-1], 1001), ax=ax, color='r', plot_data=False, add_txt=True,
                        lw=2, fitfunc_params={'phi_err': 0.0})

        plt.errorbar(x, fit_result['fitfunc'](x), y_u, fmt='x', label='fit with phase errors', color='r', zorder=10)

        p_dict = fit_result['params_dict']
        e_dict = fit_result['error_dict']

        if p_dict['A'] < 0:
            p_dict['phi'] = p_dict['phi'] + 180
            p_dict['A'] = p_dict['A'] * (-1)

        try:
            detuning = a.g.attrs['phase_detuning']
            print 'This is the phase detuning', detuning
            print 'Acquired phase per repetition (compensating for phase_detuning=) {:3.3f} +/- {:3.3f}'.format(
                round(360 * (p_dict['f']), 3) - detuning, round(360 * (e_dict['f']), 3))
            print 'phase offset {:3.3f} +/- {:3.3f}'.format(p_dict['phi'], e_dict['phi'])
        except:
            print 'no phase detuning found'
            ## save and close plot. We are done.

    save_and_close_plot(f)

    if kw.get('ret', False):
        return fit_result

    if kw.get('ret_data', False):
        return x, y, y_u

    if kw.get('ret_fit_data', False):
        return fit_result, x, y, y_u

def number_of_repetitions_stitched(contains='', do_fit=False, older_thans=None, multi_contains=None, **kw):
    '''
    gets data from a folder whose name contains the contains variable.
    Does or does not fit the data with a gaussian function
    '''

    ### kw for fitting

    g_a = kw.pop('fit_a', 0)
    g_A = kw.pop('fit_A', None)
    g_x0 = kw.pop('fit_x0', 0)
    g_T = kw.pop('fit_T', 500)
    g_n = kw.pop('fit_n', 1)
    g_f = kw.pop('fit_f', 0.0000)
    g_phi = kw.pop('fit_phi', 0)
    fixed = kw.pop('fixed', [])
    show_guess = kw.pop('show_guess', False)
    x_only = kw.pop('x_only', False)
    tomo_basis = kw.pop('tomo_basis', None)
    T2star_correction = kw.pop('T2star_correction', False)

    ### folder choice
    if contains == '':
        contains = '_sweep_number_of_reps'
    elif len(contains) == 2:
        contains = '_sweep_number_of_reps' + contains
    elif len(contains) == 1:
        contains = '_sweep_number_of_reps_' + contains

    # older_than = kw.get('older_than',None) automatically handled by kws
    ### acquire data
    multi_fs = []
    multi_as = []

    if older_thans is not None:
        for ot in older_thans:
            f = toolbox.latest_data(contains, older_than=ot, **kw)
            a = mbi.MBIAnalysis(f)

            multi_fs.append(f)
            multi_as.append(a)
    elif multi_contains is not None:
        for containy in multi_contains:
            f = toolbox.latest_data(containy, **kw)
            a = mbi.MBIAnalysis(f)

            multi_fs.append(f)
            multi_as.append(a)
    else:
        print("What do you want?")
        return

    x = np.array([])
    y = np.array([])
    y_u = np.array([])

    for a in multi_as:
        if tomo_basis is not None:
            x_new, y_new, y_u_new = get_pos_neg_data(a, adwindata_str=tomo_basis + '_', use_preceding_ssro_calib=True,
            **kw)
            ylabel = tomo_basis
        elif '_Z' in f and x_only == False:
            x_new, y_new, y_u_new = get_pos_neg_data(a, adwindata_str='Z_', use_preceding_ssro_calib=True, **kw)
            ylabel = 'Z'
        else:
            x_new, y1_new, y1_u_new = get_pos_neg_data(a, adwindata_str='X_', use_preceding_ssro_calib=True, **kw)
            if x_only:
                y_new = y1_new
                y_u_new = y1_u_new
                ylabel = 'X'
            else:
                x2_new, y2_new, y2_u_new = get_pos_neg_data(a, adwindata_str='Y_', use_preceding_ssro_calib=True, **kw)
                y_new, y_u_new = quadratic_addition(y1_new, y2_new, y1_u_new, y2_u_new)
                ylabel = 'Bloch vector length'

        x = np.append(x, x_new)
        y = np.append(y, y_new)
        y_u = np.append(y_u, y_u_new)

    if T2star_correction:
        from measurement.scripts.lt4_scripts.setup import msmt_params
        reload(msmt_params)
        import scipy.stats
        carbons = multi_as[0].g.attrs['carbons']

        sample = msmt_params.cfg['samples']['current']
        m_params = msmt_params.cfg['samples'][sample]
        e_trans = m_params['electron_transition']

        T2stars = []
        for c in carbons:
            T2star_0 = m_params['C%d_T2star_0' % c]
            T2star_1 = m_params['C%d_T2star_1%s' % (c, e_trans)]

            T2star_arr = np.array([T2star_0, T2star_1])

            # take the sqrt of the harmonic mean of the squared T2stars
            # I hope this is correct
            T2star = np.sqrt(scipy.stats.hmean(T2star_arr ** 2))
            print("avg. T2* for C%d: %.2f" % (c, T2star * 1e3))

            T2stars += [T2star]

        T2stars = np.array(T2stars)

        if kw.pop("T2star_gate_duration", True):
            print("Taking gate duration into account for T2* correction")
            T2star_envelope = np.ones_like(y)

            LDE_duration = kw.pop('LDE_element_length')
            LDE_total_duration = x*LDE_duration

            m_attrs = multi_as[0].g.attrs

            sequence_type = m_attrs['sequence_type']
            gates_per_carbon = {
                'serial_swap':  3,  # 2 swap gates and 1 tomo gate
                'MBE':          2,  # 1 MBE gate and 1 tomo gate
            }

            gates_duration = 0.0
            for i_c, c in enumerate(carbons):
                sequence_duration = LDE_total_duration + gates_duration
                T2star_envelope *= np.exp(-(sequence_duration / T2stars[i_c]) ** 2)
                gate_tau = m_attrs['C%d_Ren_tau%s' % (c, e_trans)][0]
                gate_N = m_attrs['C%d_Ren_N%s' % (c, e_trans)][0]
                gates_duration += gates_per_carbon[sequence_type] * 2.0 * gate_N * gate_tau

        else:
            combined_T2star = 1. / np.sqrt( np.sum(1./T2stars**2) )
            print("combined T2* for carbons %s: %.2f" % (str(carbons), combined_T2star * 1e3))

            def T2star_envelope(t, T2star):
                return np.exp(-(t / T2star) ** 2)

            LDE_duration = kw.pop('LDE_element_length')
            sequence_duration = x*LDE_duration

            T2star_envelope = np.exp(-(sequence_duration / T2star) ** 2)

        # let's cut out data points with T2star correction of more than 5 times
        valid_pts = T2star_envelope > 0.4

        y = y / T2star_envelope
        y_u = y_u / T2star_envelope

        x = x[valid_pts]
        y = y[valid_pts]
        y_u = y_u[valid_pts]

        print("valid pts: %d" % np.sum(valid_pts))



    ### create a plot
    xlabel = multi_as[0].g.attrs['sweep_name']
    fig, ax = create_plot(f, xlabel=xlabel, ylabel=ylabel, title='Number of repetitions')

    if g_A is None:
        min_x_pos = np.argmin(x)
        g_A = y[min_x_pos]
        print("Starting amplitude: %.3f" % g_A)
    ## plot data
    plot_data(x, y, y_u=y_u)

    ### fitting if you feel like it
    if do_fit:

        p0, fitfunc, fitfunc_str = common.fit_exp_cos(g_a, g_A, g_x0, g_T, g_n, g_f, g_phi)

        if show_guess:
            # print decay
            ax.plot(np.linspace(x[0], x[-1], 201), fitfunc(np.linspace(x[0], x[-1], 201)), ':', lw=2)

        fit_result = fit.fit1d(x, y, None, p0=p0, fitfunc=fitfunc, do_print=True, fixed=fixed, ret=True)

        if isinstance(fit_result, int):
            print "Fit failed!"
        else:
            plot.plot_fit1d(fit_result, np.linspace(x[0], x[-1], 100), ax=ax, plot_data=False)

    ## save and close plot. We are done.
    save_and_close_plot(f)

    if kw.get('ret_data', False):
        return x, y, y_u

    if kw.get('ret', False):
        return fit_result

    if kw.get('ret_data_fit', False):
        return x, y, y_u, fit_result


def analyse_delay_feedback_phase_error(contains='', name='ssro', **kw):
    if contains == '':
        contains = 'fb_delayline'

    ### acquire data
    f = toolbox.latest_data(contains, **kw)
    a = PurificationDelayFBAnalysis(f)

    x, y, y_u = a.get_phase_errors(name=name, **kw)

    label = kw.pop('label', None)

    ## plot data
    plot_data(x, y, y_u=y_u, label=label)

    plt.axhline(y=0.0)

    xlabel = a.g.attrs['sweep_name']
    ylabel = "phase error after feedback"

    plt.xlabel(xlabel)
    plt.ylabel(ylabel)

    if kw.get('ret_data', False):
        return x, y, y_u

def repump_speed(contains='repump_speed', name='adwindata', do_fit=False, **kw):
    show_guess = kw.pop('show_guess', False)
    g_a = kw.pop('offset', 0.0)
    g_A = kw.pop('A', 1.0)
    g_tau = kw.pop('tau', 100.0)
    fixed = kw.pop('fixed', [0])

    ### acquire data
    f = toolbox.latest_data(contains, **kw)
    a = mbi.MBIAnalysis(f)
    a.get_sweep_pts()
    a.get_readout_results(name=name)
    a.get_electron_ROC(**kw)

    x = a.sweep_pts.reshape(-1)
    y = 1.0 - a.p0.reshape(-1)
    y_u = a.u_p0.reshape(-1)

    xlabel = a.sweep_name
    ylabel = r'$1 - p(|0\rangle)$'

    fig, ax = create_plot(f, xlabel=xlabel, ylabel=ylabel, title='repump speed')

    plot_data(x, y, y_u=y_u)

    if do_fit:

        p0, fitfunc, fitfunc_str = common.fit_exp_decay_with_offset(g_a, g_A, g_tau)

        if show_guess:
            # print decay
            ax.plot(np.linspace(x[0], x[-1], 201), fitfunc(np.linspace(x[0], x[-1], 201)), ':', lw=2)

        fit_result = fit.fit1d(x, y, None, p0=p0, fitfunc=fitfunc, do_print=True, fixed=fixed, ret=True)

        if isinstance(fit_result, int):
            print "Fit failed!"
        else:
            plot.plot_fit1d(fit_result, np.linspace(x[0], x[-1], 100), ax=ax, plot_data=False)

    save_and_close_plot(f)

    if kw.get('ret_data', False):
        return x, y, y_u

    if kw.get('ret', False):
        return fit_result

def tomo_analysis(contains="tomo", name="", **kw):
    # older_than = kw.get('older_than',None) automatically handled by kws
    ### acquire data
    f = toolbox.latest_data(contains, **kw)
    a = mbi.MBIAnalysis(f)

    ro_array = ['positive', 'negative']
    x, y, y_u = get_pos_neg_data(a, adwindata_str=name, ro_array=ro_array, **kw)
    x = ["".join(bases) for bases in x]

    xlabel = a.g.attrs['sweep_name']
    ylabel = "Expectation value"

    fig, ax = create_plot(f, xlabel=xlabel, ylabel=ylabel, title='Delay feedback + tomography')

    fake_x = np.arange(len(x))
    rects = ax.bar(fake_x, y, yerr=y_u, align='center', ecolor='k')
    ax.set_xticks(fake_x)
    ax.set_xticklabels(x)
    ax.set_ylim(-1.0, 1.0)

    def autolabel(rects):
        for ii, rect in enumerate(rects):
            height = rect.get_height()
            plt.text(rect.get_x() + rect.get_width() / 2., 1.02 * height,
                     str(round(y[ii], 2)) + '(' + str(int(round(y_u[ii] * 100))) + ')',
                     ha='center', va='bottom')
    autolabel(rects)

    save_and_close_plot(f)

def calibrate_LDE_phase_stitched(contains='', do_fit=False, older_thans=None, multi_contains=None, multi_folders=None,
                                                                                                                 **kw):
    '''
    gets data from a folder whose name contains the contains variable.
    Does or does not fit the data with a gaussian function
    '''

    ### folder choice
    if contains == '':
        contains = 'LDE_phase_calibration'

    # tomography
    tomo = kw.pop('tomo_basis', 'X')
    # return fit
    ret = kw.get('ret', False)
    # for fitting
    freq = kw.pop('freq', None)
    decay = kw.pop('decay', 50)
    phi0 = kw.pop('phi0', 0)
    offset = kw.pop('offset', 0)
    A0 = kw.pop('A0', None)
    show_plot = kw.pop('show_plot', True)

    fixed = kw.pop('fixed', [1])
    show_guess = kw.pop('show_guess', False)

    if freq is None:
        try:
            freq = a.g.attrs['phase_detuning'] / 360.0
        except:
            freq = 1./12. # voll auf die zwoelf.

    ro_array = ['positive', 'negative']

    # older_than = kw.get('older_than',None) automatically handled by kws
    ### acquire data
    multi_fs = []
    multi_as = []

    if older_thans is not None:
        for ot in older_thans:
            f = toolbox.latest_data(contains, older_than=ot, **kw)
            multi_fs.append(f)
    elif multi_contains is not None:
        for containy in multi_contains:
            f = toolbox.latest_data(containy, **kw)
            multi_fs.append(f)
    elif multi_folders is not None:
        multi_fs.extend(multi_folders)
    else:
        print("What do you want?")
        return

    for f in multi_fs:
        a = mbi.MBIAnalysis(f)
        multi_as.append(a)

    x = np.array([])
    y = np.array([])
    y_u = np.array([])

    for a in multi_as:
        x_new, y_new, y_u_new = get_pos_neg_data(a, adwindata_str=tomo + '_', use_preceding_ssro_calib=True,
                                                 **kw)
        ylabel = tomo

        x = np.append(x, x_new)
        y = np.append(y, y_new)
        y_u = np.append(y_u, y_u_new)

    ### create a plot
    xlabel = a.g.attrs['sweep_name']
    # x = a.g.attrs['sweep_pts']  # could potentially be commented out?
    fig, ax = create_plot(f, xlabel=xlabel, ylabel=ylabel, title='Acquired phase')

    ## plot data
    plot_data(x, y, y_u=y_u)

    ### fitting if you feel like it / still needs implementation
    if do_fit:
        if A0 is None:
            A0 = 1.2*np.max(np.abs(y))

        p0, fitfunc, fitfunc_str = common.fit_decaying_cos(freq, offset, A0, phi0, decay)

        if show_guess:
            # print decay
            ax.plot(np.linspace(x[0], x[-1], 201), fitfunc(np.linspace(x[0], x[-1], 201)), ':', lw=2)

        fit_result = fit.fit1d(x, y, None, p0=p0, fitfunc=fitfunc, do_print=True, ret=True, VERBOSE=True, fixed=fixed)

        plot.plot_fit1d(fit_result, np.linspace(x[0], x[-1], 1001), ax=ax, color='r', plot_data=False, add_txt=True,
                        lw=2)

        p_dict = fit_result['params_dict']
        e_dict = fit_result['error_dict']

        fit_result['carbon_id'] = a.g.attrs['carbons'][0]

        if p_dict['A'] < 0:
            p_dict['phi'] = p_dict['phi'] + 180
            p_dict['A'] = p_dict['A'] * (-1)

        try:
            detuning = a.g.attrs['phase_detuning']
            fit_result['detuning'] = detuning
            fit_result['acq_phase_per_rep'] = 360 * (p_dict['f']) - detuning
            fit_result['u_acq_phase_per_rep'] = 360 * (e_dict['f'])
            print 'This is the phase detuning', detuning
            print 'Acquired phase per repetition (compensating for phase_detuning=) {:3.3f} +/- {:3.3f}'.format(
                round(fit_result['acq_phase_per_rep'], 3), round(fit_result['u_acq_phase_per_rep'], 3))
            print 'phase offset ', round(p_dict['phi'], 3)
        except:
            print 'no phase detuning found'
            ## save and close plot. We are done.

    save_and_close_plot(f, show=show_plot)

    if kw.get('ret', False):
        return fit_result

    if kw.get('ret_data', False):
        return x, y, y_u

    if kw.get('ret_fit_data', False):
        return fit_result, x, y, y_u

    if kw.get('ret_data_fit', False):
        return x, y, y_u, fit_result
