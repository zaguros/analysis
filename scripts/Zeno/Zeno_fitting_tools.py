"""
Fitting tools for the Zeno experiment.
"""
import numpy as np
import os
from analysis.lib.tools import toolbox; reload(toolbox)
from analysis.lib.tools import plot; reload(plot)
from analysis.lib.fitting import fit, common, ramsey;reload(common); reload(fit)
from analysis.lib.m2.ssro import mbi
from matplotlib import pyplot as plt


#######################################

#       PROCESS FIDELITY FITS

#######################################


def fit_1msmt_proc_fid(g_A0, g_offset0, g_t, g_p):
    '''

    g_offset0   -  Offset, derived from the 0 msmt case. Fixed.
    g_A0        -  Amplitude derived from the 0 msmt case. Fixed.
    g_p         -  Probability for faulty parity measurement
    g_t         -  Decay of Ramsey. Fixed.

    The function should take the data set for 1 Zeno-measurement and fit the process fidelity with a prederived function. That depends on one parameter. The fidelity of the parity measurement.
    '''

    fitfunc_str = '''analyitcal solution'''

    ### Parameters

    A0           = fit.Parameter(g_A0, 'A0')
    offset0      = fit.Parameter(g_offset0, 'offset0')

        ### Ramsey
    t   = fit.Parameter(g_t, 't')

    ### Zeno
    p   = fit.Parameter(g_p,'p')

    p0 = [A0,offset0,t,p]

    def fitfunc(x):

        decay = np.exp(- ( x**2 /(2. * (t()**2) ) ) )
        Hahn = np.exp(-(x/604.)**1.6162)
        Fx = -0.25*decay*(-1.+p())-0.25*(-1.+p())*Hahn + 0.5
        Fz0 = 2.*offset0()
        Con0 = 2.*(Fz0-0.5)
        Conp = Con0*(1-p())
        Fzp = Conp/2. + 0.5

        # am = (2*A0()+1)/2.
        # return (2.*am*Fx+Fzp-1.)/2.
        ### commented out because apparently wrong.
        Cx = 2*Fx-1.
        Cxred = 2*A0()*Cx
        Fx = (Cxred+1)/2.
        return (2.*Fx+Fzp-1.)/2.

    return p0, fitfunc, fitfunc_str



def fit_2msmt_proc_fid(g_A0, g_offset0, g_t, g_p):
    '''

    g_offset0   -  Offset, derived from the 0 msmt case. Fixed.
    g_A0        -  Amplitude derived from the 0 msmt case. Fixed.
    g_p         -  Probability for faulty parity measurement
    g_t         -  Decay of Ramsey. Fixed.

    The function should take the data set for 2 Zeno-measurements and fit the process fidelity with a prederived function. That depends on one parameter. The fidelity of the parity measurement.
    '''

    fitfunc_str = '''analyitcal solution'''

    ### Parameters

    A0           = fit.Parameter(g_A0, 'A0')
    offset0      = fit.Parameter(g_offset0, 'offset0')

        ### Ramsey
    t   = fit.Parameter(g_t, 't')

    ### Zeno
    p   = fit.Parameter(g_p,'p')

    p0 = [A0,offset0,t,p]

    def fitfunc(x):

        decay = 0.125*np.exp(- ( x**2 /(2. * (t()**2) ) ) )
        decay2 = np.exp((4./9.)*( x**2 /((t()**2) ) ) )
        Fx = 0.5+decay*(1.+3*decay2)*((p()-1.)**2)
        Fz0 = 2.*offset0()
        Con0 = 2.*(Fz0-0.5)
        Conp = Con0*(1-p())**2
        Fzp = Conp/2. + 0.5
        # am = (2*(A0()+1)/2.
        # return (2.*am*Fx+Fzp-1.)/2.
        Cx = 2*Fx-1.
        Cxred = 2*A0()*Cx
        Fx = (Cxred+1)/2.
        return (2.*Fx+Fzp-1.)/2.

    return p0, fitfunc, fitfunc_str


def fit_3msmt_proc_fid(g_A0, g_offset0, g_t, g_p):
    '''

    g_offset0   -  Offset, derived from the 0 msmt case. Fixed.
    g_A0        -  Amplitude derived from the 0 msmt case. Fixed.
    g_p         -  Probability for faulty parity measurement
    g_t         -  Decay of Ramsey. Fixed.

    The function should take the data set for 3 Zeno-measurements and fit the process fidelity with a prederived function. That depends on one parameter. The fidelity of the parity measurement.
    '''

    fitfunc_str = '''analyitcal solution'''

    ### Parameters

    A0           = fit.Parameter(g_A0, 'A0')
    offset0      = fit.Parameter(g_offset0, 'offset0')

        ### Ramsey
    t   = fit.Parameter(g_t, 't')

    ### Zeno
    p   = fit.Parameter(g_p,'p')

    p0 = [A0,offset0,t,p]

    def fitfunc(x):

        decay = -0.0625*np.exp(- ( x**2 /(2. * (t()**2) ) ) )

        decay2 = np.exp((3./8.)*( x**2 /((t()**2) ) ) )

        Hahn = np.exp(-(x/604.)**1.6162)

        Fx = -Hahn*((-11.+3*p()*(3+(3-p())*p()))/16.+0.5)+decay*((p()-1.)**3+4.*decay2*(p()-1.)**3)+0.5

        Fz0 = 2.*offset0()
        Con0 = 2.*(Fz0-0.5)
        Conp = Con0*(1-p())**3
        Fzp = Conp/2. + 0.5

        # am = (2*A0()+1)/2.
        # return (2.*am*Fx+Fzp-1.)/2.
        ### commented out because apparently wrong.
        Cx = 2*Fx-1.
        Cxred = 2*A0()*Cx
        Fx = (Cxred+1)/2.
        return (2.*Fx+Fzp-1.)/2.

    return p0, fitfunc, fitfunc_str


def fit_4msmt_proc_fid(g_A0, g_offset0, g_t, g_p):
    '''

    g_offset0   -  Offset, derived from the 0 msmt case. Fixed.
    g_A0        -  Amplitude derived from the 0 msmt case. Fixed.
    g_p         -  Probability for faulty parity measurement
    g_t         -  Decay of Ramsey. Fixed.

    The function should take the data set for 4 Zeno-measurements and fit the process fidelity with a prederived function. That depends on one parameter. The fidelity of the parity measurement.
    '''

    fitfunc_str = '''analyitcal solution'''

    ### Parameters

    A0           = fit.Parameter(g_A0, 'A0')
    offset0      = fit.Parameter(g_offset0, 'offset0')

        ### Ramsey
    t   = fit.Parameter(g_t, 't')

    ### Zeno
    p   = fit.Parameter(g_p,'p')

    p0 = [A0,offset0,t,p]

    def fitfunc(x):

        decay = np.exp(- ( x**2 /(2. * (t()**2) ) ) )/32.

        decay2 = 5.*np.exp((8./25.)*( x**2 /((t()**2) ) ) )*(p()-1.)**4

        decay3 = 10.*np.exp((12./25.)*( x**2 /((t()**2) ) ) )*(p()-1.)**4


        Fx = 0.5+decay*((p()-1)**4+decay2+decay3)

        Fz0 = 2.*offset0()
        Con0 = 2.*(Fz0-0.5)
        Conp = Con0*(1-p())**4
        Fzp = Conp/2. + 0.5

        # am = (2*A0()+1)/2.
        # return (2.*am*Fx+Fzp-1.)/2.
        ### commented out because apparently wrong.
        Cx = 2*Fx-1.
        Cxred = 2*A0()*Cx
        Fx = (Cxred+1)/2.
        return (2.*Fx+Fzp-1.)/2.

    return p0, fitfunc, fitfunc_str

def fit_5msmt_proc_fid(g_A0, g_offset0, g_t, g_p):
    '''

    g_offset0   -  Offset, derived from the 0 msmt case. Fixed.
    g_A0        -  Amplitude derived from the 0 msmt case. Fixed.
    g_p         -  Probability for faulty parity measurement
    g_t         -  Decay of Ramsey. Fixed.

    The function should take the data set for 5 Zeno-measurements and fit the process fidelity with a prederived function. That depends on one parameter. The fidelity of the parity measurement.
    '''

    fitfunc_str = '''analyitcal solution'''

    ### Parameters

    A0           = fit.Parameter(g_A0, 'A0')
    offset0      = fit.Parameter(g_offset0, 'offset0')

        ### Ramsey
    t   = fit.Parameter(g_t, 't')

    ### Zeno
    p   = fit.Parameter(g_p,'p')

    p0 = [A0,offset0,t,p]

    def fitfunc(x):

        decay = -np.exp(- ( x**2 /(2. * (t()**2) ) ) )/64.

        decay2 = 6.*np.exp((5./18.)*( x**2 /((t()**2) ) ) )

        decay3 = 15.*np.exp((4./9.)*( x**2 /((t()**2) ) ) )*(p()-1.)**5

        bracket1 = 10+(-5+p())*p()
        bracket2 = (bracket1 * p() -10)*p() + 5
        bracket3 = (bracket2*5*p()-21)*2./64.

        Hahn  = np.exp(-(x/604.)**1.6162)

        #### decaying parts + the constant part which arises from the measurement echo.

        Fx = decay*((p()-1)**5*(1+decay2)+decay3) - (bracket3+0.5)*Hahn+0.5

        Fz0 = 2.*offset0()
        Con0 = 2.*(Fz0-0.5)
        Conp = Con0*(1-p())**5
        Fzp = Conp/2. + 0.5

        # am = (2*A0()+1)/2.
        # return (2.*am*Fx+Fzp-1.)/2.
        ### commented out because apparently wrong.
        Cx = 2*Fx-1.
        Cxred = 2*A0()*Cx
        Fx = (Cxred+1)/2.
        return (2.*Fx+Fzp-1.)/2.

    return p0, fitfunc, fitfunc_str


def fit_6msmt_proc_fid(g_A0, g_offset0, g_t, g_p):
    '''

    g_offset0   -  Offset, derived from the 0 msmt case. Fixed.
    g_A0        -  Amplitude derived from the 0 msmt case. Fixed.
    g_p         -  Probability for faulty parity measurement
    g_t         -  Decay of Ramsey. Fixed.

    The function should take the data set for 6 Zeno-measurements and fit the process fidelity with a prederived function. That depends on one parameter. The fidelity of the parity measurement.
    '''

    fitfunc_str = '''analyitcal solution'''

    ### Parameters

    A0           = fit.Parameter(g_A0, 'A0')
    offset0      = fit.Parameter(g_offset0, 'offset0')

    ### Ramsey (decay time divided by sqrt(2).)
    t   = fit.Parameter(g_t, 't')

    ### Zeno
    p   = fit.Parameter(g_p,'p')

    p0 = [A0,offset0,t,p]

    def fitfunc(x):

        decay = np.exp(- ( x**2 /(2. * (t()**2) ) ) )/128.

        decay2 = 7. * np.exp((12./49.)*( x**2 /((t()**2) ) ) )*(p()-1.)**6

        decay3 = 21. * np.exp((20./49.)*( x**2 /((t()**2) ) ) )*(p()-1.)**6

        decay4 = 35. * np.exp((24./49.)*( x**2 /((t()**2) ) ) )*(p()-1.)**6

        #### decay to 0.5 only for 6 measurements no other constant parts of the function due to echo effects.
        
        Fx = 0.5 + decay*((p()-1)**6 + decay2 + decay3 + decay4)

        Fz0 = 2.*offset0()
        Con0 = 2.*(Fz0-0.5) ### contrast of the ZZ correlations for 0 measurements.
        Conp = Con0*(1-p())**6 ### after 6 measurements the contrast is reduced by (1-p)^6
        Fzp = Conp/2. + 0.5 ### calculate the fidelity.

        # am = (2*A0()+1)/2.
        # return (2.*am*Fx+Fzp-1.)/2.
        ### commented out because apparently wrong.
        Cx = 2*Fx-1.
        Cxred = 2*A0()*Cx
        Fx = (Cxred+1)/2.
        return (2.*Fx+Fzp-1.)/2.

        ### notice that the fit does not distinguish between logical Z and logical Y states. Both are assumed to decay in the same way.
        ### we also do not distinguish between orthogonal states
        ### the formula for the process fidelity is usually given by Fproc = (3*avgStateFid-1)/2
        ### here it becomes Fproc = (Fx+Fy+Fz-1)/2 = (Fx*2+Fz-1)/2


    return p0, fitfunc, fitfunc_str


    #######################################

    #       STATE FIDELITY FITS

    #######################################

"""
All of the functions below underlie the same model as the process fidelity functions above
For these functions we are only interested in the state fidelity of one state.
Namely a state with X or YZ correlations which decays according to our number of Zeno measurements.
"""
def fit_1msmt_state_fid(g_A0, g_t1,g_t2, g_p,g_repump,take_repump):
    '''


    g_A0        -  Amplitude derived from the 0 msmt case. Fixed.
    g_p         -  Probability for faulty parity measurement
    g_repump    -  Probability that a repumping event destroys the state.
    g_t1        -  Decay of Ramsey (ancilla carbon). Fixed.
    g_t2        -  decay of the data carbon. Fixed.
    take_repump -  bool which tells the fit to take detrimental repumping into account or not.

    The function should take the data set for 1 Zeno-measurement and fit the process fidelity with a prederived function. That depends on one parameter. The fidelity of the parity measurement.
    '''

    fitfunc_str = '''analyitcal solution'''

    ### Parameters

    A0           = fit.Parameter(g_A0, 'A0')

        ### Ramsey
    t2   = fit.Parameter(g_t2, 't2')

    ### Zeno
    p   = fit.Parameter(g_p,'p')

    if take_repump:
        t1   = fit.Parameter(g_t1, 't1')
        repump = fit.Parameter(g_repump,'repump')
        p0 = [A0,t1,t2,p,repump]

    else:
        p0 = [A0,t2,p]


    def fitfunc(x):

        Hahn  = np.exp(-(x/604.)**1.6162)

        ### this fit function does take detrimental repumping into account.
        ### no hahn echo time...

        if take_repump:
            decay1 = np.exp(-x**2/(2*t2()**2))*(0.125+0.125*repump())
            decay2 = np.exp((-0.125/(t2()**2)-0.125/(t1()**2))*x**2)*(1-repump())/4.

            Fx = 0.5+p()*(0.125+decay1+decay2+0.125*repump())
            C  = A0()*2*(Fx-0.5)
            Fx = C/2. + 0.5

        # """
        # Taking repumping not into account.
        # """
        else:
            
            decay = np.exp(- ( x**2 /(2. * (t2()**2) ) ) )
            Fx = -0.25*decay*(-1.+p())-0.25*(-1.+p())*Hahn + 0.5
            Cx = A0()*2.*(Fx-0.5)
            Fx = Cx/2.+0.5
            # pmix = 1-A0()
            # Fx = 0.75+pmix*(-0.25+0.25*p())-0.25*p()+decay*(0.25+pmix*(-0.25+0.25*p())-0.25*p())

        return Fx

    return p0, fitfunc, fitfunc_str



def fit_2msmt_state_fid(g_A0, g_t, g_p):
    '''


    g_A0        -  Amplitude derived from the 0 msmt case. Fixed.
    g_p         -  Probability for faulty parity measurement
    g_t         -  Decay of Ramsey. Fixed.

    The function should take the data set for 2 Zeno-measurements and fit the process fidelity with a prederived function. That depends on one parameter. The fidelity of the parity measurement.
    '''

    fitfunc_str = '''analyitcal solution'''

    ### Parameters

    A0           = fit.Parameter(g_A0, 'A0')


        ### Ramsey
    t   = fit.Parameter(g_t, 't')

    ### Zeno
    p   = fit.Parameter(g_p,'p')

    p0 = [A0,t,p]

    def fitfunc(x):

        decay = 0.125*np.exp(- ( x**2 /(2. * (t()**2) ) ) )
        decay2 = np.exp((4./9.)*( x**2 /((t()**2) ) ) )
        Fx = 0.5+decay*((p()-1.)**2+3.*decay2*(p()-1.)**2)
        Cx = A0()*2.*(Fx-0.5)
        Fx = Cx/2.+0.5
        return Fx
    return p0, fitfunc, fitfunc_str


def fit_3msmt_state_fid(g_A0, g_t, g_p):
    '''


    g_A0        -  Amplitude derived from the 0 msmt case. Fixed.
    g_p         -  Probability for faulty parity measurement
    g_t         -  Decay of Ramsey. Fixed.

    The function should take the data set for 3 Zeno-measurements and fit the process fidelity with a prederived function. That depends on one parameter. The fidelity of the parity measurement.
    '''

    fitfunc_str = '''analyitcal solution'''

    ### Parameters

    A0           = fit.Parameter(g_A0, 'A0')

        ### Ramsey
    t   = fit.Parameter(g_t, 't')

    ### Zeno
    p   = fit.Parameter(g_p,'p')

    p0 = [A0,t,p]

    def fitfunc(x):

        Hahn  = np.exp(-(x/604.)**1.6162)

        decay = -0.0625*np.exp(- ( x**2 /(2. * (t()**2) ) ) )

        decay2 = np.exp((3./8.)*( x**2 /((t()**2) ) ) )

        Fx = - Hahn*((-11.+3*p()*(3+(3-p())*p()))/16.+0.5) + 0.5 + decay*((p()-1.)**3+4.*decay2*(p()-1.)**3)

        Cx = A0()*2.*(Fx-0.5)
        Fx = Cx/2.+0.5
        return Fx

    return p0, fitfunc, fitfunc_str


def fit_4msmt_state_fid(g_A0, g_t, g_p):
    '''


    g_A0        -  Amplitude derived from the 0 msmt case. Fixed.
    g_p         -  Probability for faulty parity measurement
    g_t         -  Decay of Ramsey. Fixed.

    The function should take the data set for 4 Zeno-measurements and fit the process fidelity with a prederived function. That depends on one parameter. The fidelity of the parity measurement.
    '''

    fitfunc_str = '''analyitcal solution'''

    ### Parameters

    A0           = fit.Parameter(g_A0, 'A0')


        ### Ramsey
    t   = fit.Parameter(g_t, 't')

    ### Zeno
    p   = fit.Parameter(g_p,'p')

    p0 = [A0,t,p]

    def fitfunc(x):

        decay = np.exp(- ( x**2 /(2. * (t()**2) ) ) )/32.

        decay2 = 5.*np.exp((8./25.)*( x**2 /((t()**2) ) ) )*(p()-1.)**4

        decay3 = 10.*np.exp((12./25.)*( x**2 /((t()**2) ) ) )*(p()-1.)**4

        Fx = 0.5+decay*((p()-1)**4+decay2+decay3)

        Cx = A0()*2.*(Fx-0.5)
        Fx = Cx/2.+0.5
        return Fx

    return p0, fitfunc, fitfunc_str

def fit_5msmt_state_fid(g_A0, g_t, g_p):
    '''


    g_A0        -  Amplitude derived from the 0 msmt case. Fixed.
    g_p         -  Probability for faulty parity measurement
    g_t         -  Decay of Ramsey. Fixed.

    The function should take the data set for 5 Zeno-measurements and fit the process fidelity with a prederived function. That depends on one parameter. The fidelity of the parity measurement.
    '''

    fitfunc_str = '''analyitcal solution'''

    ### Parameters

    A0           = fit.Parameter(g_A0, 'A0')


        ### Ramsey
    t   = fit.Parameter(g_t, 't')

    ### Zeno
    p   = fit.Parameter(g_p,'p')

    p0 = [A0,t,p]

    def fitfunc(x):

        Hahn  = np.exp(-(x/604.)**1.6162)

        decay = -np.exp(- ( x**2 /(2. * (t()**2) ) ) )/64.

        decay2 = 6.*np.exp((5./18.)*( x**2 /((t()**2) ) ) )

        decay3 = 15.*np.exp((4./9.)*( x**2 /((t()**2) ) ) )*(p()-1.)**5

        bracket1 = 10+(-5+p())*p()
        bracket2 = (bracket1 * p() -10)*p() + 5
        bracket3 = (bracket2*5*p()-21)*2./64.

        #### decaying parts + the constant part which arises from the measurement echo.

        Fx = decay*((p()-1)**5*(1+decay2)+decay3) - (bracket3+0.5)*Hahn + 0.5

        Cx = A0()*2.*(Fx-0.5)
        Fx = Cx/2.+0.5
        return Fx

    return p0, fitfunc, fitfunc_str


def fit_6msmt_state_fid(g_A0, g_t, g_p):
    '''


    g_A0        -  Amplitude derived from the 0 msmt case. Fixed.
    g_p         -  Probability for faulty parity measurement
    g_t         -  Decay of Ramsey. Fixed.

    The function should take the data set for 6 Zeno-measurements and fit the process fidelity with a prederived function. That depends on one parameter. The fidelity of the parity measurement.
    '''

    fitfunc_str = '''analyitcal solution'''

    ### Parameters

    A0           = fit.Parameter(g_A0, 'A0')


    ### Ramsey (decay time divided by sqrt(2).)
    t   = fit.Parameter(g_t, 't')

    ### Zeno
    p   = fit.Parameter(g_p,'p')

    p0 = [A0,t,p]

    def fitfunc(x):

        decay = np.exp(- ( x**2 /(2. * (t()**2) ) ) )/128.

        decay2 = 7. * np.exp((12./49.)*( x**2 /((t()**2) ) ) )*(p()-1.)**6

        decay3 = 21. * np.exp((20./49.)*( x**2 /((t()**2) ) ) )*(p()-1.)**6

        decay4 = 35. * np.exp((24./49.)*( x**2 /((t()**2) ) ) )*(p()-1.)**6

        #### decay to 0.5 only for 6 measurements no other constant parts of the function due to echo effects.
        
        Fx = 0.5 + decay*((p()-1)**6 + decay2 + decay3 + decay4)

        Cx = A0()*2.*(Fx-0.5)
        Fx = Cx/2.+0.5
        return Fx


    return p0, fitfunc, fitfunc_str

def fit_8msmt_state_fid(g_A0, g_t, g_p):
    '''


    g_A0        -  Amplitude derived from the 0 msmt case. Fixed.
    g_p         -  Probability for faulty parity measurement
    g_t         -  Decay of Ramsey. Fixed.

    The function should take the data set for 8 Zeno-measurements and fit the state fidelity with a prederived function. 
    That depends on one parameter. The Probability of the parity measurement to mix the state.
    '''

    fitfunc_str = '''analyitcal solution'''

    ### Parameters

    A0           = fit.Parameter(g_A0, 'A0')


    ### Ramsey (decay time divided by sqrt(2).)
    t   = fit.Parameter(g_t, 't')

    ### Zeno
    p   = fit.Parameter(g_p,'p')

    p0 = [A0,t,p]

    def fitfunc(x):

        decay = np.exp(- ( x**2 /(2. * (t()**2) ) ) )/512.

        decay2 = 9. * np.exp((16./81.)*( x**2 /((t()**2) ) ) )*(p()-1.)**8

        decay3 = 36. * np.exp((28./81.)*( x**2 /((t()**2) ) ) )*(p()-1.)**8

        decay4 = 84. * np.exp((36./81.)*( x**2 /((t()**2) ) ) )*(p()-1.)**8

        decay5 = 126. * np.exp((40./81.)*( x**2 /((t()**2) ) ) )*(p()-1.)**8

        #### decay to 0.5 only for 6 measurements no other constant parts of the function due to echo effects.
        
        Fx = 0.5 + decay*((p()-1)**8 + decay2 + decay3 + decay4 + decay5)

        Cx = A0()*2.*(Fx-0.5)
        Fx = Cx/2.+0.5
        return Fx


    return p0, fitfunc, fitfunc_str

def fit_X_state_fid(msmts,g_A0, g_p):
    '''
    g_A0        -  Amplitude derived from the 0 msmt case. Fixed.
    g_p         -  Probability for faulty parity measurement
    msmts       -  Gives the number of measurements

    The function should take the data set for 6 Zeno-measurements and fit the process fidelity with a prederived function. That depends on one parameter. The fidelity of the parity measurement.
    '''

    fitfunc_str = '''A0*(1-p)**n+0.5'''

    ### Parameters

    A0           = fit.Parameter(g_A0, 'A0')

    ### Zeno
    p   = fit.Parameter(g_p,'p')

    p0 = [A0,p]

    def fitfunc(x):

        Fz0 = A0()
        Con0 = 2.*(Fz0-0.5) ### contrast of the ZZ correlations for 0 measurements.
        Conp = Con0*(1-p())**msmts ### after 6 measurements the contrast is reduced by (1-p)^6
        Fzp = Conp/2. + 0.5 ### calculate the fidelity.

        # Fz = A0()*(1-p())**msmts + 0.5
        return Fzp


    return p0, fitfunc, fitfunc_str