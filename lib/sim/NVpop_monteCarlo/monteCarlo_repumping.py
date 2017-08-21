### imports
import numpy as np
from matplotlib import pyplot as plt
import pop_montecarlo_c as pmc; reload(pmc)
from analysis.lib.fitting import fit, common

class repumpingMonteCarlo():

    def __init__(self,**kw):
        self.init_state = 1
        self.time_duration = 1000
        self.dt = 1.0
        self.repetitions = 1000

        self.num_states = 6
        self.init_default_branching()

        self.extra_singlet_state = False

    def run(self,**kw):

        if self.extra_singlet_state:
            self.num_states = 7

        self.rateMat = self.repumping_rateMat()
        self.rateMat = self.rateMat.copy(order='C')
        self.t_array,self.populations = pmc.monteCarlo(self.init_state,self.rateMat,self.time_duration,self.dt,self.repetitions)

        self.mean_populations = np.mean(self.populations,axis = 0)
        self.correct_for_singlet_decay()

        self.plot_populations(**kw)

    def init_default_branching(self):

        # P,M,Z,S,E1,E2
        self.S_lifetime = 300.0
        self.drive = 0.2

        self.drive_branching = normalized([1,0])

        E1_branching_raw = np.array([1.0,0.0,0.0,0.7,0.0,0.0])
        E2_branching_raw = np.array([0.0,1.0,0.0,0.7,0.0,0.0])
        self.E1_lifetime = 11.7/np.sum(E1_branching_raw)
        self.E2_lifetime = 11.7/np.sum(E2_branching_raw)
        self.E1_branching = normalized(E1_branching_raw)
        self.E2_branching = normalized(E2_branching_raw)

        self.S_branching = normalized([1.0,1.0,8.0,0.0,0.0,0.0])


    def repumping_rateMat(self):
        
        # P,M,Z,S,E1,E2

        rateMat = np.zeros([self.num_states,self.num_states])

        # P to E1/E2 under drive
        rateMat[0,4] += self.drive * self.drive_branching[0]
        rateMat[0,5] += self.drive * self.drive_branching[1]
        rateMat[4,0] += self.drive * self.drive_branching[0]
        rateMat[5,0] += self.drive * self.drive_branching[1]
        # M to E1/E2 under drive
        rateMat[1,4] += self.drive * self.drive_branching[1]
        rateMat[1,5] += self.drive * self.drive_branching[0]
        rateMat[4,1] += self.drive * self.drive_branching[1]
        rateMat[5,1] += self.drive * self.drive_branching[0]


        if self.extra_singlet_state:
            S_branching = np.append(self.S_branching,0)
            E1_branching = np.append(self.E1_branching,self.E1_branching[3]/2.0)
            E1_branching[3] = self.E1_branching[3]/2.0
            E2_branching = np.append(self.E2_branching,self.E2_branching[3]/2.0)
            E2_branching[3] = self.E2_branching[3]/2.0
            rateMat[6,:] += S_branching/self.S_lifetime
        else:
            S_branching = self.S_branching
            E1_branching = self.E1_branching
            E2_branching = self.E2_branching

        # Spontaneous emission

        rateMat[3,:] += S_branching/self.S_lifetime
        rateMat[4,:] += E1_branching/self.E1_lifetime
        rateMat[5,:] += E2_branching/self.E2_lifetime

        return rateMat

    def correct_for_singlet_decay(self):
        

        self.populations_after_decay = np.copy(self.mean_populations)
        self.populations_after_decay += np.outer(self.populations_after_decay[:,4],normalized(self.rateMat[4,:]))
        self.populations_after_decay += np.outer(self.populations_after_decay[:,5],normalized(self.rateMat[5,:]))
        self.populations_after_decay += np.outer(self.populations_after_decay[:,3],normalized(self.rateMat[3,:]))
        if self.extra_singlet_state:
            self.populations_after_decay += np.outer(self.populations_after_decay[:,6],normalized(self.rateMat[6,:]))
            self.populations_after_decay[:,6] = 0
        self.populations_after_decay[:,4] = 0
        self.populations_after_decay[:,5] = 0
        self.populations_after_decay[:,3] = 0

    def plot_populations(self,**kw):

        show_guess = kw.pop('show_guess',False)
        do_fit =  kw.pop('do_fit','none')
        log_plot = kw.pop('log_plot',True)
        invert_Z = kw.pop("invert_Z",False)
        print_end_pops = kw.pop("print_end_pops",False)
        
        fixed =  kw.pop('fixed',[])

        pops_to_plot = np.copy(self.mean_populations)
        pops_after_decay = np.copy(self.populations_after_decay)

        if print_end_pops:

            print 'Final pops before decay from singlet: ', np.around(pops_to_plot[-1],3)
            print 'Final pops after decay from singlet: ', np.around(pops_after_decay[-1],3)

        if invert_Z:
            pops_to_plot[:,2] = 1-pops_to_plot[:,2]
            pops_after_decay[:,2] = 1- pops_after_decay[:,2]
            Z_string = "1-Z"
        else:
            Z_string = "Z"


        plt.figure()
        plt.xlabel('Time (ns)')
        plt.ylabel('Population')
        plot_log_or_lin(self.t_array,pops_to_plot,log_plot)

        ax = plt.gca()
        # Shrink current axis by 20%
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
        # Put a legend to the right of the current axis
        ax.legend(["P","M",Z_string,"S","E1","E2","S2"],loc='center left', bbox_to_anchor=(1, 0.5))

        if do_fit == 'before_decay':
            x = self.t_array
            y = pops_to_plot[:,2]
            points = np.size(x)
            if points > 1000:
                reduced_points = np.linspace(0,points-1,1000).astype(int)
                x = x[reduced_points]
                y = y[reduced_points]
            param_guesses =  kw.pop('param_guesses',[0.0,1.0,300])

            p0, fitfunc, fitfunc_str =  common.fit_exp_decay_with_offset(param_guesses[0],param_guesses[1],param_guesses[2])

            if show_guess:
                plot_log_or_lin(np.linspace(x[0],x[-1],201), fitfunc(np.linspace(0,x[-1],201)),log_plot)
      
            fit_result = fit.fit1d(x,y, None, p0=p0, fitfunc=fitfunc, do_print=False, ret=True,fixed=fixed)
            ## plot fit
            plot_log_or_lin(np.linspace(x[0],x[-1],201),fit_result['fitfunc'](np.linspace(x[0],x[-1],201)),log_plot)
            
            print 'params: ', np.around(fit_result['params'],3)
            
        plt.show()
        plt.close()

        plt.figure()
        plt.xlabel('Time (ns)')
        plt.ylabel('Population')
        plot_log_or_lin(self.t_array,pops_after_decay[:,0:3],log_plot)
        
        ax = plt.gca()
        # Shrink current axis by 20%
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
        # Put a legend to the right of the current axis
        ax.legend(["P","M",Z_string],loc='center left', bbox_to_anchor=(1, 0.5))

        if do_fit == 'after_decay':
            x = self.t_array
            y = pops_after_decay[:,2]
            points = np.size(x)
            if points > 1000:
                reduced_points = np.linspace(0,points-1,1000).astype(int)
                x = x[reduced_points]
                y = y[reduced_points]
            
            param_guesses =  kw.pop('param_guesses',[0.5,0.5,10,300])

            p0, fitfunc, fitfunc_str = fit_bi_exp(param_guesses[0],param_guesses[1],param_guesses[2],param_guesses[3])

            if show_guess:
                plot_log_or_lin(np.linspace(x[0],x[-1],201), fitfunc(np.linspace(0,x[-1],201)),log_plot)
      
            fit_result = fit.fit1d(x,y, None, p0=p0, fitfunc=fitfunc, do_print=False, ret=True,fixed=fixed)
            ## plot fit
            plot_log_or_lin(np.linspace(x[0],x[-1],201),fit_result['fitfunc'](np.linspace(x[0],x[-1],201)),log_plot)
            
            print 't1,t2: ', np.around(fit_result['params'][2]),np.around(fit_result['params'][3])
            print 'relative rate of exp decays: ', np.around(fit_result['params'][0]/fit_result['params'][1],2)
        plt.show()
        plt.close()



def fit_bi_exp(g_A1, g_A2, g_t1, g_t2,*arg):
    fitfunc_str = 'A1 *exp(-x/t1) + A2 *exp(-x/t2)'

    A1 = fit.Parameter(g_A1, 'A1')
    A2 = fit.Parameter(g_A2, 'A2')
    t1   = fit.Parameter(g_t1, 't1')
    t2   = fit.Parameter(g_t2, 't2')

    p0 = [A1,A2,t1,t2]

    def fitfunc(x):
        return A1()*np.exp(-x/t1()) + A2()*np.exp(-x/t2())

    return p0, fitfunc, fitfunc_str

def plot_log_or_lin(x,y,log_plot):
    if log_plot:
        plt.semilogy(x,y)
    else:
        plt.plot(x,y)

def normalized(a, axis=-1):
    a = np.array(a)
    return a/np.linalg.norm(a, 1, axis)

