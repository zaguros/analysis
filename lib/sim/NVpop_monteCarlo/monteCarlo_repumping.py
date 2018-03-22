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

        self.monteCarlo = pmc.monteCarlo

    def run(self,**kw):

        if not isinstance(self.init_state, (np.ndarray)): # Can specify initial state using only an index, but pmc wants an array (so that can specify per rep)
            init_state = np.array([self.init_state],dtype=np.intc)

        self.rateMat = self.repumping_rateMat().copy(order='C')
        
        if self.drive_time and self.drive_time < self.time_duration: # If drive time specified

            self.no_drive_RateMat = self.repumping_rateMat(drive = False).copy(order='C')

            # First do bit with drive
            drive_t_array,drive_populations,drive_raw_pops  = self.monteCarlo(init_state,self.rateMat,self.drive_time,self.dt,self.repetitions)
            # Then add on bit with no drive (tracking where the state is)!
            no_drive_t_array,no_drive_populations,no_drive_raw_pops  = self.monteCarlo(drive_raw_pops[:,-1].astype(np.intc).copy(order='C'),self.no_drive_RateMat,(self.time_duration - self.drive_time),self.dt,self.repetitions)
            
            #Stick it all together
            self.t_array = np.concatenate((drive_t_array,drive_t_array[-1]+no_drive_t_array[1:]))
            self.populations = np.concatenate((drive_populations,no_drive_populations[:,1:,:]),axis=1)
            self.raw_pops = np.concatenate((drive_raw_pops,no_drive_raw_pops[:,1:]),axis=1)

        else:
            self.t_array,self.populations,self.raw_pops = self.monteCarlo(init_state,self.rateMat,self.time_duration,self.dt,self.repetitions)

        self.mean_populations = np.mean(self.populations,axis = 0)
        self.correct_for_singlet_decay()

        self.plot_populations(**kw)

    def init_default_branching(self):

        # P,M,Z,S,E1,E2
        self.S_lifetime = 300.0
        self.drive = 0.2
        self.drive_time = False
        self.drive_branching = normalized([1,0])

        E1_branching_raw = np.array([1.0,0.0,0.0,0.7,0.0,0.0])
        E2_branching_raw = np.array([0.0,1.0,0.0,0.7,0.0,0.0])
        self.E1_lifetime = 11.7/np.sum(E1_branching_raw)
        self.E2_lifetime = 11.7/np.sum(E2_branching_raw)
        self.E1_branching = normalized(E1_branching_raw)
        self.E2_branching = normalized(E2_branching_raw)

        self.S_branching = normalized([1.0,1.0,8.0,0.0,0.0,0.0])


    def repumping_rateMat(self,drive = True):
        
        # P,M,Z,S,E1,E2

        rateMat = np.zeros([self.num_states,self.num_states])

        if drive:
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
        self.populations_after_decay[:,4] = 0
        self.populations_after_decay[:,5] = 0
        self.populations_after_decay[:,3] = 0

    def plot_populations(self,**kw):

        show_guess = kw.pop('show_guess',False)
        do_fit =  kw.pop('do_fit','none')
        log_plot = kw.pop('log_plot',True)
        invert_Z = kw.pop("invert_Z",False)
        plot_model = kw.pop("plot_model",False)
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
        if plot_model:
            self.plot_model_for_data(log_plot)

        ax = plt.gca()
        # Shrink current axis by 20%
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
        # Put a legend to the right of the current axis
        ax.legend(["P","M",Z_string,"S","E1","E2"],loc='center left', bbox_to_anchor=(1, 0.5))

            
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

      
        plt.show()
        plt.close()

    def plot_model_for_data(self,log_plot):

        rate_into_S = 0.5*self.E2_branching[3]/self.E2_lifetime # Strong drive limit
        rate_from_S_to_Z = self.S_branching[2]/self.S_lifetime
        rate_from_S_to_PM = (1-self.S_branching[2])/self.S_lifetime
        steady_state_ratio_PM_to_S = rate_from_S_to_PM/rate_into_S
    
        t = self.t_array

        initial_decrease_M = np.exp(-rate_into_S * t)
        S_decay_during_drive = np.exp(-rate_from_S_to_Z * t) * (1- initial_decrease_M)

        Spop = S_decay_during_drive # Roughly!
        Mpop = 0.5*initial_decrease_M + 0.25*steady_state_ratio_PM_to_S*S_decay_during_drive #0.25 because can be in M,P or E1,E2
        Ppop = 0.25*steady_state_ratio_PM_to_S*S_decay_during_drive

        if self.drive_time:
            ind_drive = np.argmin(np.abs(self.t_array - self.drive_time))
            t_drive = t[(ind_drive+1):]-self.drive_time
            
            S_decay_no_drive = np.exp(-t_drive/self.S_lifetime)
            E2decay = np.exp(-t_drive/self.E2_lifetime)
            Spop[(ind_drive+1):] = Spop[ind_drive]*S_decay_no_drive # Roughly!
            Mpop[(ind_drive+1):] = (2-E2decay) * Mpop[ind_drive] + self.S_branching[1]*Spop[ind_drive]*(1-S_decay_no_drive)
            Ppop[(ind_drive+1):] = (2-E2decay) * Ppop[ind_drive] + self.S_branching[0]*Spop[ind_drive]*(1-S_decay_no_drive)

       
        plot_log_or_lin(t,np.transpose([Spop,Mpop,Ppop]),log_plot,linestyle='--',color="black")

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

def plot_log_or_lin(x,y,log_plot,**kw):
    if log_plot:
        plt.semilogy(x,y,**kw)
    else:
        plt.plot(x,y,**kw)

def normalized(a, axis=-1):
    a = np.array(a)
    return a/np.linalg.norm(a, 1, axis)



# Old fitting stuff
  #       if do_fit == 'before_decay':
  #           x = self.t_array
  #           y = pops_to_plot[:,2]
  #           points = np.size(x)
  #           if points > 1000:
  #               reduced_points = np.linspace(0,points-1,1000).astype(int)
  #               x = x[reduced_points]
  #               y = y[reduced_points]
  #           param_guesses =  kw.pop('param_guesses',[0.0,1.0,300])

  #           p0, fitfunc, fitfunc_str =  common.fit_exp_decay_with_offset(param_guesses[0],param_guesses[1],param_guesses[2])

  #           if show_guess:
  #               plot_log_or_lin(np.linspace(x[0],x[-1],201), fitfunc(np.linspace(0,x[-1],201)),log_plot)
      
  #           fit_result = fit.fit1d(x,y, None, p0=p0, fitfunc=fitfunc, do_print=False, ret=True,fixed=fixed)
  #           ## plot fit
  #           plot_log_or_lin(np.linspace(x[0],x[-1],201),fit_result['fitfunc'](np.linspace(x[0],x[-1],201)),log_plot)
            
  #           print 'params: ', np.around(fit_result['params'],3)


  # if do_fit == 'after_decay':
  #           x = self.t_array
  #           y = pops_after_decay[:,2]
  #           points = np.size(x)
  #           if points > 1000:
  #               reduced_points = np.linspace(0,points-1,1000).astype(int)
  #               x = x[reduced_points]
  #               y = y[reduced_points]
            
  #           param_guesses =  kw.pop('param_guesses',[0.5,0.5,10,300])

  #           p0, fitfunc, fitfunc_str = fit_bi_exp(param_guesses[0],param_guesses[1],param_guesses[2],param_guesses[3])

  #           if show_guess:
  #               plot_log_or_lin(np.linspace(x[0],x[-1],201), fitfunc(np.linspace(0,x[-1],201)),log_plot)
      
  #           fit_result = fit.fit1d(x,y, None, p0=p0, fitfunc=fitfunc, do_print=False, ret=True,fixed=fixed)
  #           ## plot fit
  #           plot_log_or_lin(np.linspace(x[0],x[-1],201),fit_result['fitfunc'](np.linspace(x[0],x[-1],201)),log_plot)
            
  #           print 't1,t2: ', np.around(fit_result['params'][2]),np.around(fit_result['params'][3])
  #           print 'relative rate of exp decays: ', np.around(fit_result['params'][0]/fit_result['params'][1],2)
