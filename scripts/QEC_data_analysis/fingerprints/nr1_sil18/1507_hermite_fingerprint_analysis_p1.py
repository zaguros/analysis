'''
Script to analyze the dynamical decoupling data of 21 and 22 /07
N = 4,8,16,32 and 64, contrast measurement
no Nitrogen MBI, Hermite pulses
'''
import pickle
import numpy as np
import os
from analysis.lib.tools import toolbox
from analysis.lib.m2.ssro import mbi
import analysis.lib.QEC.nuclear_spin_characterisation as SC #used for simulating FP response of spins
# import magnettools as mt # Does not work atm because of qt lab being imported in MT
import matplotlib.cm as cm
from matplotlib import pyplot as plt
from analysis.lib.tools import toolbox
reload(toolbox)

from analysis.scripts.QEC_data_analysis.fingerprints import fingerprint_funcs as fp_funcs; reload(fp_funcs)


def fingerprint_contrast_concatenate(disp_sim_spin = True,n_sim_spins = 2,pts = 51,start_1 = 3.0, start_2 = 3.0+45*50*10e-3,
        step_size = 10e-3,tau_larmor = False, name = '',
        xrange = [0,20],tag_p = '',tag_n = '', older_than = None,return_data = False, do_plot = True, load_from_data = False, Nr_of_pulses = None, save_folder = None, figsize = (20,5)):


    ###################
    # Add simulated spins #
    ###################

    if disp_sim_spin == True:
            
            HF_perp, HF_par = fp_funcs.get_hyperfine_params(ms = 'plus')
            #msmp1_f from hdf5 file
            # msm1 from hdf5 file
            # ZFG g_factor from hdf5file
            B_Field = 403.555 # use magnet tools  Bz = (msp1_f**2 - msm1_f**2)/(4.*ZFS*g_factor)
            t_L = 1/(1.0705e3*B_Field)*1e6
            tau_lst = np.linspace(0,72e-6,10000)
            Mt16 = SC.dyn_dec_signal(HF_par,HF_perp,B_Field,Nr_of_pulses,tau_lst)
            FP_signal16 = ((Mt16+1)/2)

    ## Data location ##
    
    if load_from_data == False:


      ssro_calib_folder = folder = toolbox.latest_data(contains = 'AdwinSSRO_SSROCalibration_111_1_sil18', older_than = older_than,folder = 'd:\measuring\data')
      
      a1, folder = load_mult_dat_tag(tag_p+'0',older_than, number_of_msmts = 50,pts = pts, start=start_1, step_size=step_size,ssro_calib_folder=ssro_calib_folder)
      b1, folder_b1 = load_mult_dat_tag(tag_n+'0', older_than, number_of_msmts = 50,pts = pts, start=start_1, step_size=step_size,ssro_calib_folder=ssro_calib_folder)

      a2, folder = load_mult_dat_tag(tag_p+'45',older_than, number_of_msmts = 50,pts = pts, start=start_2, step_size=step_size,ssro_calib_folder=ssro_calib_folder)
      b2, folder_b2 = load_mult_dat_tag(tag_n+'45', older_than, number_of_msmts = 50,pts = pts, start=start_2, step_size=step_size,ssro_calib_folder=ssro_calib_folder)

      print folder_b2
      print folder_b1

      ###############
      ## Plotting ###
      ###############




      y = np.concatenate(((a1.p0-b1.p0),(a2.p0-b2.p0)))
      x = np.concatenate((a1.sweep_pts, a2.sweep_pts))
      y_err = np.concatenate(((np.sqrt(a1.u_p0**2+b1.u_p0**2)),(np.sqrt(a2.u_p0**2+b2.u_p0**2))))

      data = {}
      data['x'] = x
      data['y'] = y
      data['y_err']= y_err

      print y_err



      # pickle.dump(data, open( 'sil18_fingerprint_ms_min_N'+str(Nr_of_pulses)+'.p', 'wb' ) )
      np.savetxt('sil18_fingerprint_incl_errorbar_ms_plus_N'+str(Nr_of_pulses)+'.txt',(np.c_[x],y,y_err))
      # print 'shapes'
      # print np.transpose(x[0:10])
      # print np.array(y[0:10])

      # x = np.arange(0,10,1)
      # y = x**2

      # pickle.dump(data, open( 'sil18_fingerprint_ms_plus_N'+str(Nr_of_pulses)+'.p', 'wb' ) )
      # np.savetxt('sil18_fingerprint_ms_plus_N'+str(Nr_of_pulses)+'.txt',(np.c_[x],y))

    else:
      # data = pickle.load( open( 'sil18_fingerprint_ms_plus_N'+str(Nr_of_pulses)+'.p', 'rb' ) )

      # x = data['x']
      # y = data['y']
      x, y, y_err = np.loadtxt('sil18_fingerprint_incl_errorbar_ms_plus_N'+str(Nr_of_pulses)+'.txt')
      print y_err
      folder = save_folder


    if tau_larmor == True:
        x = x/t_L

    y = (y + 1)/2.
    y_err = y_err/2.


    folder = save_folder
    if do_plot == True:
      fig,ax = plt.subplots(figsize=figsize)
      ax.set_xlim(4.9,5.1)
      ax.set_xlim(xrange)
      start, end = ax.get_xlim()
      ax.xaxis.set_ticks(np.arange(xrange[0], xrange[1]+0.1, (xrange[1]- xrange[0])/4.))

      # ax.set_ylim(-1.05,1.05)
      ax.set_ylim(-0.05,1.05)
      ax.yaxis.set_ticks([0,0.5,1])
      # ax.plot(x, y, '.-b', lw=1,label = 'data',marker = 'o') #N = 16
      # ax.plot(x, y, '.-k', lw=0.4,label = 'data') #N = 16
      # ax.plot(x, (y+1)/2., '.-b', lw=1,label = 'data', marker = 'o', ms = 8) #N = 16
      ax.plot(x, y, '.-k', lw=0.4,label = 'data', marker = 'o', ms = 3, markeredgecolor = 'm') #N = 16
      ax.set_xlabel('tau (us)')
      if tau_larmor == True:
        ax.set_xlabel('tau/tau_larmor')

      ax.set_ylabel('Fidelity')

      # ax.plot(a1.sweep_pts, (a1.p0-b1.p0), '.-b', lw=0.4,label = 'data') #N = 16
      # ax.plot(a2.sweep_pts, (a2.p0-b2.p0), '.-k', lw=0.4,label = 'data') #N = 16
      if disp_sim_spin == True:
        colors = cm.rainbow(np.linspace(0, 1, n_sim_spins))
        for tt in range(n_sim_spins):
          # print tt
          if tau_larmor == False:
            ax.plot(tau_lst*1e6, Mt16[tt,:] ,'-',lw=.8,label = 'spin' + str(tt+1), color = colors[tt])
          if tau_larmor == True:
            ax.plot(tau_lst*1e6/t_L, Mt16[tt,:] ,'-',lw=.8,label = 'spin' + str(tt+1), color = colors[tt])            
      if False:
          tot_signal = np.ones(len(tau_lst))
          for tt in range(n_sim_spins):
            tot_signal = tot_signal * Mt16[tt,:]
          fin_signal = (tot_signal+1)/2.0
          ax.plot(tau_lst*1e6, fin_signal,':g',lw=.8,label = 'tot')


      # ax.vlines([5.5,6.5],-1.1,1.1,color = '0.5',lw = 1,linestyles = '-')
      # plt.axvspan(5.5,6.5, facecolor='c', alpha=0.1)


      # lgd = plt.legend(loc=4)
      plt.show(block = False)
      # ax.vlines([5.5,6.5],-1.1,1.1,color = '0.5',lw = 1,linestyles = '-')
      # plt.axvspan(5.5,6.5, facecolor='m', alpha=0.1)
      print folder

      try: 
        # plt.savefig(os.path.join(folder, 'contrast_150924.pdf'),
        #     format='pdf',bbox_inches='tight')
        # plt.savefig(os.path.join(folder, 'contrast_150924.png'),
        #     format='png',bbox_inches='tight')
        plt.savefig(os.path.join(folder, name+'.pdf'),
            format='pdf',bbox_inches='tight')
        plt.savefig(os.path.join(folder, name+'.png'),
            format='png',bbox_inches='tight')       
      except:
        print 'Figure has not been saved'

    if return_data == True:
      return x, (y+1)/2.

def fingerprint_single(disp_sim_spin = True,n_sim_spins = 2,xrange = [0,20],tag = '', step_size = 10e-3,
                   older_than = None,return_data = False):


    ###################
    # Add simulated spins #
    ###################

    if disp_sim_spin == True:
            
            HF_perp, HF_par = fp_funcs.get_hyperfine_params(ms = 'min')
            #msmp1_f from hdf5 file
            # msm1 from hdf5 file
            # ZFG g_factor from hdf5file
            B_Field = 403.555 # use magnet tools  Bz = (msp1_f**2 - msm1_f**2)/(4.*ZFS*g_factor)

            tau_lst = np.linspace(0,72e-6,10000)
            Mt16 = SC.dyn_dec_signal(HF_par,HF_perp,B_Field,32,tau_lst)
            FP_signal16 = ((Mt16+1)/2)

    ## Data location ##
    print older_than
    timestamp, ssro_calib_folder = toolbox.latest_data(contains = 'AdwinSSRO', older_than = older_than,return_timestamp = True)
    print ssro_calib_folder
    a, folder = load_mult_dat_tag(tag,older_than, number_of_msmts = 18, ssro_calib_folder=ssro_calib_folder)

    ###############
    ## Plotting ###
    ###############

    fig = a.default_fig(figsize=(35,5))
    ax = a.default_ax(fig)
    ax.set_xlim(4.9,5.1)
    ax.set_xlim(xrange)
    start, end = ax.get_xlim()
    ax.xaxis.set_ticks(np.arange(xrange[0], xrange[1], (xrange[1]- xrange[0])/10.))

    ax.set_ylim(-0.05,1.05)

    y = a.p0

    ax.plot(a.sweep_pts, y, '.-g', lw=0.4,label = 'data', marker = 'o', ms = 3) #N = 16
    if disp_sim_spin == True:
      colors = cm.rainbow(np.linspace(0, 1, n_sim_spins))
      for tt in range(n_sim_spins):
        ax.plot(tau_lst*1e6, FP_signal16[tt,:] ,'-',lw=.8,label = 'spin' + str(tt+1), color = colors[tt])
    if False:
        tot_signal = np.ones(len(tau_lst))
        for tt in range(n_sim_spins):
          tot_signal = tot_signal * Mt16[tt,:]
        fin_signal = (tot_signal+1)/2.0
        ax.plot(tau_lst*1e6, fin_signal,':g',lw=.8,label = 'tot')


    lgd = plt.legend(loc=4)
    plt.show(block = False)

    print folder
    plt.savefig(os.path.join(folder, str(disp_sim_spin)+'fingerprint.pdf'),
        format='pdf',bbox_extra_artists = (lgd,),bbox_inches='tight')
    plt.savefig(os.path.join(folder, str(disp_sim_spin)+'fingerprint.png'),
        format='png',bbox_extra_artists = (lgd,),bbox_inches='tight')

    if return_data == True:
      return a.sweep_pts, a.p0

def load_mult_dat_tag(tag,older_than, number_of_msmts, ssro_calib_folder='',start = 3.0, pts = 51, step_size = 10e-3):
   cum_pts = 0
   for kk in range(number_of_msmts):
       folder = toolbox.latest_data(contains = tag, older_than = older_than,folder = 'd:\measuring\data')
       a = mbi.MBIAnalysis(folder)
       a.get_sweep_pts()
       a.get_readout_results(name='measurement' + str(kk))
       a.get_electron_ROC(ssro_calib_folder=ssro_calib_folder)
       cum_pts += a.pts

       # print a.reps

       if kk == 0:
           cum_sweep_pts = np.linspace(start, start+(pts-1)*step_size,pts)
           cum_p0 = a.p0
           cum_u_p0 = a.u_p0
       else:
           cum_sweep_pts = np.concatenate((cum_sweep_pts, np.linspace(start+kk*(pts-1)*step_size, start+(pts-1)*step_size+kk*(pts-1)*step_size, pts)))
           cum_p0 = np.concatenate((cum_p0, a.p0))
           cum_u_p0 = np.concatenate((cum_u_p0, a.u_p0))

   a.pts   = cum_pts
   a.sweep_pts = cum_sweep_pts
   a.p0    = cum_p0
   a.u_p0  = cum_u_p0

   return a, folder



def fingerprint(disp_sim_spin = True,n_sim_spins = 2,step_size= 10e-3, xrange = [0,20],return_data = False,    timestamp ='20141016_205842',
    ssro_calib_folder = 'D:\\measuring\\data\\20141016\\150451_AdwinSSRO_SSROCalibration_111_1_sil18'):


    ###################
    # Add simulated spins #
    ###################

    if disp_sim_spin == True:
            
            HF_perp, HF_par = fp_funcs.get_hyperfine_params(ms = 'min')
            #msmp1_f from hdf5 file
            # msm1 from hdf5 file
            # ZFG g_factor from hdf5file
            B_Field = 403.555 # use magnet tools  Bz = (msp1_f**2 - msm1_f**2)/(4.*ZFS*g_factor)

            tau_lst = np.linspace(0,72e-6,10000)
            Mt16 = SC.dyn_dec_signal(HF_par,HF_perp,B_Field,32,tau_lst)
            FP_signal16 = ((Mt16+1)/2)

    ## Data location ##

    a, folder = load_mult_dat(timestamp,step_size = step_size, number_of_msmts = 40, ssro_calib_folder=ssro_calib_folder)

    ###############
    ## Plotting ###
    ###############

    fig = a.default_fig(figsize=(35,5))
    ax = a.default_ax(fig)
    ax.set_xlim(4.9,5.1)
    ax.set_xlim(xrange)
    start, end = ax.get_xlim()
    ax.xaxis.set_ticks(np.arange(xrange[0], xrange[1], (xrange[1]- xrange[0])/10.))

    ax.set_ylim(-0.05,1.05)
    ax.plot(a.sweep_pts, a.p0, '.-k', lw=0.4,label = 'data') #N = 16
    if disp_sim_spin == True:
      colors = cm.rainbow(np.linspace(0, 1, n_sim_spins))
      for tt in range(n_sim_spins):
        ax.plot(tau_lst*1e6, FP_signal16[tt,:] ,'-',lw=.8,label = 'spin' + str(tt+1), color = colors[tt])
    if False:
        tot_signal = np.ones(len(tau_lst))
        for tt in range(n_sim_spins):
          tot_signal = tot_signal * Mt16[tt,:]
        fin_signal = (tot_signal+1)/2.0
        ax.plot(tau_lst*1e6, fin_signal,':g',lw=.8,label = 'tot')


    lgd = plt.legend(loc=4)

    print folder
    plt.savefig(os.path.join(folder, str(disp_sim_spin)+'fingerprint.pdf'),
        format='pdf',bbox_extra_artists = (lgd,),bbox_inches='tight')
    plt.savefig(os.path.join(folder, str(disp_sim_spin)+'fingerprint.png'),
        format='png',bbox_extra_artists = (lgd,),bbox_inches='tight')

    if return_data == True:
      return a.sweep_pts, a.p0

def load_mult_dat(timestamp, number_of_msmts,older_than = '', ssro_calib_folder='',start = 2.0, pts = 51, step_size = 10e-3):
   cum_pts = 0
   for kk in range(number_of_msmts):
       folder = toolbox.data_from_time(timestamp,folder = 'd:\measuring\data')
       a = mbi.MBIAnalysis(folder)
       a.get_sweep_pts()
       a.get_readout_results(name='measurement' + str(kk))
       a.get_electron_ROC(ssro_calib_folder=ssro_calib_folder)
       cum_pts += a.pts

       if kk == 0:
           cum_sweep_pts = np.linspace(start, start+(pts-1)*step_size,pts)
           cum_p0 = a.p0
           cum_u_p0 = a.u_p0
       else:
           cum_sweep_pts = np.concatenate((cum_sweep_pts, np.linspace(start+kk*(pts-1)*step_size, start+(pts-1)*step_size+kk*(pts-1)*step_size, pts)))
           cum_p0 = np.concatenate((cum_p0, a.p0))
           cum_u_p0 = np.concatenate((cum_u_p0, a.u_p0))

   a.pts   = cum_pts
   a.sweep_pts = cum_sweep_pts
   a.p0    = cum_p0
   a.u_p0  = cum_u_p0

   return a, folder


def plot_all(disp_sim_spin = True,xrange  = [2.5,52.5], n_sim_spins = 8,load_data = True):
  if load_data == False:
    data = {}
    data['x_64'] , data['y_64'] = fingerprint_contrast_concatenate(disp_sim_spin = True,n_sim_spins = 0,pts = 51,start_1 = 3.0, start_2 = 3+45*50*4e-3,
          step_size = 4e-3,
          xrange = [2.5,52.5],tag_p = 'Hermite_Fingerprint_msp1_111_1_sil18_64-x',tag_n = 'Hermite_Fingerprint_msp1_111_1_sil18_64x', 
          older_than = '20150726090000',return_data = True, do_plot= False)

    data['x_128'] ,data['y_128'] = fingerprint_contrast_concatenate(disp_sim_spin = True,n_sim_spins = 0,pts = 21,start_1 = 3.0, start_2 = 3+45*20*4e-3,
            step_size = 4e-3,
            xrange = [2.5,52.5],tag_p = 'Hermite_Fingerprint_msp1_111_1_sil18_128-x',tag_n = 'Hermite_Fingerprint_msp1_111_1_sil18_128x', 
            older_than = '20150726090000',return_data = True, do_plot= False)

    data['x_4'], data['y_4']= fingerprint_contrast_concatenate(disp_sim_spin = True,n_sim_spins = 0,pts = 51,start_1 = 3.0, start_2 = 3+45*50*10e-3,
            step_size = 10e-3,
            xrange = [2.5,52.5],tag_p = 'Hermite_Fingerprint_msp1_111_1_sil18_4-x',tag_n = 'Hermite_Fingerprint_msp1_111_1_sil18_4x', 
            older_than = '20150726090000',return_data = True, do_plot= False)


    data['x_8'], data['y_8']= fingerprint_contrast_concatenate(disp_sim_spin = True,n_sim_spins = 0,pts = 51,start_1 = 3.0, start_2 = 3+45*50*10e-3,
            step_size = 10e-3,
            xrange = [2.5,52.5],tag_p = 'Hermite_Fingerprint_msp1_111_1_sil18_8-x',tag_n = 'Hermite_Fingerprint_msp1_111_1_sil18_8x', 
            older_than = '20150726090000',return_data = True, do_plot= False)


    data['x_16'] , data['y_16'] = fingerprint_contrast_concatenate(disp_sim_spin = True,n_sim_spins = 0,pts = 51,start_1 = 3.0, start_2 = 3+45*50*10e-3,
            step_size = 10e-3,
            xrange = [2.5,52.5],tag_p = 'Hermite_Fingerprint_msp1_111_1_sil18_16-x',tag_n = 'Hermite_Fingerprint_msp1_111_1_sil18_16x', 
            older_than = '20150726090000',return_data = True, do_plot= False)


    data['x_32'] , data['y_32'] = fingerprint_contrast_concatenate(disp_sim_spin = True,n_sim_spins = 0,pts = 51,start_1 = 3.0, start_2 = 3+45*50*10e-3,
            step_size = 10e-3,
            xrange = [2.5,52.5],tag_p = 'Hermite_Fingerprint_msp1_111_1_sil18_32-x',tag_n = 'Hermite_Fingerprint_msp1_111_1_sil18_32x', 
            older_than = '20150726090000',return_data = True, do_plot= False)


    pickle.dump(data, open( "fingerprinting_msp1.p", "wb" ) )

  else:
    data = pickle.load( open( "fingerprinting_msp1.p", "rb" ) )

  ###############
  ## Plotting ###
  ###############
  xrange = xrange

  fig = plt.figure(figsize=(10,35))
  ax4 = plt.subplot(511)
  ax8 = plt.subplot(512, sharex=ax4)
  ax16 = plt.subplot(513, sharex=ax4)
  ax32 = plt.subplot(514, sharex=ax4)
  ax64 = plt.subplot(515, sharex=ax4)


  for ax in [ax4,ax8,ax16,ax32,ax64]:
    ax.set_xlim(xrange)
    start, end = ax.get_xlim()
    # ax.xaxis.set_ticks()
    ax.set_ylim(-0.05,1.05)
  ax64.xaxis.set_ticks(np.arange(xrange[0], xrange[1], (xrange[1]- xrange[0])/10.))

    


  ax4.plot(data['x_4'], data['y_4'], '.-k', lw=0.4,marker = 'o', ms = 0.5,label = 'data') #N = 16
  ax8.plot(data['x_8'], data['y_8'], '.-k', lw=0.4,marker = 'o', ms = 0.5,label = 'data8') #N = 16
  ax16.plot(data['x_16'], data['y_16'], '.-k', lw=0.4,marker = 'o', ms = 0.5,label = 'data16') #N = 16
  ax32.plot(data['x_32'], data['y_32'], '.-k', lw=0.4,marker = 'o', ms = 0.5,label = 'data32') #N = 16
  ax64.plot(data['x_64'], data['y_64'], '.-k', lw=0.4,marker = 'o', ms = 0.5,label = 'data64') #N = 16

          


  if disp_sim_spin == True:
    HF_perp, HF_par = fp_funcs.get_hyperfine_params(ms = 'plus')
    B_Field = 403.555      
    tau_lst = np.linspace(0,72e-6,10000)
    colors = cm.rainbow(np.linspace(0, 1, n_sim_spins))
    for tt in range(n_sim_spins):
      Mt4 = SC.dyn_dec_signal(HF_par,HF_perp,B_Field,4,tau_lst)
      FP_signal4 = ((Mt4+1)/2)
      ax4.plot(tau_lst*1e6, FP_signal4[tt,:] ,'-',lw=.8,label = 'spin' + str(tt+1), color = colors[tt])
      ax4.set_title('N = 4')
      Mt8 = SC.dyn_dec_signal(HF_par,HF_perp,B_Field,8,tau_lst)
      FP_signal8 = ((Mt8+1)/2)
      ax8.plot(tau_lst*1e6, FP_signal8[tt,:] ,'-',lw=.8,label = 'spin' + str(tt+1), color = colors[tt])  
      ax8.set_title('N = 8')
      Mt16 = SC.dyn_dec_signal(HF_par,HF_perp,B_Field,16,tau_lst)
      FP_signal16 = ((Mt16+1)/2)
      ax16.plot(tau_lst*1e6, FP_signal16[tt,:] ,'-',lw=.8,label = 'spin' + str(tt+1), color = colors[tt])
      ax16.set_title('N = 16')
      Mt32 = SC.dyn_dec_signal(HF_par,HF_perp,B_Field,32,tau_lst)
      FP_signal32 = ((Mt32+1)/2)
      ax32.plot(tau_lst*1e6, FP_signal32[tt,:] ,'-',lw=.8,label = 'spin' + str(tt+1), color = colors[tt])
      ax32.set_title('N = 32')
      Mt64 = SC.dyn_dec_signal(HF_par,HF_perp,B_Field,64,tau_lst)
      FP_signal64 = ((Mt64+1)/2)
      ax64.plot(tau_lst*1e6, FP_signal64[tt,:] ,'-',lw=.8,label = 'spin' + str(tt+1), color = colors[tt])
      ax64.set_title('N = 64')
                    
    plt.xlabel('tau (us)')
    plt.setp(ax4.get_xticklabels(), visible=False)
    plt.setp(ax8.get_xticklabels(), visible=False)
    plt.setp(ax16.get_xticklabels(), visible=False)
    plt.setp(ax32.get_xticklabels(), visible=False)
    lgd = ax4.legend(bbox_to_anchor=(1.01,-0.5),loc='lower left',borderaxespad = 0.)
  # lgd = plt.legend(loc=4)
  plt.show(block = False)


  

  folder = r'D:\measuring\data\LT2_Data\Fingerprinting'
  plt.savefig(os.path.join(folder, str(disp_sim_spin)+'fingerprint_all_msp1.pdf'),
      format='pdf',bbox_extra_artists = (lgd,),bbox_inches='tight')
  plt.savefig(os.path.join(folder, str(disp_sim_spin)+'fingerprint_all_msp1.png'),
      format='png',bbox_extra_artists = (lgd,),bbox_inches='tight')


if 0:
  plot_all(xrange  = [22,45])

if 1:



  # x_4 , y_4 = fingerprint_contrast_concatenate(disp_sim_spin = True,n_sim_spins = 8,pts = 51,start_1 = 3.0, start_2 = 3+45*50*10e-3,
  #         step_size = 10e-3,
  #         xrange = [2.5,22.5],tag_p = 'Hermite_Fingerprint_msp1_111_1_sil18_4-x',tag_n = 'Hermite_Fingerprint_msp1_111_1_sil18_4x', 
  #         older_than = '20150726090000',return_data = True, do_plot= False, load_from_data = False, Nr_of_pulses = 4 )


  # x_8 , y_8 = fingerprint_contrast_concatenate(disp_sim_spin = True,n_sim_spins = 8,pts = 51,start_1 = 3.0, start_2 = 3+45*50*10e-3,
  #         step_size = 10e-3,
  #         xrange = [2.5,22.5],tag_p = 'Hermite_Fingerprint_msp1_111_1_sil18_8-x',tag_n = 'Hermite_Fingerprint_msp1_111_1_sil18_8x', 
  #         older_than = '20150726090000',return_data = True, do_plot= False, load_from_data = False, Nr_of_pulses = 8 )


  # x_16 , y_16 = fingerprint_contrast_concatenate(disp_sim_spin = True,n_sim_spins = 8,pts = 51,start_1 = 3.0, start_2 = 3+45*50*10e-3,
  #         step_size = 10e-3,
  #         xrange = [2.5,22.5],tag_p = 'Hermite_Fingerprint_msp1_111_1_sil18_16-x',tag_n = 'Hermite_Fingerprint_msp1_111_1_sil18_16x', 
  #         older_than = '20150726090000',return_data = True, do_plot= False, load_from_data = False, Nr_of_pulses = 16)


  x_32 , y_32 = fingerprint_contrast_concatenate(disp_sim_spin = True,n_sim_spins = 0,pts = 51,start_1 = 3.0, start_2 = 3+45*50*10e-3,
          step_size = 10e-3,
          xrange = [9.0,9.5],tag_p = 'Hermite_Fingerprint_msp1_111_1_sil18_32-x',tag_n = 'Hermite_Fingerprint_msp1_111_1_sil18_32x',
          tau_larmor = True, figsize = (5,3),
          older_than = '20150726090000',return_data = True, do_plot= True, load_from_data = True, Nr_of_pulses = 32, name = 'fig_msp1_zoomzoom3',save_folder = r'D:\measuring\data\LT2_Data\Fingerprinting',)

  # x_64 , y_64 = fingerprint_contrast_concatenate(disp_sim_spin = True,n_sim_spins = 8,pts = 51,start_1 = 3.0, start_2 = 3+45*50*4e-3,
  #       step_size = 4e-3,
  #       xrange = [2.5,22.5],tag_p = 'Hermite_Fingerprint_msp1_111_1_sil18_64-x',tag_n = 'Hermite_Fingerprint_msp1_111_1_sil18_64x', 
  #       older_than = '20150726090000',return_data = True, do_plot= False, load_from_data = False, Nr_of_pulses = 64)


  # x_128 , y_128 = fingerprint_contrast_concatenate(disp_sim_spin = True,n_sim_spins = 8,pts = 21,start_1 = 3.0, start_2 = 3+45*20*4e-3,
  #         step_size = 4e-3,
  #         xrange = [2.5,22.5],tag_p = 'Hermite_Fingerprint_msp1_111_1_sil18_128-x',tag_n = 'Hermite_Fingerprint_msp1_111_1_sil18_128x', 
  #         older_than = '20150726090000',return_data = True, do_plot= False, load_from_data = False, Nr_of_pulses = 128)