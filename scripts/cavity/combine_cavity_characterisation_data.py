#### script to get the cavity linewidth data and combine it with the corresponding cavity length in mode number (N)
#### SvD, September 2016

import os
import numpy as np
import scipy.signal
import matplotlib.pyplot as plt
import analysis.scripts.cavity.functions_for_processing_linewidth_data as funcs

datapath = 'K:ns\qt\Diamond\Projects\Cavities\Cavity characterisation paper\data/data_for_cav_char_paper/'
datapath = 'W:Diamond\Projects\Cavities\Cavity characterisation paper\data/data_for_cav_char_paper/'




def get_RT_OFFD_POS0_data(plot_data=True):
    dir_L0 =os.path.join(datapath,'20160818/RT_OFFD_POS0_L0')
    dir_L2= os.path.join(datapath,'20160818/RT_OFFD_POS0_L2')
    dir_L3 = os.path.join(datapath,'20160819/RT_OFFD_POS0_L3')
    #I relabeled the data in the dedicated folder (my docs) from V1p129 -> V0 etc. SvD 15-9-2016
    V0s=[0,1,2]
    V2s=[0,1,2,3]
    V3s=[0,1,2,3,4,5]

    #V0s = np.array([1.129,2.473,3.482])
    #V2s = np.array([6.083,7.258,8.269,9.246])
    #V3s = np.array([1.402,2.646,3.781,4.846,5.864,6.890])
    N0s = np.array([21,22,23])#possibly the N0s should be 21,22,23->not 100%clear. -> leave out?
    N2s = np.array([13,14,15,16])
    N3s = np.array([17,18,19,20,21,22])

    folder,avgs0=funcs.get_avg_array_for_Vs(V0s,1,dir_L0,tag='')
    folder,avgs2=funcs.get_avg_array_for_Vs(V2s,2,dir_L2,tag='')
    folder,avgs3=funcs.get_avg_array_for_Vs(V3s,3,dir_L3,tag='')


    RT_OFFD_POS0_avgs = np.concatenate((avgs0,avgs2,avgs3))
    RT_OFFD_POS0_Ns = np.append(np.append(N0s,N2s),N3s)
    if plot_data:
        ax = funcs.plot_updown_lws(RT_OFFD_POS0_avgs,RT_OFFD_POS0_Ns)
        ax.set_title(folder)
        plt.savefig(os.path.join(folder,'plot_lw_vs_L.png'))
        plt.show()

        ax = funcs.plot_finesse_vs_length_from_avgs(RT_OFFD_POS0_avgs,RT_OFFD_POS0_Ns)
        ax.set_title(folder)
        plt.savefig(os.path.join(folder,'Finesse_vs_L.png'))
        plt.show()



    return RT_OFFD_POS0_avgs,RT_OFFD_POS0_Ns



def get_LT_OFFD_POS3_data(plot_data=True):
    #LT OFF DIAMOND , POS 3
    dir_L0 = os.path.join(datapath,'20160907/LT_POS3_OFFD_L0')

    N0s= [17]
    V0s = [0]

    folder,avg_array0,u_avg_array0 = funcs.get_LT_avg_array_from_json(V0s ,dir_L0)

    LT_OFFD_POS3_avg_array = avg_array0
    LT_OFFD_POS3_u_avg_array = u_avg_array0
    LT_OFFD_POS3_Ns = N0s

    if plot_data:
        ax = funcs.plot_finesse_vs_length(LT_OFFD_POS3_avg_array,LT_OFFD_POS3_u_avg_array,LT_OFFD_POS3_Ns)
        funcs.plot_mean_lws(LT_OFFD_POS3_avg_array,LT_OFFD_POS3_u_avg_array,LT_OFFD_POS3_Ns)

    return LT_OFFD_POS3_avg_array,LT_OFFD_POS3_u_avg_array,LT_OFFD_POS3_Ns


def get_RT_OND_POS0_data(plot_data=True):
    dir_L0 = os.path.join(datapath,'20160819/RT_OND_POS0_L0')
    dir_L1 = os.path.join(datapath,'20160819/RT_OND_POS0_L1')
    dir_L2 = os.path.join(datapath,'20160819/RT_OND_POS0_L2')

    V0s=np.arange(0,4,1)
    V1s=np.arange(0,8,1)
    V2s=np.arange(0,3,1)

    #V0s=np.array([2.935,4.088,5.183,6.182])
    #V1s=np.array([2.139,3.587,4.753,5.836,6.891,7.928,8.909,9.958])
    #V2s=np.array([6.957,8.738,9.866])

    N0s = np.array([54,55,56,57])
    N1s = np.arange(43,43+8,1)#np.array([44,45,46,47,48,49,50,51])
    N2s = np.array([51,52,53])#np.array([51,52,53])

    folder,avgs0=funcs.get_avg_array_for_Vs(V0s,1,dir_L0,tag='')
    folder,avgs1=funcs.get_avg_array_for_Vs(V1s,1,dir_L1,tag='')
    folder,avgs2=funcs.get_avg_array_for_Vs(V2s,1,dir_L2,tag='')

    RT_OND_POS0_avgs = np.concatenate((avgs1,avgs2,avgs0))
    RT_OND_POS0_Ns = np.concatenate((N1s,N2s,N0s))

    if plot_data:
        ax = funcs.plot_updown_lws(RT_OND_POS0_avgs,RT_OND_POS0_Ns)

        ax.set_title(folder)
        plt.savefig(os.path.join(folder,'plot_lw_vs_L.png'))
        plt.show()

        funcs.plot_finesse_vs_length_from_avgs(RT_OND_POS0_avgs,RT_OND_POS0_Ns)

    return RT_OND_POS0_avgs,RT_OND_POS0_Ns


def get_RT_OND_POS1_data(plot_data=True):
    dir_L2 = os.path.join(datapath,'20160822/RT_OND_POS1_L2')
    dir_L3 = os.path.join(datapath,'20160822/RT_OND_POS1_L3')
    dir_L4 = os.path.join(datapath,'20160822/RT_OND_POS1_L4')
    dir_L5 = os.path.join(datapath,'20160822/RT_OND_POS1_L5')

    V2s=np.arange(0,8,1)
    V3s=np.arange(0,8,1)
    V4s=np.arange(0,8,1)
    V5s=np.arange(0,8,1)

    #V0s=np.array([1.546,3.019,4.237,5.348,6.404,7.473,8.515,9.551])
    #V3s=np.array([1.373,2.859,4.062,5.152,6.218,7.224,8.198,9.268])
    #V4s=np.array([0.951,2.390,3.646,4.783,5.843,6.934,7.975,8.977])
    #V5s=np.array([1.359,3.383,4.565,5.620,6.695,7.727,8.734,9.832])

    N2s = np.arange(54,62,1)#55,63#53,61
    N3s = np.arange(49,57,1)#49,57#46,55
    N4s = np.arange(44,52,1)
    N5s = np.arange(39,47,1)

    folder,avgs2=funcs.get_avg_array_for_Vs(V2s,1,dir_L2,tag='')
    folder,avgs3=funcs.get_avg_array_for_Vs(V3s,1,dir_L3,tag='')
    folder,avgs4=funcs.get_avg_array_for_Vs(V4s,1,dir_L4,tag='')
    folder,avgs5=funcs.get_avg_array_for_Vs(V5s,1,dir_L5,tag='')

    RT_OND_POS1_avgs=np.concatenate((avgs2,avgs3,avgs4,avgs5))
    RT_OND_POS1_Ns=np.concatenate((N2s,N3s,N4s,N5s))

    if plot_data:
        ax = funcs.plot_updown_lws(RT_OND_POS1_avgs,RT_OND_POS1_Ns)
        ax.set_title(folder)
        plt.savefig(os.path.join(folder,'plot_lw_vs_L.png'))
        plt.show()

        funcs.plot_finesse_vs_length_from_avgs(RT_OND_POS1_avgs,RT_OND_POS1_Ns)
        ax.set_title(folder)
        plt.savefig(os.path.join(folder,'Finesse_vs_L.png'))
        plt.show()

    return RT_OND_POS1_avgs,RT_OND_POS1_Ns


def get_RT_OND_POS2_data(plot_data=True):
    dir_L2 = os.path.join(datapath,'20160825/OND_POS2L2')
    dir_L3 = os.path.join(datapath,'20160825/OND_POS2L3')
    dir_L4 = os.path.join(datapath,'20160825/OND_POS2L4')
    dir_L5 = os.path.join(datapath,'20160825/OND_POS2L5')
    dir_L6 = os.path.join(datapath,'20160825/OND_POS2L6')
    dir_L7 = os.path.join(datapath,'20160825/OND_POS2L7')
    dir_L8 = os.path.join(datapath,'20160825/OND_POS2L8')

    L2_N0 = 70#73
    L3_N0 = 65#68
    L4_N0 = 59#62
    L5_N0 = 53#57
    L6_N0 = 48#48#49#48
    L7_N0 = 43#44
    L8_N0 = 43#42

    nr_files_L3= 9
    V2s=np.arange(0,0+9,1)
    V3s=np.arange(0,0+9,1)
    V4s=np.arange(0,0+9,1)
    V5s=np.arange(0,0+9,1)
    V6s=np.arange(0,0+9,1)
    V7s=np.arange(0,0+9,1)
    V8s=np.arange(0,0+7,1)


    N2s = np.arange(L2_N0 ,L2_N0+9,1)#55,63
    N3s= np.arange(L3_N0 ,L3_N0+9,1)
    N4s= np.arange(L4_N0 ,L4_N0+9,1)
    N5s= np.arange(L5_N0 ,L5_N0+9,1)
    N6s= np.arange(L6_N0 ,L6_N0+9,1)
    N7s= np.arange(L7_N0 ,L7_N0+9,1)
    N8s= np.arange(L8_N0 ,L8_N0+7,1)


    folder,avgs2=funcs.get_avg_array_for_Vs(V2s,1,dir_L2,tag='')
    folder,avgs3=funcs.get_avg_array_for_Vs(V3s,1,dir_L3,tag='')
    folder,avgs4=funcs.get_avg_array_for_Vs(V4s,1,dir_L4,tag='')
    folder,avgs5=funcs.get_avg_array_for_Vs(V5s,1,dir_L5,tag='')
    folder,avgs6=funcs.get_avg_array_for_Vs(V6s,1,dir_L6,tag='')
    folder,avgs7=funcs.get_avg_array_for_Vs(V7s,1,dir_L7,tag='')
    folder,avgs8=funcs.get_avg_array_for_Vs(V8s,1,dir_L8,tag='')

    #RT_OND_POS2_avgs = np.append(np.append(np.append(np.append(np.append(np.append(avgs2,avgs3,axis=0),avgs4,axis=0),avgs5,axis=0),avgs6,axis=0),avgs7,axis=0),avgs8,axis=0)
    RT_OND_POS2_avgs = np.concatenate((avgs2,avgs3,avgs4,avgs5,avgs6,avgs7,avgs8))
    RT_OND_POS2_Ns = np.concatenate((N2s,N3s,N4s,N5s,N6s,N7s,N8s))
    if plot_data:
        ax = funcs.plot_updown_lws(RT_OND_POS2_avgs,RT_OND_POS2_Ns)
        ax.set_title(folder)
        plt.savefig(os.path.join(folder,'plot_lw_vs_L.png'))
        plt.show()
        ax = funcs.plot_finesse_vs_length_from_avgs(RT_OND_POS2_avgs,RT_OND_POS2_Ns)
        ax.set_title(folder)
        plt.savefig(os.path.join(folder,'Finesse_vs_L.png'))
        plt.show()

    return RT_OND_POS2_avgs,RT_OND_POS2_Ns


def get_LT_OND_POS4_data(plot_data=True):
    dir_L2 = os.path.join(datapath,'20160909/LT_OND_POS4_L2')
    dir_L3 = os.path.join(datapath,'20160909/LT_OND_POS4_L3')
    dir_L4 = os.path.join(datapath,'20160909/LT_OND_POS4_L4')
    dir_L5 = os.path.join(datapath,'20160909/LT_OND_POS4_L5')

    d_diamond_pos4 = 4.22e-6
    N2s= [43,44,45]#[42,43,44]#this one I still have t check. Maybe keep out Vmin=3
    N3s = [44,45,46,47]
    N4s = [46,47,48,49]#[45,46,47,48]
    N5s=[48,49,50,51]
    V2s=[0,1,2]
    V3s=[0,1,2,3]
    V4s=[0,1,2,3]
    V5s = [0,1,2,3]
    folder,avg_array2,u_avg_array2 = funcs.get_LT_avg_array_from_json(V2s ,dir_L2)
    folder,avg_array3,u_avg_array3 = funcs.get_LT_avg_array_from_json(V3s ,dir_L3)
    folder,avg_array4,u_avg_array4 = funcs.get_LT_avg_array_from_json(V4s ,dir_L4)
    folder,avg_array5,u_avg_array5 = funcs.get_LT_avg_array_from_json(V5s ,dir_L5)


    LT_OND_POS4_avg_array = np.concatenate((avg_array2,avg_array3,avg_array4,avg_array5))
    LT_OND_POS4_u_avg_array =np.concatenate((u_avg_array2,u_avg_array3,u_avg_array4,u_avg_array5))  
    LT_OND_POS4_Ns = np.concatenate((N2s,N3s,N4s,N5s))

    if plot_data:
        ax = funcs.plot_finesse_vs_length(LT_OND_POS4_avg_array,LT_OND_POS4_u_avg_array,LT_OND_POS4_Ns)
        funcs.plot_mean_lws(LT_OND_POS4_avg_array,LT_OND_POS4_u_avg_array,LT_OND_POS4_Ns)

    return LT_OND_POS4_avg_array, LT_OND_POS4_u_avg_array, LT_OND_POS4_Ns

def get_LT_OND_POS2_data(plot_data=True):

    #LT ON DIAMOND , POS 2
    (funcs)
    dir_L5 = os.path.join(datapath,'20160905/LT_OND_POS2_L5')
    dir_L6 = os.path.join(datapath,'20160905/LT_OND_POS2_L6')
    dir_L7 = os.path.join(datapath,'20160905/LT_OND_POS2_L7')

    N2s= [69]
    avg_array2 = [1.92]
    u_avg_array2 = [0.22]
    V5s= [0,1,2,3]
    N5s = [52,53,54,55]
    V6s=[0,1,2,3]
    N6s=[47,48,49,50]
    V7s = [0,1,2,3]
    N7s=[43,44,45,46]


    folder,avg_array5,u_avg_array5 = funcs.get_LT_avg_array_from_json(V5s ,dir_L5)
    folder,avg_array6,u_avg_array6 = funcs.get_LT_avg_array_from_json(V6s ,dir_L6)
    folder,avg_array7,u_avg_array7 = funcs.get_LT_avg_array_from_json(V7s ,dir_L7)


    LT_OND_POS2_avg_array = np.concatenate((avg_array5,avg_array6,avg_array7,avg_array2)) 
    LT_OND_POS2_u_avg_array =  np.concatenate((u_avg_array5,u_avg_array6,u_avg_array7,u_avg_array2)) 
    LT_OND_POS2_Ns = np.concatenate((N5s,N6s,N7s,N2s))
    
    if plot_data:
        ax = funcs.plot_finesse_vs_length(LT_OND_POS2_avg_array,LT_OND_POS2_u_avg_array,LT_OND_POS2_Ns)
        funcs.plot_mean_lws(LT_OND_POS2_avg_array,LT_OND_POS2_u_avg_array,LT_OND_POS2_Ns)

    return LT_OND_POS2_avg_array, LT_OND_POS2_u_avg_array, LT_OND_POS2_Ns


