from analysis.lib.lde import tail_cts_per_shot_v4
from analysis.lib.fitting import rabi

path = r'D:\measuring\data\LDE\analysis_data\opt_rabi_vs_CR\2012-10-24-CRvsRabi_LT1\183600_rabi_vs_cr\20nW'

files = [ f for f in os.listdir(path) if f[:1] != '.' and f[-4:] == '.txt' ]


tau_list=[]
tau_err_list=[]
fR_list=[]
fR_err_list=[]
CR=[]

for i in files:
    print i
    fR, fR_err, tau, tau_err = tail_cts_per_shot_v4.optical_rabi_resonant_fit_hannes(
            os.path.join(path,i),32,76, use_channel=0, binsize=0.128)
    tau_list.append(tau)
    tau_err_list.append(tau_err)
    fR_list.append(fR*1e3)
    fR_err_list.append(fR_err*1e3)
    CR.append(float(i[:2]))

plt.close('all')

fig=plt.figure()
plt.errorbar(CR,tau_list,yerr=tau_err_list,fmt='ro')
plt.xlabel('CR threshold')
plt.ylabel('Decay time [ns]')
plt.title(path)
fig.savefig(os.path.join(path,'CRvsTau.png'))


fig2=plt.figure()
plt.errorbar(CR,fR_list,yerr=fR_err_list,fmt='ro')
plt.xlabel('CR threshold')
plt.ylabel('Rabi frequency [MHz]')
plt.title(path)
fig2.savefig(os.path.join(path,'CRvsRabiF.png'))
