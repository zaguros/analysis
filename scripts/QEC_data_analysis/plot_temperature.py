from analysis.lib.fitting import fit,esr
from analysis.lib.tools import plot
from analysis.lib.tools import toolbox as tb
from matplotlib import pyplot as plt
reload(plot)
import os

# folder=tb.latest_data('monitor_warmup', older_than='20160723_231158')
folder=tb.latest_data('temp_monitor', older_than='20160725_170000')
fn=tb.measurement_filename(folder, ext='dat')
data=np.loadtxt(fn, skiprows=16)
x=data[:,0]
y=(data[:,1]-100)/0.385
fig=plt.figure(figsize=(20,8))
ax=plt.subplot(111)
ax.plot(x,y)
# x_f=np.linspace(x[0],x[-1],2000)*1e-9
plt.xlim((0,65./60))
# plt.ylim((1.0857e2+0.02,1.0857e2+0.04))
       
ax.set_xlabel('time (hours)')
ax.set_ylabel('Resistance')
ax.set_title(folder)
#ax.set_ylim([0,0.05])
plt.savefig(os.path.join(folder, 'darkesr_analysis.png'),
        format='png')
'''

fit_result = fit.fit1d(x, y, esr.fit_ESR_gauss, np.max(y),
        np.max(y)-np.min(y), (x[-1]-x[0])/4, x[np.argmin(y)],
        do_print=False, ret=True, do_plot = False, fixed=[4])
y_f = fit_result['fitfunc'](x_f)

ax.plot(x_f,y_f)
ax.set_title(tb.get_plot_title_from_folder(folder)+ '\n x0: {:.4f} +/- {:.4f} GHz'.format(fit_result['params_dict']['x0']/1e9, fit_result['error_dict']['x0']/1e9))
fig.savefig(os.path.splitext(fn)[0]+'.png', format='png' )
'''