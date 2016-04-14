from analysis.lib.fitting import fit,esr
from analysis.lib.tools import plot
reload(plot)
import os


guess_ctr = 2.82174#2.861

guess_splitN = 2.176e-3
guess_width = .0005
guess_splitB = 0.0008
guess_amplitude = 0.13
guess_Nsplit = guess_splitN

folder=tb.latest_data('ESR')
fn=tb.measurement_filename(folder, ext='dat')
data=np.loadtxt(fn, skiprows=10)
x=data[:,0]*1e-9
y=data[:,1]
guess_offset = max(y)
guess_amplitude = max(y)-min(y)
fig=plt.figure(figsize=(6,4))
ax=subplot(111)
#ax.plot(x,y)
x_f=np.linspace(x[0],x[-1],2000)*1e-9



fit_result = fit.fit1d(x, y, esr.fit_ESR_gauss, guess_offset,
         guess_amplitude, guess_width, guess_ctr,
                (3, guess_splitN),
                #(2, guess_splitB),
         do_print=True, ret=True, fixed=[4])
plot.plot_fit1d(fit_result, np.linspace(min(x), max(x), 1000), ax=ax, plot_data=True)
        

print' Minimum is ', np.min(y)
print 'at', x[np.argmin(y)] , 'GHz'
  
ax.set_xlabel('MW frq (GHz)')
ax.set_ylabel(r'fidelity wrt. $|0\rangle$')
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