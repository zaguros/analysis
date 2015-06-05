from analysis.lib.fitting import fit,esr
import os



folder=tb.latest_data('ESR')
fn=tb.measurement_filename(folder, ext='dat')
data=np.loadtxt(fn, skiprows=10)
x=data[:,0]
y=data[:,1]
fig=plt.figure(figsize=(6,4))
ax=plt.subplot(111)
ax.plot(x,y)
x_f=np.linspace(x[0],x[-1],2000)
fit_result = fit.fit1d(x, y, esr.fit_ESR_gauss, np.max(y),
        np.max(y)-np.min(y), (x[-1]-x[0])/4, x[np.argmin(y)],
        do_print=False, ret=True, do_plot = False, fixed=[4])
y_f = fit_result['fitfunc'](x_f)

ax.plot(x_f,y_f)
ax.set_title(tb.get_plot_title_from_folder(folder)+ '\n x0: {:.4f} +/- {:.4f} GHz'.format(fit_result['params_dict']['x0']/1e9, fit_result['error_dict']['x0']/1e9))
fig.savefig(os.path.splitext(fn)[0]+'.png', format='png' )