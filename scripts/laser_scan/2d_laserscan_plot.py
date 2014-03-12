from scipy.interpolate import interp1d
fn= r'Z:\data_000000-201306\201304-06\20130529\184353_LaserFrequencyScan_red_scan_gv_0.0\184353_LaserFrequencyScan_red_scan_gv_0.0.dat'
fn2= r'Z:\data_000000-201306\201304-06\20130529\212721_LaserFrequencyScan_red_scan__gv_0.0\212721_LaserFrequencyScan_red_scan__gv_0.0.dat'
#fn=r'Z:\data_000000-201306\201304-06\20130529\184354_LaserFrequencyScan_red_scan_gv_0.0_yellow\184354_LaserFrequencyScan_red_scan_gv_0.0_yellow.dat'
#fn2=r'Z:\data_000000-201306\201304-06\20130529\212721_LaserFrequencyScan_red_scan__gv_0.0_yellow\212721_LaserFrequencyScan_red_scan__gv_0.0_yellow.dat'
d=loadtxt(fn)
d2=loadtxt(fn2)

dd=vstack((d[:,1:],d2[:,1:]))
Y=unique(dd[:,2])
X=dd[where(dd[:,2]==Y[1]),0][0]
Z=zeros((len(X), len(Y)))

for j,y in enumerate(Y):
	fltr=dd[:,2]==y #(dd[:,2]<(y+0.05)) & (dd[:,2]>(j-0.05)) #
	xx=dd[where(fltr),0][0]
	#print len(xx), len(zz)
	zz=dd[where(fltr),1][0]
	f=interp1d(xx,zz, bounds_error=False, fill_value=0.)
	Z[:,j]=f(X)
ax=subplot(111)
ax.imshow(Z, aspect='auto', origin='lower',vmax=1000, vmin=0, cmap='binary',extent=[min(Y)*30.,max(Y)*30.,X[0],X[-1]])
ax.set_xlabel('Gate voltage [V]')
ax.set_ylabel('Laser frequency [GHz]')
savefig(r'H:\My Documents\processed data\2014-3-11 ronald\2.1.png')


