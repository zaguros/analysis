from scipy.interpolate import interp1d
fn_g= r'H:\My Documents\processed data\2014-3-11 ronald\185905_LaserFrequencyScan_red_scan_go.dat'
fn=r'H:\My Documents\processed data\2014-3-11 ronald\102317_LaserFrequencyScan_red_scan_go.dat'
d=loadtxt(fn)

dd=d[:,1:]
Y=unique(dd[:,2])
X=dd[where(dd[:,2]==Y[1]),0][0]
X=X[where((62.5<X) & (X<63.1))]
Z=zeros((len(X), len(Y)))

for j,y in enumerate(Y):
	fltr=dd[:,2]==y #(dd[:,2]<(y+0.05)) & (dd[:,2]>(j-0.05)) #
	xx=dd[where(fltr),0][0]
	#print len(xx), len(zz)
	zz=dd[where(fltr),1][0]
	f=interp1d(xx,zz, bounds_error=False, fill_value=0.)
	Z[:,j]=f(X)
ax=subplot(111)
ax.imshow(Z, aspect='auto', origin='lower',vmax=np.max(Z)*.75, vmin=0,extent=[min(Y),max(Y),X[0],X[-1]])#, cmap='binary')
ax.set_xlabel('Scan nr')
ax.set_ylabel('Laser frequency [GHz]')
savefig(r'H:\My Documents\processed data\2014-3-11 ronald\3.2.png')


