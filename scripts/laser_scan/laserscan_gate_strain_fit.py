from analysis.lib.fitting import fit,common
from analysis.lib.nv import nvlevels
from analysis.lib.tools import clipboard as cb
from matplotlib import pyplot as plt
plot_all=False
fit_all=False
plt.close('all')
d=np.loadtxt(fn)
gate=unique(d[:,3])

def plot_laserscan(d,gv):
        x=d[d[:,3]==gv,1]
        y=d[d[:,3]==gv,2]
        plt.figure()
        plt.subplot(1,1,1)
        plt.plot(x,y)
        try:
            gx= float(raw_input('Ey?'))
            gy= float(raw_input('Ex?'))
            x1,y2=nvlevels.get_ES_ExEy_plottable(gx,gy,np.max(y))
            plt.plot(x1,y2)
        except:
            print 'could not understand input'
        

print 'gate voltages:'
for i,gv in enumerate(gate):
    print i,':',gv

    if plot_all:
        plot_laserscan(d,gv)

if not(plot_all):
    i=int(raw_input('Plot which gate voltage?'))
    plot_laserscan(d,gate[i])

if fit_all:

    strain_range=(3.5,12)
    min_counts=300

    result=[]
    gate_filtered=[]
    exs=[]
    eys=[]

    for gv in gate:
        x=d[d[:,3]==gv,1]
        y=d[d[:,3]==gv,2]
        if max(y)<min_counts:
            continue
        gate_filtered.append(gv)
        (ex,ey)=nvlevels.fit_laserscan(x,y,fast=True,strain_range=strain_range)
        result.append(nvlevels.get_ES_ExEy(ex,ey,fast=True))
        exs.append(ex)
        eys.append(ey)
        
    arr=np.vstack((array(gate_filtered),array(eys),array(exs)))
    cb.copy_array_2d(arr)
