import os, sys, time
import pickle
import pprint
import numpy as np
from matplotlib import pyplot as plt
from analysis import config

config.outputdir = r'/Users/wp/Documents/TUD/LDE/analysis/output'
config.datadir = r'/Volumes/MURDERHORN/TUD/LDE/analysis/data/ssro-stats_withbasis'


class SSROStatAnalysis:

    def __init__(self):
        self.srcfolder = config.datadir
        self.lde_ssro_durations_fp = os.path.join(config.outputdir, 
                '20121105-ssro-stats', 'ssro_durations.npz')
        self.savedir = os.path.join(config.outputdir, time.strftime('%Y%m%d')+'-ssro-stats')
        self.fname0 = 'ADwin_SSRO-000_fid_vs_ROtime.npz'
        self.fname1 = 'ADwin_SSRO-001_fid_vs_ROtime.npz'
        self.lt1suffix = 'LT1'
        self.lt2suffix = 'LT2'

        self.fid = {
                'LT1' : {
                    0 : [],
                    1 : [],
                    },
                'LT2' : {
                    0 : [],
                    1 : [],
                    },
                }

        self.fiderr = {
                'LT1' : {
                    0 : [],
                    1 : [],
                    },
                'LT2' : {
                    0 : [],
                    1 : [],
                    },
                }

        self.pathtime = {
                'LT1' : [],
                'LT2' : [],
                }

        self.lt1idx = 14
        self.lt2idx = 25

    def load_data(self):
        
        print 'start reading data'

        for (path,dirs,files) in os.walk(self.srcfolder):
            for fn in files:
                if fn == self.fname0 or fn == self.fname1:
                    d = np.load(os.path.join(path, fn))
                    fid = d['fid']
                    fiderr = d['fid_error']
                    d.close()

                    b,folder = os.path.split(path)
                    timestamp = folder[:6]
                    b,_date = os.path.split(b)
                    date = _date[:8]

                    if path[-3:] == self.lt1suffix:
                        setup = 'LT1'
                    else:
                        setup = 'LT2'
               
                    datetimestr = date+timestamp
                    t = int(time.mktime(time.strptime(datetimestr, '%Y%m%d%H%M%S')))
                    
                    if fn == self.fname0:
                        self.pathtime[setup].append(t)

                    state = 0 if fn == self.fname0 else 1
                    lt1duration, lt2duration = self.get_readoutduration(t)
                    idx = lt1duration-1 if path[-3:] == self.lt1suffix else lt2duration-1

                    self.fid[setup][state].append(fid[idx])
                    self.fiderr[setup][state].append(fiderr[idx])
    
    def get_readoutduration(self, ssrotime):
        f = np.load(self.lde_ssro_durations_fp)
        t = f['times']
        lt1 = f['lt1durations']
        lt2 = f['lt2durations']
        f.close()
        
        _t = (t-ssrotime) > 0
        closestidx = argmin(abs(t[_t]-ssrotime))

        print time.strftime('%m%d / %H%M', time.localtime(ssrotime)), '--->', 
        print time.strftime('%m%d / %H%M', time.localtime(t[closestidx]))
        
        return int(lt1[closestidx]), int(lt2[closestidx])
    
    def stats(self):
        
        self.meanfid = {
                'LT1' : {
                    0 : 0.,
                    1 : 0.,
                    },
                'LT2' : {
                    0 : 0.,
                    1 : 0.,
                    },
                }

        self.fidstdev = {
                'LT1' : {
                    0 : 0.,
                    1 : 0.,
                    },
                'LT2' : {
                    0 : 0.,
                    1 : 0.,
                    },
                }

        for setup in self.fid:
            for state in self.fid[setup]:
                print setup, state
                print 'mean = %.4f' % np.mean(np.array(self.fid[setup][state]))
                print 'statistical uncertainty = %.4f' % \
                        np.std(np.array(self.fid[setup][state]))
                print ''

                self.meanfid[setup][state] = np.mean(np.array(self.fid[setup][state]))
                self.fidstdev[setup][state] = np.std(np.array(self.fid[setup][state]))
                

        fig = plt.figure(figsize=(10,10))
        
        for i, setup, state in zip(range(1,5), ['LT1', 'LT1', 'LT2', 'LT2'], [0,1,0,1]):

            ax = plt.subplot(2,2,i)
            ax.hist(self.fid[setup][state], color='k', histtype='step', hatch='/')           
            ax.set_title(setup + ', ms = ' + str(state))
            plt.text(ax.get_xlim()[0]+0.0004,4.,'%.4f +/- %.4f' % \
                    (self.meanfid[setup][state], self.fidstdev[setup][state]),
                    backgroundcolor='w')
            ax.locator_params(nbins=4)
            ax.set_xlabel('Fidelity')
            ax.set_ylabel('Occurrences')

            if i == 4:
                xvals = np.linspace(0.94, 1., 101)
                yvals = 13*np.exp(-(xvals-0.993)**2/2./(0.03**2/2./np.log(13./6.)))
                ax.plot(xvals, yvals, 'r-')
                ax.text(0.95,10., 'sigma=%.4f' % (np.sqrt(0.03**2/2./np.log(13./6.))),
                    color='r')

                ax.set_xlim(0.94, 1.)

        plt.tight_layout()

        # fig = plt.figure(figsize=(10,10))

        # for i, setup, state in zip(range(1,5), ['LT1', 'LT1', 'LT2', 'LT2'], [0,1,0,1]):

        #     ax = plt.subplot(2,2,i)
        #     ax.errorbar(np.arange(len(self.fid[setup][state])), np.array(self.fid[setup][state]),
        #             yerr=np.array(self.fiderr[setup][state]), fmt='ko')
        #     ax.axhline(y=self.meanfid[setup][state], ls=':')
        #     ax.set_title(setup + ', ms = ' + str(state))
        #     ax.set_xticks(np.arange(len(self.fid[setup][state])))
        #     ax.set_xticklabels(self.pathinfo[setup])
        #     labels = ax.get_xticklabels()
        #     plt.setp(labels, rotation=90)

        #     ax.set_ylabel('Fidelity')
        #     ax.set_xlabel('Run')

        # plt.tight_layout()


if __name__ == '__main__':
    a = SSROStatAnalysis()
    a.load_data()
    a.stats()

    if not os.path.exists(a.savedir):
        os.makedirs(a.savedir)

    f = open(os.path.join(a.savedir, 'mean_fidelities.pkl'), 'wb')
    pickle.dump(a.meanfid, f)
    f.close()

    f = open(os.path.join(a.savedir, 'mean_stdevs.pkl'), 'wb')
    pickle.dump(a.fidstdev, f)
    f.close()

    f = open(os.path.join(a.savedir, 'fidelities.pkl'), 'wb')
    pickle.dump(a.fid, f)
    f.close()

    f = open(os.path.join(a.savedir, 'stdevs.pkl'), 'wb')
    pickle.dump(a.fiderr, f)
    f.close()




