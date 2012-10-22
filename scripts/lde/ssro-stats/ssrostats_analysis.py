import os, sys, time
import pickle
import pprint
import numpy as np
from matplotlib import pyplot as plt


class SSROStatAnalysis:

    def __init__(self):
        self.srcfolder = '/Users/wp/Documents/TUD/LDE/data/ssro-stats'
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
                    }
                }

        self.fiderr = {
                'LT1' : {
                    0 : [],
                    1 : [],
                    },
                'LT2' : {
                    0 : [],
                    1 : [],
                    }
                }

        self.pathinfo = {
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
                    timestamp = folder[:4]
                    b,date = os.path.split(b)[-4:]

                    setup = 'LT1' if path [-3:] == self.lt1suffix else 'LT2'
                    state = 0 if fn == self.fname0 else 1
                    idx = self.lt1idx if path[-3:] == self.lt1suffix else self.lt2idx

                    self.fid[setup][state].append(fid[idx])
                    self.fiderr[setup][state].append(fiderr[idx])
                    self.pathinfo[setup].append(date+'/'+timestamp)

    def stats(self):
        
        self.meanfid = {
                'LT1' : {
                    0 : 0.,
                    1 : 0.,
                    },
                'LT2' : {
                    0 : 0.,
                    1 : 0.,
                    }
                }

        self.fidstdev = {
                'LT1' : {
                    0 : 0.,
                    1 : 0.,
                    },
                'LT2' : {
                    0 : 0.,
                    1 : 0.,
                    }
                }

        self.fiderrsumofsquares = {
                'LT1' : {
                    0 : 0.,
                    1 : 0.,
                    },
                'LT2' : {
                    0 : 0.,
                    1 : 0.,
                    }
                }


        for setup in self.fid:
            for state in self.fid[setup]:
                print setup, state
                print 'mean = %.4f' % np.mean(np.array(self.fid[setup][state]))
                print 'statistical uncertainty = %.4f' % \
                        np.std(np.array(self.fid[setup][state]))
                print 'error propagation uncertainty = %.4f' % \
                        np.sqrt(np.sum((np.array(self.fiderr[setup][state])/len(self.fiderr[setup][state]))**2))
                print ''

                self.meanfid[setup][state] = np.mean(np.array(self.fid[setup][state]))
                self.fidstdev[setup][state] = np.std(np.array(self.fid[setup][state]))
                self.fiderrsumofsquares[setup][state] = \
                        np.sqrt(np.sum((np.array(self.fiderr[setup][state])/len(self.fiderr[setup][state]))**2))


        fig = plt.figure(figsize=(10,10))
        
        for i, setup, state in zip(range(1,5), ['LT1', 'LT1', 'LT2', 'LT2'], [0,1,0,1]):

            ax = plt.subplot(2,2,i)
            ax.hist(self.fid[setup][state], color='k', histtype='step', hatch='/')
            ax.set_title(setup + ', ms = ' + str(state))
            plt.text(ax.get_xlim()[0]+0.0002,1.,'%.4f +/- %.4f' % \
                    (self.meanfid[setup][state], self.fidstdev[setup][state]),
                    backgroundcolor='w')

            ax.set_xlabel('Fidelity')
            ax.set_ylabel('Occurrences')

        fig = plt.figure(figsize=(10,10))

        for i, setup, state in zip(range(1,5), ['LT1', 'LT1', 'LT2', 'LT2'], [0,1,0,1]):

            ax = plt.subplot(2,2,i)
            ax.errorbar(np.arange(len(self.fid[setup][state])), np.array(self.fid[setup][state]),
                    yerr=np.array(self.fiderr[setup][state]), fmt='ko')
            ax.axhline(y=self.meanfid[setup][state], ls=':')
            ax.set_title(setup + ', ms = ' + str(state))
            ax.set_xticks(np.arange(len(self.fid[setup][state])))
            ax.set_xticklabels(self.pathinfo[setup])
            labels = ax.get_xticklabels()
            plt.setp(labels, rotation=90)

            ax.set_ylabel('Fidelity')
            ax.set_xlabel('Run')

        plt.tight_layout()


if __name__ == '__main__':
    a = SSROStatAnalysis()
    a.load_data()
    a.stats()

    f = open('fidelities.pkl', 'wb')
    pickle.dump(a.meanfid, f)
    f.close()

    f = open('stdevs.pkl', 'wb')
    pickle.dump(a.fidstdev, f)
    f.close()



