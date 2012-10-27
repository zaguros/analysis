import os, sys, time
import pickle
import pprint
import numpy as np
from matplotlib import pyplot as plt



class StarkAnalysis:

    def __init__(self):

        self.srcfolder = r'K:\ns\qt\Diamond\Data\LDE\analysis_data\stark_tuning'
        self.scanidf = r'LT25nW_mw_True_0nW_green'
        self.common_freq_axis = np.linspace(40,90,5001)

    def find_nearest(self, array,value):
        idx=(abs(array-value)).argmin()
        return idx

    def get_data(self,file_name):
        d=np.loadtxt(file_name)
        freq=d[:,1]
        counts=d[:,2]
        return freq, counts

    def change_freq_axis(self, freq_axis, counts):
        self.new_counts=np.ones(len(new_axis))*-1

        for i in np.arange(len(counts)):
            idx=self.find_nearest(self.common_freq_axis,freq_axis[i])
            try:
                if self.new_counts[idx]<0:
                    self.new_counts[idx]=counts[i]
                else:
                    self.new_counts[idx]+=counts[i]
            except:
                print 'something went wrong with index ', idx


        return self.new_counts

    def get_files(self):
        self.files=[]
        for d,sd,f in os.walk(self.srcfolder):
            self.files=np.hstack((self.files,
                [i for i in f if self.scanidf in i and i[:1] != '.' and i[-5:]== '1.dat']))

        return self.files

if __name__ == '__main__':
    a = StarkAnalysis()

