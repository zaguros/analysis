### Functions for analysing oscilloscope data

import numpy as np 
from matplotlib import pyplot as plt
import pandas as pd

def load_data(filepath):
    """
    Loading csv -files into pandas dataframe and convert to array
    Input: 
    filepath  - the filepath of csv file
    Output: 
    data  -  numpy ndarray with (time, intensity)
    """
    data = pd.read_csv(filepath, skiprows=16, names = ["none","X","Y","none"],usecols=[1,2]) #creating a dataframe in pandas and importing the data
    data = data.as_matrix()
    return data





# import os
# import sys
# import numpy as np
# sys.path.append("H:\My Documents\measuring/")

# %matplotlib inline

# import analysis.scripts.cavity.oscilloscope_analysis as oa
# import analysis.scripts.cavity.fit_oscilloscope_data as od


# data_dir = "K:/ns\qt\Diamond\Projects\Cavities\data/20160426/EOM_lw_LT/"
# file_name = "NNNNNNNLOWT007.csv"
# filename = os.path.join(data_dir,file_name)
# EOM_freq = 6 
# reload(oa)
# reload(od)
# data = oa.load_data(filename)
# od.get_linewidth(data,EOM_freq)