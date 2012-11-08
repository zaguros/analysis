#dirty hack below.
#small convenience module to quickly load or convert npz files

import numpy as np
import os, sys, time
import __main__ as mn

def load_npz(filename):
    """Loads all keys in a given npz-file into the global namespace"""
    data=np.load(filename)
    for key in data.keys():
        try:
            objtype=str(mn.__dict__[key].__class__)
            print 'key', key, 'already exists in namespace, of type',objtype
        except KeyError:
            mn.__dict__[key]=data[key]

def convert_npz(filename,keys=[],type='txt'):
    """Saves all keys in a given npz-file into separate files of given types type"""
    name,ext=os.path.splitext(filename)
    data=np.load(filename)
    if type=='txt':
        if len(keys)>0:
            for key in keys:
                np.savetxt(name+'_'+key,data[key]+'.txt')
        else:
            for key in data.keys():
                inp=raw_input('save ' + key +' ?')
                if inp=='y':
                    np.savetxt(name+'_'+key,data[key]+'.txt')
    else:
        print 'type', type, 'not supported'
    