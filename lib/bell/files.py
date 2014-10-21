import h5py
import os

import Settings

### data finding and identification

def get_msmt_name(fp):
    """
    This assumes that there is only one group, whose name is the msmt name.
    """
    _root, fn = os.path.split(fp)
    f = h5py.File(fp, 'r')
    for k in f.keys():
        if f.get(k, getclass=True) == h5py._hl.group.Group and k in fn:
            f.close()
            return k
    
    raise Exception('Cannot find the name of the measurement.')
    
def get_msmt_fp(folder, ext='hdf5'):
    dirname = os.path.split(folder)[1]
    fn = dirname+'.'+ext
    return os.path.join(folder, fn)

def get_msmt_header(fp):
    _root, fn = os.path.split(fp)
    root, folder = os.path.split(_root)
    daystamp = os.path.split(root)[1]
    return daystamp + '/' + folder


### basic functionality for analysis data
def delete_analysis_data(fp, name, subgroup=None):
    try:
        f = h5py.File(fp, 'r+')
    except:
        print "Cannot open file", fp
        raise

    try:
        agrp = f.require_group(Settings.ANALYSISGRP + ('/' + subgroup if subgroup!=None else ''))
        if name in agrp.keys():
            del agrp[name]
        f.flush()
        f.close()        
    except:
        f.close()
        raise




def has_data(fp, name, subgroup=None):
    """
    Checks if data with a specific name is in the keys of the main folder
    give filepath and name
    """
    try:
        f = h5py.File(fp, 'r')
    except:
        return False

    
    if name in f.keys():
        f.close()
        return True
    else:
        f.close()
        return False
        
def has_analysis_data(fp, name, subgroup=None):
    try:
        f = h5py.File(fp, 'r')
    except:
        print "Cannot open file", fp
        raise
    
    if Settings.ANALYSISGRP in f.keys():
        agrp = f.require_group(Settings.ANALYSISGRP + ('/' + subgroup if subgroup!=None else ''))
        if name in agrp.keys():
            f.close()
            return True
        else:
            f.close()
            return False
    else:
        f.close()
        return False
    
def get_analysis_data(fp, name, subgroup=None):
    try:
        f = h5py.File(fp, 'r')
    except:
        print "Cannot open file", fp
        raise
        
    agrp = f.require_group(Settings.ANALYSISGRP + ('/' + subgroup if subgroup!=None else ''))

    if name not in agrp.keys():
        return None
    
    dat = agrp[name].value
    attrs = {}
    for (an, av) in agrp[name].attrs.items():
        attrs[an] = av
        
    f.close()
        
    return dat, attrs
