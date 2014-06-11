# some convenience tools
#

import os
import time
import logging
import numpy as np

try:
    import qt
    datadir = qt.config['datadir']
except:
    # Added a line for Mac compatibility. Does require data to be saved in correct folder (as below)
    if os.name == 'posix':
        datadir = r'/Users/'+os.getlogin()+r'/Documents/teamdiamond/data'
    else:
        datadir = r'd:\measuring\data'

def nearest_idx(array, value):
    '''
    find the index of the value closest to the specified value.
    '''
    return np.abs(array-value).argmin()

def nearest_value(array, value):
    '''
    find the value in the array that is closest to the specified value.
    '''
    return array[nearest_idx(array,value)]

def verify_timestamp(timestamp):
    if len(timestamp) == 6:
        daystamp = time.strftime('%Y%m%d')
        tstamp = timestamp
    elif len(timestamp) == 14:
        daystamp = timestamp[:8]
        tstamp = timestamp[8:]
    elif len(timestamp) == 15: #### In case day and timestamp separted by _
        daystamp = timestamp[:8]
        tstamp = timestamp[9:]
    else:
        raise Exception("Cannot interpret timestamp '%s'" % timestamp)

    return daystamp, tstamp

def is_older(ts0, ts1):
    '''
    returns True if timestamp ts0 is an earlier data than timestamp ts1,
    False otherwise.
    '''
    if ts0 == None or ts1 == None:
        return True
    else:

        dstamp0, tstamp0 = verify_timestamp(ts0)
        dstamp1, tstamp1 = verify_timestamp(ts1)

        return (dstamp0+tstamp0) < (dstamp1+tstamp1)

def latest_data(contains='', older_than=None, newer_than=None,return_timestamp = False,raise_exc = True):
    '''
    finds the latest taken data with <contains> in its name.
    returns the full path of the data directory.

    if older_than is not None, than the latest data that fits and that
    is older than the date given by the timestamp older_than is returned.
    if newer_than is not None, than the latest data that fits and that
    is newer than the date given by the timestamp newer_than is returned

    If no fitting data is found, an exception is raised. Except when you specifically ask not to to
    this in: raise_exc = False, then a 'False' is returned.
    '''

    daydirs = os.listdir(datadir)
    if len(daydirs) == 0:
        logging.warning('No data found in datadir')
        return None

    daydirs.sort()

    measdirs = []
    i = len(daydirs)-1
    while len(measdirs) == 0 and i >= 0:
        daydir = daydirs[i]
        all_measdirs = [d for d in os.listdir(os.path.join(datadir, daydir))]
        all_measdirs.sort()

        measdirs = []

        for d in all_measdirs:

            # this routine verifies that any output directory is a 'valid' directory
            # (i.e, obeys the regular naming convention)
            _timestamp = daydir + d[:6]
            try:
                dstamp,tstamp = verify_timestamp(_timestamp)
            except:
                continue
            timestamp = dstamp+tstamp

            if contains in d:              
                
                if older_than != None:
                    if not is_older(timestamp, older_than):
                        continue
                if newer_than != None:
                    if not is_older(newer_than,timestamp):
                        continue
                measdirs.append(d)

        i -= 1

    if len(measdirs) == 0:
        if raise_exc == True:
            raise Exception('No fitting data found.')
        else:
            return False
    else:
        measdirs.sort()
        measdir = measdirs[-1]
        if return_timestamp == False:
            return os.path.join(datadir,daydir,measdir)
        else: return str(daydir)+str(measdir[:6]) , os.path.join(datadir,daydir,measdir)

# def newer_data(starttimestamp, endtimestamp=None, contains='',):
#     '''
#     finds all matching data in the given range. default end is now.
#     '''
#     if endtimestamp == None:
#         endtimestamp = time.strftime("%Y%m%d%H%M%S")

#     results = []
#     daydirs = os.listdir(datadir)
#     if len(daydirs) == 0:
#         logging.warning('No data found in datadir')
#         return None
#     daydirs.sort()

#     if len(endtimestamp) == 6:
#         endday = time.strftime('%Y%m%d')
#         endtime = endtimestamp
#     elif len(endtimestamp) == 14:
#         endday = endtimestamp[:8]
#         endtime = endtimestamp[8:]
#     else:
#         logging.warning("Cannot interpret timestamp '%s'" % endtimestamp)
#         return None

#     if len(starttimestamp) == 6:
#         startday = time.strftime('%Y%m%d')
#         starttime = starttimestamp
#     elif len(starttimestamp) == 14:
#         startday = starttimestamp[:8]
#         starttime = starttimestamp[8:]
#     else:
#         logging.warning("Cannot interpret timestamp '%s'" % starttimestamp)
#         return None

#     # TODO continue here

#     return None

def data_from_time(timestamp):
    '''
    returns the full path of the data specified by its timestamp in the
    form YYYYmmddHHMMSS.
    '''
    daydirs = os.listdir(datadir)
    if len(daydirs) == 0:
        raise Exception('No data in the data directory specified')

    daydirs.sort()
    daystamp, tstamp = verify_timestamp(timestamp)

    if not os.path.isdir(os.path.join(datadir,daystamp)):
        logging.warning("Requested day '%s' not found" % daystamp)
        return None

    measdirs = [ d for d in os.listdir(os.path.join(datadir,daystamp)) \
            if d[:6] == tstamp ]

    if len(measdirs) == 0:
        logging.warning("Requested data '%s'/'%s' not found" \
                % (daystamp, tstamp))
        return None
    elif len(measdirs) == 1:
        return os.path.join(datadir,daystamp,measdirs[0])
    else:
        logging.warning('Timestamp is not unique: ', str(measdirs))
        return None

def measurement_filename(directory=os.getcwd(), ext='hdf5'):
    dirname = os.path.split(directory)[1]
    fn = dirname+'.'+ext

    if os.path.exists(os.path.join(directory,fn)):
        return os.path.join(directory,fn)
    else:
        logging.warning("Data path '%s' does not exist" % \
                os.path.join(directory,fn))
        return None

def get_plot_title_from_folder(folder):
    measurementstring = os.path.split(folder)[1]
    timestamp = os.path.split(os.path.split(folder)[0])[1] \
            + '/' + measurementstring[:6]
    measurementstring = measurementstring[7:]
    default_plot_title = timestamp+'\n'+measurementstring
    return default_plot_title




############### 2014-06-11, Hannes: here i am putting some of wolgang's teleporation tools: file management etc.

def get_all_msmt_filepaths(folder, suffix='hdf5', pattern=''):
    filepaths = []
    suffixlen = len(suffix)
    
    for root,dirs,files in os.walk(folder):
        for f in files:
            if len(f) > suffixlen and f[-suffixlen:] == suffix and pattern in f:
                filepaths.append(os.path.join(root, f))
    
    return filepaths

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
def delete_analysis_data(fp, name, analysisgrp = 'analysis', subgroup=None):
    try:
        f = h5py.File(fp, 'r+')
    except:
        print "Cannot open file", fp
        raise

    try:
        agrp = f.require_group(analysisgrp + ('/' + subgroup if subgroup!=None else ''))
        if name in agrp.keys():
            del agrp[name]
        f.flush()
        f.close()        
    except:
        f.close()
        raise

def clear_analysis_data(fp, analysisgrp = 'analysis'):
    try:
        f = h5py.File(fp, 'r+')
    except:
        print "Cannot open file", fp
        raise

    if analysisgrp in f.keys():
        try:
            del f[analysisgrp]
            f.flush()
            f.close()
        except:
            f.close()
            raise
    else:
        f.close()

def set_analysis_data(fp, name, data, analysisgrp = 'analysis', subgroup=None, **kw):
    try:
        f = h5py.File(fp, 'r+')
    except:
        print "Cannot open file", fp
        raise
    
    try:
        agrp = f.require_group(analysisgrp + ('/' + subgroup if subgroup!=None else ''))
        if name in agrp.keys():
            del agrp[name]
        agrp[name] = data
        f.flush()
        
        for k in kw:
            agrp[name].attrs[k] = kw[k]
        
        f.flush()
        f.close()        
    except:
        f.close()
        raise
        
def has_analysis_data(fp, name, analysisgrp = 'analysis', subgroup=None):
    try:
        f = h5py.File(fp, 'r+')
    except:
        print "Cannot open file", fp
        raise
    
    agrp = f.require_group(analysisgrp + ('/' + subgroup if subgroup!=None else ''))
    if name in agrp.keys():
        f.close()
        return True
    else:
        f.close()
        return False
    
def get_analysis_data(fp, name, analysisgrp = 'analysis', subgroup=None):
    try:
        f = h5py.File(fp, 'r+')
    except:
        print "Cannot open file", fp
        raise
        
    agrp = f.require_group(analysisgrp + ('/' + subgroup if subgroup!=None else ''))

    if name not in agrp.keys():
        return None
    
    dat = agrp[name].value
    attrs = {}
    for (an, av) in agrp[name].attrs.items():
        attrs[an] = av
        
    f.close()
        
    return dat, attrs








