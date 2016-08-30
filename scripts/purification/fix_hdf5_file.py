from analysis.lib.tools import toolbox
import h5py
import os

def fix_hdf5_file(contains = '',**kw):

	folder = toolbox.latest_data(contains,**kw)
	h5filepath = toolbox.measurement_filename(folder)
	filename = os.path.split(folder)[1]

	h5mode=kw.pop('hdf5_mode', 'r+')
	
	f = h5py.File(h5filepath,h5mode)

	groups = 0
	for k in f.keys():
		if type(f[k])==h5py.Group:
			groups += 1
			group_name = k
			print k

	if groups != 1:
		print 'Ambiguous hdf5 group! Cannot automatically rename!'
		return

	f[filename] = f[k]
	#del f[k]
