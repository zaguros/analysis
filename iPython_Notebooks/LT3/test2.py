import numpy as np
import h5py as h5py

f = h5py.File("testfile.hdf5")

arr = np.ones((5,2))
f["my dataset"] = arr
dset = f["my dataset"]
dset
dset.dtype
dset.shape

print dset


