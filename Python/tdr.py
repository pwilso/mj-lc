import h5py
import numpy as np
import math
import sys
import string

class tdr_plot:
    def __init__(self, x, y):
        assert len(x) == len(y)
        self.x = x
        self.y = y
        self.max_x = np.max(x)
        self.min_x = np.min(x)
        self.max_y = np.max(y)
        self.min_y = np.min(y)

class tdr:
    def __init__(self, filename):
        f = h5py.File(filename, 'r')
        self.collection = f['collection']
        self.s0 = f['collection']['geometry_0']['state_0']
        num_plots = f['collection'].attrs['number of plots']
        print( '*Number of Plots', num_plots)
        self.plots = dict()
        i=0
        # Read datasets from TDR file and make a dict of tdr_plots
        # this only applies to 1-D (x-y) data files.
        while i < num_plots:
            ds = self.s0['dataset_' + str(i)]
            name = string.split(ds.attrs['name'], '_')[0]
            x = np.array(ds['values'])
            x_name = ds.attrs['name']
            i+= 1
            ds = self.s0['dataset_' + str(i)]
            y = np.array(ds['values'])
            y_name = ds.attrs['name']
            self.plots[name] = tdr_plot(x, y)
            i+= 1

