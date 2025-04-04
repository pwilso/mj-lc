# -*- coding: utf-8 -*-
"""
Created on Tue Nov 15 16:20:54 2022

@author: pwilson3
"""

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import csv
import numpy as np
import sys
import pandas as pd
import seaborn as sns

import matrix_plotter as mplt
import matrix_handler

# location of matrix file
folder = '/path/to/the/folder/'
filename = 'file_name_of_csv'


# parameters for scaling the data
wtot = 1
num_segs = 2

print('Reading in Matrix...')
# import the matrix with the matrix handler
m1 = matrix_handler.matrix_handler(folder + filename + '.csv')
print('Finished reading in Matrix.')

# fill this in to get custom label names
custom_labels = []
for x in m1.layers:
    if x[-2:] == 'em':
        custom_labels.append('sc' + str(num_segs))
        num_segs -= 1
    elif x == 'cap':
        custom_labels.append(x)
    elif x =='substrate':
        custom_labels.append('buffer')
    elif x == 'ARC1':
        custom_labels.append('ARC')
    else:
        custom_labels.append(' ')

sns.set_theme(font_scale = 1.5, style = 'white')

#plot data
print("Plotting data")
p1, ax1, ax2 = mplt.setup_simple_plot(1)
mplt.simple_matrix_plot(p1, ax1, ax2,
                        m1.xs + m1.xs2, 
                        m1.ys + m1.ys2, 
                        m1.zs*wtot + m1.zs2*wtot) # multiply by wtot to get reabsorption density (instead of coupling coefficient)
mplt.create_labels(ax1, custom_labels, 
                    [x.lstrip('\t').lstrip(' ') for x in m1.materials],
                    m1.thicknesses)
mplt.create_legend(p1, [x.lstrip('\t').lstrip(' ') for x in m1.materials])

# compare upper and lower halves of the matrix
# create a list of values which compares upper half of matrix to lower half
# sorted_lower = sorted(zip(m1.xs,m1.ys, m1.zs))
# sorted_upper = sorted(zip(m1.ys2, m1.xs2, m1.zs2))
# compare_z = np.array(sorted_lower) - np.array(sorted_upper)

# p2, ax3, ax4 = mplt.setup_simple_plot(2)
# mplt.simple_matrix_plot(p2, ax3, ax4,
#                         m1.xs2,
#                         m1. ys2,
#                         compare_z[:,2],
#                         tick_points = np.linspace(-10, 10, 5))

# mplt.create_labels(ax3, m1.layers, 
#                     [x.lstrip('\t').lstrip(' ') for x in m1.materials],
#                     m1.thicknesses)
# mplt.create_legend(p2, [x.lstrip('\t').lstrip(' ') for x in m1.materials])

 # compare two matrixes
# m2 = matrix_handler.matrix_handler(folder2 + filename2 + '.csv')


plt.show()
