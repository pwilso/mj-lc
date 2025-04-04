#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  4 11:38:17 2023

@author: pwilson3
"""

import operator
import tables
import numpy as np
import seaborn as sns
from cycler import cycler
import matrix_plotter as mplt
import matplotlib.pyplot as plt
import sys
import time

ops = {'==': operator.eq, 'and': operator.and_, '!=': operator.ne}

def grab_h5_data(h5_filename, table, column, conditions):
    """
    

    Parameters
    ----------
    h5_filename : name of the h5 file (ex: 'foo.h5')
    table : the specfic table to be accessed: (ex: foo.root.group.table)
    column : the column of the table you want to grab (ex: 'zA')
    conditions : list of conditions you want on the data (ex: [('pol', '==', val1]),
                                                               ('angle', '==', val2))

    Returns
    -------
    None.

    """
    
    tb = table
    # create conditions
    # for x in tb.iterrows():
    #     if all([ops[cond[1]](x[cond[0]], cond[2]) for cond in conditions]):
    #         print(x[column])
        # print('row: ', x[column], 'conditions: ', [ops[cond[1]](x[cond[0]], cond[2]) for cond in conditions])
    vals = [x[column] for x in tb.iterrows() if all([ops[cond[1]](x[cond[0]], cond[2]) for cond in conditions])]
    # h5_data.close()
    return vals

sns.set_theme(context='poster', style="whitegrid",
              rc={'lines.markeredgecolor': 'k', 'axes.edgecolor': 'k',
                  'xtick.direction': 'in', 'ytick.direction': 'in',
                  'xtick.top': True, 'ytick.right': True,
                  'lines.markersize': 10})
sns.set_theme(style = 'whitegrid',font_scale = 1.5, rc = {"font.size": 20, "lines.linewidth": 2})
colours = sns.color_palette('Set1', desat = 1)
colours2 = sns.color_palette('Set1', desat = 0.5)
cycle_colours = (cycler(color=sns.color_palette('Set1', desat = 1)))


#%% 

h5_file = 'slice.h5'

h5_data = tables.open_file(h5_file, mode = 'r')
table = h5_data.root.TMM_info.TMM_data

angles = grab_h5_data(h5_file, table, 'angle', [['pol', '==', b'TE']])
angles = list(dict.fromkeys(angles))

energies = grab_h5_data(h5_file, table, 'energy', [['pol', '==', b'TE']])
energies = list(dict.fromkeys(energies))

angle = angles[::50]
for i in range(len(angle)):
    kz = grab_h5_data(h5_file, table, 'wavevec', [['pol', '==', b'TE'],
                                                  ['energy', '==', 1.6],
                                                  ['angle', '==', angle[i]]])
    yB = grab_h5_data(h5_file, table, 'LayerB', [['pol', '==', b'TE'],
                                                 ['energy', '==', 1.6],
                                                 ['angle', '==', angle[i]]])
    plt.plot(yB, abs(np.real(kz)), label = np.degrees(angle[i]))

plt.legend()
h5_data.close()

#%%
h5_file = 'slice.h5'

h5_data = tables.open_file(h5_file, mode = 'r')
table = h5_data.root.Pos_info.Pos_data

angles = grab_h5_data(h5_file, table, 'angle', [])
angles = list(dict.fromkeys(angles))

energies = grab_h5_data(h5_file, table, 'energy', [])
energies = list(dict.fromkeys(energies))


angle = [1.57]# angles
for i in range(len(energies)):
    yA = grab_h5_data(h5_file, table, 'zA', [['energy', '==', energies[i]],
                                             ['angle', '==', angle[0]]])
    yB = grab_h5_data(h5_file, table, 'zB', [['energy', '==', energies[i]],
                                             ['angle', '==', angle[0]]])
    a_coupl= grab_h5_data(h5_file, table, 'a_TE', [['energy', '==', energies[i]],
                                                   ['angle', '==', angle[0]]])
    
    plt.plot(yB, a_coupl, marker = 'o')#, label = round(np.degrees(angle[i]), 2))
# plt.legend()

#%%
h5_file = '/path/to/LC_h5data.h5'

h5_data = tables.open_file(h5_file, mode = 'r')
table = h5_data.root.Pos_info.Pos_data

angles = grab_h5_data(h5_file, table, 'angle', [])
angles = list(dict.fromkeys(angles))

energies = grab_h5_data(h5_file, table, 'energy', [])
energies = list(dict.fromkeys(energies))

#%%
yA = []
yB = []
a_TE_list = []
# angles = angles[-3:]
for i in [0, 5, 10, 19]:
    start = time.time()
    print(i)
    A = grab_h5_data(h5_file, table, 'zA', [['angle', '==', angles[i]],
                                                       ['energy', '==', energies[18]]])
    yA.append(A)
    B= grab_h5_data(h5_file, table, 'zB', [['angle', '==', angles[i]],
                                                       ['energy', '==', energies[18]]])
    yB.append(B)
    a_TE = grab_h5_data(h5_file, table, 'a_TE', [['angle', '==', angles[i]],
                                                       ['energy', '==', energies[18]]])
    a_TE_list.append(a_TE)
    end = time.time()
    print('time:', end-start)

# for i in range(len(energies)):
#     yA = grab_h5_data(h5_file, table, 'zA', [['a_coupling', '!=', 0.0]])
#     yB = grab_h5_data(h5_file, table, 'zB', [['a_coupling', '!=', 0.0]])
#     a_coupl= grab_h5_data(h5_file, table, 'a_TE', [['a_coupling', '!=', 0.0]])
    
    # plt.plot(yB, a_coupl, marker = 'o')#, label = round(np.degrees(angle[i]), 2))
# plt.legend()
#%%
p1, ax1, ax2 = mplt.setup_simple_plot(1)

# a_coupl_nonnegative = [x if x > 1e-10 else 1e-10 for x in a_coupl]
mplt.simple_matrix_plot(p1, ax1, ax2, yA[2], yB[2], a_TE_list[2],
                        tick_points = np.linspace(-1e-12, 2, 4))

#%% 2023-08-17

h5_file = '/path/to/LC_h5data.h5'

h5_data = tables.open_file(h5_file, mode = 'r')
table = h5_data.root.Pos_info.Pos_data

angles = grab_h5_data(h5_file, table, 'angle', [])
angles = list(dict.fromkeys(angles))

energies = grab_h5_data(h5_file, table, 'energy', [])
energies = list(dict.fromkeys(energies))


# angle = [1.57]# angles
# 
# yA = grab_h5_data(h5_file, table, 'zA', [['a_coupling', '!=', 0.0]])
# yB = grab_h5_data(h5_file, table, 'zB', [['a_coupling', '!=', 0.0]])
# a_coupl= grab_h5_data(h5_file, table, 'a_coupling', [['a_coupling', '!=', 0.0]])

#%%
yA = []
yB = []
z_list = []
# angles = angles[-3:]
for i in [5]:
    start = time.time()
    print(i)
    A = grab_h5_data(h5_file, table, 'zA', [['angle', '==', angles[i]],
                                                       ['energy', '==', energies[18]]])
    yA.append(A)
    B= grab_h5_data(h5_file, table, 'zB', [['angle', '==', angles[i]],
                                                       ['energy', '==', energies[18]]])
    yB.append(B)
    # a_TE = grab_h5_data(h5_file, table, 'a_TE', [['angle', '==', angles[i]],
    #                                                    ['energy', '==', energies[18]]])
    # a_TE_list.append(a_TE)
    
    z = grab_h5_data(h5_file, table, 'a_TE', [['angle', '==', angles[i]],
                                                       ['energy', '==', energies[18]]])
    z_list.append(z)
    end = time.time()
    print('time:', end-start)

# for i in range(len(energies)):
#     yA = grab_h5_data(h5_file, table, 'zA', [['a_coupling', '!=', 0.0]])
#     yB = grab_h5_data(h5_file, table, 'zB', [['a_coupling', '!=', 0.0]])
#     a_coupl= grab_h5_data(h5_file, table, 'a_TE', [['a_coupling', '!=', 0.0]])
    
    # plt.plot(yB, a_coupl, marker = 'o')#, label = round(np.degrees(angle[i]), 2))
# plt.legend()
#%%
plt.rcParams['figure.figsize'] = 10, 10
p1, ax1, ax2 = mplt.setup_simple_plot(1)
i = 0
# a_coupl_nonnegative = [x if x > 1e-10 else 1e-10 for x in a_coupl]
mplt.simple_matrix_plot(p1, ax1, ax2, yB[i], yA[i], np.real(np.sqrt((z_list[i]*np.conj(z_list[i])))),
                        tick_points = np.linspace(-1e-12, 2, 6),
                        barlabel = r'|E| (a_TE)')
p1.savefig("/output.png")
             # )
# mplt.create_labels(ax, labels, materials, thicknesses)
#%% 2023-08-17

h5_file = '/path/to/LC_h5data.h5'

h5_data = tables.open_file(h5_file, mode = 'r')
table = h5_data.root.TMM_info.TMM_data

angles = grab_h5_data(h5_file, table, 'angle', [])
angles = list(dict.fromkeys(angles))

energies = grab_h5_data(h5_file, table, 'energy', [])
energies = list(dict.fromkeys(energies))

zAs = grab_h5_data(h5_file, table, 'zA', [])
zAs = list(dict.fromkeys(zAs))
#%%

yA = 6.27581262588501
yB = []
z_list = []
z_list2 = []
r_list = []
t_list = []
# angles = angles[-3:]
for i in range(len(angles)):
    start = time.time()
    print(i)
    B= grab_h5_data(h5_file, table, 'LayerB', [['angle', '==', angles[i]],
                                                       ['energy', '==', energies[18]],
                                                       ['pol', '==', b'TE'],
                                                       ['zA', '==', yA]])
    yB.append(B)
    z = grab_h5_data(h5_file, table, 'complex_angle', [['angle', '==', angles[i]],
                                                       ['energy', '==', energies[18]],
                                                       ['pol', '==', b'TE'],
                                                       ['zA', '==', yA]])
    z_list.append(z)
    
    z2 = grab_h5_data(h5_file, table, 'wavevec', [['angle', '==', angles[i]],
                                                       ['energy', '==', energies[18]],
                                                       ['pol', '==', b'TE'],
                                                       ['zA', '==', yA]])
    z_list2.append(z2)
    r = grab_h5_data(h5_file, table, 'r', [['angle', '==', angles[i]],
                                                       ['energy', '==', energies[18]],
                                                       ['pol', '==', b'TE'],
                                                       ['zA', '==', yA]])
    r_list.append(r)
    t = grab_h5_data(h5_file, table, 't', [['angle', '==', angles[i]],
                                                       ['energy', '==', energies[18]],
                                                       ['pol', '==', b'TE'],
                                                       ['zA', '==', yA]])
    t_list.append(t)
    end = time.time()
    print('time:', end-start)
    
#%%
plt.rcParams['figure.figsize'] = 10, 10
fig1, ax1 = plt.subplots(1,1)
# i = 1
for i in range(len(angles)):
# a_coupl_nonnegative = [x if x > 1e-10 else 1e-10 for x in a_coupl]
    ax1.plot(B, np.real(np.sqrt((r_list[i]*np.conj(r_list[i])))), marker = 'o', label = np.degrees(angles[i]))

# ax1.set_ylim(0, 1)
fig1.legend(framealpha=1).set_draggable(True)
