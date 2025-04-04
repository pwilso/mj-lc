# -*- coding: utf-8 -*-
"""
Created on Mon Nov 14 13:30:20 2022

@author: pwilson3
"""
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.tri as tri
from matplotlib.lines import Line2D 
import matplotlib.patches as patches
import seaborn as sns
import itertools

#------------------------------------------------------------------------------
def setup_simple_plot(fig_num):
    p1 = plt.figure(fig_num,figsize=(8,9))

    gs1 = matplotlib.gridspec.GridSpec(15,12)
    gs1.update(left=0.1, right=0.87, top=0.90, hspace=0.23, bottom=0.04)
    ax2 = plt.subplot( gs1[0:12,0:11])
    ax2.set_aspect('equal', 'box')
    ax3 = plt.subplot( gs1[13,3:8])

    ax2.set_ylabel(r'Emission location $\mathit{z_A}$ ($\mu$m)', labelpad=1)
    ax2.set_xlabel(r'Re-absorption location $\mathit{z_B}$ ($\mu$m)', labelpad=1)
    ax2.invert_yaxis()
    return p1, ax2, ax3

#------------------------------------------------------------------------------
def simple_matrix_plot(p1, ax2, ax3, x, y, z, tick_points = [], logmin = -4, logmax = 1,
                       barlabel = r'Re-absorbtion density ($\mu m^{-1}}$)'):
    # generate a series of color levels for the contour plot.  Add a zero level to capture everything less than the range of interest.
    zlim_min = logmin # base 10
    zlim_max = logmax #base 10
    
    if len(tick_points) == 0:
        # cmap = plt.cm.get_cmap(name='RdYlBu_r', lut=None)
        cmap = plt.cm.get_cmap(name='YlGnBu_r', lut=None)
        tick_points = [10**x for x in np.linspace(zlim_min, zlim_max, abs(zlim_min) + zlim_max + 1)]
        levels =  np.logspace(zlim_min, zlim_max, num=50)
        levels = np.concatenate((np.array([0]), levels))
        norm = matplotlib.colors.LogNorm(vmin=10**(zlim_min), vmax=10**(zlim_max), 
                                         clip=True)
        tri1 = tri.Triangulation(x, y)
        #ax2.tricontour(tri1, zs, levels=levels, cmap=cmap) 
        cs = ax2.tricontourf(tri1, np.maximum(np.array(z), np.zeros(np.array(z).shape)),
                             cmap=cmap, levels=levels, linewidth=0.25, norm=norm)
        cbar = p1.colorbar(cs, cax=ax3, orientation='horizontal', ticks=tick_points, 
                            format=matplotlib.ticker.LogFormatterMathtext(10, labelOnlyBase=True))
    else:
        cmap = plt.cm.get_cmap(name='turbo', lut=None)
        tick_points = tick_points
        levels = np.linspace(tick_points[0], tick_points[-1], 50)
        tri1 = tri.Triangulation(x, y)
        #ax2.tricontour(tri1, zs, levels=levels, cmap=cmap) 
        cs = ax2.tricontourf(tri1, np.array(z),
                             cmap=cmap, levels=levels, linewidth=0.25)#, extend = 'neither')
        ax2.set_xlim(min(x + y), max(x + y))
        ax2.set_ylim(max(x + y), min(x + y))
        cbar = p1.colorbar(cs, cax=ax3, orientation='horizontal', ticks=tick_points)
    #insert(cs)
    #ax2.triplot(tri1, linewidth=0.2)
    cbar = p1.colorbar(matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap), cax=ax3, orientation='horizontal',ticks=tick_points,format=matplotlib.ticker.LogFormatterMathtext(10, labelOnlyBase=True))
    cbar.solids.set_edgecolor("face")
    cbar.solids.set_linewidth(0.25)
    cbar.solids.set_antialiased(True)
    
    ax3.set_xlabel(r'Re-absorbtion density ($\mu m^{-1}}$)')

#------------------------------------------------------------------------------
def create_labels(ax, labels, materials, thicknesses, xbar = True, ybar = True,
                  xbar_anchor = -0.15, ybar_anchor = 0.1, bar_thickness = 0.1):
    colours = sns.color_palette('Set1')
    colours_list = {'GaAs': colours[0]} # add additional colours and materials here

    ybots = [sum(thicknesses[:i + 1]) for i in range(len(thicknesses))]
    ytops = [0] + ybots[:-1]
    ymids = [(ytops[i] + ybots[i])/2 for i in range(len(ytops))]

    for i in range(len(thicknesses)):
        if xbar == True:
            ax.add_patch(patches.Rectangle((ytops[i], xbar_anchor),
                                           thicknesses[i], bar_thickness,
                                           linewidth=0.25,
                                           facecolor=colours_list[materials[i]],
                                           edgecolor = 'k', alpha=0.6,clip_on=False))
            ax.annotate(labels[i], (ymids[i], xbar_anchor),
                        annotation_clip = False, rotation = 45).draggable()
        if ybar == True:
            ax.add_patch(patches.Rectangle((ybots[-1] + ybar_anchor, ytops[i]),
                                           bar_thickness, thicknesses[i],
                                           linewidth=0.25,
                                           facecolor=colours_list[materials[i]],
                                           edgecolor = 'k', alpha=0.6,clip_on=False))
            ax.annotate(labels[i], (ybots[-1] + 0.3, ymids[i]),
                        annotation_clip = False).draggable()

#------------------------------------------------------------------------------
def create_legend(p1, materials):
    colours = sns.color_palette('Set1')
    colours_list = {'GaAs': colours[0]} # add additional colours and materials here
    
    trimmed_list = []
    for m in materials:
        if m not in trimmed_list:
            trimmed_list.append(m)                

    legend_lines = [Line2D([0], [0], color = colours_list[m], lw = 4, alpha = 0.6) for m in trimmed_list]
    p1.legend(legend_lines, trimmed_list, loc = 'lower right', frameon = False).set_draggable(True)
    

#------------------------------------------------------------------------------

def plot_LCMatrix(folder,filename,num_segs,wtot):
	#input location of matrix file
	import matplotlib
	import matplotlib.tri as tri
	import pandas as pd
	import seaborn as sns

	import matrix_handler


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
	p1, ax1, ax2 = setup_simple_plot(1)
	simple_matrix_plot(p1, ax1, ax2,
						m1.xs + m1.xs2, 
						m1.ys + m1.ys2, 
						np.abs( m1.zs + m1.zs2)*wtot) # multiply by wtot to get reabsorption density (instead of coupling coefficient)
	create_labels(ax1, custom_labels, [x.lstrip('\t').lstrip(' ') for x in m1.materials],m1.thicknesses)
	#mplt.create_legend(p1, [x.lstrip('\t').lstrip(' ') for x in m1.materials])
    
    
#------------------------------------------------------------------------------
def plot_nkData(e):

    matplotlib.rcParams.update({'figure.autolayout': True})

    fig, (ax2_1,ax2_2) = plt.subplots(2)


    material_list = []
    colours = sns.color_palette('Set1')
    colours_list = {'GaAs': colours[0]} # add additional colours and materials here
    emission_range = np.linspace(e.wl_emission_min, e.wl_emission_max, 50)


    palette1 = itertools.cycle(colours)
    palette2 = itertools.cycle(colours)
    for li in range(len(e.layers)):
        if (e.layers[li].material).lstrip(' ') not in material_list:
          material_list.append((e.layers[li].material).lstrip(' '))
          
          n_vals = [e.layers[li].index_n(x) for x in emission_range]
          k_vals = [e.layers[li].index_k(x) for x in emission_range]
          ax2_1.plot(emission_range, n_vals, label = (e.layers[li].material).lstrip(' '),
                     color = colours_list[(e.layers[li].material).lstrip(' ')], ls = '-', marker = 'o', markeredgecolor = 'k')
          ax2_2.plot(emission_range, k_vals, label = (e.layers[li].material).lstrip(' '),
                     color = colours_list[(e.layers[li].material).lstrip(' ')], ls = '-', marker = 'o', markeredgecolor = 'k')



    ax2_1.set_xlabel('Wavelength $\mathit{\lambda}$ ($\mu m$)')
    ax2_1.set_ylabel('Re(n)')
    ax2_2.set_xlabel('Wavelength $\mathit{\lambda}$ ($\mu m$)')
    ax2_2.set_ylabel('Im(n)')


    ax2_1.legend()
    ax2_2.legend()
    ax2_1.legend()
