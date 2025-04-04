#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 27 16:38:15 2023

@author: pwilson3
"""

import sys, os
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), 
                             'tmm'))

import tdr
import h5py
import pyximport; pyximport.install()
import epi_cmd_LC
from scipy.integrate import trapz
import numpy as np
import csv
from datetime import datetime
import math
import tmm.tmm_core_mw as tmm
# cimport  tmm.tmm_core_mw as tmm
import h5py

from multiprocessing import Pool, Lock, Queue

import matplotlib.pyplot as plt
import seaborn as sns
from cycler import cycler
import time
from tables import *

def layer_setup(results_path, epi_node, sde_node, MatPar_node, angles):
    """
    Grabs the layer structure, material parameters, and mesh points from the
    Sentaurus setup files. Same as the setup from the LCMatrix_pyt.py file.

    Parameters
    ----------
    results_path: path to where the Sentaurus files can be found
    epi_node: points to the Sentaurus epi node files
    sde_node: points to the Sentaurus sdevice node files
    MatPar_node: points to the Sentaurus MatPar node files

    Returns
    -------
    e : object containing information about the layer setup

    """
    e = epi_cmd_LC.epifile(results_path + '/results/nodes/'+ str(epi_node) + '/pp'+ str(epi_node) + '_epi.cmd')
    e.bottom_material(back_mat, material_files_list, back_mat_thickness)
    # Add MatPar data to each layer
    # including n,k vs. wavelength and B_rad
    e.read_parfiles(results_path + '/results/nodes/'+ str(MatPar_node) + '/n' + str(MatPar_node) + '_mpr.par', 
                    angle_list = angles)
    # creates e.layer.yvalues, a dict indexed by yposition of vertices, by reading the mesh file.
    # e.get_layer_vertices('results/nodes/' + str(sde_node) + '/n' + str(sde_node) + '_msh.tdr')
    e.get_layer_vertices(results_path + '/results/nodes/'+ str(sde_node) + '/n' + str(sde_node) + '_msh.tdr')
    return e

class LCData_TMM(IsDescription):
    LayerA = Int32Col()
    zA = Float32Col()
    angle = Float32Col()
    energy = Float32Col()
    LayerB = Int32Col()
    pol = StringCol(2)
    wavevec = ComplexCol(16)
    complex_angle = ComplexCol(16)
    r = ComplexCol(16)
    t = ComplexCol(16)
    Ef_l1 = ComplexCol(16)
    Ef_l2 = ComplexCol(16)
    Eb_l1 = ComplexCol(16)
    Eb_l2 = ComplexCol(16)
    Ef_r1 = ComplexCol(16)
    Ef_r2 = ComplexCol(16)
    Eb_r1 = ComplexCol(16)
    Eb_r2 = ComplexCol(16)
    
class LCData_Pos(IsDescription):
    zA = Float64Col()
    zB = Float64Col()
    angle = Float64Col()
    energy = Float64Col()
    Ef1_TE = ComplexCol(16)
    Ef2_TE = ComplexCol(16)
    Eb1_TE = ComplexCol(16)
    Eb2_TE = ComplexCol(16)
    Ef1_TM = ComplexCol(16)
    Ef2_TM = ComplexCol(16)
    Eb1_TM = ComplexCol(16)
    Eb2_TM = ComplexCol(16)
    a_TE = Float64Col()
    a_TM = Float64Col()
    a_int_E = Float64Col()
    a_int_angle = Float64Col()
    a_coupling = Float64Col()

def TMM_calc(e, li, yA0, back_mat):
    n_list = np.zeros((len(e.layers)+1+len(back_mat)), dtype=complex) # +1 semi-inf front layer, + number of back layers (most bottom one will be semi-inf)
    d_list = np.zeros((len(e.layers)+1+len(back_mat)), dtype=np.float64)
    n_list[0] = 1.0+0j
    d_list[0] = tmm.inf
    for i in range(len(e.layers)):
        d_list[i+1] = e.layers[i].thickness
    for i in range(len(e.backlayers)):
        d_list[len(e.layers) +1 + i] = e.backlayers[i].thickness # changes to backmt -pwils
    d_list[-1] = tmm.inf #whatever last backlayer thickness is will be overwritten
    n_list[-1] = 1.0+0j # this will be adjusted later in new_n_list

    th_list = e.theta_list
    sin_th_list = np.sin(th_list)
    
    
    if e.layers[li].material_type == 'Semiconductor':
      yi_list = e.layers[li].yvalues
      thickness = e.layers[li].thickness
      if yA0 < thickness: # make sure the arbitrary yA0 used is less than the thickness of the layer
          E_list = e.layers[li].E_list
    
          # P_list will contain a list of tuples (yA, yB, P) for emission from li at yA and
          # absorption in each layer lk at yB
          # Pint_list contains only the P values for integration.
          # P_list = []
          # Pint_list = []
          # for lk in range(len(e.layers)):
                # P_list.append([])
                # Pint_list.append([])
    
          # P_arr = np.zeros((len(e.layers), len(yi_list)) )
          # Ps = np.zeros(len(e.layers))
          new_n_list = n_list.copy()
    
          # set up storage to solutions to TMM calculations.
          data_list_TE = []
          data_list_TM = []
          for a in E_list:
              data_list_TE.append([tmm.Tmm_data_internalsource()]*(len(th_list))) #!!! Note: this was originally len(th_list) - 1, and the 90 degree case was done separately -pwils
              data_list_TM.append([tmm.Tmm_data_internalsource()]*(len(th_list))) #!!! Note: this was originally len(th_list) - 1, and the 90 degree case was done separately -pwils
     
          #yA0 = yA - e.layers[li].ytop
     
          # Build TMM solutions at yA
          layer_flags = np.zeros(len(e.layers), dtype=int)
          # iterate over energy values near the emitting layer band gap
          for E_i in range(len(E_list)):
               E = E_list[E_i]
               wl = e.layers[li].wl_list[E_i] # in um
               for i in range(len(e.layers)):
                   new_n_list[i+1] = e.layers[i].nk(wl)
              # set semi-infinite layer at the  of the stack same as last finite layer.
        # set the semi-infinite layer at the bottom of the stack to be the same as the substrate pwils
               for i in range(len(e.backlayers)):
                   if e.backlayers[i].material != 'Air':
                       new_n_list[len(e.layers) + 1 + i] = e.backlayers[i].nk(wl) #e.bottom_mat.nk(wl)
                   else:
                       new_n_list[len(e.layers) + 1 + i] = 1. + 0j # nk for Air/Vacuum
              #print('##############new_n_list:', new_n_list, '@ wl = ', wl)
              # iterate over polar emission angles
               for th_i in range(len(th_list)): #!!! Note: this was originally len(th_list) - 1, and the 90 degree case was done separately -pwils
                   th = th_list[th_i]
                   # print(th)
                   # generate TMM solutions for TE and TM polarizations
                   data_list_TE[E_i][th_i] = tmm.coh_tmm_internalsource(tmm.TE, new_n_list.copy(), d_list, li, yA0, th, wl)
                   data_list_TM[E_i][th_i] = tmm.coh_tmm_internalsource(tmm.TM, new_n_list.copy(), d_list, li, yA0, th, wl)
                   # flag layers with non-negligible power entering from the source.
    else:
        print('Value of yA must be less than the layer thickness')
    return data_list_TE, data_list_TM

sns.set_theme(context='poster', style="whitegrid",
              rc={'lines.markeredgecolor': 'k', 'axes.edgecolor': 'k',
                  'xtick.direction': 'in', 'ytick.direction': 'in',
                  'xtick.top': True, 'ytick.right': True,
                  'lines.markersize': 10})
sns.set_theme(style = 'whitegrid',font_scale = 1.5, rc = {"font.size": 20, "lines.linewidth": 2})
colours = sns.color_palette('Set1', desat = 1)
colours2 = sns.color_palette('Set1', desat = 0.5)
cycle_colours = (cycler(color=sns.color_palette('Set1', desat = 1)))

if __name__ == '__main__':
    angles = np.linspace(0., np.pi/2, 500)
    # use_gauss_quadrature = True
    
    substrate_t = 0.5
    # sets up lists of refractive index, thickness, angles, and energy
    if substrate_t == 0.0: # when substrate_t is set to 0 the model assumes a flat gold backreflector
        back_mat = ["Gold"]
        back_mat_thickness = [1] # last layer thickness will be set to inf
    else:
        back_mat = ["InP"]
        back_mat_thickness = [1] # last layer thickness will be set to inf
    material_files = {'Gold': 'Gold.par', 'InP': 'InP.tcl', 'Air': None}
    material_files_list = [material_files[x] for x in back_mat]
    
    filepath = '/path/to/Sentaurus/project/'
    e = layer_setup(filepath, 123, 456, 789, angles)
    original_theta = e.theta_list
       
#%% PLOTS THE ABSORPTION INTEGRATED OVER EMISSION ENERGIES WITGH RESPECT TO 
### ANGLE FOR A GIVEN EMISSION LOCATION yA AND A SUBSET OF ABSORPTION LOCATIONS
### yB
    sns.set_style("white")
    fig1, ax1 = plt.subplots(1, 1)
    fig1.set_size_inches(11, 9)
    
    li = 3
    yA = e.layers[li].yvalues[3]
    yA0 = yA - e.layers[li].ytop
    
    yB_vals = e.layers[3].yvalues[::7]
    lks = [3]*33
    marker = ['o', 's', 'd', 'P', '*', 'd', '<', '+']
    no_marker = ['']*33
    for j in range(len(yB_vals)):
        lk = lks[j]
        yB = yB_vals[j]
        
        TE_data, TM_data = TMM_calc(e, li, yA0, ["InP"])
        
        P0_list_th, P1_list_th, P_th = epi_cmd_LC.calc_a2(li,lk, yA, yB, TE_data, TM_data, e, False, True)
        P0_list_E, P1_list_E, P_E = epi_cmd_LC.calc_a2(li,lk, yA, yB, TE_data, TM_data, e, True, False)
        
        
        ax1.plot(np.degrees(e.theta_list), P0_list_E + 0.6*j, marker = marker[j], c = colours[j],
                 label = 'yA = {A:.2f}, yB = {B:.2f}'.format(A = yA, B = yB))
        ax1.axhline(y = 0.6*j, c = 'lightgrey', ls = ':')
        
    ax1.set_xlabel('Emission Angle ($\degree$)')
    ax1.set_ylabel('$\int a_{AB}$($\\theta$,E)$\hat{S}(E)\sin(\\theta)dE$')
    ax1.legend(frameon = True, loc = 'best').set_draggable(True)
    
#%% SAME AS THE ABOVE PLOT EXCEPT THIS IS BEING USED TO TEST THE GUASSIAN 
### QUADRATURE INTEGRATION. THE TE_data AND TM_data NEED TO BE CALCULATED BASED
### ON THE GAUSS-LEGENDRE POSITIONS IN THE TRANSFORMED COORDINATE SPACE: 
### u = 2*cos(theta)-1

    sns.set_style("white")
    fig1, ax1 = plt.subplots(1, 1)
    fig1.set_size_inches(11, 9)
    
    li = 3
    yA = 1.8959
    yA0 = yA - e.layers[li].ytop
    
    yB_vals = [yA]
    # yB_vals.add(yA) # want the emission point to also be included
    lks = [3]*33
    marker = ['o', 's', 'd', 'P', '*', 'd', '<', '+']
    no_marker = ['']*33
    
    for use_gauss_quadrature in [False]:
        if use_gauss_quadrature == True:
            n = len(angles)
            [u, w] = np.polynomial.legendre.leggauss(n)
            theta_of_u = np.arccos(0.5*(u + 1)) # let u = 2*cos(theta) - 1
            e.theta_list = theta_of_u # when done this can be an option in the epi_cmd_LC file to use theta_of_u instead of linearly spaced theta
    
        else:
            e.theta_list = original_theta
        TE_data, TM_data = TMM_calc(e, li, yA0, ["InP"])
        
        start = time.time()
        # put all the TE and TM data into a table so it can be parsed and plotted more easily
        # LCData = open_file('LC_Data_2J_test.h5', mode = 'w', title = 'Test')
        # TMM_info = LCData.create_group('/', 'TMM_info')
        # TMM_table = LCData.create_table(TMM_info, 'TMM_data', LCData_TMM)
        # TMM_row = TMM_table.row
        
        # for i in range(np.shape(TE_data)[0]): #maybe this would be better served inside the tmm_calc function
        #     for j in range(np.shape(TE_data)[1]):
        #         for layerB in range(len(e.layers) + 1):
        #             TMM_row['LayerA'] = li
        #             TMM_row['LayerB'] = layerB
        #             TMM_row['zA'] = yA
        #             TMM_row['angle'] = TE_data[i][j].phi_i
        #             TMM_row['energy'] = TE_data[i][j].lam_vac
        #             TMM_row['pol'] = 'TE'
        #             TMM_row['wavevec'] = TE_data[i][j].kz_list[layerB]
        #             TMM_row['complex_angle'] = TE_data[i][j].th_list[layerB]
        #             TMM_row['r'] = TE_data[i][j].r_list[layerB][layerB + 1]
        #             TMM_row['t'] = TE_data[i][j].t_list[layerB][layerB + 1]
                    
        #             TMM_row.append()
                
        #             TMM_row['LayerA'] = li
        #             TMM_row['LayerB'] = layerB
        #             TMM_row['zA'] = yA
        #             TMM_row['angle'] = TM_data[i][j].phi_i
        #             TMM_row['energy'] = TM_data[i][j].lam_vac
        #             TMM_row['pol'] = 'TM'
        #             TMM_row['wavevec'] = TM_data[i][j].kz_list[layerB]
        #             TMM_row['complex_angle'] = TM_data[i][j].th_list[layerB]
        #             TMM_row['r'] = TM_data[i][j].r_list[layerB][layerB + 1]
        #             TMM_row['t'] = TM_data[i][j].t_list[layerB][layerB + 1]
                    
        #             TMM_row.append()
        
        # TMM_table.flush()
        # LCData.close()
        end = time.time()
        time_total = end - start
        print('H5 Time: ', time_total)
        
        # LayerA = Int64Col()
        # angle = Float64Col()
        # energy = Float64Col()
        # LayerB = Int64Col()
        # pol = StringCol(2)
        # wavevec = ComplexCol(16)
        # complex_angle = ComplexCol(16)
        # r = ComplexCol(16)
        # t = ComplexCol(16)
        
        for j in range(len(yB_vals)):
            lk = lks[j]
            yB = yB_vals[j]
            
            start = time.time()
            
            if use_gauss_quadrature == False:
                P0_list_E, P1_list_th, P_th = epi_cmd_LC.calc_a_ordered_integ(li,lk, yA, yB, 
                                                        TE_data, TM_data, e, 
                                                        integrate_theta_first = False,
                                                        use_gauss_quad = False)
                P0_list_E = P0_list_E/np.sin(e.theta_list)
                P0_list_E
                print(P_th)
                # P0_list_E, P1_list_E, P_E = epi_cmd_LC.calc_a_ordered_integ(li,lk, yA, yB, 
                #                                         TE_data, TM_data, e, 
                #                                         integrate_theta_first = True,
                #                                         use_gauss_quad = False)
                # print(P_E)
                ax1.plot(np.degrees(e.theta_list), P0_list_E, marker = marker[j], c = colours[j],
                           label = 'yA = {A:.2f}, yB = {B:.2f}'.format(A = yA, B = yB), ls = '-')
            else:
                P0_list_E, P1_list_th, P_th = epi_cmd_LC.calc_a_ordered_integ(li,lk, yA, yB, 
                                                        TE_data, TM_data, e, 
                                                        integrate_theta_first = False,
                                                        use_gauss_quad = True)
                P0_list_E
                print(P_th)
                # P0_list_E, P1_list_E, P_E = epi_cmd_LC.calc_a_ordered_integ(li,lk, yA, yB, 
                #                                         TE_data, TM_data, e, 
                #                                         integrate_theta_first = True,
                #                                         use_gauss_quad = True)
                # print(P_E)
                ax1.plot(np.degrees(e.theta_list), P0_list_E, marker = marker[j], c = colours2[j],
                           label = '__no_label__', ls = '--')
            end = time.time()
            time_total = end - start
            print('Time: ', time_total)
            
            ax1.axhline(y = 0, c = 'lightgrey', ls = ':')
            # ax1.set_ylim(0, 5)
    t = np.linspace(0, 89.9, 1000)
    ax1.plot(t, (1/np.cos(np.radians(t))), c = 'k', ls = '-', label = '1/cos($\\theta$)')
    # ax1.plot(t, ((np.cos(np.radians(t)))**4 - (np.cos(np.radians(t)))**2 + np.cos(np.radians(t)))/((np.cos(np.radians(t)))**3 + (np.cos(np.radians(t)))**2), c = 'k', ls = '-')
    ax1.set_xlabel('Emission Angle ($\degree$)')
    ax1.set_ylabel('$\int a_{AB}$($\\theta$,E)$\hat{S}(E)dE$')
    ax1.set_ylim(-0.01, 70)
    ax1.legend(frameon = True, loc = 'best').set_draggable(True)
    
    current_t = datetime.now().strftime("%Y-%m-%d-%H%M%S")
    fig1.savefig('/path/to/' + current_t + 'output.png')
#%% SAME AS THE ABOVE PLOT EXCEPT THIS IS BEING USED TO TEST WRITING TO
### PYTABLES. [2023-08-03]

    sns.set_style("white")
    fig1, ax1 = plt.subplots(1, 1)
    fig1.set_size_inches(11, 9)
    
    li = 3
    yA = e.layers[li].yvalues[17]
    yA0 = yA - e.layers[li].ytop
    
    # yB_vals = []
    # lks = []
    # for i in range(len(e.layers)):
    #     for pos in e.layers[i].yvalues:
    #         yB_vals.append(pos)
    #         lks.append(i)
    yB_vals = e.layers[3].yvalues[::8]
    # yB_vals.add(yA) # want the emission point to also be included
    lks = [3]*33
    marker = ['o', 's', 'd', 'P', '*', 'd', '<', '+']
    no_marker = ['']*33
    
    for use_gauss_quadrature in [False]:
        if use_gauss_quadrature == True:
            n = len(angles)
            [u, w] = np.polynomial.legendre.leggauss(n)
            theta_of_u = np.arccos(0.5*(u + 1)) # let u = 2*cos(theta) - 1
            e.theta_list = theta_of_u # when done this can be an option in the epi_cmd_LC file to use theta_of_u instead of linearly spaced theta
    
        else:
            e.theta_list = original_theta
        TE_data, TM_data = TMM_calc(e, li, yA0, ["InP"])
        
        start = time.time()
        # put all the TE and TM data into a table so it can be parsed and plotted more easily
        # create the table
        LCData = open_file('LC_Data.h5', mode = 'w', title = 'Test')
        TMM_info = LCData.create_group('/', 'TMM_info')
        TMM_table = LCData.create_table(TMM_info, 'TMM_data', LCData_TMM)
        
        Pos_info = LCData.create_group('/', 'Pos_info')
        Pos_table = LCData.create_table(Pos_info, 'Pos_data', LCData_Pos)
        #pos_row = Pos_table.row
        #LCData.close() # close the file until we need it again
        
        TMM_row = TMM_table.row
        for i in range(np.shape(TE_data)[0]): #maybe this would be better served inside the tmm_calc function
            for j in range(np.shape(TE_data)[1]):
                for layerB in range(len(e.layers) + 1):
                    TMM_row['LayerA'] = li
                    TMM_row['LayerB'] = layerB
                    TMM_row['zA'] = yA
                    TMM_row['angle'] = TE_data[i][j].phi_i
                    TMM_row['energy'] = TE_data[i][j].lam_vac
                    TMM_row['pol'] = 'TE'
                    TMM_row['wavevec'] = TE_data[i][j].kz_list[layerB]
                    TMM_row['complex_angle'] = TE_data[i][j].th_list[layerB]
                    TMM_row['r'] = TE_data[i][j].r_list[layerB][layerB + 1]
                    TMM_row['t'] = TE_data[i][j].t_list[layerB][layerB + 1]
                    if layerB <= li:
                        TMM_row['Ef_l1'] = TE_data[i][j].vw_list_l1[layerB + 1][0] #forward emission, left of source
                        TMM_row['Eb_l1'] = TE_data[i][j].vw_list_l1[layerB + 1][1]
                        TMM_row['Ef_l2'] = TE_data[i][j].vw_list_l2[layerB + 1][0]
                        TMM_row['Eb_l2'] = TE_data[i][j].vw_list_l2[layerB + 1][1]
                    if layerB >= li:
                        TMM_row['Ef_r1'] = TE_data[i][j].vw_list_r1[layerB - li][0] #forward emission, left of source
                        TMM_row['Eb_r1'] = TE_data[i][j].vw_list_r1[layerB - li][1]
                        TMM_row['Ef_r2'] = TE_data[i][j].vw_list_r2[layerB - li][0]
                        TMM_row['Eb_r2'] = TE_data[i][j].vw_list_r2[layerB - li][1]
                    
                    TMM_row.append()
                
                    TMM_row['LayerA'] = li
                    TMM_row['LayerB'] = layerB
                    TMM_row['zA'] = yA
                    TMM_row['angle'] = TM_data[i][j].phi_i
                    TMM_row['energy'] = TM_data[i][j].lam_vac
                    TMM_row['pol'] = 'TM'
                    TMM_row['wavevec'] = TM_data[i][j].kz_list[layerB]
                    TMM_row['complex_angle'] = TM_data[i][j].th_list[layerB]
                    TMM_row['r'] = TM_data[i][j].r_list[layerB][layerB + 1]
                    TMM_row['t'] = TM_data[i][j].t_list[layerB][layerB + 1]
                    if layerB <= li:
                        TMM_row['Ef_l1'] = TM_data[i][j].vw_list_l1[layerB + 1][0] #forward emission, left of source
                        TMM_row['Eb_l1'] = TM_data[i][j].vw_list_l1[layerB + 1][1]
                        TMM_row['Ef_l2'] = TM_data[i][j].vw_list_l2[layerB + 1][0]
                        TMM_row['Eb_l2'] = TM_data[i][j].vw_list_l2[layerB + 1][1]
                    if layerB >= li:
                        TMM_row['Ef_r1'] = TM_data[i][j].vw_list_r1[layerB - li][0] #forward emission, left of source
                        TMM_row['Eb_r1'] = TM_data[i][j].vw_list_r1[layerB - li][1]
                        TMM_row['Ef_r2'] = TM_data[i][j].vw_list_r2[layerB - li][0]
                        TMM_row['Eb_r2'] = TM_data[i][j].vw_list_r2[layerB - li][1]
                    # TMM_row['vw0_1'] = TM_data[i][j].vw0_1[layerB][layerB + 1]
                    # TMM_row['vw0_2'] = TM_data[i][j].vw0_2[layerB][layerB + 1]
                    # TMM_row['vw1_1'] = TM_data[i][j].vw1_1[layerB][layerB + 1]
                    # TMM_row['vw1_2'] = TM_data[i][j].vw1_2[layerB][layerB + 1]
                    
                    TMM_row.append()
        
        TMM_table.flush()
        LCData.close() # close the file until we need it again
        #end = time.time()
        #time_total = end - start
        #print('H5 Time: ', time_total)
        
        # LayerA = Int64Col()
        # angle = Float64Col()
        # energy = Float64Col()
        # LayerB = Int64Col()
        # pol = StringCol(2)
        # wavevec = ComplexCol(16)
        # complex_angle = ComplexCol(16)
        # r = ComplexCol(16)
        # t = ComplexCol(16)
        
        for j in range(len(yB_vals)):
            lk = lks[j]
            yB = yB_vals[j]
            
            start = time.time()

            if use_gauss_quadrature == False:
                P0_list_E, P1_list_th, P_th = epi_cmd_LC.calc_a_ordered_integ_writeH5_test(li,lk, yA, yB, 
                                                        TE_data, TM_data, e, 
                                                        integrate_theta_first = True,
                                                        use_gauss_quad = False,
                                                        h5_filename = '20230817_LC_Data_2J_test_more_info.h5')
                P0_list_E = P0_list_E
                P0_list_E
                print(P_th)
                # P0_list_E, P1_list_E, P_E = epi_cmd_LC.calc_a_ordered_integ_writeH5_test(li,lk, yA, yB, 
                #                                         TE_data, TM_data, e, 
                #                                         integrate_theta_first = True,
                #                                         use_gauss_quad = False)
                # print(P_E)
                # ax1.plot(np.degrees(e.theta_list), P0_list_E, marker = marker[j], c = colours[j],
                           # label = 'yA = {A:.2f}, yB = {B:.2f}'.format(A = yA, B = yB), ls = '-')
            else:
                P0_list_E, P1_list_th, P_th = epi_cmd_LC.calc_a_ordered_integ_writeH5_test(li,lk, yA, yB, 
                                                        TE_data, TM_data, e, 
                                                        integrate_theta_first = True,
                                                        use_gauss_quad = True,
                                                        h5_filename = '20230817_LC_Data_2J_test_more_info.h5')
                P0_list_E
                print(P_th)
                # P0_list_E, P1_list_E, P_E = epi_cmd_LC.calc_a_ordered_integ_writeH5_test(li,lk, yA, yB, 
                #                                         TE_data, TM_data, e, 
                #                                         integrate_theta_first = True,
                #                                         use_gauss_quad = True)
                # print(P_E)
                # ax1.plot(np.degrees(e.theta_list), P0_list_E, marker = marker[j], c = colours2[j],
                           # label = '__no_label__', ls = '--')
            end = time.time()
            time_total = end - start
            print('Time: ', time_total)
            
            ax1.axhline(y = 0, c = 'lightgrey', ls = ':')
            ax1.set_ylim(0, 5)

    ax1.set_xlabel('Emission Angle ($\degree$)')
    ax1.set_ylabel('$\int a_{AB}$($\\theta$,E)$\hat{S}(E)dE$')
    # ax1.set_ylim(-0.01, 11)
    ax1.legend(frameon = True, loc = 'best').set_draggable(True)
    
    current_t = datetime.now().strftime("%Y-%m-%d-%H%M%S")

#%%
    fig1, ax1 = plt.subplots(1, 1)
    fig1.set_size_inches(11, 9)
    
    li = 3
    yA = 1.8959
    yA0 = yA - e.layers[li].ytop
    TE_data, TM_data = TMM_calc(e, li, yA0, ["InP"])
    
    yB_vals = e.layers[10].yvalues[5]
    lks = [10]*5
    marker = ['o', 's', 'd', 'P', '*']
    no_marker = ['']*5
    absor_list = np.ones(len(e.theta_list))
    for j in range(len(yB_vals)):
        for angle_index in range(len(e.theta_list)-1):
            lk = lks[j]
            yB = yB_vals[j]
            yB0 = yB - e.layers[lks[j]].ytop
            
            TE_data, TM_data = TMM_calc(e, li, yA0, ["InP"])
            
            absor_list[angle_index] = epi_cmd_LC.calc_absor(li, lk, yA, yB, 
                                                         TE_data, TM_data, e, 
                                                         2, angle_index)
    
        ax1.plot(np.degrees(e.theta_list), absor_list, c = colours[j], marker = no_marker[j],
             label = 'yA = {A:.2f}, yB = {B:.2f}'.format(A = yA, B = yB))
    
    ax1.set_xlabel('Emission Angle ($\degree$)')
    ax1.set_ylabel('$a_{AB}$($\\theta$,E)$')
    ax1.legend(frameon = True, loc = 'best').set_draggable(True)

#%% EXTRACT DATA FROM H5 PYTABLE FILE AND PLOT

LCData = open_file('LC_Data.h5', mode = 'r')

tb = LCData.root.TMM_info.TMM_data

all_angles = []
all_wl = []
for x in tb.iterrows():
    if x['angle'] not in all_angles:
        all_angles.append(x['angle'])
    if x['energy'] not in all_wl:
        all_wl.append(x['energy'])

# grab TE data at specific angle and energies
zeta_TE_89 = [x['complex_angle'] for x in tb.iterrows() if x['pol'] == b'TE'
               and x['angle'] == all_angles[0] and x['energy'] == all_wl[1]]
zeta_TE_84 = [x['complex_angle'] for x in tb.iterrows() if x['pol'] == b'TE'
               and x['angle'] == all_angles[12] and x['energy'] == all_wl[1]]
zeta_TE_69 = [x['complex_angle'] for x in tb.iterrows() if x['pol'] == b'TE'
               and x['angle'] == all_angles[24] and x['energy'] == all_wl[1]]
zeta_TE_48 = [x['complex_angle'] for x in tb.iterrows() if x['pol'] == b'TE'
               and x['angle'] == all_angles[36] and x['energy'] == all_wl[1]]
zeta_TE_2 = [x['complex_angle'] for x in tb.iterrows() if x['pol'] == b'TE'
               and x['angle'] == all_angles[-1] and x['energy'] == all_wl[1]]

zeta_TM_89 = [x['complex_angle'] for x in tb.iterrows() if x['pol'] == b'TM'
               and x['angle'] == all_angles[0] and x['energy'] == all_wl[1]]
zeta_TM_84 = [x['complex_angle'] for x in tb.iterrows() if x['pol'] == b'TM'
               and x['angle'] == all_angles[12] and x['energy'] == all_wl[1]]
zeta_TM_69 = [x['complex_angle'] for x in tb.iterrows() if x['pol'] == b'TM'
               and x['angle'] == all_angles[24] and x['energy'] == all_wl[1]]
zeta_TM_48 = [x['complex_angle'] for x in tb.iterrows() if x['pol'] == b'TM'
               and x['angle'] == all_angles[36] and x['energy'] == all_wl[1]]
zeta_TM_2 = [x['complex_angle'] for x in tb.iterrows() if x['pol'] == b'TM'
               and x['angle'] == all_angles[-1] and x['energy'] == all_wl[1]]

kvec_TE_89 = [x['wavevec'] for x in tb.iterrows() if x['pol'] == b'TE'
               and x['angle'] == all_angles[0] and x['energy'] == all_wl[1]]
kvec_TE_84 = [x['wavevec'] for x in tb.iterrows() if x['pol'] == b'TE'
               and x['angle'] == all_angles[12] and x['energy'] == all_wl[1]]
kvec_TE_69 = [x['wavevec'] for x in tb.iterrows() if x['pol'] == b'TE'
               and x['angle'] == all_angles[24] and x['energy'] == all_wl[1]]
kvec_TE_48 = [x['wavevec'] for x in tb.iterrows() if x['pol'] == b'TE'
               and x['angle'] == all_angles[36] and x['energy'] == all_wl[1]]
kvec_TE_2 = [x['wavevec'] for x in tb.iterrows() if x['pol'] == b'TE'
               and x['angle'] == all_angles[-1] and x['energy'] == all_wl[1]]

kvec_TM_89 = [x['wavevec'] for x in tb.iterrows() if x['pol'] == b'TM'
               and x['angle'] == all_angles[0] and x['energy'] == all_wl[1]]
kvec_TM_84 = [x['wavevec'] for x in tb.iterrows() if x['pol'] == b'TM'
               and x['angle'] == all_angles[12] and x['energy'] == all_wl[1]]
kvec_TM_69 = [x['wavevec'] for x in tb.iterrows() if x['pol'] == b'TM'
               and x['angle'] == all_angles[24] and x['energy'] == all_wl[1]]
kvec_TM_48 = [x['wavevec'] for x in tb.iterrows() if x['pol'] == b'TM'
               and x['angle'] == all_angles[36] and x['energy'] == all_wl[1]]
kvec_TM_2 = [x['wavevec'] for x in tb.iterrows() if x['pol'] == b'TM'
               and x['angle'] == all_angles[-1] and x['energy'] == all_wl[1]]

# grab ydata
layer_index = [x['LayerB'] for x in tb.iterrows() if x['pol'] == b'TE'
               and x['angle'] == all_angles[0] and x['energy'] == all_wl[1]]

fig1, ax1 = plt.subplots(1, 1)
fig1.set_size_inches(11, 9)

ax1.plot(layer_index, (np.real(kvec_TE_89)))
ax1.plot(layer_index, (np.real(kvec_TE_84)))
ax1.plot(layer_index, (np.real(kvec_TE_69)))
ax1.plot(layer_index, (np.real(kvec_TE_48)))
ax1.plot(layer_index, (np.real(kvec_TE_2)))
ax1.set_xlabel('Layer Index')
ax1.set_ylabel('Complex Angle (real part)')
# ax1.plot(layer_index, np.degrees(np.imag(zeta_TE_89)))

# ax1.plot(layer_index, np.real(kvec_TE_89))
# ax1.plot(layer_index, np.imag(kvec_TE_89))
