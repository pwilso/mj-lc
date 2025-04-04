#!/usr/bin/python
#setdep @node|sde@

import os,sys
#os.environ['OPENBLAS_NUM_THREADS'] = '1'
import numpy as np
from tables import *

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
    a_TE = Float64Col()
    a_TM = Float64Col()
    Ef1_TE = ComplexCol(16)
    Eb1_TE = ComplexCol(16)
    Ef2_TE = ComplexCol(16)
    Eb2_TE = ComplexCol(16)
    Ef1_TM = ComplexCol(16)
    Eb1_TM = ComplexCol(16)
    Ef2_TM = ComplexCol(16)
    Eb2_TM = ComplexCol(16)
    power_out1_TE = Float64Col()
    power_out2_TE = Float64Col()
    power_out1_TM = Float64Col()
    power_out2_TM = Float64Col()
    a_int_E = Float64Col()
    a_int_angle = Float64Col()
    a_coupling = Float64Col()
    
class LCData_Pos_Scaled(IsDescription):
    zA = Float64Col()
    zB = Float64Col()
    angle = Float64Col()
    energy = Float64Col()
    a_TE_scaled = Float64Col()
    a_TM_scaled = Float64Col()
    a_int_E_scaled = Float64Col()
    a_int_angle_scaled = Float64Col()
    a_coupling_scaled = Float64Col()

### USER INPUT ################################################################ #!!!
run_cmdline = True # run on the command line (True) or through Sentaurus (False)
                    # running on the command line can be useful for debugging and profiling code
Pooling = True    # turn multiprocessing on or off, recommended on unless profiling or debugging =pwils
Threading = True  # turn multithreading on or off, recommended of unless profiling or debugging #TODO: does this actually work?

use_gaussian_quadrature = False # if True use guassian quadrature to perform the
                               # integral. If False use the trapezoidal rule
                               # with linearly spaced angles
test_integral_z = 'z_loop_calc_a' # this will probably get removed once I finish testing the feature
n_angles = 40 # number of angles to be used, check convergence
n_energies = 30 # number of energies to be used, check convergence
include90atend = True

# should be relative to the project directory
matrix_file = 'LC_Outputs/matrix_filename.csv'
h5_filename = 'LC_Outputs/h5_filename.h5'

# what kind of data do you want to be output?
# LC_out = True # output the LC coupling matrix
# abs_out = True # output the absorption matrix
# abs_th_E = True # output the absorption integrated in two different ways
# Efield_out = True
# TMM_out = True # output a bunch of TMM data
write_TMM = False # write the TMM data to the h5 file
write_pos = False # write the position resolved data to the h5 file

# these only need to be set if run_cmdline is true
epi_node = 123
sde_node = 456
MatPar_node = 789
wtot = 1.0 # 1.0 for 1D problems
substrate_t = 0.5 # 0.0 to include a flat Au back reflector
###############################################################################

if run_cmdline == False:
    if "@LC@" == "on": # this is always False from the command line, only relevant if Sentaurus preprocessor has run
        run_Sentaurus = True # program will run in Sentaurus if not running from command line and LC set to on
        matrix_file = 'LC_Outputs/' + "@filename_prefix@" + ".csv" #rename the file
    else:
        run_Sentaurus = False
        
if run_cmdline == True or run_Sentaurus == True:
    import os, sys

    ### path to custom python files
    ## use these paths for running outside of Sentuarus
    if run_cmdline == True:
        sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'Python'))
        sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'Python', 'tmm'))
    ## use these paths for running inside Sentaurus
    else:
        try:
            sys.path.append(os.path.join('@pwd@', '..', 'Python'))
            sys.path.append(os.path.join('@pwd@', '..', 'Python', 'tmm'))
        except SyntaxError:
            print('Ignoring Sentaurus Syntax when running from command line')

    import tdr
    import h5py
    import pyximport; pyximport.install()
    import epi_cmd_LC
    from scipy.integrate import trapz
    import csv
    import datetime
    import math
    import tmm.tmm_core_mw as tmm

    from multiprocessing import Pool, Lock, Queue


#%%
### BASIC SETUP ############################################################### #!!!
    ##suppress print statements
    #sys.stdout = open(os.devnull, 'w')
#    if Threading == False:
#        os.environ['OPENBLAS_NUM_THREADS'] = '1'
#    else:
#        pass

    if run_cmdline == True:
        # these variables are substituted by the sentaurus preprocessor when the script is run through Sentaurus Workbench.
        epi_node = epi_node
        sde_node = sde_node
        MatPar_node = MatPar_node
        wtot = wtot
        substrate_t = substrate_t

        #specify the name of the file and path to be written
        LC_filepath = './' + matrix_file

    else:
        try:
            epi_node = int("@node|epi@")
            sde_node = int("@node|sde@")
            MatPar_node = int("@node|MatPar@")
            wtot = float("@wtot@")
            substrate_t = float("@substrate_t@")

            #specify the name of the file and path to be written
            LC_filepath = '@pwd@/' + matrix_file
        except SyntaxError:
            print('Ignoring Sentaurus Syntax when running from command line')

    #back_mat = "back_mat"
    if substrate_t == 0.0: # when substrate_t is set to 0 the model assumes a flat gold backreflector
        back_mat = ["Gold"]
        back_mat_thickness = [1] # last layer thickness will be set to inf
    else:
        back_mat = ["InP"]
        back_mat_thickness = [1] # last layer thickness will be set to inf
    material_files = {'Gold': 'Gold.par', 'InP': 'InP.tcl', 'Air': None}
    material_files_list = [material_files[x] for x in back_mat]
    
    if use_gaussian_quadrature == False:
        angles = np.linspace(0, np.pi/2, n_angles) # the pi/2 case is treated separately and can be set to zero (ie. ignored) if include90==False
        #angles_above85 = np.array([np.pi/2]) #!!! this will be used to add additional angles to the main angle list, for testing
        #angles = np.concatenate((angles, angles_above85)) # + 89 deg, 89.9 deg, 89.99 deg
    else:
        [u, w] = np.polynomial.legendre.leggauss(n_angles)
        theta_of_u = np.arccos(0.5*(u + 1)) # let u = 2*cos(theta) - 1
        angles = theta_of_u

    # setup H5 file
    if write_TMM == True or write_pos == True:
        LCData = open_file(h5_filename, mode = 'w', title = '')
        if write_TMM == True:
            TMM_info = LCData.create_group('/', 'TMM_info')
            TMM_table = LCData.create_table(TMM_info, 'TMM_data', LCData_TMM)
        if write_pos == True:
            Pos_info = LCData.create_group('/', 'Pos_info')
            Pos_table = LCData.create_table(Pos_info, 'Pos_data', LCData_Pos)
            Pos_scaled_info = LCData.create_group('/', 'Pos_scaled_info')
            Pos_scaled_table = LCData.create_table(Pos_scaled_info, 'Pos_scaled_data', LCData_Pos_Scaled)
        LCData.close()
### MAIN FUNCTION#####################################################################

    # This function performs the coupling calculation for emissions from a single 
    # layer. Multiple threads are spawned, each processing one layer.
#%%
    def process_layer(li):
            # Calculate coupling for emission layer li
            print( "* Processing layer {0}. {1} y-values.".format(e.layers[li].name, len(e.layers[li].yvalues)))
            
            # Due to the threading, always flush after printing to ensure output appears in correct sequence.
            sys.stdout.flush()

            if e.layers[li].material_type == 'Semiconductor':
              yi_list = e.layers[li].yvalues
              E_list = e.layers[li].E_list

              # P_list will contain a list of tuples (yA, yB, P) for emission from li at yA and
              # absorption in each layer lk at yB
              # Pint_list contains only the P values for integration.
              P_list = []
              Pint_list = []
              for lk in range(len(e.layers)):
                    P_list.append([])
                    Pint_list.append([])

              P_arr = np.zeros((len(yi_list),len(e.layers)))
              Ps = np.zeros(len(e.layers))
              new_n_list = n_list.copy()

              # set up storage to solutions to TMM calculations.
              for z in range(len(yi_list)):
                yA = yi_list[z]
             #   P_for_given_yA = [] # want to collect yB and P for each yA in this list
                data_list_TE = []
                data_list_TM = []
                for a in E_list:
                    data_list_TE.append([tmm.Tmm_data_internalsource()]*(len(th_list))) #!!! Note: this was originally len(th_list) - 1, and the 90 degree case was done separately -pwils
                    data_list_TM.append([tmm.Tmm_data_internalsource()]*(len(th_list))) #!!! Note: this was originally len(th_list) - 1, and the 90 degree case was done separately -pwils

                yA0 = yA - e.layers[li].ytop

                # Build TMM solutions at yA
                layer_flags = np.zeros(len(e.layers), dtype=int)
                # iterate over energy values near the emitting layer band gap
                for E_i in range(len(E_list)):
                     E = E_list[E_i]
                     wl = e.layers[li].wl_list[E_i] # in um
                     for i in range(len(e.layers)):
                          new_n_list[i+1] = e.layers[i].nk(wl)
                     # set semi-infinite layer at the  of the stack same as last finite layer.
               # set the semi-infinite layer at the bottom of the stack to be the same as the substrate -pwils
                     for i in range(len(e.backlayers)):
                         if e.backlayers[i].material != 'Air':
                             new_n_list[len(e.layers) + 1 + i] = e.backlayers[i].nk(wl) 
                         else:
                            new_n_list[len(e.layers) + 1 + i] = 1. + 0j # nk for Air/Vacuum

                     # iterate over polar emission angles
                     for th_i in range(len(th_list)): #!!! Note: this was originally len(th_list) - 1, and the 90 degree case was done separately -pwils
                          th = th_list[th_i]

                          # generate TMM solutions for TE and TM polarizations
                          data_list_TE[E_i][th_i] = tmm.coh_tmm_internalsource(tmm.TE, new_n_list.copy(), d_list, li, yA0, th, wl)
                          data_list_TM[E_i][th_i] = tmm.coh_tmm_internalsource(tmm.TM, new_n_list.copy(), d_list, li, yA0, th, wl)
                          # flag layers with non-negligible power entering from the source.
                          for i in range(len(e.layers)):
                              layer_flags[i] = 1
                              # if tmm.check_layer(i, data_list_TE[E_i][th_i], data_list_TM[E_i][th_i]): # turn this back on in an attempt to correct fields blowing up -pwils
                              #    layer_flags[i] = 1
                              # else:
                              #    print('/n***Layer flag not set - calculation will be skipped') # if True, consider this layer; if False, fields entering the layer are small enough to be ignored and P is set to 0 in the block below this
############################################################################### #!!!
            # put all the TE and TM data into a table so it can be parsed and plotted more easily
                if write_TMM == True:
                    LCData = open_file(h5_filename, mode = 'r+', title = '')
                    TMM_row = LCData.root.TMM_info.TMM_data.row
                    
                    for i in range(np.shape(data_list_TE)[0]): #maybe this would be better served inside the tmm_calc function
                        for j in range(np.shape(data_list_TE)[1]):
                            for layerB in range(len(e.layers) + 1):
                                TMM_row['LayerA'] = li
                                TMM_row['LayerB'] = layerB
                                TMM_row['zA'] = yA
                                TMM_row['angle'] = data_list_TE[i][j].phi_i
                                TMM_row['energy'] = data_list_TE[i][j].lam_vac
                                TMM_row['pol'] = 'TE'
                                TMM_row['wavevec'] = data_list_TE[i][j].kz_list[layerB]
                                TMM_row['complex_angle'] = data_list_TE[i][j].th_list[layerB]
                                TMM_row['r'] = data_list_TE[i][j].r_list[layerB][layerB + 1]
                                TMM_row['t'] = data_list_TE[i][j].t_list[layerB][layerB + 1]
                                if layerB <= li:
                                    TMM_row['Ef_l1'] = data_list_TE[i][j].vw_list_l1[layerB + 1][0] #forward emission, left of source
                                    TMM_row['Eb_l1'] = data_list_TE[i][j].vw_list_l1[layerB + 1][1]
                                    TMM_row['Ef_l2'] = data_list_TE[i][j].vw_list_l2[layerB + 1][0]
                                    TMM_row['Eb_l2'] = data_list_TE[i][j].vw_list_l2[layerB + 1][1]
                                if layerB >= li:
                                    TMM_row['Ef_r1'] = data_list_TE[i][j].vw_list_r1[layerB - li][0] #forward emission, left of source
                                    TMM_row['Eb_r1'] = data_list_TE[i][j].vw_list_r1[layerB - li][1]
                                    TMM_row['Ef_r2'] = data_list_TE[i][j].vw_list_r2[layerB - li][0]
                                    TMM_row['Eb_r2'] = data_list_TE[i][j].vw_list_r2[layerB - li][1]
                                
                                TMM_row.append()
                            
                                TMM_row['LayerA'] = li
                                TMM_row['LayerB'] = layerB
                                TMM_row['zA'] = yA
                                TMM_row['angle'] = data_list_TM[i][j].phi_i
                                TMM_row['energy'] = data_list_TM[i][j].lam_vac
                                TMM_row['pol'] = 'TM'
                                TMM_row['wavevec'] = data_list_TM[i][j].kz_list[layerB]
                                TMM_row['complex_angle'] = data_list_TM[i][j].th_list[layerB]
                                TMM_row['r'] = data_list_TM[i][j].r_list[layerB][layerB + 1]
                                TMM_row['t'] = data_list_TM[i][j].t_list[layerB][layerB + 1]
                                if layerB <= li:
                                    TMM_row['Ef_l1'] = data_list_TM[i][j].vw_list_l1[layerB + 1][0] #forward emission, left of source
                                    TMM_row['Eb_l1'] = data_list_TM[i][j].vw_list_l1[layerB + 1][1]
                                    TMM_row['Ef_l2'] = data_list_TM[i][j].vw_list_l2[layerB + 1][0]
                                    TMM_row['Eb_l2'] = data_list_TM[i][j].vw_list_l2[layerB + 1][1]
                                if layerB >= li:
                                    TMM_row['Ef_r1'] = data_list_TM[i][j].vw_list_r1[layerB - li][0] #forward emission, left of source
                                    TMM_row['Eb_r1'] = data_list_TM[i][j].vw_list_r1[layerB - li][1]
                                    TMM_row['Ef_r2'] = data_list_TM[i][j].vw_list_r2[layerB - li][0]
                                    TMM_row['Eb_r2'] = data_list_TM[i][j].vw_list_r2[layerB - li][1]
                                
                                TMM_row.append()
                                
                    LCData.root.TMM_info.TMM_data.flush()
                    LCData.close()
############################################################################


                # now use TMM solutions to find absorption at specific locations yB.
                # iterate over absorption layer lk
                if test_integral_z == "SameLayer":
                    for lk in range(len(e.layers)): #[0] #!!!
                        if e.layers[lk].material_type == 'Semiconductor':
                            yk_list = e.layers[lk].yvalues
                            P_for_given_yA = [] # this will serve a similar purpose as Pint_list, clearing after each layer
                            Pint_list = np.zeros(len(yk_list))
                            for k in range(len(yk_list)):
                               yB = yk_list[k]
                               if layer_flags[lk] == 1:
                                   yB0 = yB - e.layers[lk].ytop
############################################################################ #!!!
                                   if write_pos == True:
                                       P = epi_cmd_LC.calc_a_ordered_integ_writeH5_test(li, lk, 
                                                        yA, yB,
                                                        data_list_TE, data_list_TM, e, 
                                                        integrate_theta_first = False, 
                                                        h5_filename = h5_filename, 
                                                        use_gauss_quad = use_gaussian_quadrature)[2]/(2.0*wtot)
                                   else: #!!! this needs to be updated but it is not often used so it can be done later
                                       P = epi_cmd_LC.calc_a(li, lk, yA, yB,
                                                             data_list_TE, data_list_TM, e, 
                                                             use_gauss_quad = use_gaussian_quadrature)/ (2.0*wtot)
############################################################################
                                   if P < -1e-10:
                                        print('Negative Coupling value, li=', li, ' lk=', lk)
                                        print( data_list_TE )
                                        print( data_list_TM)
                                        raise ValueError

                                   P_for_given_yA.append((yA, yB, P)) # refreshes for every new yA -pwils [2023-10-05]
                                   
                               else:
                                   P_for_given_yA.append((yA, yB, 0.)) # refreshes for every new yA -pwils [2023-10-05]
                                
                            if li == lk: # if we are in the same layer as the emission, values can get large, it is good to smear things out -pwils [2023-10-05]
                             # do the smoothing then add to P_list
                                dist_right = P_for_given_yA[1][1] - P_for_given_yA[0][1] 
                                P_right = P_for_given_yA[1][2]
                                P_at_yA = P_for_given_yA[0][2]
                                P_area = 0.25*dist_right*(1.5*P_at_yA + 0.5*P_right)
                                scaled_P = P_area/(dist_right)

                                new_tup = (yA, P_for_given_yA[0][1], scaled_P)
                                P_list[lk].append(new_tup)
                                Pint_list[k] = scaled_P
                                for index in np.arange(1, len(P_for_given_yA)-1, 1):   # deal with the 1st and last values in the list separately -pwils [2023-10-05]
                                    dist_right = P_for_given_yA[index+1][1] - P_for_given_yA[index][1]
                                    dist_left = P_for_given_yA[index][1] - P_for_given_yA[index-1][1]
                                    P_right = P_for_given_yA[index+1][2]
                                    P_at_yA = P_for_given_yA[index][2]
                                    P_left = P_for_given_yA[index-1][2]
                                    P_area = 0.25*dist_right*(1.5*P_at_yA + 0.5*P_right)
                                    P_area += 0.25*dist_left*(1.5*P_at_yA + 0.5*P_left)
                                    scaled_P = P_area/(dist_right + dist_left)

                                    new_tup = (yA, P_for_given_yA[index][1], scaled_P)
                                    P_list[lk].append(new_tup)
                                    Pint_list[k] = scaled_P
                                dist_left = P_for_given_yA[-1][1] - P_for_given_yA[-2][1]
                                P_at_yA = P_for_given_yA[-1][2]
                                P_left = P_for_given_yA[-2][2]
                                P_area = 0.25*dist_left*(1.5*P_at_yA + 0.5*P_left)
                                scaled_P = P_area/(dist_left)

                                new_tup = (yA, P_for_given_yA[-1][1], scaled_P)
                                P_list[lk].append(new_tup)
                                Pint_list[k] = scaled_P
                            else:
                            # this is probably okay for of diagonals (aka it's easier if I just leave it for now) -pwils [2023-10-05]
                                for orig_tup in P_for_given_yA:
                                    P_list[lk].append(orig_tup)
                                    Pint_list[k] = orig_tup[2]
                            
                            P = trapz(Pint_list, x=yk_list)
                            P_arr[lk,z] = P*wtot
# TEST ************************************************************************
                elif test_integral_z == 'z_loop_calc_a':
############################################################################### #!!!
                   if write_pos == True:
                       # P is now a list with length equal to the number of points in the mesh (duplicates have been removed at the boundary edges)
                       # P_list is a list of tuples with (yA, yB, P/(2*wtot)) and indexed by layer lk; this is initially a list of empty list but it will get fed back into the function and updated
                       # P_list_yA is the sum of P*dz/(2*wtot) and is indexed by lk
                       # in either case the value of P is actually scaled so that the integral of absor over yB is 1
                       P_list, P_list_int, front_loss, back_loss = epi_cmd_LC.calc_a_integrate_over_z(li, yA, data_list_TE,  data_list_TM,  e, P_list,
                                                                      integrate_theta_first = True,
                                                                      h5_filename = h5_filename, use_gauss_quad = False,
                                                                      include90=include90atend)[3:]
                   else: #!!! need to update this epi function with the H5 argument
                       P_list, P_list_int, front_loss, back_loss = epi_cmd_LC.calc_a_integrate_over_z(li, yA, data_list_TE,  data_list_TM,  e, P_list,
                                                                      integrate_theta_first = True,
                                                                      h5_filename = h5_filename, use_gauss_quad = False,
                                                                      include90=include90atend,
                                                                      write_H5 = False)[3:]
###############################################################################
                   P_arr[z] = P_list_int*wtot # P_arr is indexed by lk and z which is the index for the yA list in layer li
                
                print(' {0: <10} {1: 3d}/{2: 3d}       yA={3} Eff={4:.4f} FrontLoss={5:.4f} BackLoss={6:.4f} Err={7}'.format(e.layers[li].name, z+1, len(yi_list), yA, sum(P_arr[z]), front_loss, back_loss, 1- (sum(P_arr[z])+front_loss+back_loss)))
                sys.stdout.flush()

              # integrated coupling between 2 layers
              print('.')
              for lk in range(len(e.layers)):
                   # print(P_arr)
                   # print(P_arr[:, lk])
                   # print(yi_list)
                   Ps[lk] = trapz(P_arr[:, lk], x=yi_list)/e.layers[li].thickness ###FIXME: (actually maybe it does work) I don't think this works with the new scaled P values (sum * dz again?) 
              print('.')
              sys.stdout.flush()
              # after completing all calculations for emission from layer li, we acquire the file lock
              # and write to the output file.
              print('* Writing Output from {0}.  Total coupling efficiency: {1:.4f}'.format(e.layers[li].name, sum(Ps)))
              sys.stdout.flush()
              fileLock.acquire()
              for lk in range(len(e.layers)):
                     outwriter.writerow(['**', li, e.layers[li].name, lk, e.layers[lk].name, Ps[lk] ])
                     for i in range(len(P_list[lk])):
                          tu = P_list[lk][i]
                          outwriter.writerow([tu[0], tu[1], tu[2]])
                     outfile.flush()
              fileLock.release()
            return True
#%%
### BEGIN THE CALCULATION ############################################################################
    ### first check if an LCMatrix is already available

    if os.path.isfile(LC_filepath):
        print("*** LCMatrix file already exists, skipping calculation ***")

    else:
        starttime = datetime.datetime.now()
        print('***')
        print('*** LCMatrix - Compute optical coupling between layers in a multilayer stack ')
        print('*** started at ', starttime)
        print('***' )

        # Open an epi file and read the layer stack - materials, thicknesses etc.
        if run_cmdline == True:
            e = epi_cmd_LC.epifile('./results/nodes/'+ str(epi_node) + '/pp'+ str(epi_node) + '_epi.cmd')
            e.bottom_material(back_mat, material_files_list, back_mat_thickness)
            # Add MatPar data to each layer
            # including n,k vs. wavelength and B_rad
            e.read_parfiles('./results/nodes/'+ str(MatPar_node) + '/n' + str(MatPar_node) + '_mpr.par',
                            angle_list = angles, N_energies = n_energies)
            # creates e.layer.yvalues, a dict indexed by yposition of vertices, by reading the mesh file.
            e.get_layer_vertices('./results/nodes/'+ str(sde_node) + '/n' + str(sde_node) + '_msh.tdr')
        else:
            try:
                e = epi_cmd_LC.epifile('@pwd@' + '/results/nodes/' + str(epi_node) + '/pp'+ str(epi_node) + '_epi.cmd')
                # takes material, material file, thickness, and optional doping and molefraction; re-added this function -pwils
                e.bottom_material(back_mat, material_files_list, back_mat_thickness)
                # Add MatPar data to each layer
                # including n,k vs. wavelength and B_rad
                e.read_parfiles('@pwd@' +'/results/nodes/' + str(MatPar_node) + '/n' + str(MatPar_node) + '_mpr.par')
                # creates e.layer.yvalues, a dict indexed by yposition of vertices, by reading the mesh file.
                e.get_layer_vertices('@pwd@' + '/results/nodes/' + str(sde_node) + '/n' + str(sde_node) + '_msh.tdr')
            except SyntaxError:
                print('Ignoring Sentaurus Syntax when running from command line')
            
        sys.stdout.flush()

        print('')

        # Open output file for writing.
        #altered to write as text for the time being -pwils
        with  open(LC_filepath, 'w') as outfile:
                outwriter = csv.writer(outfile, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)

                # write the header
                outwriter.writerow(['#'])
                outwriter.writerow(['# Output from LCMatrix ' + datetime.datetime.now().strftime("%B %d, %Y at %I:%M%p")] )
                outwriter.writerow(['# device width = ' + str(wtot) + ' um'])
                outwriter.writerow(['#'])

                # Write Radiative Recombination coeffcients to the output file.
                if __name__ == '__main__':
                        for lk in range(len(e.layers)):
                                outwriter.writerow(['Brad', e.layers[lk].name, e.layers[lk].material, e.layers[lk].Brad, e.layers[lk].thickness ])

                print('\n*** Computing Optical Coupling:')
                print('\n*** Angles: {num} angles from {start} to {stop} (radians); integration method: {method}'.format(num = len(angles), 
                                                                                                                         start = angles[0], 
                                                                                                                         stop = angles[-1], 
                                                                                                                         method = 'linearly spaced, trapezoidal integration' if use_gaussian_quadrature == False else 'guassian quadrature integration'))

                # set up lists for TMM calc.  n_list will be calculated per E value.  The lists include the semi-infinite media
                # on either side of the layer stack.  Air on top, same material as last layer on the bottom.
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

                #   Now set up for multiprocessing and execute the parallel calculation.
                #   Make a lock for the output file so only one process can work on it at a time.

                # flush write buffer to file _before_ forking -- otherwise we fork the buffer contents and get repeated lines in the file.
                outfile.flush()

                # sub-processes must acquire the lock before writing to the output file, then release it.  This way only one sub-process writes
                # to file at any time.
                fileLock =  Lock()

                # Define a pool of processes for parallel processing of the LCMatrix calculation.  This will create as many
                # worker processes as there are processors in the machine.
                # For debugging, it could be useful to use Pool(1).  Then only one worker process is created
                # and the calculation is done sequentially.  The main process forks at this point, so all sub-processes start
                # with the same state as the main process at this point.
                pool = Pool()

                res = []
                # submit jobs to the process pool.  Each job calculates luminescent coupling of emissions from
                # one layer in the stack into all other layers and writes the results to file.  Results for job i will be avaialble
                # in res[i] once the job is finished.

                # we sort the jobs before submiting them to the queue so that layers with the most y-values are processed first.  This should lead
                # to more even distribution of work across the processors and shorter overall run time. (i.e. we avoid a situation where all processes have finished
                # except one thread which processes all y-values in a thick substrate.)

                # sort jobs by number of y-values
                job_list = []
                print( len(e.layers))
                for m in range(len(e.layers)):
                        job_list.append( (m, len(e.layers[m].yvalues)))
                job_list.sort(key=lambda tup: tup[1], reverse=True)

                # use the pooled jobs or skip and avoid multiprocessing
                if Pooling == True:
                    # submit jobs to queue
                    for tup in job_list:
                            # queue up calls to the process_layer() function, one for each layer.
                            res.append( pool.apply_async(process_layer, [tup[0]]) )
    
                    # Wait until all jobs have finished, then close the pool so no new jobs can be submitted, and stop all sub-processes.
                    for r in res:
                            r.wait()
                            try:
                                    r.get()
                            except IndexError:
                                    print('IndexError')
                    pool.close()
                    pool.join()
                else:
                    # process_layer(4)
                    for tup in job_list:
                        process_layer(tup[0])

                # ... and we're done.
                endtime = datetime.datetime.now()
                print( '\n***')
                print( '*** LCMatrix Finshed at', endtime)
                print( '*** Runtime: ', (endtime-starttime))
                print( '***')
