# cython: profile=True, language_level=3, linetrace=True
# distutils: define_macros=CYTHON_TRACE_NOGIL=1

# Note: This module is written in Cython, which is a variant of Python intended to be converted to
# C code and compiled.  The main difference is the use of cdef statements and defining variable types.

# make division of integers return a float as in future versions of Python.

from __future__ import division
import os,sys

import csv
import numpy as np
cimport numpy as np
import sdevice_cmd
from scipy.interpolate import interp1d, UnivariateSpline
from scipy.integrate import trapezoid as trapz
from math import cos, sin, exp, pi
import h5py
import re
import tdr
cimport  tmm.tmm_core_mw as tmm
#import blist
from tables import *


cdef double e_charge = 1.602e-19 # C
cdef double h = 6.626e-34      # J.s
cdef double c0 = 299792458.0   # m/s
cdef double kb = 1.3806488e-23 # J/K

test_var = 0            
# for collecting all the data and storing it in a pytable
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

# ***************************************************
# Class with data on a specific layer in the epi file
#    epi



cdef class epilayer:
    cdef public double doping, Eg_300, xmole, ymole
    cdef public double ytop, ybot, thickness, Brad
    cdef public name, material, material_type, materialfile, ext
    cdef public par, index_n, index_k, nk_data
    cdef public np.ndarray E_list,  alpha_list, wl_list, n_list
    cdef public double E_min, E_max, wl_min, wl_max, wl_emission_min, wl_emission_max
    cdef public double S_hat_den
    cdef public yvalues
   
    def __init__(self, name, material, materialfile, thickness, doping,
                    molefraction,  ext, double ytop):
        self.name = name
        self.material = material
        if self.material in ['GaAsOx', 'TiOx', 'MgF', 'SiN', 'SiO2', 'Si3N4', 'Oxide', 'Nitride', 'TaO']:
           self.material_type = 'Dielectric'
        elif self.material in ['Gold']:
            self.material_type = 'Metal'
        elif self.material in ['Air']:
            self.material_type = 'Air/Vacuum'
        else:
           self.material_type = 'Semiconductor'

        self.materialfile = materialfile
        self.thickness = float(thickness)
        if len(doping) > 0 and (not molefraction.isspace()):
            self.doping = float(doping) # > 0: ptype, <0: ntype
        else:
            self.doping = 0.0

        if len(molefraction) > 0 and (not molefraction.isspace()):
            if molefraction[0] == '(':
                moles = molefraction.strip('()').split(' ')
                self.xmole = float(moles[0])
                self.ymole = float(moles[1])
            else:
                self.xmole = float(molefraction)
                self.ymole = 0.0
        self.ext = ext
        self.ytop = ytop
        self.ybot = ytop + self.thickness

# numerator of the spectral distribution of PL emission. Assume 300K.
# wavelength should be in um.
# parameter i is the index to the layer's E_list property.
    def S_hat_num(self,int i):
        cdef double alpha = self.alpha_list[i] # absorption coefficient in um^-1
        cdef double n = self.n_list[i]
        cdef double T = 300.0
        cdef double  E = self.E_list[i]
        cdef double S = alpha*(n*n)*(E*E)*exp(-E/(kb*T))
        return S

# denominator.  Integrate numerator over all energies.
    def calc_S_hat_den(self):
        cdef np.ndarray[double, ndim=1, negative_indices=False] S_num = np.empty(len(self.E_list), dtype='float')
        cdef double E

        for i in range(len(self.E_list)):
             S_num[i] = self.S_hat_num(i)

        self.S_hat_den = trapz(S_num, x=self.E_list)
        assert self.S_hat_den != 0.0, 'Found S_hat_den == NaN, material {0}, k_list: {1}'.format(self.material, self.alpha_list)


# Normalized spectral emission distribution.  If denominator has not already been calculated, it will be done automatically.
# wavelength should be in um.
# parameter i is the index to the layer's E_list property.
    def S_hat(self, int i):
        if self.material_type == 'Semiconductor':
           if self.S_hat_den == 0.0:
                self.calc_S_hat_den()
           try:
               return self.S_hat_num(i)/self.S_hat_den
           except ZeroDivisionError:
               print( 'S_hat_den', self.S_hat_den)
        else:
           return 0.0

# complex index of refraction interpolated to arbitrary wavelength (um).
    def nk(self, double wl):
       return complex(self.index_n(wl), self.index_k(wl))


# **********************************************************
# class containing a list of all layers in the epi structure.
# Each layer is represented by an instance of class epilayer.
# Reads in a pre-processed epi file ( i.e. pp@node@_epi.cmd)
# 

class epifile:
    def __init__(self, filename):
        with open(filename, 'r') as csvfile:
            print( '\n*** Getting layer structure from ', filename)
            self.r = csv.reader(csvfile, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL, skipinitialspace=True)
            self.filename = filename
            self.layers = []
            self.backlayers = [] # add seperate layer list for back materials
            self.ytop = 0.0
            self.ybot = 0.0


            for row in self.r:
                if len(row) > 0:
                    if len(row[0].lstrip('\t')) >2:
                        if  row[0].lstrip('\t')[0] != '$':
                            self.add_layer(row, False)
                        elif row[0].lstrip('\t').split(' ')[0] == '$repeat':
                            N = int(row[0].split(' ')[1])
                            row = self.r.next()
                            replayers = []
                            while row[0].lstrip('\t') != '$end':
                                    replayer = epilayer(row[0].lstrip('\t'),row[1].lstrip('\t'),row[2],row[3], row[4], row[5], row[6],0)
                                    replayers.append(replayer)
                                    row = self.r.next()
                            for i in range(N):
                                for rl in replayers:
                                    self.add_layer([re.sub("\$i",str(i),rl.name),rl.material,rl.materialfile,str(rl.thickness),str(rl.doping), '(' + str(rl.xmole) +' '+str(rl.ymole)+')',rl.ext,self.ybot], False)


    # ---------------------------------------------------------------------------------

		# Sets the layers at the bottom of the stack
    def bottom_material(self, mat, mat_file, thick):
        print( '\n*** Adding additional layers to bottom of the stack')#
        for i in range(len(mat)):
            backlayer = [f'back_mat_{i}', mat[i], mat_file[i], thick[i], '\t', '\t', '\t'] #doping and moelfraction set to 0 by default (### NOTE; might want to be able to change these later), ext (mesh from Sentaurus) is left blank, ytop is ybot from whatever the previous layer was and is set in the add_layer function
            self.add_layer(backlayer, True)

		# ----------------------------------------------------------------------------------

    def add_layer(self, row, backmat):
                layer = epilayer(row[0].lstrip('\t'), row[1].lstrip('\t'), row[2].lstrip('\t'), row[3], row[4], row[5], row[6],self.ybot)
                self.ybot = layer.ybot
                if backmat == False: #distinguish between a layer in the main stack where coupling sould be considered and one where coupling can be ignored
                    self.layers.append(layer)
                else:
                    self.backlayers.append(layer)
                print( '  ', layer.name, '\t', layer.material, '\t', layer.thickness, 'um')

    # -------------------------------------------------------------------------------------
    # read par files generated by MatPar.  file argument is the main MatPar file containing 
    # names of sub-files, e.g. n137_mpr.par  Extract material parameters that we need for
    # calculations and add the data to the corresponding layer objects..

    def read_parfiles(self, file, angle_list = np.array(np.linspace(0., 1.55, 20), dtype='float64'),
                      N_energies = 30):
        cdef  epilayer i
        print( '\n*** Getting additional material parameters')
        parfiles = sdevice_cmd.sdevice_file(file)

	# list of polar angles for numerical integration. Note - Avoid using exactly pi/2.
        self.theta_list = angle_list#np.array(np.linspace(0., 1.55, 20), dtype='float64') # angles too close to pi/2 cause problems, cut off at around 77 degrees (1.57 approx pi/2)
        if len(self.theta_list) == 1:
            print('Warning: theta_list contains only 1 value')
        for i in self.layers + self.backlayers: # get material params for layers and backlayers
		# Get extra material parameters from MatPar output.
                if i.name in parfiles.region_files:
                        parfile = sdevice_cmd.sdevice_file(parfiles.region_files[i.name])
                elif i.material in parfiles.material_files:
                        print(sdevice_cmd.sdevice_file(parfiles.material_files[i.material]))
                        parfile = sdevice_cmd.sdevice_file(parfiles.material_files[i.material])
                i.par = parfile.par

		# set up n,k for tmm library
                i.nk_data = []
                for j in range(len(i.par['TableODB']['wl'])):
                        i.nk_data.append( np.complex128(i.par['TableODB']['n'][j], i.par['TableODB']['k'][j] )) # updated np.complex to np.complex128 -pwils [2026-05-19]
                if i.material == 'TaO': # want a constant interpolation for this material -pwils
                    i.index_n = interp1d(i.par['TableODB']['wl'], i.par['TableODB']['n'], kind='linear', bounds_error = False, fill_value = (np.real(i.nk_data[0]), np.real(i.nk_data[-1])))
                    i.index_k = interp1d(i.par['TableODB']['wl'], i.par['TableODB']['k'], kind='linear', bounds_error = False, fill_value = (np.imag(i.nk_data[0]), np.imag(i.nk_data[-1])))
                else: # linear interpolation with linear extrapolation, a warning will be printed if the emission wavelengths are in an extrapolated region -pwils
                    print('trying')
                    print(i.material)
                    print(i.name)
                    i.index_n = interp1d(i.par['TableODB']['wl'], i.par['TableODB']['n'], kind='linear', bounds_error = False, fill_value = 'extrapolate')
                    i.index_k = interp1d(i.par['TableODB']['wl'], i.par['TableODB']['k'], kind='linear', bounds_error = False, fill_value = 'extrapolate')

		# Track wavelength range that is available for all materials.
                i.wl_max = i.par['TableODB']['wl'][-1]
                i.wl_min = i.par['TableODB']['wl'][0]

		# for semiconductors, cache an array of energies near the bandgap, covering the range of luminescence for this layer.
                if i.material_type == 'Semiconductor':
                   i.Eg_300 = parfile.Eg_T(300.0)*e_charge   # in J
                   i.E_min = i.Eg_300 # -10*kb*300 # at some point explore reducing E_min -pwils [2023-10-10]
                   i.E_max = i.Eg_300 + 10*kb*300
                   i.E_list = np.linspace(i.E_min, i.E_max, N_energies) # number of energies is now a variable -pwils [2023-10-10]
                   i.wl_list = 1e6*h*c0/i.E_list # in um
                   i.wl_emission_min = i.wl_list[-1] # shortest emission wavelength for the material -pwils
                   i.wl_emission_max= i.wl_list[0] #longest emission wavelenght for the material
		   # also cache alpha values to minimize interpolations.
                   i.alpha_list = 4*pi*i.index_k(i.wl_list)/i.wl_list	# um^-1
                   i.alpha_list[0] = 0.
                   i.n_list = i.index_n(i.wl_list)	# um^-1
                   i.Brad = i.par['RadiativeRecombination']['C']
                   print( '  ', i.name, '\t',i.material, 'Eg=', i.Eg_300/e_charge, 'eV', 'Brad=', i.Brad, 'wl_min=', i.wl_min, 'wl_max=', i.wl_max)

# keep track of the maximum and minimum wavelengths available in each of the .par files
                if hasattr(self,'wl_min'):
                        if i.wl_min > self.wl_min:
                                self.wl_min = i.wl_min
                        if i.wl_max < self.wl_max:
                                self.wl_max = i.wl_max
                else:
                        self.wl_min = i.wl_min
                        self.wl_max = i.wl_max

# keep track of the maximum and minimum emission wavelengths required for each of the semiconductors -pwils
                if i.material_type == "Semiconductor":
                    if hasattr(self,'wl_emission_min'):
                            if i.wl_emission_min < self.wl_emission_min:
                                    self.wl_emission_min = i.wl_emission_min
                            if i.wl_emission_max > self.wl_emission_max:
                                    self.wl_emission_max = i.wl_emission_max
                    else:
                            self.wl_emission_min = i.wl_emission_min
                            self.wl_emission_max = i.wl_emission_max
        print( '\n*** All materials have nk data in the range:  [', self.wl_min, 'to', self.wl_max, '] um')
        print('\n*** Maximum emission wavelength:', self.wl_emission_max, '   Minimum emission wavelength:', self.wl_emission_min)
        if self.wl_min > self.wl_emission_min or self.wl_max < self.wl_emission_max:
            print('\n***Warning: Wavelength range of emitted photons is outside the bounds for some of the material file nk data. This may result in extrapolating the necessary data.')



    # ----------------------------------------------------------------------------------------------------
    # Open a TDR mesh file and get vertices associated with each layer (and corresponding y-values).
    #
    def get_layer_vertices(self, tdr_file):

        print( '\n*** Finding unique y-coordinates in', tdr_file)
        t = tdr.tdr(tdr_file)
        cdef np.ndarray vertices = t.collection[u'geometry_0'][u'vertex'][()]
        regions = dict()
        # Find regions in TDR file and map keys (i.e. region_0) to names (i.e. emitter)
        for i in t.collection[u'geometry_0'].keys():
                parts = i.split('_')
                if len(parts) == 2 and parts[0] == 'region':
                        r= t.collection[u'geometry_0'][i]
                        # Map region names to TDR region indices
                        regions[r.attrs[u'name']] = i

        # now iterate over layers in epi-file.  Find mesh regions matching each layer. Make a sorted set containing
        # all unique y-coodinates of mesh vertices in the region.

        # indexing with [()] forces making a copy of the whole TDR dataset as a numpy array.
        # this gives a huge speed-up over reading data from the TDR file one element at a time.

        # In the TDR file, the dataset ['elements_0'] is a list of integers.  It consists of an integer indicating the
        # element type , followed by indices of each of the vertices for that element,
        # and repeated for every element in the region.  See the PMI section of the sdevice manual for the different possible element types.
        # type=1 is a line element (2 vertices) found at contacts in a 2D problem, and type=2 is a 2D triangle element (3 vertices).

        # ['vertex'] is a list of all vertices in the mesh. For each vertex it contains the 3D coordinates as an array of double floats.
        # The vertex indices in ['elements_0'] can be used to find corresponding vertex coordinates in [u'vertex'].
        
        cdef epilayer l
        cdef np.ndarray elements
        cdef size, vert
        cdef yval

        for l in self.layers:
                if l.name.encode('ascii') in regions:
                        # build list of unique y-values in the region.
                        #l.yvalues = blist.sortedset()  # removed blist since installation with pip has compatibility issues with setuptools and python3.13
                        l.yvalues = set()                 
                        elements = t.collection[u'geometry_0'][regions[l.name.encode('ascii')]][u'elements_0'][()]
                        elem = iter(elements)
                        size = elem.__next__()
                        done = False
                        while done == False:
                                for cnt in range(size+1):
                                        vert = elem.__next__()
                                        yval = vertices[vert][1]
                                        l.yvalues.add(yval)
                                try:
                                       size = elem.__next__()
                                except StopIteration:
                                       done = True
                # overwrite calculated layer boundaries with values from mesh file
                # avoids some problems due to roundoff errors.
                l.yvalues = sorted(l.yvalues) # sorting the values since blist.sortedset was removed                                                                                    
                l.ytop = l.yvalues[0]
                l.ybot = l.yvalues[-1]					
                print( '   ', l.name, ': ', len(l.yvalues), ' unique y-coords')

# Integrate absorption over E, theta.
def calc_a(int li, int lk, double yA, double yB, data_list_TE,  data_list_TM,  e,
           use_gauss_quad = False): # added guassian quadrature integration option -pwils
                      cdef np.ndarray E_list = e.layers[li].E_list
                      cdef np.ndarray th_list = e.theta_list
                      cdef np.ndarray P0_list = np.zeros((len(th_list)), dtype=float)
                      cdef np.ndarray P1_list = np.zeros((len(E_list)), dtype=float)
                      cdef np.ndarray sin_th_list = np.sin(th_list)
                      cdef double yB0 = yB - e.layers[lk].ytop

                      if use_gauss_quad == False: # the original trapezoidal rule calculation for the angle integration
                          for E_i in range(len(E_list)):
                               for th_i in range(len(th_list)-1):
                                      P0_list[th_i]  = tmm.position_resolved_a(lk, yB0, data_list_TE[E_i][th_i]) * sin_th_list[th_i]
                                      P0_list[th_i] += tmm.position_resolved_a(lk, yB0, data_list_TM[E_i][th_i]) * sin_th_list[th_i]
    
                               # Assume all emissions at angle pi/2 are absorbed at y-position of emission.
                               if yA  == yB:
                                      P0_list[-1] = 1.
    
                               P1_list[E_i] = trapz(P0_list, x=th_list) * e.layers[li].S_hat(E_i)
                          P = trapz(P1_list, x=E_list) 
    #                      if P/60. > 5.:
    #                            print( 'P1_list\n', P1_list, 'E_list\n', E_list)
                      else: # the calculation with Gaussian quadrature for the angle integration
                          # let the number of weights be equal to the 
                          # length of the linearly spaced th_list for now
                          # so that we can compare them more easily
                          n = len(th_list) 
                          # integral has been rewritten to apply over the 
                          # domain [-1, 1], standardized Gauss-Legendre
                          # weights apply
                          [u, w] = np.polynomial.legendre.leggauss(n)
                          for E_i in range(len(E_list)):
                              # using the change of variable u = 2*cos(theta)-1
                              for u_i in range(len(u)): #!!! Note: this was originally len(th_list) - 1, and the 90 degree case was done separately -pwils
                                  P0_list[u_i] = 0.5*tmm.position_resolved_a(lk, yB0, data_list_TE[E_i][u_i])
                                  P0_list[u_i] += 0.5*tmm.position_resolved_a(lk, yB0, data_list_TM[E_i][u_i])
                            
                              # Assume all emissions at angle pi/2 (u = -1) are absorbed at y-position of emission.
                              # if yA  == yB:
                                  # P0_list[0] = 1. #TODO: feel like I should double-check this
                                
                              P1_list[E_i] = np.sum(w*P0_list) # multiply by the weights and sum
                              P1_list[E_i] *= e.layers[li].S_hat(E_i) # multiply by normalized emission for the layer
                          P = trapz(P1_list, x = E_list)
                      return P

def calc_a_integrate_over_z(int li, double yA, data_list_TE,  data_list_TM,  e, P_list, integrate_theta_first,
            h5_filename, use_gauss_quad = False, include90 = True, wtot = 1, write_H5 = True):
    # """
    # Updated version of calc_a function. Currently being used to test bringing
    # the absor function integral (sum) to 1 over all given zB. Instead of a 
    # single lk and yB value, the loop over all yB is done inside this function
    # rather than outside. The function returns a list of P values rather than a
    # single value.
    #
    # Parameters
    # ----------
    # int li : Layer i (where emission occurs)
    # int lk : Layer k (where absorption occurs)
    # double yA : Position where emission occurs (yA psotion is with respect to 
    #             the entire device, yA0 is with respect to the layer itself)
    # double yB : Position where re-absorption occurs (yB psotion is with respect 
    #             to the entire device, yA0 is with respect to the layer itself)
    # data_list_TE : TMM calculation data produced from TE polarized emissions
    # data_list_TM : TMM calculation data produced from TM polarized emissions
    # e : class instance produced from epi_cmd_LC, contains information about
    #     layer structures
    # integrate_theta_first : what order the integration should be performed in, 
    #                         if True the integral over theta is P0_list, else
    #                         P0_list is the integral over E
    # use_guass_quad : TYPE, optional
    #     Whether to us guassian quadrature to perform the integral over angle.
    #     If True gaussian quadrature is used. If false trapezoidal integration
    #     is used. The default is False.

    # Returns
    # -------
    # P0_list : The theta part
    # P1_list : The E part
    # P : The result of the full integration. P should be the same regardless of 
    #     whether integrate_theta_first is true. This is a good check to make sure
    #     the integration has been performed correctly.

    # """
    cdef np.ndarray E_list = e.layers[li].E_list
    cdef np.ndarray th_list = e.theta_list
    # cdef np.ndarray P0_list = np.zeros((len(th_list)), dtype=float)
    # cdef np.ndarray P1_list = np.zeros((len(E_list)), dtype=float)
    # cdef list list_out = []
    cdef np.ndarray sin_th_list = np.sin(th_list)
    # cdef double yB0 = yB - e.layers[lk].ytop
    
    # original method used the trapezoidal method to perform
    # the integral over theta and E
    if use_gauss_quad == False:
        # P0_list is calculated over thetam P1_list is calculated over E
        # if integrate_E is true the E integral is performed first and P0_list is with respect to theta
        # if integrate_th is true the theta integral is performed first and P1_list is with respect to E
        # if write_H5 is True then write all data to H5 files, THIS CAN RESULT IN A VERY LARGE FILE
        if write_H5 == True:
            LCData = open_file(h5_filename, mode = 'r+', title = 'Test')
            Pos_table = LCData.root.Pos_info.Pos_data
            Pos_table_scaled = LCData.root.Pos_scaled_info.Pos_scaled_data
            pos_row = Pos_table.row
            pos_row_scaled = Pos_table_scaled.row
        
        full_lk_list = []
        full_yk_list = []
        x_list = [[] for i in range(len(e.layers))] #another list because I didn't plan things out right the first time
        P_list_yA = [[] for i in range(len(e.layers))]#np.zeros(len(e.layers)) # this is now inside the function, but it will sum up the contributions from yA throughout the device
        P_list_int = np.zeros(len(e.layers)) # this will hold integral of scaled values for each layer
        #collects all the layer and position data
        for lk in range(len(e.layers)):
            #create a list of pairs (layers, yB) absorption layer and absorption position
            lk_list = [lk for position in e.layers[lk].yvalues if e.layers[lk].material_type == 'Semiconductor']
            yk_list = [position for position in e.layers[lk].yvalues if e.layers[lk].material_type == 'Semiconductor']
            try:
                if lk != range(len(e.layers))[-1]:
                    lk_list.pop() # remove duplicates
                    yk_list.pop() # the last yvalue in the layer list is the same as the first value in the next layer, remove it from the full list to prevent duplicates
            except IndexError: # if the material is not a semiconductor the
                pass           # list will be empty and these values will not
                               # be calculated anyway
                               # TODO: This could cause problems if the non-
                               # semiconductor material is in the middle
                               # of the slab

            full_lk_list += lk_list
            full_yk_list += yk_list

        # P0 and P1 are now a matrix with the yk values as columns and either theta or E as rows
        P0_list = np.zeros((len(th_list), len(full_yk_list)), dtype=float)
        front_loss0_list = np.zeros(len(th_list), dtype=float)
        back_loss0_list = np.zeros(len(th_list), dtype=float)
        P1_list = np.zeros((len(E_list), len(full_yk_list)), dtype=float)
        front_loss1_list = np.zeros(len(E_list), dtype=float)
        back_loss1_list = np.zeros(len(E_list), dtype=float)
#!!!***************************************************************************
        if integrate_theta_first == True:
            for E_i in range(len(E_list)):
                for th_i in range(len(th_list)-1): # -1 assumes the last angle in the list is 90 degrees and will be dealt with separately
                    # take the outer integral and put it in this function
                    # for lk in range(len(e.layers)):
                    #     #create a list of pairs (layers, yB) absorption layer and absorption position
                    #     full_lk_list = [lk for position in e.layers[lk].yvalues if e.layers[lk].material_type == 'Semiconductor']
                    #     full_yk_list = [position for position in e.layers[lk].yvalues if e.layers[lk].material_type == 'Semiconductor']
                    all_a_TE = np.zeros(len(full_yk_list))
                    all_a_TM = np.zeros(len(full_yk_list))
                    
                    for index in range(len(full_yk_list)):
                        lk = full_lk_list[index]
                        yB = full_yk_list[index]
                        yB0 = yB - e.layers[full_lk_list[index]].ytop # pair[0] = lk

                        a_TE_abs, Ef1_TE, Eb1_TE, Ef2_TE, Eb2_TE, power_output1_TE, power_output2_TE, f_loss_TE, b_loss_TE = tmm.position_resolved_a_test_h5(lk, yB0, data_list_TE[E_i][th_i])
                        a_TM_abs, Ef1_TM, Eb1_TM, Ef2_TM, Eb2_TM, power_output1_TM, power_output2_TM, f_loss_TM, b_loss_TM = tmm.position_resolved_a_test_h5(lk, yB0, data_list_TM[E_i][th_i])
                        
                        all_a_TE[index] = a_TE_abs
                        all_a_TM[index] = a_TM_abs
                        if write_H5 == True:
                            pos_row['zA'] = yA
                            pos_row['zB'] = yB
                            pos_row['angle'] = th_list[th_i]
                            pos_row['energy'] = E_list[E_i]
                            pos_row['a_TE'] = a_TE_abs
                            pos_row['a_TM'] = a_TM_abs
                            pos_row['Ef1_TE'] = Ef1_TE
                            pos_row['Ef2_TE'] = Ef2_TE
                            pos_row['Eb1_TE'] = Eb1_TE
                            pos_row['Eb2_TE'] = Eb2_TE
                            pos_row['Ef1_TM'] = Ef1_TM
                            pos_row['Ef2_TM'] = Ef2_TM
                            pos_row['Eb1_TM'] = Eb1_TM
                            pos_row['Eb2_TM'] = Eb2_TM
                            pos_row.append()

                    # scale the values of a so that the area under the integral is <= 1 when the discrete mesh points are used
                    a_TE_interp = UnivariateSpline(full_yk_list, all_a_TE, k = 1, s = 0)
                    a_TM_interp = UnivariateSpline(full_yk_list, all_a_TM, k = 1, s = 0)

                    # take integral around each point, starting and ending midway between the points on either side
                    # endpoints are treated separately
                    scaled_a_TE = np.zeros(len(full_yk_list))
                    scaled_a_TM = np.zeros(len(full_yk_list))
                    dz_list = np.zeros(len(full_yk_list))
                    for index in np.arange(1, len(full_yk_list)-1):
                        x1 = 0.5*(full_yk_list[index-1] + full_yk_list[index])
                        x2 = 0.5*(full_yk_list[index] + full_yk_list[index+1])
                        dz = x2-x1
                        dz_list[index] = dz
                        scaled_a_TE[index] = a_TE_interp.integral(x1, x2)/dz
                        scaled_a_TM[index] = a_TM_interp.integral(x1, x2)/dz
                        #######################################################
                        if write_H5 == True:
                            pos_row_scaled['zA'] = yA
                            pos_row_scaled['zB'] = full_yk_list[index]
                            pos_row_scaled['angle'] = th_list[th_i]
                            pos_row_scaled['energy'] = E_list[E_i]
                            pos_row_scaled['a_TE_scaled'] = scaled_a_TE[index]
                            pos_row_scaled['a_TM_scaled'] = scaled_a_TM[index]
                            pos_row_scaled.append()
                            ################################################
                        
                    # now do the endpoints
                    x1 = full_yk_list[0]
                    x2 = 0.5*(full_yk_list[0] + full_yk_list[1])
                    dz = x2-x1
                    dz_list[0] = dz
                    scaled_a_TE[0] = a_TE_interp.integral(x1, x2)/dz
                    scaled_a_TM[0] = a_TM_interp.integral(x1, x2)/dz
                    ########################################################
                    if write_H5 == True:
                        pos_row_scaled['zA'] = yA
                        pos_row_scaled['zB'] = full_yk_list[0]
                        pos_row_scaled['angle'] = th_list[th_i]
                        pos_row_scaled['energy'] = E_list[E_i]
                        pos_row_scaled['a_TE_scaled'] = scaled_a_TE[0]
                        pos_row_scaled['a_TM_scaled'] = scaled_a_TM[0]
                        pos_row_scaled.append()
                    ########################################################
                    x1 = 0.5*(full_yk_list[-2] + full_yk_list[-1])
                    x2 = full_yk_list[-1]
                    dz=x2-x1
                    dz_list[-1] = dz
                    scaled_a_TE[-1] = a_TE_interp.integral(x1, x2)/dz
                    scaled_a_TM[-1] = a_TM_interp.integral(x1, x2)/dz
                    ########################################################
                    if write_H5 == True:
                        pos_row_scaled['zA'] = yA
                        pos_row_scaled['zB'] = full_yk_list[-1]
                        pos_row_scaled['angle'] = th_list[th_i]
                        pos_row_scaled['energy'] = E_list[E_i]
                        pos_row_scaled['a_TE_scaled'] = scaled_a_TE[-1]
                        pos_row_scaled['a_TM_scaled'] = scaled_a_TM[-1]
                        pos_row_scaled.append()
                    ########################################################

                    ph_absorbed = 0.5*(sum(scaled_a_TE*dz_list) + sum(scaled_a_TM*dz_list)) #!!! sum of the integrated portions
                    front_loss = 0.5*(f_loss_TE + f_loss_TM)
                    back_loss = 0.5*(b_loss_TE + b_loss_TM)
                    ph_loss = front_loss + back_loss #!!!

                    if abs(1 - (ph_loss + ph_absorbed)) >= 0.3:
                    #print('divided list 2, x:', x_vals)
                # print(len(x_vals))
                #print('divided list 2, P', P_vals)
                # print(len(P_vals))
                # P_list_int[index] = trapz(P_vals, x = x_vals)                    print('*** Warning: sum of transmitted and absorbed photons != 1.')
                        print('    yA = {0} theta = {1} E = {2}'.format(yA, th_list[th_i], E_list[E_i]))
                        #print('    all_a_TE:')
                        #print(list(all_a_TE))
                        #print('    all_a_TM:')
                        #print(list(all_a_TM))
                        #print('    all_yk:')
                        #print(list(full_yk_list))
                        print('    Absorbed: {0} Front loss: {1} Back loss: {2}'.format(ph_absorbed, front_loss, back_loss))
                    P0_list[th_i]  = scaled_a_TE * sin_th_list[th_i]
                    P0_list[th_i] += scaled_a_TM * sin_th_list[th_i]
                    front_loss0_list[th_i] = front_loss* sin_th_list[th_i]
                    back_loss0_list[th_i] = back_loss* sin_th_list[th_i]

                    #pos_row.append()

                # Assume all emissions at angle pi/2 are absorbed at y-position of emission
                # If a takes a delta function shape at the point scale the value by dz where dz is the width of the element.
                deg90_case = np.zeros(len(full_yk_list))
                for index in range(len(full_yk_list)):
                    if yA == full_yk_list[index]:
                        if index == 0: # if its the first point in the list
                            x1 = full_yk_list[0]
                            x2 = 0.5*(full_yk_list[0] + full_yk_list[1])
                            dz = x2-x1
                        elif index == len(full_yk_list) - 1: # if its the last point in the list
                            x1 = 0.5*(full_yk_list[-2] + full_yk_list[-1])
                            x2 = full_yk_list[-1]
                            dz=x2-x1
                        else:
                            x1 = 0.5*(full_yk_list[index-1] + full_yk_list[index])
                            x2 = 0.5*(full_yk_list[index] + full_yk_list[index + 1])
                            dz = x2 - x1
                        deg90_case[index] = 1/dz
                        ####################################################
                        if write_H5 == True:
                            pos_row_scaled['zA'] = yA
                            pos_row_scaled['zB'] = full_yk_list[index]
                            pos_row_scaled['angle'] = th_list[th_i]
                            pos_row_scaled['energy'] = E_list[E_i]
                            pos_row_scaled['a_TE_scaled'] = 1/dz # I think this will get scaled by 2 later but might be good to doublecheck at some point -pwils
                            pos_row_scaled['a_TM_scaled'] = 1/dz
                            pos_row_scaled.append()
                        ####################################################
                    else:
                        ####################################################
                        if write_H5 == True:
                            pos_row_scaled['zA'] = yA
                            pos_row_scaled['zB'] = full_yk_list[index]
                            pos_row_scaled['angle'] = th_list[th_i]
                            pos_row_scaled['energy'] = E_list[E_i]
                            pos_row_scaled['a_TE_scaled'] = 0.
                            pos_row_scaled['a_TM_scaled'] = 0.
                            pos_row_scaled.append()
                        ####################################################
                if include90 == False:
                    P0_list[-1] = np.zeros(len(deg90_case))
                else:
                    P0_list[-1] = np.array(deg90_case)

                # axis = 0 will integrate along a column, return integral over all theta for given yB
                P1_list[E_i] = trapz(P0_list, x=th_list, axis = 0) * e.layers[li].S_hat(E_i)
                front_loss1_list[E_i] = trapz(front_loss0_list, x=th_list) * e.layers[li].S_hat(E_i)
                back_loss1_list[E_i] = trapz(back_loss0_list, x=th_list) * e.layers[li].S_hat(E_i)
                ############################################################
                for index in range(len(full_yk_list)):
                    if write_H5 == True:
                        pos_row_scaled['zA'] = yA
                        pos_row_scaled['zB'] = full_yk_list[index]
                        pos_row_scaled['energy'] = E_list[E_i]
                        pos_row_scaled['a_int_angle_scaled'] = P1_list[E_i][index]
                        pos_row_scaled.append()
                ############################################################
            # axis = 0 will integrate along a column, return integral over all energies for given yB
            P = trapz(P1_list, x=E_list, axis = 0)
            FrontLoss = trapz(front_loss1_list, x=E_list)
            BackLoss = trapz(back_loss1_list, x=E_list)
            ################################################################
            for index in range(len(full_yk_list)):
                if write_H5 == True:
                    pos_row_scaled['zA'] = yA
                    pos_row_scaled['zB'] = full_yk_list[index]
                    pos_row_scaled['a_coupling_scaled'] = P[index]
                    pos_row_scaled.append()

                new_tup = (yA, full_yk_list[index], P[index]/(2*wtot))
                P_list[full_lk_list[index]].append(new_tup)
                
                if index == 0: # if its the first point in the list
                    x1 = full_yk_list[0]
                    x2 = 0.5*(full_yk_list[0] + full_yk_list[1])
                    dz = x2-x1
                elif index == len(full_yk_list) - 1: # if its the last point in the list
                    x1 = 0.5*(full_yk_list[-2] + full_yk_list[-1])
                    x2 = full_yk_list[-1]
                    dz=x2-x1
                else:
                    x1 = 0.5*(full_yk_list[index-1] + full_yk_list[index])
                    x2 = 0.5*(full_yk_list[index] + full_yk_list[index + 1])
                    dz = x2 - x1

                x_list[full_lk_list[index]].append(full_yk_list[index])
                P_list_yA[full_lk_list[index]].append(P[index]*dz/(2*wtot))

            for index in range(len(e.layers)):
                if index != range(len(e.layers))[-1]:
                    if len(x_list[index+1]) != 0: # if the material is not a
                                                  # semiconductor the x_list
                                                  # for that index will be
                                                  # empty
                        x_vals = x_list[index]
                        x_vals.append(x_list[index+1][0])
                        P_vals = P_list_yA[index]
                        P_vals.append(P_list_yA[index+1][0])
                    else:
                        pass
                else:
                    x_vals = x_list[index]
                    P_vals = P_list_yA[index]
                P_list_int[index] = sum(P_list_yA[index])
            #if sum(P_list_int) >= 1:
            #    print("*** Warning: efficiency larger than 1***")
            #    print('P_list_int:', P_list_int)
            #    print(list(all_a_TE))
            #    print(list(all_a_TM))
            #    print(list(scaled_a_TE))
            #    print(list(scaled_a_TM))
                # P_list_yA[full_lk_list[index]] = trapz()
            ###########################################################
#!!!***************************************************************************
        else:
            print("'Integrate energy first' option option doesn't work yet")
        #     for th_i in range(len(th_list)): #!!! Note: this was originally len(th_list) - 1, and the 90 degree case was done separately -pwils
        #         for E_i in range(len(E_list)):
        #             pos_row['zA'] = yA
        #             pos_row['zB'] = yB
        #             pos_row['angle'] = th_list[th_i]
        #             pos_row['energy'] = E_list[E_i]
                    
        #             # print(tmm.position_resolved_a(lk, yB0, data_list_TE[E_i][th_i]))
        #             # a_TE_abs = tmm.position_resolved_a_test_h5(lk, yB0, data_list_TE[E_i][th_i])[0]
        #             # a_TM_abs = tmm.position_resolved_a_test_h5(lk, yB0, data_list_TM[E_i][th_i])[0]
        #             # pos_row['a_TE'] = a_TE_abs
        #             # pos_row['a_TM'] = a_TM_abs
        #             a_TE_abs, Ef1_TE, Eb1_TE, Ef2_TE, Eb2_TE, power_output1_TE, power_output2_TE = tmm.position_resolved_a_test_h5(lk, yB0, data_list_TE[E_i][th_i])
        #             a_TM_abs, Ef1_TM, Eb1_TM, Ef2_TM, Eb2_TM, power_output1_TM, power_output2_TM = tmm.position_resolved_a_test_h5(lk, yB0, data_list_TM[E_i][th_i])
        #             pos_row['a_TE'] = a_TE_abs
        #             pos_row['a_TM'] = a_TM_abs
        #             pos_row['Ef1_TE'] = Ef1_TE
        #             pos_row['Ef2_TE'] = Ef2_TE
        #             pos_row['Eb1_TE'] = Eb1_TE
        #             pos_row['Eb2_TE'] = Eb2_TE
        #             pos_row['Ef1_TM'] = Ef1_TM
        #             pos_row['Ef2_TM'] = Ef2_TM
        #             pos_row['Eb1_TM'] = Eb1_TM
        #             pos_row['Eb2_TM'] = Eb2_TM
                    
        #             P1_list[E_i]  = a_TE_abs * e.layers[li].S_hat(E_i)
        #             P1_list[E_i] += a_TM_abs * e.layers[li].S_hat(E_i)
                    
        #             pos_row.append()
        #         P0_list[th_i] = trapz(P1_list, x=E_list)* sin_th_list[th_i] # * e.layers[li].S_hat(E_i)
        #         pos_row['zA'] = yA
        #         pos_row['zB'] = yB
        #         pos_row['angle'] = th_list[th_i]
        #         pos_row['a_int_E'] = P0_list[th_i]
        #         pos_row.append()
        #         # Assume all emissions at angle pi/2 are absorbed at y-position of emission.
        #     if yA == yB:
        #         P0_list[-1] = 1.#* e.layers[li].S_hat(E_i)
        #     P = trapz(P0_list, x=th_list) 
        #     pos_row['zA'] = yA
        #     pos_row['zB'] = yB
        #     pos_row['a_coupling'] = P0_list[th_i]
        #     pos_row.append()
        # Pos_table.flush()
        # LCData.close()
        
    #                      if P/60. > 5.:
    #                            print( 'P1_list\n', P1_list, 'E_list\n', E_list
    else:
        print("'Use gaussian quadrature' option doesn't work yet")
        # # let the number of weights be equal to the 
        # # length of the linearly spaced th_list for now
        # # so that we can compare them more easily
        # n = len(th_list) 
        # # integral has been rewritten to apply over the 
        # # domain [-1, 1], standardized Gauss-Legendre
        # # weights apply
        # [u, w] = np.polynomial.legendre.leggauss(n)
        # if integrate_theta_first == True:
        #     for E_i in range(len(E_list)):
        #         # using the change of variable u = 2*cos(theta)-1
        #         for u_i in range(len(u)): #TODO: check the pi/2 case TE_data is one angle short #!!! Note: this was originally len(th_list) - 1, and the 90 degree case was done separately -pwils
        #             P0_list[u_i] = 0.5*tmm.position_resolved_a(lk, yB0, data_list_TE[E_i][u_i])
        #             P0_list[u_i] += 0.5*tmm.position_resolved_a(lk, yB0, data_list_TM[E_i][u_i])
                
                    
                
        #         # Assume all emissions at angle pi/2 (u = -1) are absorbed at y-position of emission.
        #         # if yA  == yB:
        #             # P0_list[0] = 1. #TODO: feel like I should double-check this
                    
        #         P1_list[E_i] = np.sum(w*P0_list) # multiply by the weights and sum
        #         P1_list[E_i] *= e.layers[li].S_hat(E_i) # multiply by normalized emission for the layer
        #     P = trapz(P1_list, x = E_list)
        # else:
        #     for u_i in range(len(u)): #!!! Note: this was originally len(th_list) - 1, and the 90 degree case was done separately -pwils
        #         for E_i in range(len(E_list)):
        #             P1_list[E_i]  = tmm.position_resolved_a(lk, yB0, data_list_TE[E_i][u_i]) * e.layers[li].S_hat(E_i)
        #             P1_list[E_i] += tmm.position_resolved_a(lk, yB0, data_list_TM[E_i][u_i]) * e.layers[li].S_hat(E_i)
                
        #         P0_list[u_i] = trapz(P1_list, x = E_list)
                
        #         # Assume all emissions at angle pi/2 (u = -1) are absorbed at y-position of emission.
        #         # if yA  == yB:
        #             # P0_list[0] = 1. #TODO: feel like I should double-check this
                
        #         P = np.sum(w*0.5*P0_list) # multiply by the weights and sum
    #print(full_yk_list)
    return P0_list, P1_list, P, P_list, P_list_int, FrontLoss, BackLoss


# Integrate absorption over either E , theta, or not at all depending on what extra data is reqiured
def calc_a_ordered_integ_writeH5_test(int li, int lk, double yA, double yB, data_list_TE,  data_list_TM,  e, integrate_theta_first,
            h5_filename, use_gauss_quad = False):
    # """
    # Updated version of calc_a function. Currently being used in get_TMM_data to 
    # examine the dependence of emission angle and test the implementation of
    # Gaussian quadrature integration [2023-07-19]

    # Parameters
    # ----------
    # int li : Layer i (where emission occurs)
    # int lk : Layer k (where absorption occurs)
    # double yA : Position where emission occurs (yA psotion is with respect to 
    #             the entire device, yA0 is with respect to the layer itself)
    # double yB : Position where re-absorption occurs (yB psotion is with respect 
    #             to the entire device, yA0 is with respect to the layer itself)
    # data_list_TE : TMM calculation data produced from TE polarized emissions
    # data_list_TM : TMM calculation data produced from TM polarized emissions
    # e : class instance produced from epi_cmd_LC, contains information about
    #     layer structures
    # integrate_theta_first : what order the integration should be performed in, 
    #                         if True the integral over theta is P0_list, else
    #                         P0_list is the integral over E
    # use_guass_quad : TYPE, optional
    #     Whether to us guassian quadrature to perform the integral over angle.
    #     If True gaussian quadrature is used. If false trapezoidal integration
    #     is used. The default is False.

    # Returns
    # -------
    # P0_list : The theta part
    # P1_list : The E part
    # P : The result of the full integration. P should be the same regardless of 
    #     whether integrate_theta_first is true. This is a good check to make sure
    #     the integration has been performed correctly.

    # """
    cdef np.ndarray E_list = e.layers[li].E_list
    cdef np.ndarray th_list = e.theta_list
    cdef np.ndarray P0_list = np.zeros((len(th_list)), dtype=float)
    cdef np.ndarray P1_list = np.zeros((len(E_list)), dtype=float)
    # cdef list list_out = []
    cdef np.ndarray sin_th_list = np.sin(th_list)
    cdef double yB0 = yB - e.layers[lk].ytop
    
    # original method used the trapezoidal method to perform
    # the integral over theta and E
    if use_gauss_quad == False:
        # P0_list is calculated over thetam P1_list is calculated over E
        # if integrate_E is true the E integral is performed first and P0_list is with respect to theta
        # if integrate_th is true the theta integral is performed first and P1_list is with respect to E
        LCData = open_file(h5_filename, mode = 'r+', title = 'Test')
        #Pos_info = LCData.create_group('/', 'Pos_info')
        #Pos_table = LCData.create_table(Pos_info, 'Pos_data', LCData_Pos)
        Pos_table = LCData.root.Pos_info.Pos_data
        pos_row = Pos_table.row
        
        
        if integrate_theta_first == True:
            for E_i in range(len(E_list)):
                for th_i in range(len(th_list)): #!!! Note: this was originally len(th_list) - 1, and the 90 degree case was done separately -pwils
                    pos_row['zA'] = yA
                    pos_row['zB'] = yB
                    pos_row['angle'] = th_list[th_i]
                    pos_row['energy'] = E_list[E_i]
                    
                    # print(tmm.position_resolved_a(lk, yB0, data_list_TE[E_i][th_i]))
                    # a_TE_abs = tmm.position_resolved_a(lk, yB0, data_list_TE[E_i][th_i])
                    # a_TM_abs = tmm.position_resolved_a(lk, yB0, data_list_TM[E_i][th_i])
                    # list_out = tmm.position_resolved_a_test_h5(lk, yB0, data_list_TE[E_i][th_i])
                    # [absor_TE, Ef1_TE, Eb1_TE, Ef2_TE, Eb2_TE, power_output1_TE, power_output2_TE]
                    a_TE_abs, Ef1_TE, Eb1_TE, Ef2_TE, Eb2_TE, power_output1_TE, power_output2_TE = tmm.position_resolved_a_test_h5(lk, yB0, data_list_TE[E_i][th_i])
                    a_TM_abs, Ef1_TM, Eb1_TM, Ef2_TM, Eb2_TM, power_output1_TM, power_output2_TM = tmm.position_resolved_a_test_h5(lk, yB0, data_list_TM[E_i][th_i])
                    pos_row['a_TE'] = a_TE_abs
                    pos_row['a_TM'] = a_TM_abs
                    pos_row['Ef1_TE'] = Ef1_TE
                    pos_row['Ef2_TE'] = Ef2_TE
                    pos_row['Eb1_TE'] = Eb1_TE
                    pos_row['Eb2_TE'] = Eb2_TE
                    pos_row['Ef1_TM'] = Ef1_TM
                    pos_row['Ef2_TM'] = Ef2_TM
                    pos_row['Eb1_TM'] = Eb1_TM
                    pos_row['Eb2_TM'] = Eb2_TM
                    
                    P0_list[th_i]  = a_TE_abs * sin_th_list[th_i]
                    P0_list[th_i] += a_TM_abs * sin_th_list[th_i]

                    pos_row.append()
                # Assume all emissions at angle pi/2 are absorbed at y-position of emission.
                if yA  == yB:
                    P0_list[-1] = 1.

                P1_list[E_i] = trapz(P0_list, x=th_list) * e.layers[li].S_hat(E_i)
                pos_row['zA'] = yA
                pos_row['zB'] = yB
                pos_row['energy'] = E_list[E_i]
                pos_row['a_int_angle'] = P1_list[E_i]
                pos_row.append()
            P = trapz(P1_list, x=E_list) 
            pos_row['zA'] = yA
            pos_row['zB'] = yB
            pos_row['a_coupling'] = P
            pos_row.append()
    #                      if P/60. > 5.:
    #                            print( 'P1_list\n', P1_list, 'E_list\n', E_list)
        else:
            for th_i in range(len(th_list)): #!!! Note: this was originally len(th_list) - 1, and the 90 degree case was done separately -pwils
                for E_i in range(len(E_list)):
                    pos_row['zA'] = yA
                    pos_row['zB'] = yB
                    pos_row['angle'] = th_list[th_i]
                    pos_row['energy'] = E_list[E_i]
                    
                    # print(tmm.position_resolved_a(lk, yB0, data_list_TE[E_i][th_i]))
                    # a_TE_abs = tmm.position_resolved_a_test_h5(lk, yB0, data_list_TE[E_i][th_i])[0]
                    # a_TM_abs = tmm.position_resolved_a_test_h5(lk, yB0, data_list_TM[E_i][th_i])[0]
                    # pos_row['a_TE'] = a_TE_abs
                    # pos_row['a_TM'] = a_TM_abs
                    a_TE_abs, Ef1_TE, Eb1_TE, Ef2_TE, Eb2_TE, power_output1_TE, power_output2_TE = tmm.position_resolved_a_test_h5(lk, yB0, data_list_TE[E_i][th_i])
                    a_TM_abs, Ef1_TM, Eb1_TM, Ef2_TM, Eb2_TM, power_output1_TM, power_output2_TM = tmm.position_resolved_a_test_h5(lk, yB0, data_list_TM[E_i][th_i])
                    pos_row['a_TE'] = a_TE_abs
                    pos_row['a_TM'] = a_TM_abs
                    pos_row['Ef1_TE'] = Ef1_TE
                    pos_row['Ef2_TE'] = Ef2_TE
                    pos_row['Eb1_TE'] = Eb1_TE
                    pos_row['Eb2_TE'] = Eb2_TE
                    pos_row['Ef1_TM'] = Ef1_TM
                    pos_row['Ef2_TM'] = Ef2_TM
                    pos_row['Eb1_TM'] = Eb1_TM
                    pos_row['Eb2_TM'] = Eb2_TM
                    
                    P1_list[E_i]  = a_TE_abs * e.layers[li].S_hat(E_i)
                    P1_list[E_i] += a_TM_abs * e.layers[li].S_hat(E_i)
                    
                    pos_row.append()
                P0_list[th_i] = trapz(P1_list, x=E_list)* sin_th_list[th_i] # * e.layers[li].S_hat(E_i)
                pos_row['zA'] = yA
                pos_row['zB'] = yB
                pos_row['angle'] = th_list[th_i]
                pos_row['a_int_E'] = P0_list[th_i]
                pos_row.append()
                # Assume all emissions at angle pi/2 are absorbed at y-position of emission.
            if yA == yB:
                P0_list[-1] = 1.#* e.layers[li].S_hat(E_i)
            P = trapz(P0_list, x=th_list) 
            pos_row['zA'] = yA
            pos_row['zB'] = yB
            pos_row['a_coupling'] = P0_list[th_i]
            pos_row.append()
        Pos_table.flush()
        LCData.close()
        
    #                      if P/60. > 5.:
    #                            print( 'P1_list\n', P1_list, 'E_list\n', E_list
    else:
        # let the number of weights be equal to the 
        # length of the linearly spaced th_list for now
        # so that we can compare them more easily
        n = len(th_list) 
        # integral has been rewritten to apply over the 
        # domain [-1, 1], standardized Gauss-Legendre
        # weights apply
        [u, w] = np.polynomial.legendre.leggauss(n)
        if integrate_theta_first == True:
            for E_i in range(len(E_list)):
                # using the change of variable u = 2*cos(theta)-1
                for u_i in range(len(u)): #TODO: check the pi/2 case TE_data is one angle short #!!! Note: this was originally len(th_list) - 1, and the 90 degree case was done separately -pwils
                    P0_list[u_i] = 0.5*tmm.position_resolved_a(lk, yB0, data_list_TE[E_i][u_i])
                    P0_list[u_i] += 0.5*tmm.position_resolved_a(lk, yB0, data_list_TM[E_i][u_i])
                
                    
                
                # Assume all emissions at angle pi/2 (u = -1) are absorbed at y-position of emission.
                # if yA  == yB:
                    # P0_list[0] = 1. #TODO: feel like I should double-check this
                    
                P1_list[E_i] = np.sum(w*P0_list) # multiply by the weights and sum
                P1_list[E_i] *= e.layers[li].S_hat(E_i) # multiply by normalized emission for the layer
            P = trapz(P1_list, x = E_list)
        else:
            for u_i in range(len(u)): #!!! Note: this was originally len(th_list) - 1, and the 90 degree case was done separately -pwils
                for E_i in range(len(E_list)):
                    P1_list[E_i]  = tmm.position_resolved_a(lk, yB0, data_list_TE[E_i][u_i]) * e.layers[li].S_hat(E_i)
                    P1_list[E_i] += tmm.position_resolved_a(lk, yB0, data_list_TM[E_i][u_i]) * e.layers[li].S_hat(E_i)
                
                P0_list[u_i] = trapz(P1_list, x = E_list)
                
                # Assume all emissions at angle pi/2 (u = -1) are absorbed at y-position of emission.
                # if yA  == yB:
                    # P0_list[0] = 1. #TODO: feel like I should double-check this
                
                P = np.sum(w*0.5*P0_list) # multiply by the weights and sum

    return P0_list, P1_list, P

def calc_a_def_int(lk, z1, z2, pol, power_output):
    # calculate the definite integral of a under its corresponding mesh element
    layer = lk
    power_output = data.power_output1
    w = vw_list_l1[layer][0]
    v = vw_list_l1[layer][1]
    n = n_list[layer+1]
    th = th_list[layer+1]
    kz = kz_list[layer+1]


    if pol == 'TE':
        A1 = np.imag(n*cos(th)*kz)*w*w.conjugate()/power_output
        A2 =np.imag(n*cos(th)*kz)*v*v.conjugate()/power_output
        A3 = np.imag(n*cos(th)*kz)*v*w.conjugate()/power_output
        A3_star = A3.conjugate()

    elif pol == 'TM':
        A1 = 2*np.imag(kz)*np.real(n*cos(th).conjugate())*w*w.conjugate()/power_output
        A2 =2*np.imag(kz)*np.real(n*cos(th).conjugate())*v*v.conjugate()/power_output
        A3 = -2*np.real(kz)*np.imag(n*cos(th).conjugate())*v*w.conjugate()/power_output
        A3_star = A3.conjugate()

    def_int = 0.5*(A1*np.exp(2*z1*np.imag(kz)) - A2*exp(-2*z1*np.imag(kz)) - 1j*A3*np.exp(2*1j*z1*np.real(kz)) + 1j*A3_star*exp(-2*1j*z1*np.real(kz)))/z1
    def_int -= 0.5*(A1*np.exp(2*z2*np.imag(kz)) - A2*exp(-2*z2*np.imag(kz)) - 1j*A3*np.exp(2*1j*z2*np.real(kz)) + 1j*A3_star*exp(-2*1j*z1*np.real(kz)))/z2

    return def_int


# Integrate absorption over either E , theta, or not at all depending on what extra data is reqiured
def calc_a_ordered_integ(int li, int lk, double yA, double yB, data_list_TE,  data_list_TM,  e, integrate_theta_first,
            use_gauss_quad = False):
    # """
    # Updated version of calc_a function. Currently being used in get_TMM_data to 
    # examine the dependence of emission angle and test the implementation of
    # Gaussian quadrature integration [2023-07-19]

    # Parameters
    # ----------
    # int li : Layer i (where emission occurs)
    # int lk : Layer k (where absorption occurs)
    # double yA : Position where emission occurs (yA psotion is with respect to 
    #             the entire device, yA0 is with respect to the layer itself)
    # double yB : Position where re-absorption occurs (yB psotion is with respect 
    #             to the entire device, yA0 is with respect to the layer itself)
    # data_list_TE : TMM calculation data produced from TE polarized emissions
    # data_list_TM : TMM calculation data produced from TM polarized emissions
    # e : class instance produced from epi_cmd_LC, contains information about
    #     layer structures
    # integrate_theta_first : what order the integration should be performed in, 
    #                         if True the integral over theta is P0_list, else
    #                         P0_list is the integral over E
    # use_guass_quad : TYPE, optional
    #     Whether to us guassian quadrature to perform the integral over angle.
    #     If True gaussian quadrature is used. If false trapezoidal integration
    #     is used. The default is False.

    # Returns
    # -------
    # P0_list : The theta part
    # P1_list : The E part
    # P : The result of the full integration. P should be the same regardless of 
    #     whether integrate_theta_first is true. This is a good check to make sure
    #     the integration has been performed correctly.

    # """
    cdef np.ndarray E_list = e.layers[li].E_list
    cdef np.ndarray th_list = e.theta_list
    cdef np.ndarray P0_list = np.zeros((len(th_list)), dtype=float)
    cdef np.ndarray P1_list = np.zeros((len(E_list)), dtype=float)
    cdef np.ndarray sin_th_list = np.sin(th_list)
    cdef double yB0 = yB - e.layers[lk].ytop
    
    # original method used the trapezoidal method to perform
    # the integral over theta and E
    if use_gauss_quad == False:
        # P0_list is calculated over thetam P1_list is calculated over E
        # if integrate_E is true the E integral is performed first and P0_list is with respect to theta
        # if integrate_th is true the theta integral is performed first and P1_list is with respect to E
        
        if integrate_theta_first == True:
            for E_i in range(len(E_list)):
                for th_i in range(len(th_list)): #!!! Note: this was originally len(th_list) - 1, and the 90 degree case was done separately -pwils
                    P0_list[th_i]  = tmm.position_resolved_a(lk, yB0, data_list_TE[E_i][th_i]) * sin_th_list[th_i]
                    P0_list[th_i] += tmm.position_resolved_a(lk, yB0, data_list_TM[E_i][th_i]) * sin_th_list[th_i]

                # Assume all emissions at angle pi/2 are absorbed at y-position of emission.
                if yA  == yB:
                    P0_list[-1] = 1.

                P1_list[E_i] = trapz(P0_list, x=th_list) * e.layers[li].S_hat(E_i)
            P = trapz(P1_list, x=E_list) 
    #                      if P/60. > 5.:
    #                            print( 'P1_list\n', P1_list, 'E_list\n', E_list)
        else:
            for th_i in range(len(th_list)): #!!! Note: this was originally len(th_list) - 1, and the 90 degree case was done separately -pwils
                for E_i in range(len(E_list)):
                    P1_list[E_i]  = tmm.position_resolved_a(lk, yB0, data_list_TE[E_i][th_i]) * e.layers[li].S_hat(E_i)
                    P1_list[E_i] += tmm.position_resolved_a(lk, yB0, data_list_TM[E_i][th_i]) * e.layers[li].S_hat(E_i)
    
                P0_list[th_i] = trapz(P1_list, x=E_list)* sin_th_list[th_i] # * e.layers[li].S_hat(E_i)
                # Assume all emissions at angle pi/2 are absorbed at y-position of emission.
            if yA == yB:
                P0_list[-1] = 1.#* e.layers[li].S_hat(E_i)
            P = trapz(P0_list, x=th_list) 
    #                      if P/60. > 5.:
    #                            print( 'P1_list\n', P1_list, 'E_list\n', E_list
    else:
        # let the number of weights be equal to the 
        # length of the linearly spaced th_list for now
        # so that we can compare them more easily
        n = len(th_list) 
        # integral has been rewritten to apply over the 
        # domain [-1, 1], standardized Gauss-Legendre
        # weights apply
        [u, w] = np.polynomial.legendre.leggauss(n)
        if integrate_theta_first == True:
            for E_i in range(len(E_list)):
                # using the change of variable u = 2*cos(theta)-1
                for u_i in range(len(u)): #TODO: check the pi/2 case TE_data is one angle short #!!! Note: this was originally len(th_list) - 1, and the 90 degree case was done separately -pwils
                    P0_list[u_i] = 0.5*tmm.position_resolved_a(lk, yB0, data_list_TE[E_i][u_i])
                    P0_list[u_i] += 0.5*tmm.position_resolved_a(lk, yB0, data_list_TM[E_i][u_i])
                
                # Assume all emissions at angle pi/2 (u = -1) are absorbed at y-position of emission.
                # if yA  == yB:
                    # P0_list[0] = 1. #TODO: feel like I should double-check this
                    
                P1_list[E_i] = np.sum(w*P0_list) # multiply by the weights and sum
                P1_list[E_i] *= e.layers[li].S_hat(E_i) # multiply by normalized emission for the layer
            P = trapz(P1_list, x = E_list)
        else:
            for u_i in range(len(u)): #!!! Note: this was originally len(th_list) - 1, and the 90 degree case was done separately -pwils
                for E_i in range(len(E_list)):
                    P1_list[E_i]  = tmm.position_resolved_a(lk, yB0, data_list_TE[E_i][u_i]) * e.layers[li].S_hat(E_i)
                    P1_list[E_i] += tmm.position_resolved_a(lk, yB0, data_list_TM[E_i][u_i]) * e.layers[li].S_hat(E_i)
                
                P0_list[u_i] = trapz(P1_list, x = E_list)
                
                # Assume all emissions at angle pi/2 (u = -1) are absorbed at y-position of emission.
                # if yA  == yB:
                    # P0_list[0] = 1. #TODO: feel like I should double-check this
                
                P = np.sum(w*0.5*P0_list) # multiply by the weights and sum

    return P0_list, P1_list, P


def calc_absor(int li, int lk, double yA, double yB, data_list_TE,  data_list_TM,  e, E_index, th_index):
                      cdef np.ndarray E_list = e.layers[li].E_list
                      cdef np.ndarray th_list = e.theta_list
                      cdef np.ndarray P0_list = np.zeros((len(th_list)), dtype=float)
                      cdef np.ndarray P1_list = np.zeros((len(E_list)), dtype=float)
                      cdef np.ndarray sin_th_list = np.sin(th_list)
                      cdef double yB0 = yB - e.layers[lk].ytop

                      # P0_list is calculated over thetam P1_list is calculated over E
                      # if integrate_E is true the E integral is performed first and P0_list is with respect to theta
                      # if integrate_th is true the theta integral is performed first and P1_list is with respect to E

                      absor  = tmm.position_resolved_a(lk, yB0, data_list_TE[E_index][th_index])
                      absor += tmm.position_resolved_a(lk, yB0, data_list_TM[E_index][th_index])
                      
                      return absor