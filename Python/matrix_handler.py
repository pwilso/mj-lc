# -*- coding: utf-8 -*-
"""
Created on Fri Nov 11 14:11:44 2022

@author: pwilson3
"""

import csv
import numpy as np
import sys
import pandas as pd

#------------------------------------------------------------------------------
# class for manipulating the luminescent coupling matrix files
class matrix_handler:
    """
    # Reads in the matrix file and stores information on the layer names,
    # materials, thicknesses, (x, y, z) coordinates
    # (x = absorption point, y = emission point, z = coupling coefficient)
    # NOTE: (xs, ys, zs) is the lower matrix (backward travelling emission, yA >= yB)
    #       (xs2, ys2, zs2) is the upper matrix (forward travelling emission, yB > yA)
    # yA refers to an emission point, yB is an absorption point
    """
    def __init__(self, filepath):
        # coordinate lists
        self.xs = []
        self.ys = []
        self.zs = []
        self.xs2 = []
        self.ys2 = []
        self.zs2 = []
        self.a = [] # this contains only the value of z elements as a 1D list
    
        self.large_values = []
        self.neg_values = []
        
        with open(filepath) as infile:
            inreader = csv.reader(infile)
    
            self.layer_combos = []
            self.layers = []
            self.thicknesses = []
            self.materials= []
            yB1 = []
            P1 = []
            yB2 = []
            P2 = []
            yB3 = []
            P3 = []
        
            for r in inreader:
                
                # skip beginning of file
                if len(r) == 0 or r[0][0] == '#':
                    pass

                # get information on layer names, thicknesses, and materials
                elif r[0] == 'Brad':
                    # build a list of layers
                    self.layers.append(r[1])
                    self.thicknesses.append(float(r[4]))
                    self.materials.append(r[2])

                # pulls the actual data from the file
                elif len(r)==3 and float(r[2]) != np.inf and np.isnan( float(r[2]) ) == False:
                    yA = float(r[0])
                    yB = float(r[1])
                    P = float(r[2])
                    # we have a valid coupling value.
                    if yA >= yB:
                        self.ys.append(yA)
                        self.xs.append(yB)
                        if P < 0 and P > -1e-10 or P == 0: 
                            self.zs.append(sys.float_info.epsilon) # if 0 set to a small value to make the log scale happy
                        else:
                            self.zs.append(P)
                    if yA <= yB:
                        self.ys2.append(yA)
                        self.xs2.append(yB)
                        if P < 0 and P > -1e-10 or P == 0:
                            self.zs2.append(sys.float_info.epsilon) # if 0 set to a small value to make the log scale happy
                        else:
                            self.zs2.append(P)
                    self.a.append(P)
                    # Data for profile
                    if yA == 0.5045:
                        yB1.append(yB)
                        P1.append(P)
                    # Data for profile
                    if yA == 3.592:
                        yB2.append(yB)
                        P2.append(P)
                    if yA == 3.5925917968749994:
                        yB3.append(yB)
                        P3.append(P)
                # track whether we have extremely large or small values
                    if P < 0:
                        self.neg_values.append(P)
                    if P > 10:
                        self.large_values.append(P)
                elif len(r) == 6 and r[0] == '**':
                    # we have a new layer combo.
                    self.layer_combos.append( (int(r[1]), r[2], int(r[3]), r[4], float(r[5]) ))

#------------------------------------------------------------------------------
    # Returns indices and values around (xval, yval) within +/- (xint, yint)
    # Also returns the distance between (xval, yval) and all values within the
    # interval
    def find_index(self, xint, yint, xval, yval): # since it searches the lower half first xval < yval
        ind_list_upper = []
        ind_list_lower = []
        dist = []
        for (i, (x, y, z)) in enumerate(zip(self.xs, self.ys, self.zs)):
            if abs(x - xval) < xint and abs(y - yval) < yint:
                ind_list_lower.append((i, x, y, z))
                dist.append(np.sqrt((x - xval)*(x - xval) + (y - yval)*(y - yval)))
        for (i, (x, y, z)) in enumerate(zip(self.xs2, self.ys2, self.zs2)):
            if abs(x - yval) < yint and abs(y - xval) < xint:
                ind_list_upper.append((i, x, y, z))
        return pd.DataFrame({'LowerVals': ind_list_lower, 
                             'UpperVals': ind_list_upper, 
                             'Distances': dist})

#------------------------------------------------------------------------------
    # Given an approximate value of (x, y) returns the coordinates that are closest
    def nearest_xy(self, xval, yval):
        subset = self.find_index(1, 1, xval, yval)
        min_ind = subset['Distances'].idxmin()
        #min_dist = subset[0][subset[2]==min(subset[2])]
        return subset['LowerVals'][min_ind][1], subset['LowerVals'][min_ind][2]

#------------------------------------------------------------------------------
    # Take a slice across the matrix. Given an emission point, yA (in ys or ys2),
    # return all absorption points, yB (xs or xs2) and coupling coefficients (z)
    def get_slice(self, yA):
        lower = zip(self.xs, self.ys, self.zs)
        upper = zip(self.xs2, self.ys2, self.zs2)
        emission_slice = [(x, z) for (x, y, z) in lower if y == yA]
        emission_slice += [(x, z) for (x, y, z) in upper if y == yA]
        return pd.DataFrame({'yB (Absorption)': [x[0] for x in emission_slice],
                             'z (Coupling)': [x[1] for x in emission_slice]})
    
#------------------------------------------------------------------------------
    def efficiency_matrix(self, matrix_name):
        """
        Writes the coupling efficiency between layers as a matrix file.
        Adapted from plot_matrix4.py file.
        """
        num_layers = len(self.layers)
        coupling_matrix = np.zeros((num_layers, num_layers))
        thicknesses = np.array(self.thicknesses)
        for combo in self.layer_combos:
             coupling_matrix[combo[0], combo[2]] = combo[4]
        
        with open(matrix_name, 'wb') as matrixfile:
           writer = csv.writer(matrixfile)
           print('Writing Coupling Matrix to ', matrix_name)
           writer.writerow(['Coupling Matrix'])
           writer.writerow([' '] + self.layers)
           writer.writerow(['thickness'] + thicknesses.tolist() )
           for i in range(num_layers):
              writer.writerow([self.layers[i]] + coupling_matrix[i,:].tolist())

#------------------------------------------------------------------------------
    def yA_efficiency(self, logfile_name): #TODO: move this to its own class or file
        """
        Coupling efficiencies (% of emissions at yA that get reabsobed 
        somewhere else in the device) are written in logfile_name
        """
        yA = []
        Eff = []
        BackLoss = []
        FrontLoss = []
        Err = []
        x = 0
        with open(logfile_name) as logfile:
            for line in logfile:
                if len(line) > 3 and line[0] != '*' and line[0] != '(' \
                and line[0] != '' and line[0] != '<' and line[0] != '.' and line[:3] != '   '\
                and line[0]!='l' and line[0]!='[' and line[0]!=' [':
                    print(len(line))
                    print(line[0])
                    # the file is written with spaces to separate info (could
                    # be changed to tab delimited #TODO)
                    # [1]: layer name, [5] & [6], simulation progress
                    # [-2]: yA coord in layer, [-1] Coupling efficiency
                    # [2024-01-03] file has been updated to include front and back losses
                    #              [1]: layer name, [5] & [6], simulation progress
                    #              [-5]: yA coord in layer, [-4] Coupling efficiency, [-3] Front loss, [-2] Back loss, [-1] Err (1-Eff-FL-BL)
                    split_line = line.split(' ')
                    try:
                        if split_line[-5][0] == 'y':
                            print('Line Number:', x)
                            print(line)
                            print(split_line[-5])
                            print(split_line[-4])
                            print(split_line[-3])
                            print(split_line[-2])
                            print(split_line[-1])
                            yA.append(float(split_line[-5].lstrip('yA=')))
                            Eff.append(float(split_line[-4].lstrip('Eff=')))
                            FrontLoss.append(float(split_line[-3].lstrip('FrontLoss=')))
                            BackLoss.append(float(split_line[-2].lstrip('BackLoss=')))
                            Err.append(float(split_line[-1].rstrip('/n').lstrip('Err=')))
                    except IndexError:
                        pass

                x += 1
        return pd.DataFrame({'yA (Emission)': yA,
                             'Efficiency': Eff,
                             'FrontLoss': FrontLoss,
                             'BackLoss': BackLoss,
                             'Err': Err})


