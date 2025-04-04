#!/usr1/bin/python
import csv
import string
import re
import math

e_charge = 1.602e-19 # C
h = 6.626e-34      # J.s
c0 = 299792458.0   # m/s
kb = 1.3806488e-23 # J/K

# Class to read in a .cmd file (pre-processed epi command file)
#    epi

class epilayer:
    def __init__(self, name, material, materialfile, thickness, doping, molefraction,ext, ytop):
        self.name = name
        # Used to refine labels in plots.
        self.label = name
        self.offset = 0.

        self.material = material
        if self.material in ['GaAsOx', 'TiOx', 'MgF', 'SiN', 'SiO2', 'Si3N4']:
           self.material_type = 'Dielectric'	
        else:
           self.material_type = 'Semiconductor'

        self.materialfile = materialfile
        self.thickness = float(thickness)
        if len(doping) > 0:
            self.doping = float(doping) # > 0: ptype, <0: ntype
        else:
            self.doping = 0.0
        if len(molefraction) > 0:
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



# class containing a list of all layers in the epi structure.
class epifile:
    def __init__(self, filename, layers=''):
        with open(filename, 'rb') as csvfile:
            self.r = csv.reader(csvfile, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL, skipinitialspace=True)
            self.filename = filename
            self.layers = []
            self.ytop = 0.0
            self.ybot = 0.0

            for row in self.r:
                if len(row) > 0:
                    if len(row[0].lstrip('\t')) >2:
                        if  row[0].lstrip('\t')[0] != '$':
                            self.add_layer(row)
                        elif row[0].lstrip('\t').split(' ')[0] == '$repeat':
                            N = int(row[0].split(' ')[1])
                            row = self.r.next()
                            replayers = []
                            while row[0].lstrip('\t') != '$end':
                                    replayer = epilayer(row[0].lstrip('\t'),row[1],row[2],row[3], row[4], row[5], row[6],0)
                                    replayers.append(replayer)
                                    row = self.r.next()
                            for i in range(N):
                                for rl in replayers:
                                    self.add_layer([re.sub("\$i",str(i),rl.name),rl.material,rl.materialfile,str(rl.thickness),str(rl.doping), '(' + str(rl.xmole) +' '+str(rl.ymole)+')',rl.ext,self.ybot])

        with open(layers, 'rb') as layersfile:
              reader = csv.reader(layersfile)
              while True:
                      try:
                             row = reader.next()
                             if len(row) == 3:
                                   name = row[0]
                                   label = row[1]
                                   offset = float(row[2])
                                   print( name, label, offset)
                                   for l in self.layers:
                                          if l.name == name:
                                                 print( 'matched', l.name)
                                                 l.label = label
                                                 l.offset = offset
                      except StopIteration:
                            break
                                
    def add_layer(self, row):
                layer = epilayer(row[0].lstrip('\t'), row[1], row[2], row[3], row[4], row[5], row[6],self.ybot)
                self.ybot = layer.ybot
                self.layers.append(layer)


