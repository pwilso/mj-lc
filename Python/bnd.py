#!/usr1/bin/python
import csv
import numpy as np
import glob
import matplotlib.pyplot as plt
import string

# Class to read in a .bnd file produced by Sentaurus epi

class region:
    def __init__(self, name):
        self.name = name
        self.elements = []
        self.material = ''


class bnd:
    def __init__(self, filename):
        print 'Reading DF-ISE file', filename
        with open(filename, 'rb') as csvfile:
            self.r = csv.reader(csvfile, delimiter=' ', quotechar='"', quoting=csv.QUOTE_MINIMAL, skipinitialspace=True)
            row = self.r.next()
            assert row[0] == "DF-ISE"
            self.filename = filename
            if len(row) > 0:
                for row in self.r:
                    if len(row) > 0:
                        if row[0] == 'Info':
                            self.readInfo()
                        if row[0] == 'Vertices':
                            self.readVertices()
                        if row[0] == 'Edges':
                            self.readEdges()
                        if row[0] == 'Elements':
                            self.readElements()
                        if row[0] == 'Region':
                            self.readRegion(row[1].strip('()"'))

    def readInfo(self):
#        print 'starting readInfo'
        self.regions = dict([])
        self.materials = []
        for row in self.r:
            if row[0] == '}':
                break
            elif row[0] == 'version':
                self.dfise_version = row[2]
            elif row[0] == 'type':
                self.dfise_type = row[2]
                assert self.dfise_type == 'boundary'
            elif row[0] == 'dimension':
                self.dimension = int(row[2])
            elif row[0] == 'nb_vertices':
                self.vertices = int(row[2])
            elif row[0] == 'regions':
                for i in row[3:-1]:
                    self.regions[i] = region(i)
            elif row[0] == 'materials':
                for i in row[3:-1]:
                    self.materials.append(i)
            
    def readVertices(self):
        self.vertices = []
        row = self.r.next()
        while len(row) > 0 and row[0] != '}':
            self.vertices.append([float(row[0]), float(row[1])])
            row = self.r.next()

    def readEdges(self):
        self.edges = []
        row = self.r.next()
        while len(row) ==2:
            self.edges.append([int(row[0]), int(row[1])])
            row = self.r.next()
            
    def readElements(self):
        self.elements = []
        row = self.r.next()
        while len(row) > 1 and row[0] != '}':
            self.elements.append([int(row[2]), int(row[3]),int(row[4]), int(row[5])])
            row = self.r.next()

    def readRegion(self, rname):
        row = self.r.next()
        while row[0] != '}':
            if row[0] == 'material':
                self.regions[rname].material = row[2]                  
            row = self.r.next()

