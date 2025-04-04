This is a collection of Python modules, originally put together by Matt Wilkins.

The modules should work under both Windows and Linux. They were developed using Python 2.7
which is more widely used than the newer Python 3.x. Expect some minor incompatibilities if using Python 3.x.

The tools may require any of these packages:

Python 2.7	Core Python interpreter and base libraries
numpy		Numerical library (matrices, linear algebra, etc)
scipy		Scientific math library (ODEs, numerical integration, optimization, root solving, etc.)
matplotlib	High quality plotting, capabilities similar to Matlab

(and less frequently used)
pyodbc		Read/write to SQL databases
h5py		Read/write HDF5 data (i.e. Sentaurus TDR file format)
PySide		bindings to Qt application programming framework

Transfer Matrix Method calculations can be done with the TMM library at
https://pypi.python.org/pypi/tmm.

Cython is a variant of Python where code is converted to C and pre-compiled rather than run under an interpreter.
This can provide large performance improvements in some cases.

Also, ipython is a more advanced python shell with graphical, text and web-based interfaces.

Any of the files below can be imported to make the classes and functions avaialble to your script.



iv_analysis.py:
_______________

process arrays of current and voltage.  Fit a bicubic spline to the area around maximum power point,
and extract maximum power and other parameters.


dfise.py:
_________
Reads in a file in DF-ISE format.  Sentaurus CurrentPlot (.plt) files are in this format, and so are
some XY-data files exported from Tecplot.

Makes all available datasets available in a dict indexed by dataset name.  Also includes functions to
calculate EQE and integrate EQE with ASTM spectra to get Jsc.

requires: ASTMG173.py

example:

import dfise  # loads the module dfise.py.
import matplotlib.pyplot as plt


# create an instance of class dfise from the given file.
d = dfise.dfise('n1781_des.plt')

# print all available dataset names
print 'Available datasets:'
for i in d.data.keys():
    print i

# plot I-V curve (assuming this file contains I-V data).
plt.plot(d.data['anode OuterVoltage'], d.data['cathode TotalCurrent'])

# label axes
plt.xlabel('Bias Voltage [$\mathrm{V}$]')
plt.ylabel('Current [$\mathrm{A}$]')
plt.show()



epi_cmd.py:
___________

Reads in a pre-processed epi command file, ie. pp132_epi.csv.  Makes a list of all layers in the structure
with name, thickness, material and doping, plus top and bottom y-coordinate of each layer.


epi_cmd_LC.pyx:
___________

More advanced version of epi_cmd.py used only for luminescent coupling calculations.   Implemented as
Cython code.  
- Makes a list of all layers
- Reads in MatPar files for each layer and extracts material parameters (n,k data, radiative recombination coefficient,
band gap, etc).  
- Reads in .msh file (TDR format) and finds all unique y-coordinates in each layer. 
- Also contains some code for luminescent coupling calculations.

requires: sdevice_cmd.py, tmm_core_mw.pyx


sdevice_cmd.py:
_______________

Parse files in MatPar .par format and extract parameters.

tdr.py:
_______

Read in a TDR file and make all datasets available indexed by name.


svisual_1d_csv.py:
__________________

Read in a csv file output by svisual containing a 1D cut of Sentaurus simulation results.


read_qe.py:
___________

Read in a .log file produced by Oriel IQE-200 machine.  All columns are avaialble indexed by name.
Can also calculate integrated Jsc based on ASTM G173 spectra.

requires: ASTMG173.py


Brad.py and Brad_PR.py:
__________

Calculates radiative recombination rate coefficient from nk data.
..PR.py accounts for photon recycling, when calculating Brad.


Curve_compare.py:
_________________

Functions to compare 2 datasets. Useful for fitting to measured data.


k_ndoping_dependent.py:
_______________________

Calculates the n-doping dependent extinction coefficient data assuming the Anderson model. Useful for close to band edge + small bandgap materials.


py2prepper.py:
______________

Sentaurus preprocessor runs this file, which writes the pyCollector tool cmd file. Required to amass database.


StitchMesh.py:
______________

v2 of importing externaly calculated carrier generation rates.
Deprecated.
For converting TDR to grd and dfsi files, then writing generation rates, then converting back to TDR files.

StitchMesh_v3.py:
_________________

v3 of importing externaly calculated carrier generation rates.
This will directly write the generation rate within the TDR file.




