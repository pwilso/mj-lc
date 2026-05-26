Installation
===================================

The project can be run as a standalone project or in conjunction with a Synopsys Sentaurus simulation. Python is required in both cases.

Compatibility with different Python versions and different Sentaurus versions is important. We have successfully implemented the code with Python v3.6.8 and with Python v3.13.5. Significant changes affecting installation and package management occured between these two versions so we provide two versions of the code with minor changes to each. Compatibility between different Sentaurus versions can also impact the way the code is implemented. We provide some troubleshooting advice for these concerns. The code has been successfully implemented with Sentaurus versions S-2021.06 and U-2022.12-SP1.

Python Packages
----------------------------------

The Python environment requires the Python.h header files in order to compile the Cython extensions. Ensure that these are installed on your system (using the command *sudo apt install python3-dev* or similar) or by building your Python environment from source.

We recommend the user create a virtual environment and install the necessary python packages via the requirements.txt file. For example:

.. code-block:: console

  mkdir venvs
  python3 -m venv ./venvs/LC_env
  source ./venvs/LC_env/bin/activate

  pip install -r requirements.txt

Compiling with Cython
----------------------------------

The *epi_cmd_LC.pyx* and the *tmm_core_mw.pyx* files must be compiled with Cython before the program can be successfully run.

.. code-block:: console

  cythonize /path/to/Python/tmm/tmm_core_mw.pyx
  cythonize /path/to/Python/epi_cmd_LC.pyx

Integration with Sentuarus
----------------------------------

Since this code was originally written for exclusive use with Sentaurus even running the program from the terminal or a Python IDE requires some initial Sentaurus integration as described in :doc:`basic_usage`.

The Sentaurus mesh output (usually named *n##_msh.tdr*) may be generated as version 16 or version 17 depending on the Sentaurus version. Version 16 is required in order for the file format to be recognized as HDF5 and parsable by the h5py functions. To check which version has been produced you can run the following in the terminal:

.. code-block:: console

  file /path/to/n##_msh.tdr

The output should return *Hierarchical Data Format (version 5) data*. If it does not then you may need to convert the file to the appropriate version. From within Sentaurus open a terminal at the mesh node location. The Sentaurus Data Explorer commands can be used to inspect and act on the tdr files. To display information about the tdr file and convert the version:

.. code-block:: console

  tdx --info n##_mesh.tdr
  tdx --tdr2tdr -v16 n##_msh.tdr n##_msh.tdr

To have the program run this command automatically, the following line can be added to the file *gtooldb.tcl*, which has commands for controlling the prologue and epilogue commands for the various Sentuarus tools.

.. code-block:: tcl

  set WB_tool(sde, epilogue) {exec tdx --tdr2tdr -v16 n@node@_msh.tdr n@node@_msh.tdr}
