SHARPy
======

Simulation of High-Aspect-Ratio Planes in Python (SHARPy).

Nonlinear and dynamically linearized models of very flexible aircraft dynamics
for design, analysis, and control law synthesis.

Install
-------

1. Clone this git repository onto your local machine at location
[SharPyProjectDir].
2. Edit the file [SharPyProjectDir]/src/Main/SharPySettings.py so that

	SharPyProjectDir = [SharPyProjectDir]

Note: The BeamLib shared library, BeamLib.so, is compiled automatically by
unit tests in the PyBeam package.
However, the UVLMLib.so library must be compiled manually at present;
to do this import the UVLMLib project in eclipse (with the C/C++ SDK installed)
by navigating to 'File->Import->General->Existing Projects into Workspace',
and set [SharPyProjectDir]/UVLMLib as the root
directory.
Then build a Debug configuration of the project.

### Dependencies ###

* Python 3.2
* numpy (1.7.x+)
* scipy
* h5py
* slycot (routines for preditive controller design)
* muAO-MPC (MPC controller design and C-code generation)
* Eigen (C++ libraries)
* Boost (C++ libraries)
* LAPACK (Fortran libraries)
* BLAS (Fortran libraries)

The fortran and C/C++ projects must be reconfigured to reflect the locations
of these libraries on your system.

Run
---

The first thing you should do is run some python unittests
which can be found in the test/ directory of PyAero, PyBeam and PyCoupled
solver scripts.


### Installation Troubleshot ###

- Anaconda3-2.3.0-Linux-x86-64 (https://repo.continuum.io/archive/) contains all the required python packages

- ensure paths to all python files are set-up correctly. This is achieved adding this path:
    ./SHARPy/src/Main/
  to sys.path list.

- The BeamLib.so can be compiled manually:
    1. go to:
        /home/sm6110/git/SHARPy/BeamLib/src/wrapper
    2. type:
        make all
  
- Run the shell script:
    /home/sm6110/git/SHARPy/BeamLib/src/wrapper/f2py/runf2py.sh
  The script relies on the F2PY, which is installed with numpy. If multiple python installations are icluded, ensure 
  that the correct script is called by modifying the call to 'f2py3.2' to 'f2py3' or the full path of the script.

- Eigen and boost installation:
    1. You can try to download from: 
        http://eigen.tuxfamily.org/index.php?title=Main_Page#Download 
        http://sourceforge.net/projects/boost/?source=typ_redirect
    and follow the instructions. Usually no compilation is required.
    2. In the eclipse project/makefile include the path to X, where X is the folder containing Eigen/Dense. If no 
    compilation was required, these are the folders you extracted
