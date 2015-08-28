SHARPy
======

Simulation of High-Aspect-Ratio Planes in Python (SHARPy).

Nonlinear and dynamically linearized models of very flexible aircraft dynamics
for design, analysis, and control law synthesis.

Install
-------

Clone this git repository onto your local machine at location
[SharPyProjectDir].
Edit the file [SharPyProjectDir]/src/Main/SharPySettings.py so that

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
  The script relies on the F2PY executable f2py3. This is shipped with numpy but may require a manual installation. On
  Ubuntu the executable comes with the standard numpy package from the Ubuntu disctribution. If a manual installation is 
  required:
    - find location where numpy is installed
    - go to the f2py directory
    - python3 setup.py install 
