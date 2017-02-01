SHARPy
======

Simulation of High-Aspect-Ratio Planes in Python (SHARPy).

Nonlinear and consitently linearized models of very flexible aircraft dynamics
for design, analysis, and control law synthesis.



### Install ###

1. Make sure Anaconda is installed (see below)
2. From terminal, run 'python install.py'. See comments in install.py for more
details and installation options.



### Dependencies ###

- Anaconda3-2.3.0-Linux-x86-64 or higher: this contains all the python packages
required to run SHARPy (see detailed list below)
- gfortran compiler
- C++ libraries Eigen and Boost. If not availale, these can be automatically 
downloaded using install.py

- Details:
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

During the installation, the fortran and C/C++ projects Makefiles are 
reconfigured to reflect the locations of these libraries on your system.



### Run ###

The first thing you should do is run some python unittests
which can be found in the test/ directory of PyAero, PyBeam and PyCoupled
solver scripts.



### Installation Troubleshot ###

- see comments in install.py



### Uninstall ###
- Remove the SHARPy directory from your system. 

