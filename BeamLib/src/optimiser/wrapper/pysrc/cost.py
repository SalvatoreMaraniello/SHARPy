'''
# Salvatore Maraniello 7 Aug 2014

This module contains methods to evaluate the cost functions related to the 
beam solver.

The routines computing the global cost and the constraints cost are not 
extracted. These can be rebuilt in python environment.

Remark: there is no particular advantage to writen further 

'''

import ctypes as ct
import numpy as np

import shared
import beamvar


# Load Dynamic Library
xb = ct.cdll.LoadLibrary(shared.wrapso_abspath)

# ------------------------------------------------------------ Fortran routines
# The following routine are not extracted - not worth to built an interface
#fcost  =xb.__opt_cost_MOD_cost_global
#fconstr=xb.__opt_cost_MOD_cost_constraints

# Total structural mass
f_total_mass = xb.__opt_cost_MOD_cost_total_mass_wrap

# Nodal Dfisplacement (absolute value)
f_node_disp  = xb.__opt_cost_MOD_cost_node_disp_wrap













