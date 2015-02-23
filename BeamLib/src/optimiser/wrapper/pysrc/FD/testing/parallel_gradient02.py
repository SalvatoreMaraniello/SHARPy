'''
Created on 19 feb 2015

Compute numerical FD gradient in parallel using XBeamsolver class

Use the Pool class to distribute the jobs. Two optional codes can be built,
depending on how the XBeamSolver class object is passed.

    a. if the object is passed by deepcopy, only fwd_run needs to be deleted and
    rebuilt. this is extremely easy

    b. if the object is passed without deepcopy, the issue of not being able to 
    pickle the ctypes and the pointer fwd_run, is overcome by deleting all 
    pointers (ctypes and fwd_run) from the class copied and rebuilding them 
    inside each process.
    
'''

import numpy as np
import multiprocessing as mpr


codedir = '/home/sm6110/git/SHARPy/BeamLib/src/optimiser/wrapper'
import sys
sys.path.append(codedir)
sys.path.append(codedir+'/pysrc')

import shared
from xbcomponent import XBeamSolver
import lib.read

import time

import ctypes as ct
wrapso = ct.cdll.LoadLibrary(shared.wrapso_abspath)


method='b'




def perturb_a(xbpert,attr, ii):
    ''' Requires deepcopy of the class! '''
    
    xbpert.fwd_run=getattr(wrapso, "__opt_routine_MOD_opt_main") 

    setattr(xbpert, attr, 1.1*getattr(xbpert,attr) ) # 10% step change    
              
    xbpert._savedir='/home/sm6110/git/SHARPy_studies/GEBM/20141217_fdparallel/20150219_multipr/fdparallel/'
    xbpert.TestCase='fdnx%.2d.h5' %(ii)
    xbpert.execute()
    
    jv=0

    return jv


def perturb_b(xbpert,attr, ii):
    ''' Does not require a deepcopy '''
    
    xbpert.set_ctypes()
    xbpert.fwd_run=getattr(wrapso, "__opt_routine_MOD_opt_main") 
  
    setattr(xbpert, attr, 1.1*getattr(xbpert,attr) ) # give a 10% step change    
              
    xbpert._savedir='/home/sm6110/git/SHARPy_studies/GEBM/20141217_fdparallel/20150219_multipr/fdparallel/'
    xbpert.TestCase='fdnx%.2d.h5' %(ii)
    xbpert.execute()
    
    jv=0

    return jv




# ------------------------------------------------------------------------------
# reading here: the file is read from all the components. In an FD application 
# this will be used as perturbed
xborigin=lib.read.h5comp('/home/sm6110/git/SHARPy_studies/GEBM/20141217_fdparallel/20150219_multipr/NC9p3sol932NE10_dt0.00100_damp0.010_Fmax2.000_001.h5')
xborigin.PrintInfo=False

T0 = time.time()

# List of Scalar Design Variables
LxScalar = ['cs_l2','cs_l3','BeamLength1','NewmarkDamp','MinDelta']
NxScalar = len(LxScalar)
NxTot = NxScalar # + Vector quantities

# create pool
PROCESSES = 5
print 'Creating pool with %d processes\n' % PROCESSES
pool = mpr.Pool(processes=PROCESSES)
# Launch processes in parallel
results = []
xborigin.fwd_run=None


if method=='a':
    # parallel code with async
    for ii in range(NxScalar):
        results.append( pool.apply_async(perturb_a,args=(xborigin.copy(),LxScalar[ii],ii) ))   
    # retrieve the output  
    output = [p.get() for p in results]
    print(output)

if method=='b':
    # add ctypes to make it realistic
    xborigin.set_ctypes()
    # delete c_types
    xborigin.del_ctypes(printinfo=True)
    # parallel code with async
    for ii in range(NxScalar):
        results.append( pool.apply_async(perturb_b,args=(xborigin,LxScalar[ii],ii) ))   
    # retrieve the output  
    output = [p.get() for p in results]
    print(output)


'''
#serial code
results = []
xbpert=xborigin
for ii in range(NxScalar):
    xbpert=xborigin
    #xbpert=xborigin.copy() ## <-- this would mess up the code, as some methods 
    # don't seen to be copied. It is unnecessary as xbpert is passed into function, 
    # and changes to it are not see outside.
    results.append( perturb(xbpert,LxScalar[ii],ii) )
'''

# and build the jacobian
Ttot = time.time() - T0
print '------------------------ Jacobian computed in %f seconds' %(Ttot)

