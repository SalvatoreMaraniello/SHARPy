'''

Salvatore Maraniello. 18/Aug/2014

Module to define nodal pretwist properties of a beam given:
    a. a method to derive how these change along the span

The module is analogous to design.beamelems, but operates on the beam nodes.

The only output currently supported is a 1D array (PhiNodes in Fortran code) 
containing the pretwist at each node of the beam

    pretwist = (NumNodes)
    
Remark:
    - arguments in input to the methods constant and spline are kept general and
    written as argslist in case future developments may require new input/output
    to be defined.

'''

import numpy as np
import scipy.interpolate as scint
import matplotlib.pylab as mp

import sys
sys.path.append('..')
import beamelem



def constant(NumNodes,phi):
    ''' 
    define constant nodal properties along the beam:
        NumNodes: number of elements in the beam
        phi: pre-twist value to assign to all the nodes.
    '''
    
    # define output
    PreTwist = phi * np.ones((NumNodes),dtype=float,order='F')

    return PreTwist
    
  
def spline(NumNodes,PhiControl,ControlNodes,SplineOrder):
    '''
    Define beam nodal with properties varying along the span.
    The properties are defined on the nodes defined in ControlNodes 
    and the global reconstruction is done according to a spline of order 
    SplineOrder.
    
    NumNodes: number of nodes in the beam
    PhiControl: pre-twist value at the control nodes
    ControlNodes: integer list with the number of the nodes to be used as
    control points
    SplineOrder: order of spline as per scint.interpolate.spline
    
    Remark:
    all cross-sectional variables are reconstructed using the spline scheme 
    defined, even those which are not chosen to be design variables in the 
    optimisation problem. However, at the level of the optimiser is it possible
    to keep fixed the value of some of them, such that these will not change
    during the optimisation
    '''
    
    NC = len(ControlNodes)
    # check input:
    if len(PhiControl) != NC:
        raise NameError('PhiControl should have the same length as ControlNodes!')

    if ControlNodes[0] != 0 or ControlNodes[-1] != NumNodes-1:
        raise NameError('Root and Tip nodes must be included in the ControlNodes')
    
    # define design of each element based on control elements
    xv = np.arange(NumNodes)
    xc = np.array(ControlNodes)
    
    # reconstruct pre-twist distribution
    yc = PhiControl    
    PreTwist = scint.interpolate.spline(xc,yc,xv,order=SplineOrder)
    
    
    #-------------------------------------------------------------------- check:
    mp.plot(xc,yc,'kd')
    mp.plot(xv,PreTwist,'r+--')
    mp.show()
    #---------------------------------------------------------------------------
   
    return PreTwist
 
    

def spline_derivative(NumNodes,ControlNodes,SplineOrder):
    '''
    Returns the derivative of the spline curves, evaluated on each node, in
    respect to the design variable.
    
    Being the spline reconstruction based on a linear relation, the sought 
    derivatives are the spline bases function themselves.
    
    '''
    
    Base = beamelem.spline_derivative(NumNodes,ControlNodes,SplineOrder)
    
    return Base
    

'''--------------------------------------------------------------------------''' 

if __name__=='__main__':
    
    import numpy as np
    
    #--------------------------------------------- Constant Properties Beam Test
    print 'Testing constant beam building'
    NumNodes = 5
    phi = 0.5
    
    PreTwist = constant(NumNodes,phi)
    print 'Pretwist:', PreTwist
    
    #------------------------------------------- Spline Beam Reconstruction Test
    print 'Testing spline beam building'
    NumNodes=10
    ControlNodes = [0, 2, 6, 10-1]
    SplineOrder = 3
    
    # constant beam by spline
    PhiControl = phi*np.ones(len(ControlNodes))
    PreTwist = spline(NumNodes,PhiControl,ControlNodes,SplineOrder)
    print 'Pretwist:', PreTwist 

    # variable span properties
    PhiControl = np.zeros(len(ControlNodes))
    for ii in range(len(ControlNodes)):
        PhiControl[ii] = phi * (1.0 + 100.0*float(ii)**5)
    print PhiControl

        
    # Variable Spline
    PreTwist = spline(NumNodes,PhiControl,ControlNodes,SplineOrder)
    print 'Pretwist:', PreTwist
    
    mp.plot(ControlNodes,PhiControl,'o')
    mp.plot(np.arange(NumNodes),PreTwist,'r')
    mp.title('Pretwist')
    mp.show() 
       
    
    #------------------------------------------------------- Splines derivatives
    Base = spline_derivative(NumNodes,ControlNodes,SplineOrder)
    for ii in range(len(Base[0,:])):
        mp.plot(np.arange(NumNodes),Base[:,ii],'o-')
    mp.title('Analytically Computed Basis (derivatives)')
    mp.show()    
    
    

    
