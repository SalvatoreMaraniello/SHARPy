'''

Salvatore Maraniello. 12/Aug/2014

Module to define elements stiffness and mass properties of a beam given:
    a. the properties of one/more cross-section
    b. a method to derive how these change along the span

The output two are 3D array of size:
    Kbeam = (6,6,NumElems)
    Mbeam = (6,6,NumElems)
where NumElems are the number of elements in the beam.

'''
import numpy as np
import scipy.interpolate as scint
import matplotlib.pylab as mp

import sys
sys.path.append('..')
import design.isosec


def constant(NumElems,cross_section_type,argslist):
    ''' 
    define beam with constant properties:
        NumElems: number of elements in the beam
        cross_section_type: see extract_cross_section_method
        argslist: arguments for the cross-section calculation 
                  (see design.isosec)
    '''
        
    # define method:
    csfun = extract_cross_section_method(cross_section_type)
    
    # evaluate cross-section
    M, K = csfun(*argslist)
    
    # define output
    Kbeam = np.zeros((NumElems,6,6),dtype=float,order='F')
    Mbeam = np.zeros((NumElems,6,6),dtype=float,order='F')
    
    # allocate variables
    for ii in range( NumElems ): 
        M, K = csfun(*argslist)
        Kbeam[ii,:,:]=K
        Mbeam[ii,:,:]=M
    
    return Mbeam, Kbeam
    
  
def spline(NumElems,cross_section_type,ArgsControlList,ControlElems,SplineOrder):
    '''
    Define beam with properties varying along the span.
    The properties are defined on the elements defined in ControlElems 
    and the global reconstruction is done according to a spline of order 
    SplineOrder.
    
    Input:
    NumElems: number of elements in the beam
    cross_section_type: see extract_cross_section_method
    ArgsControlList = [argslist]: each item of ArgList is a tuple with the input
    argument for building the cross-sectional properties (see design.isosec).
    ControlElems: integer list with the number of the elements to be used as
    control points
    SplineOrder: order of spline as per scint.interpolate.spline
    
    Output:
    Mbeam, Kbeam: 3d arrays of shape (ii,6,6), containing the mass and stiffness
    matrices of the ii-th element
    ArgsBeamList: list containing the arguments to build the ii-th cross-section
    using the cross-section_type method.
    
    Remark:
    all cross-sectional variables are reconstructed using the spline scheme 
    defined, even those which are not chosen to be design variables in the 
    optimisation problem. However, at the level of the optimiser is it possible
    to keep fixed the value of some of them, such that these will not change
    during the optimisation
    '''
    
    NC = len(ControlElems)
    # check input:
    if len(ArgsControlList) != NC:
        raise NameError('ArgsControl should have the same length as ControlElems!')
    for ii in range(NC-1):
        if ArgsControlList[ii][-1] != ArgsControlList[-1][-1]:
            raise NameError('Not all ControlCrossSections have the same material!') 
     
     
    if ControlElems[0] != 0 or ControlElems[-1] != NumElems-1:
        raise NameError('Root and Tip cross sections must be included in the ControlElems')
    

    # get number of numerical input:
    if cross_section_type == 'isocirc':
        Nargs=1        
    elif cross_section_type == 'isorect' or cross_section_type == 'isoellip' or cross_section_type == 'isohollowcirc':
        Nargs=2 #(l2, l3, material) - (r,t,material)
    elif cross_section_type == 'isorect_fact_torsion':
        Nargs=3
    elif cross_section_type == 'isohollowrect' or cross_section_type == 'isohollowellip':
        Nargs=4
    else:
        raise NameError('Cross Section Type "%s" not found!' %cross_section_type )    


    # allocate list of input for all cross sections
    ArgsBeamArray = np.zeros((NumElems,Nargs)) # [None]*NumElems
    
    # define design of each element based on control elements
    xc = np.array(ControlElems)
    xv = np.arange(NumElems)
    
    # spline interpolation
    for jj in range(Nargs):
        yc=np.zeros((NC))
        # extract values
        for ii in range(NC):
            yc[ii]=ArgsControlList[ii][jj]    
        ArgsBeamArray[:,jj] = scint.interpolate.spline(xc,yc,xv,order=SplineOrder)
    ArgsBeamList = ArgsBeamArray.tolist()    
    
    # add material information
    for ii in range(NumElems):
        ArgsBeamList[ii].append(ArgsControlList[0][-1])
    
    #-------------------------------------------------------------------- check:
    for jj in range(Nargs):
        for ii in range(NC):
            mp.plot(xc[ii],ArgsControlList[ii][jj],'kd')
        mp.plot(xv,ArgsBeamArray[:,jj],'r+--')
        mp.show()
    #---------------------------------------------------------------------------

    # define method
    csfun = extract_cross_section_method(cross_section_type)
   
    # define output
    Kbeam = np.zeros((NumElems,6,6),dtype=float,order='F')
    Mbeam = np.zeros((NumElems,6,6),dtype=float,order='F')    

    # build output
    for ii in range(NumElems):
        Mbeam[ii,:,:], Kbeam[ii,:,:] = csfun(*ArgsBeamList[ii])
        
    return Mbeam, Kbeam, ArgsBeamList
 

def constant_hardcoded(NumElems,Mroot,Kroot):
    ''' 
    Define beam with constant properties given root Mass and stiffness matrices
    The method can be used to reproduce existing results.
    
    Input:
        NumElems: number of elements in the beam
        Mroot, Kroot: mass and stiffness matrices of one element      
    '''
    
    # define output
    Kbeam = np.zeros((NumElems,6,6),dtype=float,order='F')
    Mbeam = np.zeros((NumElems,6,6),dtype=float,order='F')
    
    # allocate variables
    for ii in range( NumElems ): 
        Kbeam[ii,:,:]=Kroot
        Mbeam[ii,:,:]=Mroot
    
    return Mbeam, Kbeam    
        
   
def spline_derivative(NumElems,ControlElems,SplineOrder):
    '''
    Return the derivative of the spline curves, evaluated on each element, in
    respect to the design variable.
    
    Being the spline reconstruction based on a linear relation, the sought 
    derivatives are the spline bases function themselves.
    
    '''
    NC = len(ControlElems)
    Base = np.zeros((NumElems,NC))
    
    xc = np.array(ControlElems)
    x  = np.arange(NumElems)
    
    for ii in range(NC):
        yc=np.zeros((NC))
        yc[ii]=1.0
        Base[:,ii] = scint.interpolate.spline(xc,yc,x,order=SplineOrder)
    
    return Base
    
    
def extract_cross_section_method(cross_section_type):
    ''' returns the method to compute mass and stiffness properties of the beam
        cross section '''
    
    if   cross_section_type == 'isorect':
        csfun = design.isosec.rect #(l2, l3, material)
    elif cross_section_type == 'isorect_fact_torsion':
        csfun = design.isosec.rect_fact_torsion
    elif cross_section_type == 'isoellip':
        csfun = design.isosec.ellip #(l2, l3, material)
    elif cross_section_type == 'isocirc':
        csfun = design.isosec.circ #(r, material)
    elif cross_section_type == 'isohollowrect':
        csfun = design.isosec.hollowrect #(l2, l3, t2, t3, material) 
    elif cross_section_type == 'isohollowellip':
        csfun = design.isosec.hollowellip #(l2, l3, t2, t3, material)
    elif cross_section_type == 'isohollowcirc':
        csfun = design.isosec.hollowcirc #(r, t, material)
    else:
        raise NameError('Cross Section Type "%s" not found!' %cross_section_type )  
    
    return csfun    
    
    
'''--------------------------------------------------------------------------''' 


if __name__=='__main__':
    
    import numpy as np
    
    
    #------------------------------------ Cross Sectional Method extraction Test
    print 'Testing cross-section extraction method'
    print extract_cross_section_method('isorect')
    print extract_cross_section_method('isoellip')
    print extract_cross_section_method('isocirc')
    print extract_cross_section_method('isohollowrect')
    print extract_cross_section_method('isohollowellip')
    print extract_cross_section_method('isohollowcirc')
    #extract_cross_section_method('give_error!')
    
    
    #--------------------------------------------- Constant Properties Beam Test
    print 'Testing constant beam building'
    material='titnA'
    NumElems = 3
    
    l2=0.2; l3=0.5;
    t2=0.1; t3=0.25; # full beam
    
    argslist=[l2,l3,t2,t3,material]
    M, K = design.isosec.hollowrect(*argslist)
    Mbeam, Kbeam = constant(NumElems,'isohollowrect',argslist)
    
    print 'Expected:'
    print 'diag(K)', np.diag(K)
    print 'diag(M)', np.diag(M)
    
    print 'Span Variation'
    for ii in range(NumElems):
        print 'element %i' %ii
        print 'diag(K)', np.diag(Kbeam[ii,:,:])
        print 'diag(M)', np.diag(Mbeam[ii,:,:]) 
    
    
    #------------------------------------------- Spline Beam Reconstruction Test
    print 'Testing spline beam building'
    NumElems=10
    ControlElems = [0, 2, 6, 10-1]
    SplineOrder = 1
    
    # constant beam by spline
    ArgsControlList = [argslist]*len(ControlElems)  ##### warning: this creates links to same object!!! 
    Mbeam,Kbeam, ArgsList = spline(NumElems,'isohollowrect',ArgsControlList,ControlElems,SplineOrder)
    for ii in range(NumElems):
        print 'element %i' %ii
        print 'diag(K)', np.diag(Kbeam[ii,:,:])
        print 'diag(M)', np.diag(Mbeam[ii,:,:])         

    """
    # Variation of the JJarg-th argument via power law
    JJarg = 1
    ArgsControlList=[]
    for ii in range(len(ControlElems)):
        ArgsControlList.append(list(argslist))  ##### warning : the list command is not to convert into list (argslist is already a list) but to deepcopy!!!
      
    for ii in range(len(ArgsControlList)):
        ArgsControlList[ii][JJarg] = ArgsControlList[0][JJarg] * (1.0 + 100.0*float(ii)**5)
    print ArgsControlList

        
    # Variable Spline
    Mbeam,Kbeam, ArgsList = spline(NumElems,'isohollowrect',ArgsControlList,ControlElems,SplineOrder)
    for ii in range(NumElems):
        print 'element %i' %ii
        print 'diag(K)', np.diag(Kbeam[ii,:,:])
        print 'diag(M)', np.diag(Mbeam[ii,:,:])
    for ii in range(NumElems):    
        mp.plot(float(ii),Mbeam[ii,0,0],'o')
    mp.title('Cross Sectional Mass')
    mp.show() 
    for ii in range(NumElems):    
        mp.plot(float(ii),Kbeam[ii,0,0],'o')
    mp.title('Cross Sectional EA')
    mp.show() 
       
    """
    #------------------------------------------------------- Splines derivatives
    Base = spline_derivative(NumElems,ControlElems,SplineOrder)
    
    for ii in range(len(Base[0,:])):
        mp.plot(np.arange(NumElems),Base[:,ii],'o-')
    mp.title('Analytically Computed Basis (derivatives)')
    mp.show()    
    
    

    #-------------------------------------------------------- Constant Hardcoded
    Mroot = np.random.random((6,6))
    Kroot = np.random.random((6,6))
    NumElems = 3
    Mbeam, Kbeam = constant_hardcoded(NumElems,Mroot,Kroot)
    for ii in range(len(Mbeam[:,0,0])):
        print Mbeam[ii,:,:]-Mroot
        print Kbeam[ii,:,:]-Kroot
    
        

    

