'''

Salvatore Maraniello. 12/Aug/2014

Module to define stiffness and mass properties of a beam given:
    a. the properties of one/more cross-section
    b. a method to derive how these change along the span

The output two are 3D array of size:
    Kbeam = (6,6,NumElems)
    Mbeam = (6,6,NumElems)
where NumElems are the number of elements in the beam.

'''

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
        K, M = csfun(*argslist)
        Kbeam[ii,:,:]=K
        Mbeam[ii,:,:]=M
    
    return Mbeam, Kbeam
    
  
#----------------------------------------------------------------------------- -
def extract_cross_section_method(cross_section_type):
    ''' returns the method to compute mass and stiffness properties of the beam
        cross section '''
    
    if   cross_section_type == 'isorect':
        csfun = design.isosec.rect #(l2, l3, material)
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
    
    print 'Testing cross-section extraction method'
    print extract_cross_section_method('isorect')
    print extract_cross_section_method('isoellip')
    print extract_cross_section_method('isocirc')
    print extract_cross_section_method('isohollowrect')
    print extract_cross_section_method('isohollowellip')
    print extract_cross_section_method('isohollowcirc')
    #extract_cross_section_method('give_error!')
    
    
    print 'Testing constant beam building'
    material='allum'
    NumElems = 3
    
    l2=0.2; l3=0.5;
    t2=0.1; t3=0.25; # full beam
    
    argslist=(l2,l3,t2,t3,material)
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
    
    
    
    
    

