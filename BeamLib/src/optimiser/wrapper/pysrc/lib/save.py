'''
Created on 20 Aug 2014

@author: sm6110

Tools for saving. 

List of routines:
  - h5sol (superseded): saves results from a simulation into a HDF5 binary data 
    format.
  - h5comp: saves all attributes of the class XBeamSolver into a hdf5 file. the 
    routine can be called before/after execution. If an attribute is not found
    or has no value assigned to it, a 'not found' or 'no value' string will be
    assigned.
  - conditional_saving: utility routine to check if an attribute is present or
    not and, if yes, to check if it has a value assigned or not.

'''

import h5py


def h5sol(PosIni,PsiIni,PosDef,PsiDef,filename='./sol.h5'):
    ''' 
    Saves nodes displacements and rotations in h5 format.
    
    Pos,Psi: position and rotation arrays in initial (Ini) or deformed (Def)
    configuration
    
    
    Remark: the format differs from the standard format of the Fortran code
    
    Fortran code format:
    Num.Elem   Num.Node.local   Position x,y,z   Rotation x,y,z
    
    xml format:
    Position x,y,z   Rotation x,y,z
    (where row ii-th corresponds to the ii-th node in the global ordering)
    '''

    hdfile = h5py.File(filename,'w') 
    
    hdfile['PosIni']=PosIni
    hdfile['PsiIni']=PsiIni
    hdfile['PosDef']=PosDef
    hdfile['PsiDef']=PsiDef
    
    hdfile.close()    
    
    return None
   
   
def h5comp(XBinst,filename='./fwd_model.h5'):
    ''' 
    Saves All variables of an instance of class XBeamSolver into a hdf5 file.
    
    For a details of the variable see: xbcomponent.XBeamSolver class.
    
    Last update: 21 Aug 2014
    
    Input:
    - XBinst: instance of XBeamSolver class
    - filename: saving name
    
    Remarks:
    the function can be called to save the component before or after the 
    solution is run. 
    '''   
    
    hdfile = h5py.File(filename,'w')
    
    #-------------------------------------------------------- input with default
    
    # options:
    hdfile['FollowerForce']=XBinst.FollowerForce
    hdfile['FollowerForceRig']=XBinst.FollowerForceRig
    hdfile['PrintInfo']=XBinst.PrintInfo
    hdfile['OutInBframe']=XBinst.OutInBframe
    hdfile['OutInaframe']=XBinst.OutInaframe 
    hdfile['ElemProj']=XBinst.ElemProj
    hdfile['Solution']=XBinst.Solution    
    hdfile['MaxIterations']=XBinst.MaxIterations
    hdfile['NumLoadSteps']=XBinst.NumLoadSteps
    hdfile['DeltaCurved']=XBinst.DeltaCurved
    hdfile['MinDelta']=XBinst.MinDelta
    hdfile['NewmarkDamp']=XBinst.NewmarkDamp

    # input shared
    hdfile['NumElems']=XBinst.NumElems
    hdfile['NumNodesElem']=XBinst.NumNodesElem
    hdfile['ElemType']=XBinst.ElemType
    hdfile['BConds']=XBinst.BConds
    hdfile['TestCase']=XBinst.TestCase
    
    # design
    hdfile['BeamLength1']=XBinst.BeamLength1
    hdfile['TipMass']=XBinst.TipMass
    hdfile['TipMassY']=XBinst.TipMassY
    hdfile['TipMassZ']=XBinst.TipMassZ
    
    hdfile['BeamAxis']=XBinst.BeamAxis
    hdfile['beam_geometry']=XBinst.beam_geometry    
    
    hdfile['beam_shape']=XBinst.beam_shape        
    hdfile['cross_section_type']=XBinst.cross_section_type
    hdfile['material']=XBinst.material

    # static loads
    hdfile['AppStaLoadType']=XBinst.AppDynLoadType
    hdfile['NodeAppStaForce']=XBinst.NodeAppStaForce

    # dynamics loads
    hdfile['NumSteps']=XBinst.NumSteps
    hdfile['AppDynLoadType']=XBinst.AppDynLoadType
    hdfile['AppDynLoadVariation']=XBinst.AppDynLoadVariation
    hdfile['NodeAppDynForce']=XBinst.NodeAppDynForce
    hdfile['TimeRamp']=XBinst.TimeRamp
    hdfile['Omega']=XBinst.Omega

    # Prescribed velocities
    hdfile['PrVelType']=XBinst.PrVelType

    # ---------------------------------------------------- input without default
    # these values will appear as 'not found' if not allocated
    
    # parameters for constant beam span reconstruction 
    hdfile = conditional_saving(hdfile,XBinst,'cs_l2')
    hdfile = conditional_saving(hdfile,XBinst,'cs_l3')
    hdfile = conditional_saving(hdfile,XBinst,'cs_t2')
    hdfile = conditional_saving(hdfile,XBinst,'cs_t3')
    hdfile = conditional_saving(hdfile,XBinst,'cs_r')
    hdfile = conditional_saving(hdfile,XBinst,'cs_t')
    hdfile = conditional_saving(hdfile,XBinst,'cs_pretwist')
    
    # parameters for spline reconstruction
    hdfile = conditional_saving(hdfile,XBinst,'ControlElems')        
    hdfile = conditional_saving(hdfile,XBinst,'SplineOrder')
    hdfile = conditional_saving(hdfile,XBinst,'CS_l2')
    hdfile = conditional_saving(hdfile,XBinst,'CS_l3')
    hdfile = conditional_saving(hdfile,XBinst,'CS_t2')
    hdfile = conditional_saving(hdfile,XBinst,'CS_t3')
    hdfile = conditional_saving(hdfile,XBinst,'CS_r')
    hdfile = conditional_saving(hdfile,XBinst,'CS_t')
    
    
    # ------------------------------------------------------------ preprocessing
    # these values will appear as 'not found' if not allocated

    hdfile = conditional_saving(hdfile,XBinst,'NumNodes')

    # spline control Elements/Nodes
    hdfile = conditional_saving(hdfile,XBinst,'ControlNodes')
    hdfile = conditional_saving(hdfile,XBinst,'CS_pretwist')
    
    # cost/constrains
    hdfile = conditional_saving(hdfile,XBinst,'TotalMass')
    hdfile = conditional_saving(hdfile,XBinst,'NodeDisp')
    hdfile = conditional_saving(hdfile,XBinst,'ZDisp')
    hdfile = conditional_saving(hdfile,XBinst,'YDisp')
    hdfile = conditional_saving(hdfile,XBinst,'XDisp')
    hdfile = conditional_saving(hdfile,XBinst,'NNode')
    
    # __init__
    hdfile = conditional_saving(hdfile,XBinst,'ExtForce')
    hdfile = conditional_saving(hdfile,XBinst,'ExtMomnt')
    hdfile = conditional_saving(hdfile,XBinst,'ExtForceDyn')
    hdfile = conditional_saving(hdfile,XBinst,'ExtMomntDyn')
    hdfile = conditional_saving(hdfile,XBinst,'PrTrVelAmpl')
    hdfile = conditional_saving(hdfile,XBinst,'PrRotVelAmpl')
    
    ###hdfile = conditional_saving(hdfile,XBinst,'BeamSpanStiffness')
    ###hdfile = conditional_saving(hdfile,XBinst,'BeamSpanMass')
    ###hdfile = conditional_saving(hdfile,XBinst,'PhiNodes')
    hdfile = conditional_saving(hdfile,XBinst,'Mroot')
    hdfile = conditional_saving(hdfile,XBinst,'Kroot')
    
    # output
    hdfile = conditional_saving(hdfile,XBinst,'DensityVector')
    hdfile = conditional_saving(hdfile,XBinst,'LengthVector')
    hdfile = conditional_saving(hdfile,XBinst,'PosIni')
    hdfile = conditional_saving(hdfile,XBinst,'PosDef')
    hdfile = conditional_saving(hdfile,XBinst,'PsiIni')
    hdfile = conditional_saving(hdfile,XBinst,'PsiDef')
    hdfile = conditional_saving(hdfile,XBinst,'InternalForces')
  
    #------------------------------------------------------------ Static Related
    hdfile = conditional_saving(hdfile,XBinst,'ForceStatic')  

    #---------------------------------------------------------- Dynamics Related
    # input
    hdfile = conditional_saving(hdfile,XBinst,'Time')   
    hdfile = conditional_saving(hdfile,XBinst,'ForceDynAmp')   
    hdfile = conditional_saving(hdfile,XBinst,'ForceTime')   
    hdfile = conditional_saving(hdfile,XBinst,'ForcedVel')   
    hdfile = conditional_saving(hdfile,XBinst,'ForcedVelDot') 
    # output
    hdfile = conditional_saving(hdfile,XBinst,'PosDotDef')   
    hdfile = conditional_saving(hdfile,XBinst,'PsiDotDef')   
    hdfile = conditional_saving(hdfile,XBinst,'PosPsiTime')   
    hdfile = conditional_saving(hdfile,XBinst,'VelocTime')   
    hdfile = conditional_saving(hdfile,XBinst,'DynOut')   

    #---------------------------------------------------------------- Rigid Body
    # output
    hdfile = conditional_saving(hdfile,XBinst,'RefVel')  
    hdfile = conditional_saving(hdfile,XBinst,'RefVelDot')  
    hdfile = conditional_saving(hdfile,XBinst,'Quat')  
    
    
    #-------------------------------------------------------------- Optimisation
    hdfile = conditional_saving(hdfile,XBinst,'fval')  
    hdfile = conditional_saving(hdfile,XBinst,'fname')  
    hdfile = conditional_saving(hdfile,XBinst,'fargs')  
    
    
    '''
    # to add future fields...
    hdfile['']=XBinst.
    hdfile = conditional_saving(hdfile,XBinst,'')    
    '''
    
    hdfile.close() 
     
    return None
 
 
def conditional_saving(hdfile,obj,attrname):
    ''' 
    Given the object obje and the attribute 'attrname', the routine:
        a. if the attribute is found:
            1. if obj[attrname] has a value assigned, saves the value under the 
            field attrname of the h5 file hdfile 
            2. saves 'no value'   
        b. saves 'not found' otherwise    
    '''

    try:
        val = getattr(obj,attrname)
        #print 'attribute found'
    except:
        #print 'attribute not found'
        val = 'not found'

    try:
        hdfile[attrname]=val
    except:
        #print 'no value for this attribute'
        hdfile[attrname]='no value'
    
    return hdfile
    
    
    
     
 
 
 
    
    
if __name__=='__main__':
    
    import numpy as np
    
    import sys
    sys.path.append('..')
    sys.path.append('../..')
    
    from xbcomponent import XBeamSolver
    
    #---------------------------------------------------------------- Test h5sol
    PosIni=np.array([ [0.0, 10.0, 0.0 ],
                      [1.0, 11.0, 5.0 ] ])
    PosDef=np.array([ [1.0, 90.0, 4.0 ],
                      [5.0,  9.0, 2.0 ] ])
    
    PsiIni, PsiDef = 2.0*PosIni, 10.0*PosDef
    filename ='save_test.h5'
    h5sol(PosIni,PsiIni,PosDef,PsiDef,filename)
    
    
    #--------------------------------------------------------------- Test h5comp
    XBinst=XBeamSolver()
    
    
    XBinst.PsiIni=np.array([1,2,3])
    
    h5comp(XBinst)
    #h5comp(XBinst,'lalalacomp.h5')
    
    
        
    

