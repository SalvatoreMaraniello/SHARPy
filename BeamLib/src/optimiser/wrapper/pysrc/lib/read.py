'''
Created on 20 Aug 2014

@author: sm6110

Tools for reading.
-h5sol (superseded): reads displacements and rotations from a solution h5 file.
-h5comp: reads all variables of a *h5 solution and creates a component from them
- h5series: given a list of attributes, reads all the fields and puts them into an
array.


'''

import h5py
from warnings import warn

import sys
sys.path.append('..')
sys.path.append('../..')
from xbcomponent import XBeamSolver



def h5comp(filename):
    ''' 
    Reads All variables of an instance of class XBeamSolver from a hdf5 file.
    
    For a details of the variable see: xbcomponent.XBeamSolver class.
    
    Last update: 21 Aug 2014
    
    Input:
    - filename: name of file
    
    Remarks:
    - If filename contains the solution of a simulation, the arrays read will have
    a Fortran ordering. This is not strictly required though, as for restart, 
    only the inputs are required - the rest is reallocated during the XBeamSolver
    execution
    - fields expected to be found will lead the execution to an error (inputs)
    - fields that are created during the execution, if not found, will only lead
    produce a warning message
    '''    
    
    
    # allocate instance
    XBinst = XBeamSolver()
    
    # open file
    hdfile = h5py.File(filename,'r')
    
    
    #-------------------------------------------------------- input with default
    # options:
    setattr(XBinst,'FollowerForce',hdfile['FollowerForce'].value)
    setattr(XBinst,'FollowerForceRig',hdfile['FollowerForceRig'].value)
    setattr(XBinst,'PrintInfo',hdfile['PrintInfo'].value)
    setattr(XBinst,'OutInBframe',hdfile['OutInBframe'].value)
    setattr(XBinst,'OutInaframe',hdfile['OutInaframe'].value)
    setattr(XBinst,'ElemProj',hdfile['ElemProj'].value)
    setattr(XBinst,'Solution',hdfile['Solution'].value)
    setattr(XBinst,'MaxIterations',hdfile['MaxIterations'].value)
    setattr(XBinst,'NumLoadSteps',hdfile['NumLoadSteps'].value)
    setattr(XBinst,'DeltaCurved',hdfile['DeltaCurved'].value)
    setattr(XBinst,'MinDelta',hdfile['MinDelta'].value)
    setattr(XBinst,'NewmarkDamp',hdfile['NewmarkDamp'].value)

    # input shared
    setattr(XBinst,'NumElems',hdfile['NumElems'].value)
    setattr(XBinst,'NumNodesElem',hdfile['NumNodesElem'].value)
    setattr(XBinst,'ElemType',hdfile['ElemType'].value)
    setattr(XBinst,'BConds',hdfile['BConds'].value)
    setattr(XBinst,'TestCase',hdfile['TestCase'].value)    
    
    # design
    setattr(XBinst,'BeamLength1',hdfile['BeamLength1'].value)
    setattr(XBinst,'TipMass',hdfile['TipMass'].value)
    setattr(XBinst,'TipMassY',hdfile['TipMassY'].value)    
    setattr(XBinst,'TipMassZ',hdfile['TipMassZ'].value)

    setattr(XBinst,'BeamAxis',hdfile['BeamAxis'].value)
    setattr(XBinst,'beam_geometry',hdfile['beam_geometry'].value)


    setattr(XBinst,'beam_shape',hdfile['beam_shape'].value)
    setattr(XBinst,'cross_section_type',hdfile['cross_section_type'].value)
    setattr(XBinst,'material',hdfile['material'].value)

    # static loads
    setattr(XBinst,'AppStaLoadType',hdfile['AppStaLoadType'].value)
    setattr(XBinst,'NodeAppStaForce',hdfile['NodeAppStaForce'].value) 
 
    # dynamics loads
    setattr(XBinst,'NumSteps',hdfile['NumSteps'].value)
    setattr(XBinst,'AppDynLoadType',hdfile['AppDynLoadType'].value)
    setattr(XBinst,'AppDynLoadVariation',hdfile['AppDynLoadVariation'].value)
    setattr(XBinst,'NodeAppDynForce',hdfile['NodeAppDynForce'].value)
    setattr(XBinst,'TimeRamp',hdfile['TimeRamp'].value)
    setattr(XBinst,'Omega',hdfile['Omega'].value)

    # Prescribed velocities
    setattr(XBinst,'PrVelType',hdfile['PrVelType'].value) 
    
    # ---------------------------------------------------- input without default
    # these values will appear as 'not found' if not allocated
    
    # parameters for constant beam span reconstruction 
    XBinst = conditional_reading(hdfile,XBinst,'cs_l2')
    XBinst = conditional_reading(hdfile,XBinst,'cs_l3')
    XBinst = conditional_reading(hdfile,XBinst,'cs_t2')
    XBinst = conditional_reading(hdfile,XBinst,'cs_t3')
    XBinst = conditional_reading(hdfile,XBinst,'cs_r')
    XBinst = conditional_reading(hdfile,XBinst,'cs_t')
    XBinst = conditional_reading(hdfile,XBinst,'cs_pretwist')
    
    # parameters for spline reconstruction
    XBinst = conditional_reading(hdfile,XBinst,'ControlElems')        
    XBinst = conditional_reading(hdfile,XBinst,'SplineOrder')
    XBinst = conditional_reading(hdfile,XBinst,'CS_l2')
    XBinst = conditional_reading(hdfile,XBinst,'CS_l3')
    XBinst = conditional_reading(hdfile,XBinst,'CS_t2')
    XBinst = conditional_reading(hdfile,XBinst,'CS_t3')
    XBinst = conditional_reading(hdfile,XBinst,'CS_r')
    XBinst = conditional_reading(hdfile,XBinst,'CS_t')
    
    
    # ------------------------------------------------------------ preprocessing
    # these values will appear as 'not found' if not allocated

    XBinst = conditional_reading(hdfile,XBinst,'NumNodes')

    # spline control Elements/Nodes
    XBinst = conditional_reading(hdfile,XBinst,'ControlNodes')
    XBinst = conditional_reading(hdfile,XBinst,'CS_pretwist')
    
    # cost/constrains
    XBinst = conditional_reading(hdfile,XBinst,'TotalMass')
    XBinst = conditional_reading(hdfile,XBinst,'NodeDisp')
    XBinst = conditional_reading(hdfile,XBinst,'ZDisp')
    XBinst = conditional_reading(hdfile,XBinst,'YDisp')
    XBinst = conditional_reading(hdfile,XBinst,'XDisp')
    XBinst = conditional_reading(hdfile,XBinst,'NNode')
    
    # __init__
    XBinst = conditional_reading(hdfile,XBinst,'ExtForce')
    XBinst = conditional_reading(hdfile,XBinst,'ExtMomnt')
    XBinst = conditional_reading(hdfile,XBinst,'ExtForceDyn')
    XBinst = conditional_reading(hdfile,XBinst,'ExtMomntDyn')  
    XBinst = conditional_reading(hdfile,XBinst,'PrTrVelAmpl')
    XBinst = conditional_reading(hdfile,XBinst,'PrRotVelAmpl')  
  
    XBinst = conditional_reading(hdfile,XBinst,'Mroot')
    XBinst = conditional_reading(hdfile,XBinst,'Kroot')
    ###XBinst = conditional_reading(hdfile,XBinst,'BeamSpanStiffness')
    ###XBinst = conditional_reading(hdfile,XBinst,'BeamSpanMass')
    ###XBinst = conditional_reading(hdfile,XBinst,'PhiNodes')
    
    # output
    XBinst = conditional_reading(hdfile,XBinst,'DensityVector')
    XBinst = conditional_reading(hdfile,XBinst,'LengthVector')
    XBinst = conditional_reading(hdfile,XBinst,'PosIni')
    XBinst = conditional_reading(hdfile,XBinst,'PosDef')
    XBinst = conditional_reading(hdfile,XBinst,'PsiIni')
    XBinst = conditional_reading(hdfile,XBinst,'PsiDef')
    XBinst = conditional_reading(hdfile,XBinst,'InternalForces')    
    
    #------------------------------------------------------------ Static Related
    XBinst = conditional_reading(hdfile,XBinst,'ForceStatic')     
    
    #---------------------------------------------------------- Dynamics Related
    # input
    XBinst = conditional_reading(hdfile,XBinst,'Time')   
    XBinst = conditional_reading(hdfile,XBinst,'ForceDynAmp')   
    XBinst = conditional_reading(hdfile,XBinst,'ForceTime')   
    XBinst = conditional_reading(hdfile,XBinst,'ForcedVel')   
    XBinst = conditional_reading(hdfile,XBinst,'ForcedVelDot')
    XBinst = conditional_reading(hdfile,XBinst,'ForceDynamic')
    
    XBinst = conditional_reading(hdfile,XBinst,'Acf') # fourier
    XBinst = conditional_reading(hdfile,XBinst,'Bcf')
    XBinst = conditional_reading(hdfile,XBinst,'Fcf')
    XBinst = conditional_reading(hdfile,XBinst,'Scf') # spline
    XBinst = conditional_reading(hdfile,XBinst,'TCint')
    XBinst = conditional_reading(hdfile,XBinst,'DynFrc_spline_order')
    
    
    # output
    XBinst = conditional_reading(hdfile,XBinst,'PosDotDef')   
    XBinst = conditional_reading(hdfile,XBinst,'PsiDotDef')   
    XBinst = conditional_reading(hdfile,XBinst,'PosPsiTime')   
    XBinst = conditional_reading(hdfile,XBinst,'VelocTime')   
    XBinst = conditional_reading(hdfile,XBinst,'DynOut')     
    
    #---------------------------------------------------------------- Rigid Body
    # output
    XBinst = conditional_reading(hdfile,XBinst,'RefVel')  
    XBinst = conditional_reading(hdfile,XBinst,'RefVelDot')  
    XBinst = conditional_reading(hdfile,XBinst,'Quat')
    
    
    #----------------------------------------------------- ics (restart sol 932)
    XBinst = conditional_reading(hdfile,XBinst,'Quat0')
    XBinst = conditional_reading(hdfile,XBinst,'RefVel0') 
    XBinst = conditional_reading(hdfile,XBinst,'RefVelDot0') 
    XBinst = conditional_reading(hdfile,XBinst,'PosDotDef0')   
    XBinst = conditional_reading(hdfile,XBinst,'PsiDotDef0')       
    XBinst = conditional_reading(hdfile,XBinst,'PosDef0')
    XBinst = conditional_reading(hdfile,XBinst,'PsiDef0')        
    
    #-------------------------------------------------------------- Optimisation
    XBinst = conditional_reading(hdfile,XBinst,'fval') 
    XBinst = conditional_reading(hdfile,XBinst,'fname') 
    XBinst = conditional_reading(hdfile,XBinst,'fargs') 
    
    XBinst = conditional_reading(hdfile,XBinst,'geqval') 
    XBinst = conditional_reading(hdfile,XBinst,'geqname') 
    XBinst = conditional_reading(hdfile,XBinst,'geqargs') 

    XBinst = conditional_reading(hdfile,XBinst,'gdisval') 
    XBinst = conditional_reading(hdfile,XBinst,'gdisname') 
    XBinst = conditional_reading(hdfile,XBinst,'gdisargs') 

    
    '''
    setattr(XBinst,'',hdfile[''].value)
    XBinst = conditional_reading(hdfile,XBinst,'')
    '''
    
    # close and return
    hdfile.close()
    
    return XBinst
    
    
def h5series(rootname,attrlist):
    ''' 
    Given a list of attributes, creates a list of lists for all the solutions
    run for a DOE or optimisation. the output is in a
     list format to allow also 
    non numerical variables to change.
    '''
    
    outlist=[]
    
    cc=0
    go_on=True
    
    while go_on is True:
        cc_str =  '%.3d' % (cc)
        filename = rootname + cc_str + '.h5'
        
        try:
            print 'Reading: %s' %(filename)
            hdfile = h5py.File(filename,'r') 
            inlist=[]
            for attr in attrlist:
                try:
                    inlist.append(hdfile[attr].value)
                except:
                    warn('%s attribute not found!!!' %(attr)) 
                    inlist.append('not found!!!')      
            outlist.append(inlist)
            hdfile.close()
            cc=cc+1
        except:
            print '%s not found. %s files read in total!' %(filename,cc_str)
            go_on=False
    
    return outlist
     
  
def conditional_reading(hdfile,obj,attrname): 
    ''' 
    Given a hdf5 file 'hdfile' and the object obj, the routine:
        a. if the attribute attrname is found and has a value, assignes it to the
           attribute obj.attrname.  
        b. does nothing otherwise    
    '''
    
    try:
        val = hdfile[attrname].value
        if val!='not found' and val!='no value':
            setattr(obj,attrname,val)
    except:
        warn('Attribute "%s" not found!!!' %attrname)
            
    return obj      
     


if __name__=='__main__':
    
    import numpy as np

    
    #--------------------------------------------------------------- h5comp test
    #filename = 'fwd_model.h5' 
    #filename='save_test.h5'
    XBinst = h5comp(filename)
    # print all the attributes
    #print XBinst.get_attributes()
    for tt in XBinst.items():
        print tt
    
    
    
    '''
    #------------------------------------------------------------- h5series test
    rootname = '/home/sm6110/git/SHARPy_studies/OPT/20140820_validation_static_opt/res_opt45/isorec_thick_45deg_6kN_sol102_'
    attrlist = ['cs_l2', 'cs_l3','TotalMass','XDisp','YDisp','ZDist','NodeDisp']
    outlist = h5series(rootname,attrlist)
    outvals = np.array(outlist)
    print outvals
    '''
 
    