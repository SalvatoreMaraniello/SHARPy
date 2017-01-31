'''
@author: salvatore maraniello
@contact: salvatore.maraniello10@imperial.ac.uk
@date: 25 Jan 2017
@brief: unittest class to test the flight-dynamics solution with constraints 
        kinematics. A HALE wing response to a prescribed aileron input is 
        computed. Wing orientation and deflections are verified at both the
        initial and final state. 
@warning: A very low fidelity model is used.
@quote: You will pay for your sins. If you have already paid, please disregard 
        this message.
'''

import os
import sys
import numpy as np
import ctypes as ct
import unittest

sys.path.append( os.environ["SHARPYDIR"]+'/src' )
sys.path.append( os.environ["SHARPYDIR"]+'/src/Main' )

import SharPySettings as Settings
import DerivedTypes, DerivedTypesAero
from PyCoupled.Utils.DerivedTypesAeroelastic import AeroelasticOps
from PyCoupled.Coupled_NlnFlightDynamic_lag import Solve_Py

from DerivedTypesAero import ControlSurf
import PyLibs.CVP.spline

TOL=1e-8


class TestConstrainedFlightDynamics(unittest.TestCase):
    '''
    Each method defined in this class contains a test case.
    @warning: by default, only functions whose name starts with 'test' will be 
    run during testing.
    '''


    def setUp(self):
        '''
        Common piece of code run by each test
        '''        

        # wing model: common specifications
        c = 1. # chord
        Umag = 0. # speed of FoR A
        Umag_flow = 30. # flow speed
        fi = 5.0*np.pi/180.0 # attitude of FoR A
        alpha = 0.0*np.pi/180.0 # flow angle
        NumElems = 2 * 8  # finite elements (multiple of 8 to match ailerons)
        M = 1 * 4  # chord-wise panels (multiple of 4 to match ailerons)
        NumNodesElem = 3 
        NumElemsHere = int(0.5*NumElems)
        Ueff = Umag_flow # magnitude of relative speed (flow vs. FoR A)
        delTime = c/( M*(Ueff) ) # time-step
        Ttot= 2.0 # simulation time

        # saving options
        self.savedict=Settings.SaveDict
        self.savedict['OutputDir'] = Settings.SharPyProjectDir + \
                                      'output/tests/PyCoupled/NlnFlightDynamic/'
        self.savedict['OutputFileRoot'] = 'HALE_hinged_%0.1fms%0.2fdegM%0.2dNE%0.2d' \
                                          %(Umag_flow,fi*180.0/np.pi,M,NumElems,)
        self.savedict['Format']='h5'
        self.savedict['NumSavePoints']=20
        self.savedict['SaveWake']=False
        self.savedict['WaveSaveFreq']=50


        # ------------------------------------------------ Beam solution options
        self.xbopts = DerivedTypes.Xbopts(FollowerForce = ct.c_bool(False),
                                     MaxIterations = ct.c_int(25),
                                     PrintInfo = ct.c_bool(True),
                                     OutInaframe = ct.c_bool(True),
                                     NumLoadSteps = ct.c_int(15),
                                     Solution = ct.c_int(912),
                                     MinDelta = ct.c_double(1e-5),
                                     NewmarkDamp = ct.c_double( 5e-2))


        # ----------------------------------------------------------- Beam input
        self.xbinput = DerivedTypes.Xbinput(NumNodesElem,NumElemsHere,
                                                            BConds='MS',g=9.754)
        self.xbinput.BeamLength = 2.0*16.0
        self.xbinput.addRBMass(node=NumElemsHere,  Mmat=np.diag(3*[0.0] + 3*[0.0]))
   
        # get mass for sig1=1 and stiffness for sig1=sig1
        mvec=np.array([.75,.75,.75,.1,.001,.001])
        for ii in range(6): self.xbinput.BeamMass[ii,ii]=mvec[ii]

        # pitch-plunge coupling term (b-frame coordinates)
        distIA = 0.5 # perc of chord from LE    
        distEA = 0.5 # perc of chord from LE  
        cg = c*np.array([ 0.0, distEA-distIA, 0.0  ]) # positive if CG ahead
        cgSkew = np.array([[   0.0, -cg[2], cg[1] ],\
                           [ cg[2],    0.0, -cg[0]],\
                           [-cg[1],  cg[0],   0.0] ])
        self.xbinput.BeamMass[:3,3:] = -self.xbinput.BeamMass[0,0]*cgSkew
        self.xbinput.BeamMass[3:,:3] = self.xbinput.BeamMass[:3,3:].T
        self.xbinput.BeamMass[0,5] = 0.0
        self.xbinput.BeamMass[5,0] = 0.0
        
        # damping model
        self.xbinput.sph_joint_damping=(1e-1)*2000/6.7
        self.xbinput.str_damping_model = 'prop'
        self.xbinput.str_damping_param={'alpha': 0.005, 'beta': 0.0007}
        
        # Unsteady parameters
        self.xbinput.dt = delTime
        self.xbinput.t0 = 0.0
        self.xbinput.tfin = Ttot
        self.xbinput.PsiA_G=fi*np.array([1,0,0])
        self.xbinput.EnforceAngVel_FoRA=[True, False, True]
        self.xbinput.EnforceTraVel_FoRA=3*[True]
        
        # Set motion of wing.
        NumSteps = np.ceil( (self.xbinput.tfin + self.xbinput.dt - 
                                            self.xbinput.t0) / self.xbinput.dt)
        self.xbinput.ForcedVel = np.zeros((NumSteps,6),ct.c_double,'F')
        AOA=fi # because input velocity is in A FoR!!!
        for i in range(self.xbinput.ForcedVel.shape[0]):
            self.xbinput.ForcedVel[i,:] = [0.0, Umag*np.cos(AOA), 
                                              -Umag*np.sin(AOA), 0.0, 0.0, 0.0]
        self.xbinput.ForcedVelDot = np.zeros((NumSteps,6),ct.c_double,'F')


        #--------------------------------------------------------- # VLM options
        N=self.xbinput.NumNodesTot - 1 # span-wise UVLM panels
        WakeLength = 20.0*c
        self.vmopts = DerivedTypesAero.VMopts(M = M, 
                                         N = N,
                                         ImageMethod = False,
                                         Mstar = int(WakeLength/(delTime*Ueff)),
                                         Steady = False,
                                         KJMeth = True,
                                         NewAIC = True,
                                         DelTime = delTime,
                                         NumCores = 4)


        #------------------------------------------------ # HALE ailerons Common
        typeMotion='asInput'
        # control parametrisation
        NC = 8 # to get a spline impulse around 20 Hz
        p=3
        NS=NC-1+p
        iMin = np.int(np.round( 0.75*M ))  #M - M/4
        iMax = M
        Non = 1.5/Ttot*NC
        iistart=3
        iion = [ii for ii in range(iistart,iistart+int(Non))]
        

        #------------------------------------------------- # HALE ailerons right
        # indices given in python numbering (from zero)
        jMinR = np.int(np.round( 0.75*N )) 
        jMaxR = N
        
        def asInput_fun_R():
            # note: some input arg defined out of this function
            tv = np.arange(self.xbinput.t0,self.xbinput.tfin  + 
                                              self.xbinput.dt, self.xbinput.dt)
            p=3
            tcint = np.linspace(tv[0],tv[-1],NC)
            scfv=np.zeros(NS)
            scfv[iion]=6.5*np.pi/180.0
            (y,dy)=PyLibs.CVP.spline.build_vec(tv, tcint,scfv,p,EvalDeriv=True)   
            return tv,y,dy
        
        ctrlSurf_R = ControlSurf(iMin,     # index at which control surface starts (in chord-wise).
                               iMax,       #                ....             ends      .... 
                               jMinR,      #                ....            starts (in span-wise).
                               jMaxR,      #                ....             ends      ....
                               typeMotion,
                               0.0, # amplitude of control surface angle.
                               0.0, # control surface rates (backward difference)
                               asInput_fun = asInput_fun_R,
                               asInput_args = () )    
        
        
        #------------------------------------------------- # HALE ailerons right
        # indices given in python numbering (from zero)
        jMinL = 0 #np.int(np.round( 0.75*N )) 
        jMaxL = np.int(np.round( 0.25*N ))
        
        # build sample input
        def asInput_fun_L():
            # note: some input arg defined out of this function
            tv = np.arange(self.xbinput.t0,self.xbinput.tfin  + 
                                            self.xbinput.dt, self.xbinput.dt)
            p=3
            tcint = np.linspace(tv[0],tv[-1],NC)
            scfv=np.zeros(NS)
            scfv[iion]=-6.5*np.pi/180.0
            (y,dy)=PyLibs.CVP.spline.build_vec(tv, tcint,scfv,p,EvalDeriv=True)   
            return tv,y,dy
        
        ctrlSurf_L = ControlSurf(iMin,     # index at which control surface starts (in chord-wise).
                               iMax,       #                ....             ends      .... 
                               jMinL,      #                ....            starts (in span-wise).
                               jMaxL,      #                ....             ends      ....
                               typeMotion, # string describing type of prescribed, or other, motion.
                               0.0, # amplitude of control surface angle.
                               0.0, # control surface rates (backward difference)
                               asInput_fun = asInput_fun_L,
                               asInput_args = () )    
        
        
        #---------------------------------------------------------- # UVLM input
        self.vminput = DerivedTypesAero.VMinput(c = c,
                                           b = self.xbinput.BeamLength,
                                           U_mag = Umag_flow,
                                           alpha = alpha,
                                           theta = 0.0,
                                           WakeLength = WakeLength,
                                           ctrlSurf = [ctrlSurf_R, ctrlSurf_L],
                                           gust = None)
         

        #---------------------------------------- # aeroelastic solution options    
        self.aelaopts = AeroelasticOps(ElasticAxis = 0.0,
                                  InertialAxis = 0.0,
                                  AirDensity = 0.08891,
                                  Tight = False,
                                  ImpStart = False,
                                  MinRes = 1e10)
        self.aelaopts.MaxRes=1e10
        self.aelaopts.LimRes=1e10



    def testVeryFlexibleWing(self):

        # wing specification
        sig1=1.1
        kvec=np.array([1e7,1e7,1e7,sig1*1e4,sig1*2e4,5.e6])
        for ii in range(6): self.xbinput.BeamStiffness[ii,ii]=kvec[ii]

        self.savedict['OutputFileRoot'] = self.savedict['OutputFileRoot'] \
                                                          + '_sig%0.1f'%(sig1,)

        xbout=Solve_Py(self.xbinput, self.xbopts, self.vmopts, self.vminput,
                       self.aelaopts, SaveDict=self.savedict)


        # ------------------------------------------------------------- Testing

        # Tip initial/final position (FoR A)
        Tip0 = np.array([-14.074999557383318, -0.10975597617063762, 
                                                             7.063461979098654])
        TipEnd=np.array([-14.178217935221744, -0.09523382983611592,	
                                                             6.891872959990548])
        self.assertTrue(np.linalg.norm(xbout.DynOut[0,:]-Tip0) < TOL,
                        msg='Node 0 initial position: relative difference '  
                                                   'above tolerance %.2e' %TOL )
        self.assertTrue(np.linalg.norm(xbout.DynOut[-self.xbinput.NumNodesTot,:]
                        -TipEnd) < TOL, msg='Node 0 final position: relative '
                                        'difference above tolerance %.2e' %TOL )


        # Cartesian rotation vectors
        PsiTip0x = 0.059785862902468422
        PsiTipEndx = 0.076647741437134667
        self.assertTrue( np.abs(PsiTip0x-xbout.PsiDeforStatic[0,0,0]) < TOL,
                         msg='Node 0 initial x rotation: relative difference '
                                                    'above tolerance %.2e'%TOL )
        self.assertTrue( np.abs(PsiTipEndx-xbout.PsiList[-1][0,0,0]) < TOL,
                         msg='End nNode final x rotation: relative difference '
                                                   'above tolerance %.2e'%TOL  )

        # Quaternions final
        QuatEnd=np.array([ 0.9911874,  0.04327618, -0.12507956, -0.00546109])
        self.assertTrue( np.linalg.norm(xbout.QuatList[-1]-QuatEnd) < TOL,
                         msg='Final quaternion: relative difference above '
                                                          'tolerance %.2e'%TOL )


        return xbout



    def testVeryStiffWing(self):

        # wing specification
        sig1=50.0
        kvec=np.array([1e7,1e7,1e7,sig1*1e4,sig1*2e4,5.e6])
        for ii in range(6): self.xbinput.BeamStiffness[ii,ii]=kvec[ii]

        self.savedict['OutputFileRoot'] = self.savedict['OutputFileRoot'] + \
                                                             '_sig%0.1f'%(sig1,)

        xbout=Solve_Py(self.xbinput,self.xbopts,self.vmopts,self.vminput,self.aelaopts,
                       SaveDict=self.savedict)


        # ---------------------------------------------------------------Testing

        # Tip initial/final position (FoR A)
        Tip0 = np.array([-15.999638684431325, 0.0014610025322748543, 
                                                           0.10079232151504619])
        TipEnd=np.array([-15.999717588696065, 0.003799444708618084, 
                                                           0.09428770140427575])
        self.assertTrue(np.linalg.norm(xbout.DynOut[0,:]-Tip0) < TOL,
                        msg='Node 0 initial position: relative difference '  
                                                   'above tolerance %.2e' %TOL )
        self.assertTrue(np.linalg.norm(xbout.DynOut[-self.xbinput.NumNodesTot,:]
                        -TipEnd) < TOL, msg='Node 0 final position: relative '
                                        'difference above tolerance %.2e' %TOL )

        # Cartesian rotation vectors
        PsiTip0x = 0.0012716072782286982
        PsiTipEndx = 0.0020322238658857287
        self.assertTrue( np.abs(PsiTip0x-xbout.PsiDeforStatic[0,0,0]) < TOL,
                         msg='Node 0 initial x rotation: relative difference '
                                                    'above tolerance %.2e'%TOL )
        self.assertTrue( np.abs(PsiTipEndx-xbout.PsiList[-1][0,0,0]) < TOL,
                         msg='End nNode final x rotation: relative difference '
                                                   'above tolerance %.2e'%TOL  )

        # Quaternions final
        QuatEnd=np.array([ 0.9872585711139502, 0.04310464010934549, 
                                     -0.1530289606836371, -0.00668138871575723])
        self.assertTrue( np.linalg.norm(xbout.QuatList[-1]-QuatEnd) < TOL,
                         msg='Final quaternion: relative difference above '
                                                          'tolerance %.2e'%TOL )

        return xbout



if __name__=='__main__':

    unittest.main()

    #T=TestConstrainedFlightDynamics()
    #T.setUp()
    #xbout=T.testVeryFlexibleWing()
    #xbout=T.testVeryStiffWing()
    
 

    