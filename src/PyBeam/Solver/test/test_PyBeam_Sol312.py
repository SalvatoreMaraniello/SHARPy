'''
@author: salvatore maraniello
@contact: salvatore.maraniello10@imperial.ac.uk
@date: 16 Feb 2017
@brief: unittest class to test the gravity loads implementation in the PyBeam
        solution methods.
@warning: - Very low fidelity model is used for sol 312
          - Sol112 and Sol112F90 produce different results (within a 0.1% 
          relative error)
          - a higher tolerance (1e-3) is used for cases with rigid follower force
@quote: 
'''

import os
import sys
import numpy as np
import scipy.optimize as scopt
import scipy.integrate as scint

import ctypes as ct
import matplotlib.pyplot as plt
import unittest
from IPython import embed

sys.path.append( os.environ["SHARPYDIR"]+'/src' )
sys.path.append( os.environ["SHARPYDIR"]+'/src/Main' )

import SharPySettings as Settings
import DerivedTypes
from PyBeam.Solver.NonlinearStatic import Solve_Py as Sol112
#from PyBeam.Solver.NonlinearStatic import Solve_F90 as Sol112F90
from PyBeam.Solver.NonlinearDynamic import Solve_Py as Sol312

import PyLibs.numerics.integr
import PyLibs.plot.dyn
from PyLibs.plot.shared import fontlabel, params
import PyBeam.Utils.PostPr
import XbeamLib
import PyLibs.CVP.fourier


TOL=1e-8
TOL_rig_false=1e-3

# ------------------------------------------------------------------------------


class Test312(unittest.TestCase):
    '''
    Each method defined in this class contains a test case.
    @warning: by default, only functions whose name starts with 'test' will be 
    run during testing.

    Reproduce free falling clamped beam very flexible case from:
    Wang, Q. & Yu, W., 2011. Sensitivity Analysis of Geometrically Exact Beam 
    Theory (GEBT) Using the Adjoint Method with Hydra. 
    In 52nd AIAA/ASME/ASCE/AHS/ASC Structures, Structural Dynamics and Materials 
    Conference. Denver, Colorado, pp. 1–15.Wang
    '''



    def setUp(self):
        '''
        Common piece of code run by each test
        '''        
        self.PLOT=True

        # Beam solution options
        self.xbopts = DerivedTypes.Xbopts(
                             FollowerForce = ct.c_bool(False),
                             FollowerForceRig=ct.c_bool(True),
                             MaxIterations = ct.c_int(20),
                             PrintInfo = ct.c_bool(False),
                             OutInaframe = ct.c_bool(True),
                             NumLoadSteps = ct.c_int(15),
                             Solution = ct.c_int(312),
                             MinDelta = ct.c_double(1e-6),
                             NewmarkDamp = ct.c_double(1e-2),
                             ImpStart=True)

        # Beam input
        self.xbinput = DerivedTypes.Xbinput(NumNodesElem=3,NumElems=4,
                                                             BConds='CF',g=9.80)
        self.xbinput.BeamLength = 1.0

        # Mass/Stiffness
        mvec=np.array([.15,.15,.15,.1,.001,.001])
        for ii in range(6): self.xbinput.BeamMass[ii,ii]=mvec[ii]
        kvec=np.array([1e7,1e7,1e7,.15,.15,.15])
        for ii in range(6): self.xbinput.BeamStiffness[ii,ii]=kvec[ii]

        # Unsteady parameters
        self.xbinput.dt = 0.001
        self.xbinput.t0 = 0.0
        self.xbinput.tfin = 1.0
        self.xbinput.PsiA_G=30.*np.pi/180.*np.array([1,0,0]) # rotation to FoR A

        # saving options
        self.savedict=Settings.SaveDict
        self.savedict['OutputDir'] = Settings.SharPyProjectDir + \
                                      'output/tests/PyBeam/NonlinearDynamic/'

        return self



    def test312_ImpTrue_FollRigTrue(self):

        self.xbopts.ImpStart=True
        self.xbopts.FollowerForceRig.value=True
        self.savedict['OutputFileRoot'] = 'beam%s_Imp%s_FollRig%s_NE%0.2d' \
                 %( self.xbinput.BConds, self.xbopts.ImpStart, 
                    self.xbopts.FollowerForceRig.value, 
                    self.xbinput.NumElems)
        xbout=Sol312(self.xbinput, self.xbopts, SaveDict=self.savedict)


        # Ref solution run with sol 312 and no rotation of FoR A 
        RefTipFinal_G=np.array([0.7587529028679988,0.0,-0.6062938032858286])
        Cao=XbeamLib.Rot(xbout.QuatList[-1]).transpose()
        
        # Project tip position at the end of simulation in FoR G
        TipG = np.dot(Cao,xbout.DynOut[-1,:])
        Error = np.linalg.norm(TipG-RefTipFinal_G)  

        self.assertTrue(Error<TOL, msg='Tip final position error of %.3e above '
                                                'tolerance %.2e!' %(Error,TOL) )

        return xbout


    def test312_ImpTrue_FollRigFalse(self):

        self.xbopts.ImpStart=True
        self.xbopts.FollowerForceRig.value=False
        self.savedict['OutputFileRoot'] = 'beam%s_Imp%s_FollRig%s_NE%0.2d' \
                 %( self.xbinput.BConds, self.xbopts.ImpStart, 
                    self.xbopts.FollowerForceRig.value, 
                    self.xbinput.NumElems)
        xbout=Sol312(self.xbinput, self.xbopts, SaveDict=self.savedict)


        # Ref solution run with sol 312 and no rotation of FoR A 
        RefTipFinal_G=np.array([0.7587529028679988,0.0,-0.6062938032858286])
        Cao=XbeamLib.Rot(xbout.QuatList[-1]).transpose()
        #Cao_check=XbeamLib.RotCRV(self.xbinput.PsiA_G)
        
        # Project tip position at the end of simulation in FoR G
        TipG = np.dot(Cao,xbout.DynOut[-1,:])
        Error = np.linalg.norm(TipG-RefTipFinal_G)  

        # check
        self.assertTrue(Error<TOL_rig_false, msg='Tip final position error of '
                           '%.3e above tolerance %.2e!' %(Error,TOL_rig_false) )

        return xbout


    def test312_ImpFalse_FollRigTrue(self):

        self.xbopts.ImpStart=False
        self.xbopts.FollowerForceRig.value=True
        self.savedict['OutputFileRoot'] = 'beam%s_Imp%s_FollRig%s_NE%0.2d' \
                 %( self.xbinput.BConds, self.xbopts.ImpStart, 
                    self.xbopts.FollowerForceRig.value, 
                    self.xbinput.NumElems)
        xbout=Sol312(self.xbinput, self.xbopts, SaveDict=self.savedict)

        # Ref solution run with sol 312 and no rotation of FoR A 
        RefTipInitial_G=np.array([0.6614130665000533,0.0,-0.6958598401759009])
        RefTipFinal_G=np.array([  0.6614130665000557,0.0,-0.6958598401758554])
        CaoIni=XbeamLib.RotCRV(self.xbinput.PsiA_G)
        CaoFin=XbeamLib.Rot(xbout.QuatList[-1]).transpose()
        
        # Project tip position at the end of simulation in FoR G
        TipGIni = np.dot(CaoIni,xbout.DynOut[-1,:])
        ErIni = np.linalg.norm(TipGIni-RefTipInitial_G)  
        TipGFin = np.dot(CaoFin,xbout.DynOut[-1,:])
        ErFin = np.linalg.norm(TipGFin-RefTipFinal_G)  

        # checks
        self.assertTrue(ErIni<TOL,msg='Tip initial position error of %.3e above '
                                                'tolerance %.2e!' %(ErIni,TOL) )
        self.assertTrue(ErFin<TOL, msg='Tip final position error of %.3e above '
                                                'tolerance %.2e!' %(ErFin,TOL) )

        return xbout


    def test312_ImpFalse_FollRigFalse(self):

        self.xbopts.ImpStart=False
        self.xbopts.FollowerForceRig.value=False
        self.savedict['OutputFileRoot'] = 'beam%s_Imp%s_FollRig%s_NE%0.2d' \
                 %( self.xbinput.BConds, self.xbopts.ImpStart, 
                    self.xbopts.FollowerForceRig.value, 
                    self.xbinput.NumElems)
        xbout=Sol312(self.xbinput, self.xbopts, SaveDict=self.savedict)

        # Ref solution run with sol 312 and no rotation of FoR A 
        RefTipInitial_G=np.array([0.6614130665000533,0.0,-0.6958598401759009])
        RefTipFinal_G=np.array([  0.6614130665000557,0.0,-0.6958598401758554])
        CaoIni=XbeamLib.RotCRV(self.xbinput.PsiA_G)
        CaoFin=XbeamLib.Rot(xbout.QuatList[-1]).transpose()
                
        # Project tip position at the end of simulation in FoR G
        TipGIni = np.dot(CaoIni,xbout.DynOut[-1,:])
        ErIni = np.linalg.norm(TipGIni-RefTipInitial_G)  
        TipGFin = np.dot(CaoFin,xbout.DynOut[-1,:])
        ErFin = np.linalg.norm(TipGFin-RefTipFinal_G)  

        # checks
        self.assertTrue(ErIni<TOL_rig_false,msg='Tip initial position error of '
                           '%.3e above tolerance %.2e!' %(ErIni,TOL_rig_false) )
        self.assertTrue(ErFin<TOL_rig_false, msg='Tip final position error of '
                           '%.3e above tolerance %.2e!' %(ErFin,TOL_rig_false) )

        return xbout



    def test312_SS_ImpTrue_FollRigTrue(self):
        '''
        Dynamics of a simply supported beam on both sides. The beam is initially
        straight (undeformed) and vibrates due to its own weight.
        Warning: this tests compares the average mid-point position against a 
        static solution (the system is conservative) and the vibration frequency
        against linear theory (Euler-Bernoulli beam)
        '''

        self.xbinput.dt = 0.002
        self.xbinput.tfin = 2.0
        self.xbinput.PsiA_G=30.*np.pi/180.*np.array([1,0,0]) # rotation to FoR A
        self.xbinput.BConds='CC'
        self.xbopts.ImpStart=True
        self.xbopts.FollowerForceRig.value=True
        self.xbopts.PrintInfo = ct.c_bool(True)
        self.savedict['OutputFileRoot'] = 'beam%s_Imp%s_FollRig%s_NE%0.2d' \
                 %( self.xbinput.BConds, self.xbopts.ImpStart, 
                    self.xbopts.FollowerForceRig.value, 
                    self.xbinput.NumElems)
        xbout=Sol312(self.xbinput, self.xbopts, SaveDict=self.savedict)

        self.savedict['OutputFileRoot']=self.savedict['OutputFileRoot']+'_static'
        self.xbopts.Solution = ct.c_int(112)
        self.xbopts.MinDelta = ct.c_double(1e-5)
        xbout_static=Sol112(self.xbinput, self.xbopts, SaveDict=self.savedict)
        PosDefSta=xbout_static.PosDeforStatic

        # Post-processing
        dt = self.xbinput.dt
        tv = xbout.Time
        Vrel = xbout.Vrel
        #F90case:
        #Xi = PyLibs.numerics.integr.rotations(Vrel[:,3:],tv,
        #                                   xi0=np.array([1.0,0.0,0.0,0.0]))    
        # SolvePy
        Xi=np.array( xbout.QuatList )
        # Convert Translational Velocities from A to G frame
        if self.xbopts.OutInaframe == True:
            for ii in range(len(Xi)):
                Cga = PyBeam.Utils.XbeamLib.Rot(Xi[ii,:]).transpose()
                Vrel[ii,:3] = np.dot(Cga,Vrel[ii,:3])    
                Vrel[ii,3:] = np.dot(Cga,Vrel[ii,3:]) 
        ### computed Deformed Beam at each time-step in global FoR G
        DynOut=xbout.DynOut
        THPosDefGlobal = PyBeam.Utils.PostPr.THPosDefGlobal(DynOut,tv,Vrel,
                                                       set_origin='G',Xi=Xi)
        # compute Nodal velocities at each time-step in global FoR
        THVel=PyBeam.Utils.PostPr.compute_velocity(THPosDefGlobal,tv)

        # check average z position of mid-node
        midnode=int((self.xbinput.NumNodesTot-1)/2)
        ttvec=tv<1.5*self.xbinput.tfin #(tv>0.2)*(tv<0.8*self.xbinput.tfin)
        MidGz_avg=np.average(THPosDefGlobal[ttvec,midnode,2])

        # delta freq of 1/tv[-1]
        frv, cfv = PyLibs.CVP.fourier.fft(tv,THPosDefGlobal[:,midnode,2])
        cfabs = np.abs(cfv)
        iimax = np.argmax(cfabs[1:]) # exclude freq zero
        frv_res = frv[1+iimax]

        # Analytical frequency (Vibration and Shock Handbook, Da Silva)
        lamb=np.pi/self.xbinput.BeamLength
        fn = 0.5/np.pi*lamb**2*np.sqrt( 
                    self.xbinput.BeamStiffness[4,4]/self.xbinput.BeamMass[4,4] )

        # Post-process static
        Cao=XbeamLib.RotCRV(self.xbinput.PsiA_G)
        for ii in range(self.xbinput.NumNodesTot):
            PosDefSta[ii,:]=np.dot(Cao,PosDefSta[ii,:])
        MidGz_static=PosDefSta[midnode,2]

        if self.PLOT:
            ## Plot Velocities
            fig = plt.figure('Mid Velocities')
            ax = fig.add_subplot(111)
            hz, = ax.plot(tv,THVel[:,midnode,2],'b')
            ax.set_xlabel(r'$t$ [s]')
            ax.set_ylabel(r'[m s$^{-1}$]')
            ax.set_xlim(tv[0],tv[-1])
            plt.legend((hz,),(r'$v_Z$',),loc=2)
            fig.savefig(self.savedict['OutputDir']+
                                  self.savedict['OutputFileRoot']+'_tipvel.png')
            ## Plot Tip Position
            fig = plt.figure('Tip Position')
            ax = fig.add_subplot(111)
            hz, = ax.plot(tv,THPosDefGlobal[:,midnode,2],'b')
            ax.set_xlim(tv[0],tv[-1])
            ax.set_xlabel(r'$t$ [s]')
            ax.set_ylabel(r'[m]')
            plt.legend((hz,),(r'$Z$',),loc=2)
            fig.savefig(self.savedict['OutputDir']+'tipopt'+
                                         self.savedict['OutputFileRoot']+'.png')
            ### Snapshots:
            mask = np.in1d(tv, np.arange(0,1.0,0.04))# 25 fpt
            mask[-1]=True
            ttvec = np.where(mask)[0]
            SnapPosDefGlobal_XZ = THPosDefGlobal[ttvec,:,:][:,:,[0,2]]
            fig,ax=PyLibs.plot.dyn.beam_snapshots(SnapPosDefGlobal_XZ,
                                                                 asp_ratio=None)
            ax.set_xlabel(r'$X$ [m]')
            ax.set_ylabel(r'$Z$ [m]')
            fig.savefig(self.savedict['OutputDir']+
                                   self.savedict['OutputFileRoot']+'_snaps.png')
            #plt.show()
            plt.close()

        # asset vibration frequency (tolerance measured on FFT uncertainty and 
        # previous numerical simulations on the same low fidelity time grid)
        Error_num = np.abs(frv_res-18.5) # Hz
        self.assertTrue(Error_num<TOL,msg='Frequency of vibration changed of more'
            'then  %.3e  Hz with respect to previous simulations!'%Error_num)
        #Error_an = np.abs(frv_res/fn-1.0) # Hz
        #self.assertTrue(Error_an<0.05,msg='Frequency of vibration %.3e 5perc '
        #    'away from first natural frequency of vibration!'%Error_an)

        # assert average tip position - this is done only against numerical, as
        # a more careful average would be required (from peak to peak) to match
        # the static displacement
        MidGz_avg_num=-0.0013718020206028 # Cao unit diagonal
        Error_num = np.abs(MidGz_avg/MidGz_avg_num-1.0) 
        self.assertTrue(Error_num<TOL,msg='Average mid position %.3e above the'
            '%.3e tolerance'%(Error_num,TOL))

        return xbout




if __name__=='__main__':

    unittest.main()

    # still in progress
    #T=Test312()
    #T.setUp()
    #xbout=T.test312_SS_ImpTrue_FollRigTrue()





    
