'''
@author: salvatore maraniello
@contact: salvatore.maraniello10@imperial.ac.uk
@date: 16 Feb 2017
@brief: unittest class to test the gravity loads implementation in the PyBeam
        solution methods.
@warning: - Very low fidelity model is used for sol 312
          - Sol112 and Sol112F90 produce different results (within a 0.1% 
          relative error)
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
from PyBeam.Solver.NonlinearStatic import Solve_F90 as Sol112F90
#from PyBeam.Solver.NonlinearDynamic import Solve_Py as Sol312

import PyLibs.numerics.integr
import PyLibs.numerics.diff
import PyLibs.plot.dyn
from PyLibs.plot.shared import fontlabel, params
import PyBeam.Utils.PostPr
import XbeamLib
import PyLibs.CVP.fourier
import lib_fem


TOL=1e-8


# ------------------------------------------------------------------------------


class Test112F90(unittest.TestCase):
    '''
    Each method defined in this class contains a test case.
    @warning: by default, only functions whose name starts with 'test' will be 
    run during testing.

    This test looks predominantly at boundary conditions. The beam deflections
    are computed against reference analytical solutions (self.linbeamunif) from 
    linear theory for small displacements.

    A higher tolerance is used in this test as, especially in the hinged case,
    the beam extensional stiffness affects the results also for low 
    displacements

    A higher tolerance is used in these tests because the nonlinear solution 
    strongly depends on the extensional stiffness EA when the beam is supported
    on both sides while the linear solution neglects this effect.

    A further test with low tolerance is used to check whether the solution
    changes from one implementation to another.

    '''


    def setUp(self):
        ''' Common piece of code run by each test '''        

        self.TOL112=1e-2
        self.PLOT=True

        # Beam solution options
        self.xbopts = DerivedTypes.Xbopts(FollowerForce = ct.c_bool(False),
                                     FollowerForceRig=ct.c_bool(False),
                                     MaxIterations = ct.c_int(20),
                                     PrintInfo = ct.c_bool(True),
                                     OutInaframe = ct.c_bool(True),
                                     NumLoadSteps = ct.c_int(15),
                                     Solution = ct.c_int(112),
                                     MinDelta = ct.c_double(1e-6),
                                     NewmarkDamp = ct.c_double(1e-2))

        # Beam input
        self.xbinput = DerivedTypes.Xbinput(NumNodesElem=3,NumElems=10,g=0.980)
        self.xbinput.BeamLength = 1.0
        self.xbinput.PsiA_G=0.0*np.array([1,0,0])
        # Mass/Stiffness
        mvec=np.array([.15,.15,.15,.1,.001,.001])
        for ii in range(6): self.xbinput.BeamMass[ii,ii]=mvec[ii]
        kvec=np.array([1e5,1e5,1e5,.15,15,15])
        for ii in range(6): self.xbinput.BeamStiffness[ii,ii]=kvec[ii]

        # saving options
        self.savedict=Settings.SaveDict
        self.savedict['OutputDir'] = Settings.SharPyProjectDir + \
                                          'output/tests/PyBeam/NonlinearStatic/'

        return self


    def linbeamunif(self,xv,EI,q,BCs):
        '''
        Linear solution under uniform load
        '''
        L=xv[-1]-xv[0]

        if BCs=='CF': 
            zv = q/(24.*EI) * xv**2 * ( 6.*L**2 - 4.*L*xv + xv**2 )
        elif BCs=='ST' or BCs=='TS': 
            zv = q/(24.*EI) * xv * ( L**3 - 2.*L*xv**2 + xv**3 )
        elif BCs=='CC': 
            zv = q/(24.*EI) * xv**2 * (L-xv)**2
        else: 
            print('Analytical solution not available') 
            zv=xv
        return zv


    def test112_F90_CFlinunif(self):

        self.xbinput.BConds='CF'
        self.savedict['OutputFileRoot'] = 'F90beam%s_NE%0.2d' \
                                 %( self.xbinput.BConds,  self.xbinput.NumElems)
        xbout=Sol112F90(self.xbinput, self.xbopts, SaveDict=self.savedict)

        # Reference solution
        PosDef=xbout.PosDeforStatic
        PosIni=xbout.PosIni
        EI=self.xbinput.BeamStiffness[-1,-1]
        L=self.xbinput.BeamLength
        q=-self.xbinput.g*self.xbinput.BeamMass[0,0]
        if PosDef[0,0]>-1e-15: xvref=PosIni[:,0]
        else: xvref=PosIni[:,0]-PosIni[0,0]
        zvref = self.linbeamunif(xvref,EI,q,self.xbinput.BConds)
        del xvref
        TipGz_ref = zvref[-1]

        # Extract SHARPy solution
        Cao=XbeamLib.RotCRV(self.xbinput.PsiA_G)
        PosDefG = 0.0*PosDef
        for ii in range(self.xbinput.NumNodesTot):
            PosDefG[ii,:]=np.dot(Cao,PosDef[ii,:])
        #midnode=int((self.xbinput.NumNodesTot-1)/2)

        # Extract tip position
        TipGz = PosDefG[-1,2]
        
        # check accuracy against analytical
        Error = np.abs((TipGz-TipGz_ref)/TipGz_ref)  
        self.assertTrue(Error<self.TOL112, msg='Tip position error of %.3e '
                                   'above tolerance %.2e!' %(Error,self.TOL112))

        # check accuracy against previous numerical solution
        TipGz_num = -0.0012267543687054731
        Error_num = np.abs((TipGz-TipGz_num)/TipGz_num)  
        self.assertTrue(Error_num<TOL, msg='Tip position error of %.3e '
                                       'above tolerance %.2e!' %(Error_num,TOL))

        ## Plot Tip Position
        if self.PLOT:
            fig = plt.figure('Deformed Beam')
            ax = fig.add_subplot(111)
            hsharpy, = ax.plot( PosDefG[:,0], PosDefG[:,2] ,'r+')
            href, = ax.plot(PosIni[:,0], zvref,'k')
            ax.set_xlabel(r'$x$ [m]')
            ax.set_ylabel(r'$z$ [m]')
            plt.legend((hsharpy,href,),(r'SHARPy',r'Analytical',),loc=2)
            plt.savefig( self.savedict['OutputDir'] +
                                         self.savedict['OutputFileRoot']+'.png')
            #plt.show()
            plt.close()
 
        return xbout


    def test112_F90_CClinunif(self):

        self.xbinput.BConds='CC'
        self.savedict['OutputFileRoot'] = 'F90beam%s_NE%0.2d' \
                                 %( self.xbinput.BConds,  self.xbinput.NumElems)
        xbout=Sol112F90(self.xbinput, self.xbopts, SaveDict=self.savedict)

        # Reference solution
        midnode=int((self.xbinput.NumNodesTot-1)/2)
        PosDef=xbout.PosDeforStatic
        PosIni=xbout.PosIni
        EI=self.xbinput.BeamStiffness[-1,-1]
        L=self.xbinput.BeamLength
        q=-self.xbinput.g*self.xbinput.BeamMass[0,0]
        if PosDef[0,0]>-1e-15: xvref=PosIni[:,0]
        else: xvref=PosIni[:,0]-PosIni[0,0]
        zvref = self.linbeamunif(xvref,EI,q,self.xbinput.BConds)
        del xvref
        MidGz_ref = zvref[midnode]

        # Extract SHARPy solution
        Cao=XbeamLib.RotCRV(self.xbinput.PsiA_G)
        PosDefG = 0.0*PosDef
        for ii in range(self.xbinput.NumNodesTot):
            PosDefG[ii,:]=np.dot(Cao,PosDef[ii,:])
        
        # Extract tip position
        MidGz = PosDefG[midnode,2]
        
        # check accuracy against analytical
        Error = np.abs((MidGz-MidGz_ref)/MidGz_ref)  
        self.assertTrue(Error<self.TOL112, msg='Mid position error of %.3e '
                                   'above tolerance %.2e!' %(Error,self.TOL112))

        # check accuracy against previous numerical solution
        MidGz_num = -2.5704576527357301e-05
        Error_num = np.abs((MidGz-MidGz_num)/MidGz_num)  
        self.assertTrue(Error_num<TOL, msg='Mid position error of %.3e '
                                       'above tolerance %.2e!' %(Error_num,TOL))

        ## Plot Tip Position
        if self.PLOT:
            fig = plt.figure('Deformed Beam')
            ax = fig.add_subplot(111)
            hsharpy, = ax.plot( PosDefG[:,0], PosDefG[:,2] ,'r+')
            href, = ax.plot(PosIni[:,0], zvref,'k')
            ax.set_xlabel(r'$x$ [m]')
            ax.set_ylabel(r'$z$ [m]')
            plt.legend((hsharpy,href,),(r'SHARPy',r'Analytical',),loc=2)
            plt.savefig( self.savedict['OutputDir'] +
                                         self.savedict['OutputFileRoot']+'.png')
            plt.close()
 
        return xbout


    def test112_F90_TSlinunif(self):

        self.xbinput.BConds='TS'
        self.savedict['OutputFileRoot'] = 'F90beam%s_NE%0.2d' \
                                 %( self.xbinput.BConds,  self.xbinput.NumElems)

        xbout=Sol112F90(self.xbinput, self.xbopts, SaveDict=self.savedict)

        # Reference solution
        midnode=int((self.xbinput.NumNodesTot-1)/2)
        PosDef=xbout.PosDeforStatic
        PosIni=xbout.PosIni
        EI=self.xbinput.BeamStiffness[-1,-1]
        L=self.xbinput.BeamLength
        q=-self.xbinput.g*self.xbinput.BeamMass[0,0]
        if PosDef[0,0]>-1e-15: xvref=PosIni[:,0]
        else: xvref=PosIni[:,0]-PosIni[0,0]
        zvref = self.linbeamunif(xvref,EI,q,'TS')#self.xbinput.BConds)
        del xvref
        MidGz_ref = zvref[midnode]

        # Extract SHARPy solution
        Cao=XbeamLib.RotCRV(self.xbinput.PsiA_G)
        PosDefG = 0.0*PosDef
        for ii in range(self.xbinput.NumNodesTot):
            PosDefG[ii,:]=np.dot(Cao,PosDef[ii,:])
        
        # Extract tip position
        MidGz = PosDefG[midnode,2]
        
        ## Plot Tip Position
        if self.PLOT:
            fig = plt.figure('Deformed Beam')
            ax = fig.add_subplot(111)
            hsharpy, = ax.plot( PosDefG[:,0], PosDefG[:,2] ,'r+')
            href, = ax.plot(PosIni[:,0], zvref,'k')
            ax.set_xlabel(r'$x$ [m]')
            ax.set_ylabel(r'$z$ [m]')
            plt.legend((hsharpy,href,),(r'SHARPy',r'Analytical',),loc=2)
            plt.savefig( self.savedict['OutputDir'] +
                                         self.savedict['OutputFileRoot']+'.png')
            #plt.show()
            plt.close()

        # check accuracy against analytical
        Error = np.abs((MidGz-MidGz_ref)/MidGz_ref)  
        self.assertTrue(Error<self.TOL112, msg='Mid position error of %.3e '
                                   'above tolerance %.2e!' %(Error,self.TOL112))
 
        # check accuracy against previous numerical solution
        MidGz_num = -0.00012752922029878939
        Error_num = np.abs((MidGz-MidGz_num)/MidGz_num)  
        self.assertTrue(Error_num<TOL, msg='Mid position error of %.3e '
                                       'above tolerance %.2e!' %(Error_num,TOL))

        return xbout



    def test112_F90_TS_forced_displacements(self):
        '''
        Simply supported beam on both sides under gravity and forced 
        displacements. The beam properties are modified to obtain better 
        condition numbers while ensuring that shear and extentional deflections 
        are contained
        '''

        self.xbopts.NumLoadSteps=ct.c_int(20)

        kvec=np.array([1e3,5e2,5e2,.3,.015,.015])
        for ii in range(6): self.xbinput.BeamStiffness[ii,ii]=kvec[ii]

        self.xbinput.PsiA_G=0.0*np.pi*30.0/180.*np.array([0,1,0])

        self.xbinput.BConds='TS'
        self.xbinput.g=9.81

        # equivalent load to gravity
        #ds=self.xbinput.BeamLength/(self.xbinput.NumNodesTot-1)
        #dFg=-9.81*self.xbinput.BeamMass[0,0]*ds
        #self.xbinput.ForceStatic[:,2]=dFg
        #self.xbinput.ForceStatic[[0,-1],2]=0.5*dFg  

        # add forced displacement at the tip
        compr_fact=0.2
        self.xbinput.addForcedDisp(
            node=-1, 
            pos=np.array([(1.-compr_fact)*self.xbinput.BeamLength,0.,0.]),
            FoR='A')

        self.savedict['OutputFileRoot']='catenary_beam%s_NE%0.2d_compr%.0fperc'\
                    %(self.xbinput.BConds,self.xbinput.NumElems,100.*compr_fact)

        #xbout=Sol112(self.xbinput, self.xbopts, SaveDict=self.savedict)
        xbout=Sol112F90(self.xbinput, self.xbopts, SaveDict=self.savedict)

        # Extract SHARPy solution
        Cao=XbeamLib.RotCRV(self.xbinput.PsiA_G)
        PosDef=xbout.PosDeforStatic
        PosDefG = 0.0*PosDef
        for ii in range(self.xbinput.NumNodesTot):
            PosDefG[ii,:]=np.dot(Cao,PosDef[ii,:])
        # Extract tip position
        midnode=int((self.xbinput.NumNodesTot-1)/2)
        MidGz = PosDefG[midnode,2]
        # Estimate length of deformed cable (check extensional deform are low)

        #SegLength=np.diff(PosDef.T).T
        Lcable=0.0
        for nn in range(1,self.xbinput.NumNodesTot):
            Lcable+=np.linalg.norm(PosDef[nn,:]-PosDef[nn-1,:])
        
        # check numerics
        shear_param = self.xbinput.BeamStiffness[1,1]/\
                     self.xbinput.BeamStiffness[4,4]* self.xbinput.BeamLength**2


        # ------------------------- derive reference solution (catenary problem)

        # utility functions
        def funycable(xv,C,C1,C2):
            '''Compute cable coordinates'''
            return (C*np.cosh( (xv+C1)/C ) + C2)
        def fundycable(xv,C,C1,C2):
            '''Derivative of cable shape'''
            return (np.sinh( (xv+C1)/C ))
        def funLcable(xv,C,C1,C2):
            '''Integrates for the cable length'''
            dyv=fundycable(xv,C,C1,C2)
            Iv = np.sqrt(1.+dyv**2)
            return scint.trapz(Iv,x=xv)

        # position in FoR G
        posB=xbout.PosIni[0,:]
        posT = self.xbinput.ForcedDisp[-1]['pos']

        if self.xbinput.ForcedDisp[-1]['FoR']=='A':
            posT=np.dot(Cao,posT)
        zB,zT=posB[2],posT[2]

        # cable length
        Ltarget=self.xbinput.BeamLength
        assert Ltarget>np.linalg.norm(posT-posB), 'Cable too short!'

        # set-up nonlinear system of equations
        Xv=np.linspace(posB[0],posT[0],600)
        def ResEval(cv):
            ''' Evaluate residual associated to the constraints that the 
            coefficient of the catenary solution need to satisfy '''
            Res=[ funycable(Xv[0],*cv)-zB,    # left BC
                  funycable(Xv[-1],*cv)-zT,   # right BC
                  funLcable(Xv,*cv)-Ltarget,] # Length constraint
            return Res
        cv0=np.array([0.3382, -0.4, -0.6036]) 
        Sol=scopt.root(fun=ResEval,x0=cv0,jac=False,tol=1e-6)#method='SLSQP',
        cvsol=Sol['x']
        Yvsol=funycable(Xv,*cvsol)

        ## Plot Tip Position
        if self.PLOT:
            fig = plt.figure('Deformed Beam - F90_TS_forced_displacements')
            ax = fig.add_subplot(111)
            hsharpy, = ax.plot( PosDefG[:,0], PosDefG[:,2] ,'r',marker='o',lw=2,
                                                           markevery=(0.05,0.1))
            href, = ax.plot( Xv, Yvsol ,'b',marker='s',markevery=(0.0,0.1))
            #href,=
            ax.set_xlabel(r'$x$ [m]')
            ax.set_ylabel(r'$z$ [m]')
            plt.legend((hsharpy,href,),(r'SHARPy',r'Analytical',),loc=2,)
            plt.savefig( self.savedict['OutputDir'] +
                                         self.savedict['OutputFileRoot']+'.png')
            plt.close()

        # Comparisons
        # Admissible stretch of cable (estimate)
        dLtension_max=0.5*self.xbinput.BeamLength*self.xbinput.g*\
                      self.xbinput.BeamMass[0,0]/self.xbinput.BeamStiffness[0,0]
             
        self.assertTrue(Lcable<self.xbinput.BeamLength+dLtension_max, 
            msg='The numerically estimated length of the cable should be shorter'
                 ' of the original length is extensional deflections are small')
        self.assertTrue( np.abs(Lcable/self.xbinput.BeamLength-1.)<\
                         dLtension_max/self.xbinput.BeamLength, 
                                 msg='The length of the cable changed too much')
        # numerical vs analytical
        if Sol['success']:
            self.assertTrue(
                np.abs(np.min(Yvsol)-MidGz)/self.xbinput.BeamLength<self.TOL112,
                msg='Mid displacements not matching analytical solution')

        return xbout



# ------------------------------------------------------------------------------


class Test112(unittest.TestCase):
    '''
    Each method defined in this class contains a test case.
    @warning: by default, only functions whose name starts with 'test' will be 
    run during testing.

    This test looks predominantly at boundary conditions. The beam deflections
    are computed against reference analytical solutions (self.linbeamunif) from 
    linear theory for small displacements.

    A higher tolerance is used in these tests because the nonlinear solution 
    strongly depends on the extensional stiffness EA when the beam is supported
    on both sides while the linear solution neglects this effect.

    A further test with low tolerance is used to check whether the solution
    changes from one implementation to another.

    '''


    def setUp(self):
        '''
        Common piece of code run by each test
        '''        

        self.TOL112=1e-2
        self.TOL112_num=1e-3 # numerical comparisons against Sol112F90 solution
        self.PLOT=True

        # Beam solution options
        self.xbopts = DerivedTypes.Xbopts(FollowerForce = ct.c_bool(False),
                                     FollowerForceRig=ct.c_bool(False),
                                     MaxIterations = ct.c_int(20),
                                     PrintInfo = ct.c_bool(True),
                                     OutInaframe = ct.c_bool(True),
                                     NumLoadSteps = ct.c_int(10),
                                     Solution = ct.c_int(112),
                                     MinDelta = ct.c_double(1e-6),
                                     NewmarkDamp = ct.c_double(1e-2))

        # Beam input
        self.xbinput = DerivedTypes.Xbinput(NumNodesElem=3,NumElems=10,g=0.980)
        self.xbinput.BeamLength = 1.0
        self.xbinput.PsiA_G=30.0*np.array([1,0,0])

        # Mass/Stiffness
        mvec=np.array([.15,.15,.15,.1,.001,.001])
        for ii in range(6): self.xbinput.BeamMass[ii,ii]=mvec[ii]
        kvec=np.array([1e6,1e6,1e6,15,15,15])
        for ii in range(6): self.xbinput.BeamStiffness[ii,ii]=kvec[ii]

        # saving options
        self.savedict=Settings.SaveDict
        self.savedict['OutputDir'] = Settings.SharPyProjectDir + \
                                          'output/tests/PyBeam/NonlinearStatic/'

        return self


    def linbeamunif(self,xv,EI,q,BCs):
        '''
        Linear solution under uniform load
        '''
        L=xv[-1]-xv[0]

        if BCs=='CF': 
            zv = q/(24.*EI) * xv**2 * ( 6.*L**2 - 4.*L*xv + xv**2 )
        elif BCs=='TS' or BCs=='ST': 
            zv = q/(24.*EI) * xv * ( L**3 - 2.*L*xv**2 + xv**3 )
        elif BCs=='CC': 
            zv = q/(24.*EI) * xv**2 * (L-xv)**2
        else: 
            print('Analytical solution not available') 
            zv=xv
        return zv


    def test112_CFlinunif(self):

        self.xbinput.BConds='CF'
        self.savedict['OutputFileRoot'] = 'beam%s_NE%0.2d' \
                                 %( self.xbinput.BConds,  self.xbinput.NumElems)
        xbout=Sol112(self.xbinput, self.xbopts, SaveDict=self.savedict)

        # Reference solution
        PosDef=xbout.PosDeforStatic
        PosIni=xbout.PosIni
        EI=self.xbinput.BeamStiffness[-1,-1]
        L=self.xbinput.BeamLength
        q=-self.xbinput.g*self.xbinput.BeamMass[0,0]
        if PosDef[0,0]>-1e-15: xvref=PosIni[:,0]
        else: xvref=PosIni[:,0]-PosIni[0,0]
        zvref = self.linbeamunif(xvref,EI,q,self.xbinput.BConds)
        del xvref
        TipGz_ref = zvref[-1]

        # Extract SHARPy solution
        Cao=XbeamLib.RotCRV(self.xbinput.PsiA_G)
        PosDefG = 0.0*PosDef
        for ii in range(self.xbinput.NumNodesTot):
            PosDefG[ii,:]=np.dot(Cao,PosDef[ii,:])
        #midnode=int((self.xbinput.NumNodesTot-1)/2)

        # Extract tip position
        TipGz = PosDefG[-1,2]
        
        # check accuracy against analytical
        Error = np.abs((TipGz-TipGz_ref)/TipGz_ref)  
        self.assertTrue(Error<self.TOL112, msg='Tip position error of %.3e '
                                   'above tolerance %.2e!' %(Error,self.TOL112))

        # check accuracy against previous numerical solution 
        TipGz_num = -0.0012260928687054726
        Error_num = np.abs((TipGz-TipGz_num)/TipGz_num)  
        self.assertTrue(Error_num<TOL, msg='Tip position error of %.3e '
                                       'above tolerance %.2e!' %(Error_num,TOL))

        ## Plot Tip Position
        if self.PLOT:
            fig = plt.figure('Deformed Beam')
            ax = fig.add_subplot(111)
            hsharpy, = ax.plot( PosDefG[:,0], PosDefG[:,2] ,'r+')
            href, = ax.plot(PosIni[:,0], zvref,'k')
            ax.set_xlabel(r'$x$ [m]')
            ax.set_ylabel(r'$z$ [m]')
            plt.legend((hsharpy,href,),(r'SHARPy',r'Analytical',),loc=2)
            plt.savefig( self.savedict['OutputDir'] +
                                         self.savedict['OutputFileRoot']+'.png')
            #plt.show()
            plt.close()
 
        return xbout


    def test112_CClinunif(self):

        self.xbinput.BConds='CC'
        self.savedict['OutputFileRoot'] = 'beam%s_NE%0.2d' \
                                 %( self.xbinput.BConds,  self.xbinput.NumElems)
        xbout=Sol112(self.xbinput, self.xbopts, SaveDict=self.savedict)

        # Reference solution
        midnode=int((self.xbinput.NumNodesTot-1)/2)
        PosDef=xbout.PosDeforStatic
        PosIni=xbout.PosIni
        EI=self.xbinput.BeamStiffness[-1,-1]
        L=self.xbinput.BeamLength
        q=-self.xbinput.g*self.xbinput.BeamMass[0,0]
        if PosDef[0,0]>-1e-15: xvref=PosIni[:,0]
        else: xvref=PosIni[:,0]-PosIni[0,0]
        zvref = self.linbeamunif(xvref,EI,q,self.xbinput.BConds)
        del xvref
        MidGz_ref = zvref[midnode]

        # Extract SHARPy solution
        Cao=XbeamLib.RotCRV(self.xbinput.PsiA_G)
        PosDefG = 0.0*PosDef
        for ii in range(self.xbinput.NumNodesTot):
            PosDefG[ii,:]=np.dot(Cao,PosDef[ii,:])
        
        # Extract tip position
        MidGz = PosDefG[midnode,2]
        
        # check accuracy against analytical
        Error = np.abs((MidGz-MidGz_ref)/MidGz_ref)  
        self.assertTrue(Error<self.TOL112, msg='Mid position error of %.3e '
                                   'above tolerance %.2e!' %(Error,self.TOL112))

        # check accuracy against previous numerical solution 
        MidGz_num = -2.5539140776011136e-05
        Error_num = np.abs((MidGz-MidGz_num)/MidGz_num)  
        self.assertTrue(Error_num<TOL, msg='Mid position error of %.3e '
                                       'above tolerance %.2e!' %(Error_num,TOL))

        ## Plot Tip Position
        if self.PLOT:
            fig = plt.figure('Deformed Beam')
            ax = fig.add_subplot(111)
            hsharpy, = ax.plot( PosDefG[:,0], PosDefG[:,2] ,'r+')
            href, = ax.plot(PosIni[:,0], zvref,'k')
            ax.set_xlabel(r'$x$ [m]')
            ax.set_ylabel(r'$z$ [m]')
            plt.legend((hsharpy,href,),(r'SHARPy',r'Analytical',),loc=2)
            plt.savefig( self.savedict['OutputDir'] +
                                         self.savedict['OutputFileRoot']+'.png')
            plt.close()
 
        return xbout


    def test112_TSlinunif(self):

        self.xbinput.BConds='ST'
        self.savedict['OutputFileRoot'] = 'beam%s_NE%0.2d' \
                                 %( self.xbinput.BConds,  self.xbinput.NumElems)
        xbout=Sol112(self.xbinput, self.xbopts, SaveDict=self.savedict)

        # Reference solution
        midnode=int((self.xbinput.NumNodesTot-1)/2)
        PosDef=xbout.PosDeforStatic
        PosIni=xbout.PosIni
        EI=self.xbinput.BeamStiffness[-1,-1]
        L=self.xbinput.BeamLength
        q=-self.xbinput.g*self.xbinput.BeamMass[0,0]
        if PosDef[0,0]>-1e-15: xvref=PosIni[:,0]
        else: xvref=PosIni[:,0]-PosIni[0,0]
        zvref = self.linbeamunif(xvref,EI,q,'TS')#self.xbinput.BConds)
        del xvref
        MidGz_ref = zvref[midnode]

        # Extract SHARPy solution
        Cao=XbeamLib.RotCRV(self.xbinput.PsiA_G)
        PosDefG = 0.0*PosDef
        for ii in range(self.xbinput.NumNodesTot):
            PosDefG[ii,:]=np.dot(Cao,PosDef[ii,:])
        
        # Extract tip position
        MidGz = PosDefG[midnode,2]
        
        ## Plot Tip Position
        self.PLOT=True 
        if self.PLOT:
            fig = plt.figure('Deformed Beam')
            ax = fig.add_subplot(111)
            hsharpy, = ax.plot( PosDefG[:,0], PosDefG[:,2] ,'r+')
            href, = ax.plot(PosIni[:,0], zvref,'k')
            ax.set_xlabel(r'$x$ [m]')
            ax.set_ylabel(r'$z$ [m]')
            plt.legend((hsharpy,href,),(r'SHARPy',r'Analytical',),loc=2)
            plt.savefig( self.savedict['OutputDir'] +
                                         self.savedict['OutputFileRoot']+'.png')
            plt.close()

        # check accuracy against analytical
        Error = np.abs((MidGz-MidGz_ref)/MidGz_ref)  
        self.assertTrue(Error<self.TOL112, msg='Mid position error of %.3e '
                                   'above tolerance %.2e!' %(Error,self.TOL112))

        # check accuracy against previous numerical solution 
        MidGz_num = -0.00012733255#4801028591
        Error_num = np.abs((MidGz-MidGz_num)/MidGz_num)  
        self.assertTrue(Error_num<1e1*TOL, msg='Mid position error of %.3e '
                                   'above tolerance %.2e!' %(Error_num,1e1*TOL))
        
        return xbout


    def test112_TSlinunif_linearity(self):
        ''' Increase load of factor 10 w.r.t. test112_SSlinunif and displacements
        scale with loads '''

        self.xbinput.PsiA_G=np.zeros((3,))
        LoadFactor=1.
        self.xbinput.BConds='ST'
        self.xbinput.g=LoadFactor*self.xbinput.g
        self.savedict['OutputFileRoot'] = 'beam%s_NE%0.2d_Loadfactor%.3d' \
                     %( self.xbinput.BConds,  self.xbinput.NumElems, LoadFactor)
        xbout=Sol112(self.xbinput, self.xbopts, SaveDict=self.savedict)
        # Extract SHARPy solution
        Cao=XbeamLib.RotCRV(self.xbinput.PsiA_G)
        PosDef=xbout.PosDeforStatic
        PosDefG = 0.0*PosDef
        for ii in range(self.xbinput.NumNodesTot):
            PosDefG[ii,:]=np.dot(Cao,PosDef[ii,:])
        # Extract tip position
        midnode=int((self.xbinput.NumNodesTot-1)/2)
        MidGz = PosDefG[midnode,2]

        # increase load factor
        LoadFactor=10.
        self.xbinput.BConds='ST'
        self.xbinput.g=LoadFactor*self.xbinput.g
        self.savedict['OutputFileRoot'] = 'beam%s_NE%0.2d_Loadfactor%.3d' \
                     %( self.xbinput.BConds,  self.xbinput.NumElems, LoadFactor)
        xbout=Sol112(self.xbinput, self.xbopts, SaveDict=self.savedict)
        # Extract SHARPy solution
        Cao=XbeamLib.RotCRV(self.xbinput.PsiA_G)
        PosDef=xbout.PosDeforStatic
        PosDefG_fact = 0.0*PosDef
        for ii in range(self.xbinput.NumNodesTot):
            PosDefG_fact[ii,:]=np.dot(Cao,PosDef[ii,:])
        # Extract mid position
        MidGz_fact = PosDefG_fact[midnode,2]
        
        ErrorRel = np.abs(MidGz_fact/MidGz/LoadFactor-1.0) 
        TolFact=0.03

        self.assertTrue(ErrorRel<TolFact, 
            msg='Factor of displacements increase %.3f is not within %.2f perc '
            'of the prescribed load factor %.1f!' %(ErrorRel,TolFact,LoadFactor))

        return xbout



    def test112_TS_forced_displacements(self):
        '''
        Simply supported beam on both sides under gravity and forced 
        displacements. The beam properties are modified to obtain better 
        condition numbers while ensuring that shear and extentional deflections 
        are contained
        '''

        self.xbopts.NumLoadSteps=ct.c_int(20)

        kvec=np.array([1e3,5e2,5e2,.3,.015,.015])
        for ii in range(6): self.xbinput.BeamStiffness[ii,ii]=kvec[ii]
        self.xbinput.PsiA_G=np.zeros((3,))
        self.xbinput.BConds='TS'
        self.xbinput.g=9.81

        # add forced displacement at the tip
        compr_fact=0.2
        self.xbinput.addForcedDisp(
            node=-1, 
            pos=np.array([(1.-compr_fact)*self.xbinput.BeamLength,0.,0.]),
            FoR='A')

        self.savedict['OutputFileRoot']='catenary_beam%s_NE%0.2d_compr%.0fperc'\
                    %(self.xbinput.BConds,self.xbinput.NumElems,100.*compr_fact)


        xbout=Sol112(self.xbinput, self.xbopts, SaveDict=self.savedict)
        # Extract SHARPy solution
        Cao=XbeamLib.RotCRV(self.xbinput.PsiA_G)
        PosDef=xbout.PosDeforStatic
        PosDefG = 0.0*PosDef
        for ii in range(self.xbinput.NumNodesTot):
            PosDefG[ii,:]=np.dot(Cao,PosDef[ii,:])
        # Extract tip position
        midnode=int((self.xbinput.NumNodesTot-1)/2)
        MidGz = PosDefG[midnode,2]
        # Estimate length of deformed cable (check extensional deform are low)

        #SegLength=np.diff(PosDef.T).T
        Lcable=0.0
        for nn in range(1,self.xbinput.NumNodesTot):
            Lcable+=np.linalg.norm(PosDef[nn,:]-PosDef[nn-1,:])
        
        # check numerics
        shear_param = self.xbinput.BeamStiffness[1,1]/\
                     self.xbinput.BeamStiffness[4,4]* self.xbinput.BeamLength**2
        codK = np.linalg.cond(xbout.K)

        # ------------------------- derive reference solution (catenary problem)

        # utility functions
        def funycable(xv,C,C1,C2):
            '''Compute cable coordinates'''
            return (C*np.cosh( (xv+C1)/C ) + C2)
        def fundycable(xv,C,C1,C2):
            '''Derivative of cable shape'''
            return (np.sinh( (xv+C1)/C ))
        def funLcable(xv,C,C1,C2):
            '''Integrates for the cable length'''
            dyv=fundycable(xv,C,C1,C2)
            Iv = np.sqrt(1.+dyv**2)
            return scint.trapz(Iv,x=xv)

        # position in FoR G
        posB=xbout.PosIni[0,:]
        posT = self.xbinput.ForcedDisp[-1]['pos']
        if self.xbinput.ForcedDisp[-1]['FoR']=='A':
            posT=np.dot(Cao,posT)
        zB,zT=posB[2],posT[2]

        # cable length
        Ltarget=self.xbinput.BeamLength
        assert Ltarget>np.linalg.norm(posT-posB), 'Cable too short!'

        # set-up nonlinear system of equations
        Xv=np.linspace(posB[0],posT[0],600)
        def ResEval(cv):
            ''' Evaluate residual associated to the constraints that the 
            coefficient of the catenary solution need to satisfy '''
            Res=[ funycable(Xv[0],*cv)-zB,    # left BC
                  funycable(Xv[-1],*cv)-zT,   # right BC
                  funLcable(Xv,*cv)-Ltarget,] # Length constraint
            return Res
        cv0=np.array([0.3382, -0.4, -0.6036]) 
        Sol=scopt.root(fun=ResEval,x0=cv0,jac=False,tol=1e-6)#method='SLSQP',
        cvsol=Sol['x']
        Yvsol=funycable(Xv,*cvsol)

        ## Plot Tip Position
        if self.PLOT:
            fig = plt.figure('Deformed Beam - TS_forced_displacements')
            ax = fig.add_subplot(111)
            hsharpy, = ax.plot( PosDefG[:,0], PosDefG[:,2] ,'r',marker='o',lw=2,
                                                           markevery=(0.05,0.1))
            href, = ax.plot( Xv, Yvsol ,'b',marker='s',markevery=(0.0,0.1))
            #href,=
            ax.set_xlabel(r'$x$ [m]')
            ax.set_ylabel(r'$z$ [m]')
            plt.legend((hsharpy,href,),(r'SHARPy',r'Analytical',),loc=2,)
            plt.savefig( self.savedict['OutputDir'] +
                                         self.savedict['OutputFileRoot']+'.png')
            #plt.show()
            plt.close()

        # Comparisons
        # Admissible stretch of cable (estimate)
        dLtension_max=0.5*self.xbinput.BeamLength*self.xbinput.g*\
                      self.xbinput.BeamMass[0,0]/self.xbinput.BeamStiffness[0,0]
        self.assertTrue(Lcable<self.xbinput.BeamLength+dLtension_max, 
            msg='The numerically estimated length of the cable should be shorter'
                 ' of the original length is extensional deflections are small')
        self.assertTrue( np.abs(Lcable/self.xbinput.BeamLength-1.)<\
                         dLtension_max/self.xbinput.BeamLength, 
                                 msg='The length of the cable changed too much')
        # numerical vs analytical
        if Sol['success']:
            self.assertTrue(
                np.abs(np.min(Yvsol)-MidGz)/self.xbinput.BeamLength<self.TOL112,
                msg='Mid displacements not matching analytical solution')

        return xbout




    def test112_TS_forced_displacements_Gdef(self):
        '''
        Simply supported beam on both sides under gravity and forced 
        displacements. The beam properties are modified to obtain better 
        condition numbers while ensuring that shear and extentional deflections 
        are contained

        Checks for general problem is position of boundaries defined in FoR G
        '''

        self.xbopts.NumLoadSteps=ct.c_int(40)

        kvec=np.array([1e3,5e2,5e2,.3,.0015,.0015])
        for ii in range(6): self.xbinput.BeamStiffness[ii,ii]=kvec[ii]
        mvec=np.array([0.15,0.15,0.15,.1,1e-3,1e-3])

        self.xbinput.PsiA_G=np.zeros((3,))
        self.xbinput.BConds='TS'
        self.xbinput.g=9.81

        # add forced displacement at the tip
        self.xbinput.addForcedDisp(
            node=-1, 
            pos=np.array([0.8*self.xbinput.BeamLength,
                          0.0,
                          0.2*self.xbinput.BeamLength]),
            FoR='A')

        self.savedict['OutputFileRoot']='catenary_beam%s_NE%0.2d_random_pos'\
                                    %(self.xbinput.BConds,self.xbinput.NumElems)


        xbout=Sol112(self.xbinput, self.xbopts, SaveDict=self.savedict)
        # Extract SHARPy solution
        Cao=XbeamLib.RotCRV(self.xbinput.PsiA_G)
        PosDef=xbout.PosDeforStatic
        PosDefG = 0.0*PosDef
        for ii in range(self.xbinput.NumNodesTot):
            PosDefG[ii,:]=np.dot(Cao,PosDef[ii,:])
        # Extract tip position
        midnode=int((self.xbinput.NumNodesTot-1)/2)
        MidGz = PosDefG[midnode,2]
        # Estimate length of deformed cable (check extensional deform are low)
        Lcable=0.0
        for nn in range(1,self.xbinput.NumNodesTot):
            Lcable+=np.linalg.norm(PosDef[nn,:]-PosDef[nn-1,:])
        
        # check numerics
        shear_param = self.xbinput.BeamStiffness[1,1]/\
                     self.xbinput.BeamStiffness[4,4]* self.xbinput.BeamLength**2
        codK = np.linalg.cond(xbout.K)
        print('Shear parameter: %.2e'%shear_param)
        print('Stiffness condition number: %.2e'%codK)
        print('MidGz=%f'%MidGz)
        print('Estimated length of cable: %f'%Lcable)

        # ------------------------- derive reference solution (catenary problem)

        # utility functions
        def funycable(xv,C,C1,C2):
            '''Compute cable coordinates'''
            return (C*np.cosh( (xv+C1)/C ) + C2)
        def fundycable(xv,C,C1,C2):
            '''Derivative of cable shape'''
            return (np.sinh( (xv+C1)/C ))
        def funLcable(xv,C,C1,C2):
            '''Integrates for the cable length'''
            dyv=fundycable(xv,C,C1,C2)
            Iv = np.sqrt(1.+dyv**2)
            return scint.trapz(Iv,x=xv)

        # position in FoR G
        posB=xbout.PosIni[0,:]
        posT = self.xbinput.ForcedDisp[-1]['pos']
        if self.xbinput.ForcedDisp[-1]['FoR']=='A':
            posT=np.dot(Cao,posT)
        zB,zT=posB[2],posT[2]

        # cable length
        Ltarget=self.xbinput.BeamLength
        assert Ltarget>np.linalg.norm(posT-posB), 'Cable too short!'

        # set-up nonlinear system of equations
        Xv=np.linspace(posB[0],posT[0],600)
        def ResEval(cv):
            ''' Evaluate residual associated to the constraints that the 
            coefficient of the catenary solution need to satisfy '''
            Res=[ funycable(Xv[0],*cv)-zB,    # left BC
                  funycable(Xv[-1],*cv)-zT,   # right BC
                  funLcable(Xv,*cv)-Ltarget,] # Length constraint
            return Res
        cv0=np.array([0.3382, -0.4, -0.6036])
        Sol=scopt.root(fun=ResEval,x0=cv0,jac=False,tol=1e-6)#method='SLSQP',
        cvsol=Sol['x']
        Yvsol=funycable(Xv,*cvsol)

        ## Plot Tip Position
        if self.PLOT:
            fig = plt.figure('Deformed Beam - TS_forced_displacements_Gdef')
            ax = fig.add_subplot(111)
            hsharpy, = ax.plot( PosDefG[:,0], PosDefG[:,2] ,'r',marker='o',lw=2,
                                                           markevery=(0.05,0.1))
            href, = ax.plot( Xv, Yvsol ,'b',marker='s',markevery=(0.0,0.1))
            #href,=
            ax.set_xlabel(r'$x$ [m]')
            ax.set_ylabel(r'$z$ [m]')
            ax.set_xlim((-.1,.9))
            ax.set_ylim((-.5,.5))
            plt.legend((hsharpy,href,),(r'SHARPy',r'Analytical',),loc=2,)
            plt.savefig( self.savedict['OutputDir'] +
                                         self.savedict['OutputFileRoot']+'.png')
            #plt.show()
            plt.close()


        # ----------------------------------
        # reshape CRV matrix
        PsiDef=PyBeam.Utils.PostPr.reshape_PsiMat(xbout.PsiDeforStatic)
        Fval=np.zeros((3,))
        Dval=np.zeros((3,))
        #embed()
        # weigths to be multiplied by nodal values
        Fval,Dval=lib_fem.fem_1d_shapefun(numnodeelem=3, z=0.7)
        ####################################

        # Comparisons
        # Admissible stretch of cable (estimate)
        dLtension_max=0.5*self.xbinput.BeamLength*self.xbinput.g*\
                      self.xbinput.BeamMass[0,0]/self.xbinput.BeamStiffness[0,0]
        self.assertTrue(Lcable<self.xbinput.BeamLength+dLtension_max, 
            msg='The numerically estimated length of the cable should be shorter'
                 ' of the original length is extensional deflections are small')
        self.assertTrue( np.abs(Lcable/self.xbinput.BeamLength-1.)<\
                         dLtension_max/self.xbinput.BeamLength, 
                                 msg='The length of the cable changed too much')
        # numerical vs analytical
        #embed()
        if Sol['success']:
            self.assertTrue(
                np.abs(np.min(Yvsol)-np.min(PosDefG[:,2]))/\
                self.xbinput.BeamLength<self.TOL112,
                msg='Mid displacements not matching analytical solution')

        return xbout



class Test112_Cardona(unittest.TestCase):
    '''
    This test looks reproduces the example is sec 6.9 of Cardona and looks at
    shear locking. The beam tip deflections and rotations are computed
    '''


    def setUp(self):
        '''
        Common piece of code run by each test
        '''        
        self.PLOT=True
        self.tol=5e-3

        # Beam solution options
        self.xbopts = DerivedTypes.Xbopts(
                                     FollowerForce = ct.c_bool(False),
                                     FollowerForceRig=ct.c_bool(False),
                                     MaxIterations = ct.c_int(20),
                                     PrintInfo = ct.c_bool(True),
                                     OutInaframe = ct.c_bool(True),
                                     NumLoadSteps = ct.c_int(10),
                                     Solution = ct.c_int(112),
                                     MinDelta = ct.c_double(1e-6),
                                     NewmarkDamp = ct.c_double(1e-2))
        # saving options
        self.savedict=Settings.SaveDict
        self.savedict['OutputDir'] = Settings.SharPyProjectDir + \
                                          'output/tests/PyBeam/NonlinearStatic/'

        return self



    def test112_Cardona_shear_locking(self):
        '''
        For physical beams (in which the shear factor is not too high), using
        quadratic elements avoids shear locking even with very low 
        discretisations. Residual bending flexibility corrections are also not 
        required.
        With reference to Gerardin and Cardona, Tab.1 to 4, results computed in 
        this test should therefore show better convergence properties then those 
        from Gerardin and Cardona. 
        '''
        Nelvec=[1,2,4,8,16,]
        Ntot=len(Nelvec)
        zvec=np.zeros((Ntot,2))
        fivec=np.zeros((Ntot,2))
        Cvec = np.zeros((Ntot,2))
        Ez, Efi=np.zeros((Ntot,2)),np.zeros((Ntot,2))

        # Mass/Stiffness (Cardona)
        mvec=np.array([1.,1.,1.,1.,1.,1.])
        kvec=np.array([5e8,3.231e8,3.231e8,1e7,9.345e6,9.345e6])
        Ftip=[600.,600000.0]
        tip_ref = [2.681900e-3,2.159081]
        fi_ref =  [-8.025680e-4,-6.720002e-1]

        for ff in range(2):
            for nn in range(Ntot):
                Nel=Nelvec[nn]
                # Beam input
                self.xbinput = DerivedTypes.Xbinput(NumNodesElem=3,NumElems=Nel,
                                                    g=0.0, BConds='CF')
                self.xbinput.BeamLength = 5.0
                self.xbinput.PsiA_G=0.0*np.array([1,0,0]) # rotations will be wrong 
                #                                           if this is non-zero
                for ii in range(6): self.xbinput.BeamMass[ii,ii]=mvec[ii]
                for ii in range(6): self.xbinput.BeamStiffness[ii,ii]=kvec[ii]
                self.xbinput.ForceStatic[-1,2]=Ftip[ff]

                self.savedict['OutputFileRoot'] = 'Cardona_beam%s_NE%0.2d' \
                                     %( self.xbinput.BConds,  self.xbinput.NumElems)
                xbout=Sol112(self.xbinput, self.xbopts, SaveDict=self.savedict)

                # Extract SHARPy solution
                Cao=XbeamLib.RotCRV(self.xbinput.PsiA_G)
                PosDefG = 0.0*xbout.PosDeforStatic
                for ii in range(self.xbinput.NumNodesTot):
                    PosDefG[ii,:]=np.dot(Cao,xbout.PosDeforStatic[ii,:])
                # Extract tip position
                zvec[nn,ff] = PosDefG[-1,2]
                fivec[nn,ff] = xbout.PsiDeforStatic[-1,1,1]
                Cvec[nn,ff] = np.linalg.cond(xbout.K)

            Ez[:,ff]=np.abs(zvec[:,ff]/tip_ref[ff]-1.0)
            Efi[:,ff]=np.abs(fivec[:,ff]/fi_ref[ff]-1.0)

        shear_param = kvec[1]/kvec[4]*self.xbinput.BeamLength**2
        #embed()
        self.assertTrue(np.any(Ez>TOL),msg='The maximum tip displacement relative '
            'error, %.3e, is above the %.3e tolerance'%(np.max(Ez),self.tol))
        self.assertTrue(np.max(Efi)>TOL,msg='The maximum tip rotation relative '
            'error, %.3e, is above the %.3e tolerance'%(np.max(Efi),self.tol))

        return xbout



if __name__=='__main__':


    
    

    
    #T=Test112()
    #T.setUp()
    #xbout=T.test112_CFlinunif()
    #xbout=T.test112_CClinunif()
    #xbout=T.test112_TSlinunif()
    #xbout=T.test112_TSlinunif_linearity()
    #T.test112_TS_forced_displacements()
    #T.test112_TS_forced_displacements_Gdef()

    T=Test112F90()
    T.setUp()
    #xbout=T.test112_F90_CFlinunif()
    #xbout=T.test112_F90_CClinunif()
    #xbout=T.test112_F90_TSlinunif()
    T.test112_F90_TS_forced_displacements()
    #T=Test112_Cardona()
    #T.setUp()
    #xbout=T.test112_Cardona_shear_locking()

    # still in progress
    #T=Test312()
    #T.setUp()
    #xbout=T.test312_TS_ImpTrue_FollRigTrue()




    unittest.main()
    
