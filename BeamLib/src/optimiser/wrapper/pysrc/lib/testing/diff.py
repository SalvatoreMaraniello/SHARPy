'''
Created on 30 Oct 2014

@author: sm6110

Test Post-Process methods related to lib.diff package.

The test is based on Qiqi validation case for hinge model.

'''


import numpy as np
import matplotlib.pyplot as plt


codedir = '/home/sm6110/git/SHARPy/BeamLib/src/optimiser/wrapper'
import sys
sys.path.append(codedir)
sys.path.append(codedir+'/pysrc')

import shared

from xbcomponent import XBeamSolver
import lib.postpr, lib.plot, lib.read
import PyBeam.Utils.XbeamLib

import lib.plot.dyn

#----------------------------------------------------------------- Set Testcase
root='/home/sm6110/git/SHARPy_studies/GEBM/20141024_qiqi_validation'
h5file=root+'/res_hinge/Sol912_NE10_NN3_Ttot_ 1s_dt0.010s_uniform_Newmark_1.0e-01_Mnode_0.0e+00_Iyy_1.0e-07_000.h5'
#h5file='./res_hinge/Sol912_NE30_NN3_Ttot_ 1s_dt0.001s_uniform_Newmark_1.0e-03_Mnode_0.0e+00_Iyy_1.0e-07_000.h5'

# --------------------------------------------------------------------------- #

xbsol = lib.read.h5comp(h5file)

print 'Few Options Summary:'
print 'beam initial axis: ', xbsol.BeamAxis


# a frame origin
aOrigin = lib.postpr.integrate_function(xbsol.RefVel[:,:3],xbsol.Time)
# compute quaternions associated with a frame rotation
Xi = lib.postpr.integrate_rotations(xbsol.RefVel[:,3:],xbsol.Time,method='2tay',dOmega=xbsol.RefVelDot[:,3:] )
# compute position of x axis
xth = lib.postpr.vector_time_history(Xi, vec0=np.array([1,0,0]))
#compute position of beam mean axis
baxth = xbsol.BeamLength1*lib.postpr.vector_time_history(Xi, vec0=xbsol.BeamAxis)
#computed Deformed Beam at each time-step in global FoR
THPosDefGlobal = lib.postpr.THPosDefGlobal(xbsol.DynOut,xbsol.RefVel,xbsol.Time,set_origin='G')


# compute Nodal velocities at each time-step in global FoR
#TipVel=f_xyz_vel(DynOut,RefVel,Time,NNode=-1,TTime=-1,dir='x',FoR='G'):
THVelold=lib.postpr.compute_velocity_Iord(THPosDefGlobal,xbsol.Time)
THVelnew=lib.postpr.compute_velocity(THPosDefGlobal,xbsol.Time)

## Plot Velocities
# this is a preparation of optimisation exercise. We compute the velocity 
# of the whole beam with old and new method
fig = plt.figure('tip velocities - old/new method')
ax = fig.add_subplot(111)
hxold, = ax.plot(xbsol.Time[1:],THVelold[:,-1,0],'0.6')
hyold, = ax.plot(xbsol.Time[1:],THVelold[:,-1,2],'k')
hxnew, = ax.plot(xbsol.Time,THVelnew[:,-1,0],'b')
hynew, = ax.plot(xbsol.Time,THVelnew[:,-1,2],'r')
#ax.set_xlim(0,0.7)
#x.set_ylim(-1.1,1.1)
ax.set_xlabel('t')
ax.set_ylabel('velocity')
plt.legend([hxold,hyold,hxnew,hynew],('x old','z old','x new','z new'),2)
#fig.savefig(figfold+'tipvel'+figname+'.png')
#fig.savefig(figfold+'tipvel'+figname+'.pdf')
plt.show()




## Plot Velocities
# this is a preparation of optimisation exercise. We compute the velocity 
# of the whole beam with old and new method
fig = plt.figure('node 15 velocities - old/new method')
ax = fig.add_subplot(111)
hxold, = ax.plot(xbsol.Time[1:],THVelold[:,15,0],'0.6')
hyold, = ax.plot(xbsol.Time[1:],THVelold[:,15,2],'k')
hxnew, = ax.plot(xbsol.Time,THVelnew[:,15,0],'b')
hynew, = ax.plot(xbsol.Time,THVelnew[:,15,2],'r')
#ax.set_xlim(0,0.7)
#x.set_ylim(-1.1,1.1)
ax.set_xlabel('t')
ax.set_ylabel('velocity')
plt.legend([hxold,hyold,hxnew,hynew],('x old','z old','x new','z new'),2)
#fig.savefig(figfold+'n15vel'+figname+'.png')
#fig.savefig(figfold+'n15vel'+figname+'.pdf')
plt.show()


