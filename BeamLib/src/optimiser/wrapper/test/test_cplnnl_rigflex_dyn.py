'''
27 Aug 2015
@author: S. Maraniello
Flexible pendulum test run
'''


import os
codedir = os.environ["SHARPYDIR"]+'/BeamLib/src/optimiser/wrapper'
import sys
sys.path.append(codedir)
sys.path.append(codedir+'/pysrc')
import shared
from xbcomponent import XBeamSolver
import lib.postpr, lib.read, lib.plot.dyn, cost
import PyBeam.Utils.XbeamLib

import numpy as np
import matplotlib.pyplot as plt

from lib.plot.shared import fontlabel, params
plt.rcParams.update(params)


''' --------------------------------------------------------------------- '''


figfold='./output/fig/'
try:
    os.mkdir(figfold)
except:
    pass

xb = lib.read.h5comp('./input/flex_pendulum_DSS.h5')
xb._savedir='./output'
xb.TestCase='flex_pendulum_new'
xb.NewmarkDamp=0.1
xb.Time=np.linspace(0.0,2.0,201)
xb.PrintInfo=True

# define a cost function
def vx_cost(DynOut,RefVel,Time,ForceDynamic):
    vx = cost.f_xyz_vel(DynOut,RefVel,Time,NNode=-1,TTime=-1,dir='x',FoR='G')
    return vx
xb.fargs=['DynOut','RefVel','Time','ForceDynamic']
xb.ffun=vx_cost

# run solver
xb.execute()



''' ------------------------------------------------------------------- '''
# Torque
Torque = xb.ForceDynamic[0,4,:]
print 'Final Tip Leftward Speed (computed): %f' %xb.fval     
print 'Final Tip Leftward Speed (expected): %f' %(-8.562338) # damping 1e-1
#print 'Final Tip Leftward Speed (expected): %f' %(-14.138675141940212) # damping 1e-2



# a frame origin
aOrigin = lib.integr.function(xb.RefVel[:,:3],xb.Time)
# compute quaternions associated with a frame rotation
Xi = lib.integr.rotations(xb.RefVel[:,3:],xb.Time,method='2tay',dOmega=xb.RefVelDot[:,3:] )
# compute position of x axis
xth = lib.postpr.vector_time_history(Xi, vec0=np.array([1,0,0]))
#compute position of beam mean axis
baxth = xb.BeamLength1*lib.postpr.vector_time_history(Xi, vec0=xb.BeamAxis)
#computed Deformed Beam at each time-step in global FoR
THPosDefGlobal = lib.postpr.THPosDefGlobal(xb.DynOut,xb.RefVel,xb.Time,set_origin='G')
# compute Nodal velocities at each time-step in global FoR
THVel=lib.postpr.compute_velocity(THPosDefGlobal,xb.Time)



## Plot Torque
fig=plt.figure('Applied Torque')
ax=fig.add_subplot(111)
ht,=ax.plot(xb.Time,Torque,'0.6')
ax.set_xlabel(r'$t$ [s]',fontsize=fontlabel)
ax.set_ylabel(r'$M_y$ [Nm]',fontsize=fontlabel)
ax.set_xlim(0,2)
ax.set_ylim(-4.0,4.0)
fig.savefig(figfold+'Trq.png')
  

## Plot Tip Position
fig = plt.figure('Tip coordinates (Global FoR)')
ax = fig.add_subplot(111)
hxpos, = ax.plot(xb.Time[1:],THPosDefGlobal[1:,-1,0],'0.6')
hypos, = ax.plot(xb.Time[1:],THPosDefGlobal[1:,-1,2],'k')
ax.set_xlim(0,2.0)
ax.set_ylim(-1.1,1.1)
ax.set_xlabel(r'$t$ [s]',fontsize=fontlabel)
ax.set_ylabel(r'[m]',fontsize=fontlabel)
plt.legend([hxpos,hypos],(r'$x$',r'$z$'),2)
fig.savefig(figfold+'tippos.png')



# Extract Specific time-steps for beam deformed
mask = np.in1d(xb.Time, np.arange(0,2,0.04)) # 25 fps frequency
mask[-1]=True
ttvec = np.where(mask)[0]
jjvec = [0,2]

SnapPosDefGlobal = THPosDefGlobal[ttvec,:,:][:,:,jjvec]
SnapTime = xb.Time[ttvec]
t1, t2 = 0.7, 1.6
SSlist=[]
SSlist.append( np.where(SnapTime<=t1)[0] )
SSlist.append( np.where( (SnapTime>t1) * (SnapTime<=t2))[0] )
SSlist.append( np.where(SnapTime>t2)[0] )


cc=1
for ssvec in SSlist:
    fig,ax=lib.plot.dyn.beam_snapshots(SnapPosDefGlobal[ssvec])
    ax.set_xlim(-1.0,1.1)
    ax.set_ylim(-1.1,0.3)
    ax.set_xlabel(r'$X$ [m]',fontsize=fontlabel)
    ax.set_ylabel(r'$Z$ [m]',fontsize=fontlabel)
    fig.savefig(figfold+'THpart%.1d' %(cc) +'.png')
    cc=cc+1  


os.system('rm openmdao_log.txt')


