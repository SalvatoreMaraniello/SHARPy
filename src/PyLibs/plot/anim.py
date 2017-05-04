'''
Created on 17 Oct 2015

@author: Salvatore Maraniello
@summary: tool for quick results animations.

'''

import os

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation








def aero2dcomp(tv, span, Faero01, Faero02, filename='./foo', savepng=True, 
           fps=30, ND=1.0, Nloops=1, Tdelay=1.0, 
           lineprop = {'color'    : ['b','r'],
                       'linewidth': [2, 2]   }  ):
    '''
@warning: this function requires an update of the input (pass markes and colors for the lines
          you want to plot and generalise)
    '''

    (NT,Nnodes,six)=Faero01.shape
    assert Faero02.shape[0] == NT, ('The number of time-steps in the time histories is not the same!')
    assert Faero02.shape[1] == Nnodes, ('The number of nodes is not the same!')
    assert len(span)==Nnodes, ('The number of nodes is not the same!')


    #-------------------------------------------------------------------- Set-up animation variables

    #total frames required for final video
    NF = fps*ND*tv[-1] + 1
    
    # compute step for mask
    step_ex = float(NT)/float(NF)
    print('theoretical step: %f' %step_ex)
    step = int(np.round(step_ex))
    print('rounded step %f' %step)
    ttvec = range(0,NT,step)
    # (to ensure we always get the last time-step)
    ttvec = np.concatenate((ttvec,np.array([NT-1]))) 
    
    Faero01Maskcheap=Faero01[ttvec,:,:]
    Faero02Maskcheap=Faero02[ttvec,:,:]
    TimeMaskcheap=tv[ttvec]
    NTMaskcheap=len(Faero01Maskcheap)
    
    #interval on video mask
    interv = 1e3*(TimeMaskcheap[1] - TimeMaskcheap[0])*float(ND) #ms
    
    #real slow down
    NDreal = float(ND)*step_ex/float(step)
    print('Input slow down: %f' %ND)
    print('Real slow down: %f' %NDreal)
    
   
    # ------------------------------------------------------------------------------- Prepare figure
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    # add final configuration
    #line01x,       = ax.plot([], [], lw=4, color=lineprop['color'][0], marker='o', linewidth=1)
    #line02x,       = ax.plot([], [], lw=4, color=lineprop['color'][1], marker='o', linewidth=1)
    #line01z,       = ax.plot([], [], lw=2, color=lineprop['color'][0], marker='s', linewidth=3)
    #line02z,       = ax.plot([], [], lw=2, color=lineprop['color'][1], marker='s', linewidth=3)
    
    line01x,       = ax.plot([], [], lw=2, color=lineprop['color'][0],  linewidth=1)
    line02x,       = ax.plot([], [], lw=2, color=lineprop['color'][1],  linewidth=1)
    line01z,       = ax.plot([], [], lw=4, color=lineprop['color'][0],  linewidth=3)
    line02z,       = ax.plot([], [], lw=4, color=lineprop['color'][1],  linewidth=3)
    
    
    time_text   = ax.text(0.02, 0.95, '', transform=ax.transAxes)
    vel_text = ax.text(0.02, 0.90, '', transform=ax.transAxes)
    
    
    #----------------------------------------------------------------------- Create required methods
   
    # initialisation function (plot the background of each frame)
    def initfun():
        line01x.set_data([], [])
        line02x.set_data([], [])
        line01z.set_data([], [])
        line02z.set_data([], [])
        time_text.set_text('')
        vel_text.set_text('')
        return line01x
    
    # animate video
    def animfun(tt):
        
        # extract data
        fx01=Faero01Maskcheap[tt,:,0]
        fz01=Faero01Maskcheap[tt,:,2]
        fx02=Faero02Maskcheap[tt,:,0]
        fz02=Faero02Maskcheap[tt,:,2]         
        
        line01x.set_data(span, fx01)
        line02x.set_data(span, fx02)
        line01z.set_data(span, fz01)
        line02z.set_data(span, fz02)
        
        time_text.set_text(r'$t = %.2f  \ s$' % TimeMaskcheap[tt])

        if savepng:
            os.system('mkdir -p %s_pngdir' %filename )
            fig.savefig('%s_pngdir/snap%.4d.png'%(filename,tt),format='png')
            fig.savefig('%s_pngdir/snap%.4d.pdf'%(filename,tt),format='pdf')
        return line01x,
    
    # create a generator to save a file with Nloops
    def GenFun(NT,NL,fps,Tdel=1.0):
        snaps=[]
        ll=0
        Ndelay=np.round(Tdel*float(fps))
        while ll<NL:
            nn=0
            while nn<NT:
                snaps.append(nn)
                nn+=1
            ll+=1
            if NL>1 and ll<NL: # add delay snaps
                dd=0
                while dd<Ndelay:
                    snaps.append(NT-1)
                    dd+=1
        return snaps 
    
    frames = GenFun(NTMaskcheap,Nloops,fps,Tdelay)
  
    '''
    # -------------------------------------------------------------------------------------- Animate
    ### one loop
    #animcheap = animation.FuncAnimation(fig, animfun, init_func=initfun,
    #                                    frames=NTMaskcheap, interval=interv, 
    #                                    blit=True, repeat=True, repeat_delay=1500)
    ### more loops
    animcheap = animation.FuncAnimation(fig, animfun, 
                                        frames=GenFun(NTMaskcheap,Nloops,fps,Tdelay), 
                                        init_func=initfun,
                                        interval=interv, blit=True, repeat=True)
    
    #plt.show()
    codecname='libx264' # good and small
    #codecname='mpeg4'  # grains are likely to appear
    #codecname='gif'    # fails
    animcheap.save(filename, codec=codecname)#,fps=fps3x) #equivalent, bitrate=148kbps 
    '''

    return ax, fig, animfun, frames, initfun, interv






def beam2d(tv, THPosDef, filename='./foo', savepng=True, 
           fps=30, ND=1.0, Nloops=1, Tdelay=1.0, fig=None, ax=None):
    '''
    @summary: Creates objects required for a 2D animation of the beam deformation
    
    @param tv: time vector used for dynamic simulation of length NT
    @param THPosDef: 3D array containing beam deformed time history in the format:
                         THPosDef[time_step,node_number,(x,y)coordinates]
                     and shape (NT,NumNodesTotal,2)
    @param filename: file name without extension
    @param savepng: if True, a folder containing the png of the snapshots is created
    @param fps: frames per second
    @param ND: slow time factor
    @param Nloops: number of repetitions
    @param time between consecutive loops
    '''

    #-------------------------------------------------------------------- Set-up animation variables

    #total frames required for final video
    NF = fps*ND*tv[-1] + 1
    
    # compute step for mask
    NT=len(THPosDef)
    step_ex = float(NT)/float(NF)
    print('theoretical step: %f' %step_ex)
    step = int(np.round(step_ex))
    print('rounded step %f' %step)
    ttvec = range(0,NT,step)
    # (to ensure we always get the last time-step)
    ttvec = np.concatenate((ttvec,np.array([NT-1]))) 
    
    THPosDefMaskcheap=THPosDef[ttvec,:,:]
    #THVelMaskcheap=THVel[ttvec,:,:]
    TimeMaskcheap=tv[ttvec]
    NTMaskcheap=len(THPosDefMaskcheap)
    
    #interval on video mask
    interv = 1e3*(TimeMaskcheap[1] - TimeMaskcheap[0])*float(ND) #ms
    
    #real slow down
    NDreal = float(ND)*step_ex/float(step)
    print('Input slow down: %f' %ND)
    print('Real slow down: %f' %NDreal)
    
   
    # ------------------------------------------------------------------------------- Prepare figure
    if fig==None and ax==None:
        fig = plt.figure()
        ax = fig.add_subplot(111)

    # add final configuration
    xend=THPosDef[-1,:,0]
    yend=THPosDef[-1,:,1]
    ax.plot(xend,yend,color='0.8',linewidth=2)
    line,       = ax.plot([], [], lw=4, color='0.6')
    time_text   = ax.text(0.02, 0.95, '', transform=ax.transAxes)
    vel_text = ax.text(0.02, 0.90, '', transform=ax.transAxes)
    
    #----------------------------------------------------------------------- Create required methods
   
    # initialisation function (plot the background of each frame)
    def initfun():
        line.set_data([], [])
        time_text.set_text('')
        vel_text.set_text('')
        return line,
    
    # animate video
    def animfun(tt):
        x=THPosDefMaskcheap[tt,:,0]
        y=THPosDefMaskcheap[tt,:,1]
        line.set_data(x, y)
        time_text.set_text(r'$t = %.2f  \ s$' % TimeMaskcheap[tt])
        #vel_text.set_text(r'$v_X = %.1f \ ms^{-1}$' % THVelMaskcheap[tt,-1,0])
        # save all the png
        if savepng:
            os.system('mkdir -p %s_pngdir' %filename )
            fig.savefig('%s_pngdir/snap%.4d.png'%(filename,tt),format='png')
            fig.savefig('%s_pngdir/snap%.4d.pdf'%(filename,tt),format='pdf')
        return line,
    
    # create a generator to save a file with Nloops
    def GenFun(NT,NL,fps,Tdel=1.0):
        snaps=[]
        ll=0
        Ndelay=np.round(Tdel*float(fps))
        while ll<NL:
            nn=0
            while nn<NT:
                snaps.append(nn)
                nn+=1
            ll+=1
            if NL>1 and ll<NL: # add delay snaps
                dd=0
                while dd<Ndelay:
                    snaps.append(NT-1)
                    dd+=1
        return snaps 
    
    frames = GenFun(NTMaskcheap,Nloops,fps,Tdelay)
  
    '''
    # -------------------------------------------------------------------------------------- Animate
    ### one loop
    #animcheap = animation.FuncAnimation(fig, animfun, init_func=initfun,
    #                                    frames=NTMaskcheap, interval=interv, 
    #                                    blit=True, repeat=True, repeat_delay=1500)
    ### more loops
    animcheap = animation.FuncAnimation(fig, animfun, 
                                        frames=GenFun(NTMaskcheap,Nloops,fps,Tdelay), 
                                        init_func=initfun,
                                        interval=interv, blit=True, repeat=True)
    
    #plt.show()
    codecname='libx264' # good and small
    #codecname='mpeg4'  # grains are likely to appear
    #codecname='gif'    # fails
    animcheap.save(filename, codec=codecname)#,fps=fps3x) #equivalent, bitrate=148kbps 
    '''

    return ax, fig, animfun, frames, initfun, interv
    





def beam2dsub(tv, THPosDefList, filename='./foo', savepng=True, 
           fps=30, ND=1.0, Nloops=1, Tdelay=1.0, axList=None, fig=None ):
    '''
    @summary: As per beam2d, but created several subplots or view of the beam

    @param tv: time vector used for dynamic simulation of length NT
    @param THPosDefList: list of 3D array containing beam deformed time history 
        in the format:
            THPosDef[case,time_step,node_number,(x,y)coordinates]
            and shape (Ncases,NT,NumNodesTotal,2)
    @param filename: file name without extension
    @param savepng: if True, a folder containing the png of the snapshots is 
        created
    @param fps: frames per second
    @param ND: slow time factor
    @param Nloops: number of repetitions
    @param time between consecutive loops
    @param fig: figure where to plot
    @param axList: List of axes objects belonging to fig
    '''

    # Define number of plots
    Ncases=len(THPosDefList)
    if axList is not None:
        assert Ncases==len(axList),'The number of axis does not match the '\
                                                        'length of THPosDefList'
    # set standard format
    if fig==None:
        fig = plt.figure('Animation Figure',(Ncases*5,5))
    if axList==None:
        axList=[]
        code=int('%.1d%.1d%.1d' %(1,Ncases,1))
        for ii in range(Ncases):
            axList.append(fig.add_subplot(code))
            code+=1

    #----------------------------------------------- Set-up animation variables

    #total frames required for final video
    NF = fps*ND*tv[-1]+1
    # compute step for mask
    NT=len(tv)
    step_ex = float(NT)/float(NF)
    print('theoretical step: %f' %step_ex)
    step = int(np.round(step_ex))
    print('rounded step %f' %step)
    ttvec = range(0,NT,step)
    # ensure we get the last time-step
    ttvec = np.concatenate((ttvec,np.array([NT-1]))) 
    # reduce arrays size
    THPosDefMaskcheap_List=[]
    for nn in range(Ncases):
        THPosDefMaskcheap_List.append(THPosDefList[nn][ttvec,:,:])
    TimeMaskcheap=tv[ttvec]
    NTMaskcheap=len(THPosDefMaskcheap_List[0])   
    #interval on video mask
    interv = 1e3*(TimeMaskcheap[1] - TimeMaskcheap[0])*float(ND) #ms
    #real slow down
    NDreal = float(ND)*step_ex/float(step)
    print('Input slow down: %f' %ND)
    print('Real slow down: %f' %NDreal)
    
    ##### Prepare figure
    LineList=[]
    for nn in range(Ncases):
        # add final configuration
        #axList[nn].plot(THPosDefList[nn][-1,:,0],THPosDefList[nn][-1,:,1],
        #                                               color='0.8',linewidth=2)
        line, = axList[nn].plot([], [], lw=3, color='k')
        LineList.append(line) 
    time_text=axList[0].text(0.02, 0.95, '', transform=axList[0].transAxes)
    vel_text=axList[0].text(0.02, 0.85, '', transform=axList[0].transAxes)
    
    #### Create required methods
    # initialisation function (plot the background of each frame)
    def initfun():
        for nn in range(Ncases):
            LineList[nn].set_data([], [])
        time_text.set_text('')
        return LineList
    
    # animate video
    def animfun(tt):
        for nn in range(Ncases):
            x=THPosDefMaskcheap_List[nn][tt,:,0]
            y=THPosDefMaskcheap_List[nn][tt,:,1]
            LineList[nn].set_data(x, y)
        time_text.set_text(r'$t = %.2f  \ s$' % TimeMaskcheap[tt])
        # save all the png
        if savepng:
            os.system('mkdir -p %s_pngdir' %filename )
            fig.savefig('%s_pngdir/snap%.4d.png'%(filename,tt),format='png')
            #fig.savefig('%s_pngdir/snap%.4d.pdf'%(filename,tt),format='pdf')
        return LineList
    
    # create a generator to save a file with Nloops
    def GenFun(NT,NL,fps,Tdel=1.0):
        snaps=[]
        ll=0
        Ndelay=np.round(Tdel*float(fps))
        while ll<NL:
            nn=0
            while nn<NT:
                snaps.append(nn)
                nn+=1
            ll+=1
            if NL>1 and ll<NL: # add delay snaps
                dd=0
                while dd<Ndelay:
                    snaps.append(NT-1)
                    dd+=1
        return snaps 
    
    frames = GenFun(NTMaskcheap,Nloops,fps,Tdelay)
  
    return axList, fig, animfun, frames, initfun, interv



 
 
def beam2dcomp(tv, THPosDef01, THPosDef02, filename='./foo', savepng=True, 
           fps=30, ND=1.0, Nloops=1, Tdelay=1.0, 
           lineprop = {'color'    : ['b','r'],
                       'linewidth': [2, 2]   }  ):
    '''
    @summary: Creates objects required for a 2D animation of the beam deformation
    
    @param tv: time vector used for dynamic simulation of length NT
    @param THPosDef01/02: 3D array containing beam deformed time history in the format:
                         THPosDef01[time_step,node_number,(x,y)coordinates]
                     and shape (NT,NumNodesTotal,2)
    @param filename: file name without extension
    @param savepng: if True, a folder containing the png of the snapshots is created
    @param fps: frames per second
    @param ND: slow time factor
    @param Nloops: number of repetitions
    @param time between consecutive loops
    '''

    NT=len(THPosDef01)
    assert len(THPosDef02) == NT, ('The number of time-steps in the time histories is not the same!')

    #-------------------------------------------------------------------- Set-up animation variables

    #total frames required for final video
    NF = fps*ND*tv[-1] + 1
    
    # compute step for mask
    step_ex = float(NT)/float(NF)
    print('theoretical step: %f' %step_ex)
    step = int(np.round(step_ex))
    print('rounded step %f' %step)
    ttvec = range(0,NT,step)
    # (to ensure we always get the last time-step)
    ttvec = np.concatenate((ttvec,np.array([NT-1]))) 
    
    THPosDef01Maskcheap=THPosDef01[ttvec,:,:]
    THPosDef02Maskcheap=THPosDef02[ttvec,:,:]
    TimeMaskcheap=tv[ttvec]
    NTMaskcheap=len(THPosDef01Maskcheap)
    
    #interval on video mask
    interv = 1e3*(TimeMaskcheap[1] - TimeMaskcheap[0])*float(ND) #ms
    
    #real slow down
    NDreal = float(ND)*step_ex/float(step)
    print('Input slow down: %f' %ND)
    print('Real slow down: %f' %NDreal)
    
   
    # ------------------------------------------------------------------------------- Prepare figure
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    # add final configuration
    #xend=THPosDef01[-1,:,0]
    #yend=THPosDef01[-1,:,1]
    #ax.plot(xend,yend,color='0.8',linewidth=2)
    line01,       = ax.plot([], [], lw=4, color=lineprop['color'][0])
    line02,       = ax.plot([], [], lw=4, color=lineprop['color'][1])
    time_text   = ax.text(0.02, 0.95, '', transform=ax.transAxes)
    vel_text = ax.text(0.02, 0.90, '', transform=ax.transAxes)
    
    
    #----------------------------------------------------------------------- Create required methods
   
    # initialisation function (plot the background of each frame)
    def initfun():
        line01.set_data([], [])
        line02.set_data([], [])
        time_text.set_text('')
        vel_text.set_text('')
        return line01
    
    # animate video
    def animfun(tt):
        
        # extract data
        x01=THPosDef01Maskcheap[tt,:,0]
        y01=THPosDef01Maskcheap[tt,:,1]
        x02=THPosDef02Maskcheap[tt,:,0]
        y02=THPosDef02Maskcheap[tt,:,1]         
        
        line01.set_data(x01, y01)
        line02.set_data(x02, y02)
        
        time_text.set_text(r'$t = %.2f  \ s$' % TimeMaskcheap[tt])

        if savepng:
            os.system('mkdir -p %s_pngdir' %filename )
            fig.savefig('%s_pngdir/snap%.4d.png'%(filename,tt),format='png')
            fig.savefig('%s_pngdir/snap%.4d.pdf'%(filename,tt),format='pdf')
        return line01,
    
    # create a generator to save a file with Nloops
    def GenFun(NT,NL,fps,Tdel=1.0):
        snaps=[]
        ll=0
        Ndelay=np.round(Tdel*float(fps))
        while ll<NL:
            nn=0
            while nn<NT:
                snaps.append(nn)
                nn+=1
            ll+=1
            if NL>1 and ll<NL: # add delay snaps
                dd=0
                while dd<Ndelay:
                    snaps.append(NT-1)
                    dd+=1
        return snaps 
    
    frames = GenFun(NTMaskcheap,Nloops,fps,Tdelay)
  
    '''
    # -------------------------------------------------------------------------------------- Animate
    ### one loop
    #animcheap = animation.FuncAnimation(fig, animfun, init_func=initfun,
    #                                    frames=NTMaskcheap, interval=interv, 
    #                                    blit=True, repeat=True, repeat_delay=1500)
    ### more loops
    animcheap = animation.FuncAnimation(fig, animfun, 
                                        frames=GenFun(NTMaskcheap,Nloops,fps,Tdelay), 
                                        init_func=initfun,
                                        interval=interv, blit=True, repeat=True)
    
    #plt.show()
    codecname='libx264' # good and small
    #codecname='mpeg4'  # grains are likely to appear
    #codecname='gif'    # fails
    animcheap.save(filename, codec=codecname)#,fps=fps3x) #equivalent, bitrate=148kbps 
    '''

    return ax, fig, animfun, frames, initfun, interv
    

def beam_proj(tv, THPosDef, filename='./foo', savepng=True, 
              fps=30, ND=1.0, Nloops=1, Tdelay=1.0):
    '''
    @summary: Creates objects required for a 2D animation of the beam deformation. 
              Function similar to beam2d, but creates two sub-plots with xy and yz
              projections of beam deformed
    
    @param tv: time vector used for dynamic simulation of length NT
    @param THPosDef: 3D array containing beam deformed time history in the format:
                         THPosDef[time_step,node_number,coordinates]
                     and shape (NT,NumNodesTotal,3)
    @param filename: file name without extension
    @param savepng: if True, a folder containing the png of the snapshots is created
    @param fps: frames per second
    @param ND: slow time factor
    @param Nloops: number of repetitions
    @param time between consecutive loops
    '''


    #-------------------------------------------------------------------- Set-up animation variables

    #total frames required for final video
    NF = fps*ND*tv[-1] + 1
    
    # compute step for mask
    NT=len(THPosDef)
    step_ex = float(NT)/float(NF)
    print('theoretical step: %f' %step_ex)
    step = int(np.round(step_ex))
    print('rounded step %f' %step)
    ttvec = range(0,NT,step)
    # (to ensure we always get the last time-step)
    ttvec = np.concatenate((ttvec,np.array([NT-1]))) 
    
    THPosDefMaskcheap=THPosDef[ttvec,:,:]
    #THVelMaskcheap=THVel[ttvec,:,:]
    TimeMaskcheap=tv[ttvec]
    NTMaskcheap=len(THPosDefMaskcheap)
    
    #interval on video mask
    interv = 1e3*(TimeMaskcheap[1] - TimeMaskcheap[0])*float(ND) #ms
    
    #real slow down
    NDreal = float(ND)*step_ex/float(step)
    print('Input slow down: %f' %ND)
    print('Real slow down: %f' %NDreal)
    
   
    # ------------------------------------------------------------------------------- Prepare figure
    
    fig = plt.figure()
    
    # add final 0xz configuration
    ax_xz = fig.add_subplot(211)
    xend=THPosDef[-1,:,0]
    zend=THPosDef[-1,:,2]
    ax_xz.plot(xend,zend,color='0.8',linewidth=2)
    line_xz,       = ax_xz.plot([], [], lw=4, color='0.6')
    time_text   = ax_xz.text(0.02, 0.95, '', transform=ax_xz.transAxes)
    #vel_text = ax_xz.text(0.02, 0.90, '', transform=ax_xz.transAxes)

    # add final 0yz configuration
    ax_yz = fig.add_subplot(212)
    yend=THPosDef[-1,:,1]
    zend=THPosDef[-1,:,2]
    ax_yz.plot(yend,zend,color='0.8',linewidth=2)
    line_yz,       = ax_yz.plot([], [], lw=4, color='0.6')
    time_text   = ax_yz.text(0.02, 0.95, '', transform=ax_yz.transAxes)
    #vel_text = ax_xz.text(0.02, 0.90, '', transform=ax_xz.transAxes)    
    
    
    
    #----------------------------------------------------------------------- Create required methods
   
    # initialisation function (plot the background of each frame)
    def initfun():
        line_xz.set_data([], [])
        line_yz.set_data([], [])
        time_text.set_text('')
        #vel_text.set_text('')
        return line_xz, line_yz
    
    # animate video
    def animfun(tt):
        x=THPosDefMaskcheap[tt,:,0]
        y=THPosDefMaskcheap[tt,:,1]
        z=THPosDefMaskcheap[tt,:,2]
        
        line_xz.set_data(x, z)
        line_yz.set_data(y, z)
        
        time_text.set_text(r'$t = %.2f  \ s$' % TimeMaskcheap[tt])
        #vel_text.set_text(r'$v_X = %.1f \ ms^{-1}$' % THVelMaskcheap[tt,-1,0])
        # save all the png
        if savepng:
            os.system('mkdir -p %s_pngdir' %filename )
            fig.savefig('%s_pngdir/snap%.4d.png'%(filename,tt),format='png')
        return line_xz, line_yz
    
    # create a generator to save a file with Nloops
    def GenFun(NT,NL,fps,Tdel=1.0):
        snaps=[]
        ll=0
        Ndelay=np.round(Tdel*float(fps))
        while ll<NL:
            nn=0
            while nn<NT:
                snaps.append(nn)
                nn+=1
            ll+=1
            if NL>1 and ll<NL: # add delay snaps
                dd=0
                while dd<Ndelay:
                    snaps.append(NT-1)
                    dd+=1
        return snaps 
    
    frames = GenFun(NTMaskcheap,Nloops,fps,Tdelay)
  
    '''
    # -------------------------------------------------------------------------------------- Animate
    ### one loop
    #animcheap = animation.FuncAnimation(fig, animfun, init_func=initfun,
    #                                    frames=NTMaskcheap, interval=interv, 
    #                                    blit=True, repeat=True, repeat_delay=1500)
    ### more loops
    animcheap = animation.FuncAnimation(fig, animfun, 
                                        frames=GenFun(NTMaskcheap,Nloops,fps,Tdelay), 
                                        init_func=initfun,
                                        interval=interv, blit=True, repeat=True)
    
    #plt.show()
    codecname='libx264' # good and small
    #codecname='mpeg4'  # grains are likely to appear
    #codecname='gif'    # fails
    animcheap.save(filename, codec=codecname)#,fps=fps3x) #equivalent, bitrate=148kbps 
    '''

    return [ax_xz, ax_yz], fig, animfun, frames, initfun, interv



def GenFun(NT,NL,fps,Tdel=1.0):
    '''
    Create a generator when saving an animation with Nloops.
    
    This is an utility function
    '''
    snaps=[]
    ll=0
    Ndelay=np.round(Tdel*float(fps))
    while ll<NL:
        nn=0
        while nn<NT:
            snaps.append(nn)
            nn+=1
        ll+=1
        if NL>1 and ll<NL: # add delay snaps
            dd=0
            while dd<Ndelay:
                snaps.append(NT-1)
                dd+=1
    return snaps 




if __name__=='__main__':
    

    fps=30 # fps of final video
    ND=2.5   # slow time 3 -> 3 times slower
    Nloops=3
    Tdelay=2.0 # seconds
    
    
    '''
    ax, figxz, animfun, frames, initfun, interv = PyLibs.plot.anim.beam2d(tv, 
                                                    THPosDefGlobal[:,:,[0,2]], 
                                                    filename=saverootanim+'/THxzbeam',
                                                    ND=0.4)
    ax.set_xlim(-17,17)
    ax.set_ylim(-6,6)
    ax.set_xlabel(r'$x \ [m]$',fontsize=fontlabel)
    ax.set_ylabel(r'$z \ [m]$',fontsize=fontlabel)
    ax.set_aspect('equal')
    animcheap = animation.FuncAnimation( figxz, animfun, frames=frames, 
                                         init_func=initfun, interval=interv, 
                                         blit=True, repeat=True )
    animcheap.save(saverootanim+'/THxz.mp4', codec='libx264')
    '''
