'''
Created on 6 Feb 2015

@author: sm6110

Test spline module for spline basis generation

'''

import matplotlib.pyplot as plt
import numpy as np

import sys
sys.path.append('../../..')

import input.load.spline as sp

# ------------------------------------------------------------------------------

colors=['0.6', 'b', 'k', 'r', 'y', 'g', 
               'b+','k+','r+','y+','g+']


def test_zero_order(tv,tc):
    
    NI=len(tc)-1 # number of intervals and number of zero order splines
    NT=len(tv)
    NS=NI
    Smat = np.empty((NS,NT))
   # 
    # compute bases
    for ss in range(NS):
        Smat[ss,:]=sp.zero_base(tv, tc, ss)
    
    # plot results
    leglist=[]
    fig = plt.figure('Bases zero-order')
    ax = fig.add_subplot(111)
    for ss in range(NS):
        leglist.append('%d' %(ss))
        ax.plot(tv,Smat[ss,:],colors[ss])
    ax.set_xlim([tv[0]-0.25*tv[-1], tv[-1]+0.25*tv[-1] ])
    ax.set_ylim([-0.5, 1.5])
    plt.legend(tuple(leglist))
    plt.show()
        

def test_p_order(tv,tc,p):

    tc0 = tc[ (tc>=tv[0]) * (tc<=tv[-1]) ]
    NK = len( tc0 )
    NI = NK-1

    NT=len(tv)
    NS=NI+p
    #Smat = np.empty((NS,NT))
    
    print 'NK = %d' %(NK)
    print 'NI = %d' %(NI)
    print 'NS = %d' %(NS)
    
    Smat=sp.p_basis(tv,tc,p)
    Smat = sp.normalise_uniform_base(p, Smat)
    
    # plot results
    leglist=[]
    fig = plt.figure('Bases %d-th order' %(p))
    ax = fig.add_subplot(111)
    for ss in range(NS):
        leglist.append('%d' %(ss))
        ax.plot(tv,Smat[ss,:],colors[ss])
    ax.set_xlim([tv[0]-0.25*tv[-1], tv[-1]+0.25*tv[-1] ])
    ax.set_ylim([-0.5, 1.5])
    plt.legend(tuple(leglist))
    plt.show()
    
    # check single base knots points
    if p==2 or p==3:
        if p==2:
            fknots=np.array([0, 2.0/3.0, 2.0/3.0, 0])
            fothers=np.array([1./6, 1.0, 1./6.])
            tothers=tc[p:p+p+1]+(tc[1]-tc[0])/2.0
        else: # p=3
            fknots=np.array([0.0, 0.25, 1.0, 0.25, 0.0])
            fothers=1.0
            tothers=tc[p+2]
        ###### evaluate at closes relevant points
        ######tcmid = tc[p:p+p+1]+(tc[1]-tc[0])/2.0   <--- continue from ehre!
 
        fig = plt.figure('Base knots points - after scaling')
        ax = fig.add_subplot(111)
        ax.plot(tv,Smat[p,:],'0.6')
        ax.plot(tc[p:p+p+2],fknots,'ro')
        ax.plot(tothers,fothers,'b+')
        ax.set_xlim([tv[0]-0.25*tv[-1], tv[-1]+0.25*tv[-1] ])
        ax.set_ylim([-0.5, 1.5])
        plt.legend(('%d order base'%(p), 'knots'))
        plt.show()
        
    print 'Maximum value of the base: %f' %(np.max(Smat[p,:]))
    
    return Smat
      
    
       
# ------------------------------------------------------------------------------    
    

def test_spline_series_vec(tv, tcint, p):
    
    NI=len(tcint)-1
    NS=p+NI
    
    '''
    print 'reconstructing all the bases' # ok
    for ss in range(NS):
        Scf=np.zeros((NS,))
        Scf[ss]=1.0
        fv=sp.spline_series_vec(tv, tcint, Scf, p)
        fig = plt.figure('bases function')
        ax = fig.add_subplot(111)
        ax.plot(tv,fv,'0.6',linewidth=3)
        ax.set_xlim([tv[0]-0.25*tv[-1], tv[-1]+0.25*tv[-1] ])
        ax.set_ylim([-0.5, 1.5])
        plt.legend(('base no. %d' %(ss)))
        plt.show()
    '''
    
    '''
    print 'testing random bases'
    Scf=np.ones((NS,))
    Scf[0]=-1.0
    Scf[1]=-1.0
    Scf[2]= 1.0
    Scf[3] = -1.0
    Scf[-1]=5.0
    fv=sp.spline_series_vec(tv, tcint, Scf, p)
    fig = plt.figure('random')
    ax = fig.add_subplot(111)
    ax.plot(tv,fv,'0.6',linewidth=3)
    ax.set_xlim([tv[0]-0.25*tv[-1], tv[-1]+0.25*tv[-1] ])
    plt.show()    
    '''

    print 'testing maximum value'
    Scf=np.zeros((NS,))
    #Scf[2]=0.9
    Scf[p]=1.0
    Scf[p+1]=1.0
    fv=sp.spline_series_vec(tv, tcint, Scf, p)    
    fig = plt.figure('random')
    ax = fig.add_subplot(111)
    ax.plot(tv,fv,'0.6',linewidth=3)
    ax.set_xlim([tv[0]-0.25*tv[-1], tv[-1]+0.25*tv[-1] ])
    plt.show() 
    print 'Max value: %f' %(np.max(fv))
    
    return    
 

def test_change_of_base():
 
    # old base
    NIold=6
    pold=3
    
    # new base
    pnew=3
    NInew=23
    
    # domain in high resolution
    T=1.0
    NT=501
    tv=np.linspace(0,T,NT)    
    
    # set random coefficients for old base
    NSold=pold+NIold 
    tcint_old = np.linspace(0,T,NIold+1)
    Scf_old = np.random.rand(NSold)
    
    # compute coefficients new base
    tcint_new=np.linspace(0,T,NInew)
    Scf_new=sp.change_series_base(tcint_new, pnew, tcint_old, Scf_old, pold, tc_uniform_flag=True)
    
    
    # representation of function
    fvold=sp.spline_series_vec(tv, tcint_old, Scf_old, pold, tc_uniform_flag=True) 
    fvnew=sp.spline_series_vec(tv, tcint_new, Scf_new, pnew, tc_uniform_flag=True) 
        
    # plotting
    fig = plt.figure('change of base')
    ax = fig.add_subplot(111)
    ax.plot(tv,fvold,'r',linewidth=2)
    ax.plot(tv,fvnew,'0.6',linewidth=2)
    plt.legend(('old','new'))
    ax.set_xlim([tv[0]-0.25*tv[-1], tv[-1]+0.25*tv[-1] ])
    plt.show()     
    
     
def return_maxima_2_combinations(tv,tcint,p,pincr):
    NI=len(tcint)-1
    NS=p+NI
    
    print 'testing maximum value'
    Scf=np.zeros((NS,))
    N=51
    ssvals = np.linspace(0,1,N)
    MaxMat=np.empty((N,N))

    val=0.0
    pos=[]
    for ii in range(N):
        for jj in range(N):
            Scf[p]=ssvals[ii]
            Scf[p+pincr]=ssvals[jj]
            fv=sp.spline_series_vec(tv, tcint, Scf, p)
            maxlocal=np.max(fv)
            print 'Max(%f,%f)=%f' %(Scf[p],Scf[p+pincr],maxlocal)
            MaxMat[ii,jj]=maxlocal
            
            if maxlocal>val:
                val=maxlocal
                pos=[]
                pos.append([ii,jj])
            if maxlocal==val:
                val=maxlocal
                pos.append([ii,jj])
     
    # visualise worse case    
    Scf[p]=ssvals[pos[0][0]]
    Scf[p+pincr]=ssvals[pos[0][1]]
    fv=sp.spline_series_vec(tv, tcint, Scf, p)       
    fig = plt.figure('worst case')
    ax = fig.add_subplot(111)
    ax.plot(tv,fv,'0.6',linewidth=3)
    ax.set_xlim([tv[0]-0.25*tv[-1], tv[-1]+0.25*tv[-1] ])
    plt.show() 
    print 'Max value: %f' %(np.max(fv))
    return val, pos
    

def return_maxima_3_combinations(tv,tcint,p,pincr1,pincr2):
    NI=len(tcint)-1
    NS=p+NI
    
    print 'testing maximum value'
    Scf=np.zeros((NS,))
    N=51
    ssvals = np.linspace(0,1,N)
    MaxMat=np.empty((N,N))

    val=0.0
    pos=[]
    for ii in range(N):
        for jj in range(N):
            for kk in range(N):
                Scf[p]=ssvals[ii]
                Scf[p+pincr1]=ssvals[jj]
                Scf[p+pincr2]=ssvals[kk]
                fv=sp.spline_series_vec(tv, tcint, Scf, p)
                maxlocal=np.max(fv)
                print 'Max(%f,%f,%f)=%f' %(Scf[p],Scf[p+pincr1],Scf[p+pincr2],maxlocal)
                MaxMat[ii,jj]=maxlocal
                
                if maxlocal>val:
                    val=maxlocal
                    pos=[]
                    pos.append([ii,jj,kk])
                if maxlocal==val:
                    val=maxlocal
                    pos.append([ii,jj,kk])
     
    # visualise worse case    
    Scf[p]=ssvals[pos[0][0]]
    Scf[p+pincr1]=ssvals[pos[0][1]]
    Scf[p+pincr2]=ssvals[pos[0][2]]
    fv=sp.spline_series_vec(tv, tcint, Scf, p)       
    fig = plt.figure('worst case')
    ax = fig.add_subplot(111)
    ax.plot(tv,fv,'0.6',linewidth=3)
    ax.set_xlim([tv[0]-0.25*tv[-1], tv[-1]+0.25*tv[-1] ])
    plt.show() 
    print 'Max value: %f' %(np.max(fv))
    return val, pos


def spline_relevant_points(p):
    ''' Evaluate spline at knots and at intermediate knots points.  '''
    
    T=p+1
    NI=p+1
    NT=1e2*NI + 1.0
    NS=p+NI
    tv=np.linspace(0,T,NT) # high definition vector for plotting
    tcint = np.linspace(0,T,NI+1)   
    tveval = np.linspace(0,T,2*NI+1)  
    
    # 1 spline only
    Scf=np.zeros((NS,))
    Scf[p]=1.0
    
    fv=sp.spline_series_vec(tv, tcint, Scf, p)  
    fveval=sp.spline_series_vec(tveval, tcint, Scf, p)  
    
    fig = plt.figure('base function')
    ax = fig.add_subplot(111)
    ax.plot(tv,fv,'0.6',linewidth=3)
    ax.plot(tveval,fveval,'r+')
    ax.set_xlim([tv[0]-0.25*tv[-1], tv[-1]+0.25*tv[-1] ])
    plt.show() 
    
    print 'Intermediate points value:'
    print 'knot    value'
    for kk in range(2*NI+1):
        print '%f    %f' %(tveval[kk],fveval[kk])
    
    return None
     

def check_bound_conditions_02spl(tv,tcint,p):
    pincr=1
    NI=len(tcint)-1
    NS=p+NI
    
    print 'testing maximum value'
    Scf=np.zeros((NS,))
    N=51
    ssvals = np.linspace(0,1,N)
    MaxMat=np.empty((N,N))

    val=0.0
    pos=[]
    for ii in range(N):
        for jj in range(N):
            Scf[p]=ssvals[ii]
            Scf[p+pincr]=ssvals[jj]
            if p==3:
                #b2 only
                #if 0.71875*(Scf[p]+Scf[p+1])<=1.0:
                # b1 and b2 conditions
                if 0.71875*(Scf[p]+Scf[p+1])<=1.0 and (Scf[p]+0.25*Scf[p+1])<=1.0 and (Scf[p+1]+0.25*Scf[p])<=1.0:
                #b1 only
                #if (Scf[p]+0.25*Scf[p+1])<=1.0 and (Scf[p+1]+0.25*Scf[p])<=1.0:
                # b1, b2, b3 conditions
                #if 0.71875*(Scf[p]+Scf[p+1])<=1.0 and (Scf[p]+0.25*Scf[p+1])<=1.0 and (Scf[p+1]+0.25*Scf[p])<=1.0 and (0.71875*Scf[p]+0.03125*Scf[p+1])<=1.0 and (0.71875*Scf[p+1]+0.03125*Scf[p])<=1.0 :
                # b1+b2 synthetic:
                #if (1.71875*Scf[p]+1.25*Scf[p+1])<=2.0 and (1.71875*Scf[p+1]+1.25*Scf[p])<=2.0:
                    fv=sp.spline_series_vec(tv, tcint, Scf, p)
                    maxlocal=np.max(fv)
                else:
                    maxlocal=0.0
            elif p==2:
                # only b
                #if (Scf[p]+Scf[p+1])<=4./3.:
                # only c
                #if (Scf[p-1]+6.0*Scf[p]+Scf[p+1])<6.0 and (Scf[p]+6.0*Scf[p+1]+Scf[p+2])<6.0 and (Scf[p+1]+6.0*Scf[p+2]+Scf[p+3])<6.0:
                # b and c
                if (Scf[p]+Scf[p+1])<=4./3.  and (Scf[p-1]+6.0*Scf[p]+Scf[p+1])<6.0 and (Scf[p]+6.0*Scf[p+1]+Scf[p+2])<6.0 and (Scf[p+1]+6.0*Scf[p+2]+Scf[p+3])<6.0:
                    fv=sp.spline_series_vec(tv, tcint, Scf, p)
                    maxlocal=np.max(fv) 
                else:
                    maxlocal=0.0           
            
            print 'Max(%f,%f)=%f' %(Scf[p],Scf[p+pincr],maxlocal)
            MaxMat[ii,jj]=maxlocal
            
            if maxlocal>val:
                val=maxlocal
                pos=[]
                pos.append([ii,jj])
            if maxlocal==val:
                val=maxlocal
                pos.append([ii,jj])
     
    # visualise worse case    
    Scf[p]=ssvals[pos[0][0]]
    Scf[p+pincr]=ssvals[pos[0][1]]
    fv=sp.spline_series_vec(tv, tcint, Scf, p)       
    fig = plt.figure('worst case')
    ax = fig.add_subplot(111)
    ax.plot(tv,fv,'0.6',linewidth=3)
    ax.set_xlim([tv[0]-0.25*tv[-1], tv[-1]+0.25*tv[-1] ])
    plt.show() 
    print 'Max value: %f' %(np.max(fv))
    print 'Limit exceeded by: %f percent' %(1e2*val-1e2)
    return val, pos


def check_bound_conditions_03spl(tv,tcint,p):
    
    pincr1=1
    pincr2=2
    NI=len(tcint)-1
    NS=p+NI
    
    print 'testing maximum value'
    Scf=np.zeros((NS,))
    N=51
    ssvals = np.linspace(0,1,N)
    MaxMat=np.empty((N,N))

    val=0.0
    pos=[]
    for ii in range(N):
        for jj in range(N):
            for kk in range(N):
                Scf[p]=ssvals[ii]
                Scf[p+pincr1]=ssvals[jj]
                Scf[p+pincr2]=ssvals[kk]
                
                if p==2:
                    # check boundary
                    check=(Scf[p]+Scf[p+1])<=4./3. and (Scf[p+1]+Scf[p+2])<=4./3. and \
                          (Scf[p-1]+6.0*Scf[p]+Scf[p+1])<6.0 and (Scf[p]+6.0*Scf[p+1]+Scf[p+2])<6.0 and (Scf[p+1]+6.0*Scf[p+2]+Scf[p+3])<6.0
                    if check==True:
                        fv=sp.spline_series_vec(tv, tcint, Scf, p)
                        maxlocal=np.max(fv)
                    else:
                        maxlocal=0.0
                if p==3:
                    check = (0.25*Scf[p] + Scf[p+1] + 0.25*Scf[p+2])<=1 and 0.71875*(Scf[p]+Scf[p+1])<=1.0 and  0.71875*(Scf[p+1]+Scf[p+2])<=1.0
                    #check = 0.71875*(Scf[p]+Scf[p+1])<=1.0 and  0.71875*(Scf[p+1]+Scf[p+2])<=1.0 and \
                    #        (Scf[p]+0.25*Scf[p+1])<=1.0 and (Scf[p+1]+0.25*Scf[p])<=1.0 and (Scf[p+1]+0.25*Scf[p+2])<=1.0 and (Scf[p+2]+0.25*Scf[p+1])<=1.0 and \
                    #        (0.25*Scf[p] + Scf[p+1] + 0.25*Scf[p+2])<=1 #and \
                            # conditions  c1 and c3 uneffective
                            #(0.71875*(Scf[p]+Scf[p+1])+0.03125*Scf[p+2])<=1.0 and  (0.03125*Scf[p]+0.71875*(Scf[p+1]+Scf[p+2]))<=1.0 and \
                            #(0.71875*(Scf[p-1]+Scf[p])+0.03125*Scf[p+1])<=1.0 and  (0.03125*Scf[p-1]+0.71875*(Scf[p]+Scf[p+1]))<=1.0 and \
                            #(0.71875*(Scf[p+1]+Scf[p+2])+0.03125*Scf[p+3])<=1.0 and  (0.03125*Scf[p+1]+0.71875*(Scf[p+2]+Scf[p+3]))<=1.0 and \
                            #(0.71875*(Scf[p+2]+Scf[p+3])+0.03125*Scf[p+4])<=1.0 and  (0.03125*Scf[p+2]+0.71875*(Scf[p+3]+Scf[p+4]))<=1.0
                    
                    if check==True:
                        fv=sp.spline_series_vec(tv, tcint, Scf, p)
                        maxlocal=np.max(fv)
                    else:
                        maxlocal=0.0        
                
                print 'Max(%f,%f,%f)=%f' %(Scf[p],Scf[p+pincr1],Scf[p+pincr2],maxlocal)
                MaxMat[ii,jj]=maxlocal
                
                if maxlocal>val:
                    val=maxlocal
                    pos=[]
                    pos.append([ii,jj,kk])
                if maxlocal==val:
                    val=maxlocal
                    pos.append([ii,jj,kk])
     
    # visualise worse case    
    Scf[p]=ssvals[pos[0][0]]
    Scf[p+pincr1]=ssvals[pos[0][1]]
    Scf[p+pincr2]=ssvals[pos[0][2]]
    fv=sp.spline_series_vec(tv, tcint, Scf, p)       
    fig = plt.figure('worst case')
    ax = fig.add_subplot(111)
    ax.plot(tv,fv,'0.6',linewidth=3)
    ax.set_xlim([tv[0]-0.25*tv[-1], tv[-1]+0.25*tv[-1] ])
    plt.show() 
    print 'Max value: %f' %(np.max(fv))
    print 'Limit exceeded by: %f percent' %(1e2*val-1e2)
    return val, pos



    
    
    
if __name__=='__main__':
    
    T=1.0
    NT=101
    p=3
    NI=10
    NS=p+NI
    tv=np.linspace(0,T,NT)

    tcint = np.linspace(0,T,NI+1)
    #tcint=np.array([0,0.2,0.5,0.8,0.9,1.0])
    
    ### test base generator
    #tc = sp.extended_knots_vector(tcint, p, tv) 
    #test_zero_order(tv,tc)
    #S=test_p_order(tv,tc,p)
    
    
    #test_spline_series_vec(tv, tcint, p)
    # test_change_of_base()
    
    #pincr=2
    #val, pos = return_maxima_2_combinations(tv,tcint,p,pincr)

    #pincr1=1
    #pincr2=2
    #val, pos = return_maxima_3_combinations(tv,tcint,p,pincr1,pincr2)

    #spline_relevant_points(p)
    #check_bound_conditions_03spl(tv,tcint,p)
    



