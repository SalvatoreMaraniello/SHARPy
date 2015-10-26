'''
Created on 19 Oct 2015

@author: sm6110

Testing a class for optimisation based on SLSQP

with in-house jacobian

http://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.minimize.html#scipy.optimize.minimize

Based on xbcompFD in optimiser/wrapper


@todo:  - make choice of solver an input (current only flight dynamic solution)
        - run method has to be changed into a pool, so the FD are computed as well
        - add FD based Jacobian
            - parallel FD: distinction between total number of processor and those used
            in one fwd run.
        - eval_functional: method implemented here (and not in DerivedTypes.Functional)
          as the evaluation requires knowing the values of all the design/parameters/state.
          The code can be made more modular by:
              - moving the method in DerivedTypes.Functionals
              - adding as in put to the method the numerical values
              of the variables required to evaluate the functionals 
        - improve step choice for FD  



@note: Original methods xbcompFD

@warning: bug found when running without Jacobian for the equality 
        inequality constraints. In this case, crash is avoided by
        assigning an initial value to the cost function (and not None).
        The optimisation, however, does not run. The issue was found also
        for simple test cases.
        

            -   
'''

import copy
import ctypes as ct
import numpy as np
import scipy.optimize as scopt
import multiprocessing as mpr

import Main.SharPySettings as Settings
import PyBeam.Utils.DerivedTypes
#import PyOpt.Utils.DerivedTypes
import PyLibs.io.save






class OptComp:
    

    
    def __init__(self, DESIGN, FUNC, PARAM,                   # optimisation
                 XBINPUT, XBOPTS, VMOPTS, VMINPUT, AELAOPTS,  # solution
                 **kwargs):                                   # other
        
        # attach input classes
        self.DESIGN = DESIGN
        self.FUNC   = FUNC
        self.PARAM  = PARAM
        
        self.XBINPUT  = XBINPUT
        self.XBOPTS   = XBOPTS
        self.VMOPTS   = VMOPTS
        self.VMINPUT  = VMINPUT
        self.AELAOPTS = AELAOPTS
        
        self.XBOUTPUT = PyBeam.Utils.DerivedTypes.Xboutput()
        
        # saving
        self.SaveDict=Settings.SaveDict
        if 'SaveDict' in kwargs: self.SaveDict=kwargs['SaveDict']
        self._counter=-1
        #self.SaveDict['OutputFileRoot']=self.SaveDict['OutputFileRoot']+'%.3d'%self._counter 
        
        # fwd run
        if XBOPTS.Solution.value == 912:
            from PyCoupled.Coupled_NlnFlightDynamic import Solve_Py
            self.fwd_run=Solve_Py
        else:
            raise NameError('Solution not supported!')
        
        #optimiser options:
        self.driver_options={'SLSQP': {'disp': False, # True to print convergence message
                                       'iprint': 1  , 
                                       'eps': 1e-04 , # step size Jacobian
                                       'maxiter': 30, 
                                       'ftol': 1e-04} }
        
        #build constraint dictionary
        self.build_constr_dict()
        
        # initialise functionals values
        self.FUNC.gdis.val = np.zeros(sum(self.FUNC.gdis.len))
        self.FUNC.geq.val = np.zeros(sum(self.FUNC.geq.len))
        
        # and related jacobian matrices
        self.FUNC.gdis.jac = np.zeros( (sum(self.FUNC.gdis.len), self.DESIGN.Nx) )
        self.FUNC.geq.jac = np.zeros( (sum(self.FUNC.geq.len), self.DESIGN.Nx) )
        self.FUNC.cost.jac = np.zeros( (self.DESIGN.Nx, ) )
        
        #self.XBOPTS.del_ctypes()
        
        # FD parallel
        self.PROCESSORS=1
        self.parallelFDs = False
        # allocate 1 process for each fwd/perturbed solution
        #self.NumCoresVM = (self.DESIGN.Nx + 1) * [1]
        

     
    def optimise(self):
        
        # Manage performance
        self.set_number_processes()
        
        scopt.minimize(self.run_all, 
                       self.DESIGN.x, 
                       method='SLSQP', 
                       jac=self.get_cost_jac, hess=None, hessp=None, 
                       bounds=self.DESIGN.bounds, 
                       constraints=self.constr_dict, 
                       tol=None, callback=None, 
                       options=self.driver_options['SLSQP'])
        
        
    def run_wrap(self, xin):
          
        # xin input from optimiser
        self.DESIGN.x = xin
        
        # unpack x0
        self.DESIGN.unpack()
        
        # Extract Design and set-up the problem
        self.update_solver_input()
        
        # solve aeroelastic analysis
        
        #self.XBOPTS.set_ctypes()
        self.XBOUTPUT=self.fwd_run(self.XBINPUT, self.XBOPTS, 
                           self.VMOPTS, self.VMINPUT, self.AELAOPTS,
                           SaveDict=self.SaveDict)
        #self.XBOPTS.del_ctypes()
        
        # evaluate functionals
        self.eval_functionals()
        
        return self.FUNC.cost.val
        

    def run_all(self,xin):
        '''
        calls wrapper but also computes the jacobian
        '''
        
        self._counter += 1       
        print('starting iteration %.3d'%self._counter)
        print('state: ', self.DESIGN.scfv) 
        self.save()  
        
        # I/O
        self.SaveDict['OutputFileRoot'] = \
            self.SaveDict['OutputFileRoot'] + '_%.3d'%self._counter    

        # fwd run
        self.run_wrap(xin)
        
        print('fwd run terminated')
        print('cost: ', self.FUNC.cost.val)
        # save state and functional values
        # jacobian not available yet
        self.save()

        #-------------------------------------------------------------------------- compute Jacobian
 
        print('Computing FD based Jacobian')  
        self.delta = self.DESIGN.Nx * [ self.driver_options['SLSQP']['eps'] ]
        
        # ctypes set again in pertub method
        #self.XBOPTS.del_ctypes()
        
        if self.parallelFDs is False:
            
            # debugging: prepare dummy input
            #a = np.array([dd for dd in range(self.DESIGN.Nx)])
            
            for dd in range(self.DESIGN.Nx):
                
                ( self.FUNC.gdis.jac[:,dd], 
                  self.FUNC.geq.jac[:,dd] , 
                  self.FUNC.cost.jac[dd]  ) = self.perturb( copy.deepcopy( self ), dd) 
                '''
                ( self.FUNC.gdis.jac[:,dd], 
                  self.FUNC.geq.jac[:,dd] , 
                  self.FUNC.cost.jac[dd]  ) = self.dummy_perturb( a, dd, 
                                                         sum(self.FUNC.geq.len), 
                                                         sum(self.FUNC.gdis.len),
                                                         copy.deepcopy(self.XBOPTS),
                                                         copy.deepcopy(self) )
                '''
        else:
            
            # debugging: prepare dummy input
            #a = np.array([dd for dd in range(self.DESIGN.Nx)])
            
            # create pool of processes
            pool = mpr.Pool(processes=self.PROCESSORS)  # @UndefinedVariable
            results=[]
            for dd in range(self.DESIGN.Nx):     
                results.append( pool.apply_async(
                                     self.perturb,
                                     args=( copy.deepcopy( self ), dd ) ))   
                '''
                results.append( pool.apply_async( 
                                      self.dummy_perturb, args=(a, dd, 
                                                         sum(self.FUNC.geq.len), 
                                                         sum(self.FUNC.gdis.len),
                                                         copy.deepcopy(self.XBOPTS),
                                                         copy.deepcopy(self)
                                                         #copy.deepcopy(self.XBINPUT) 
                                                         )))           
                '''
            # retrieve results
            jac_list = [p.get() for p in results]
            
            for dd in range(self.DESIGN.Nx):
                ( self.FUNC.gdis.jac[:,dd], 
                  self.FUNC.geq.jac[:,dd] , 
                  self.FUNC.cost.jac[dd]  ) = jac_list[dd]
            # - 1. close the pool (memory in workers goes to zero) 
            # - 2. exit the worker processes (processes are killed)
            pool.close()
            pool.join() 
            
            #self.XBOPTS.set_ctypes()

        print('completed iteration %.3d'%self._counter)
        # re-save to include jacobian
        self.save()
        self.SaveDict['OutputFileRoot'] = self.SaveDict['OutputFileRoot'][:-4]
        
        return self.FUNC.cost.val 
                
 

    def dummy_perturb(self,a,dd,Neq,Ndis,cpXBOPTS,cpself):
        '''
        Dummy perturb function for debugging
        '''
        
        cpXBOPTS.set_ctypes()
        
        #jeq_dd  = a[dd]*np.ones(Neq)
        jeq_dd = cpself.XBINPUT.BeamMass[dd,dd]*np.ones(Neq)
        jdis_dd = a[dd]*np.ones(Ndis)
        jcost_dd = cpXBOPTS.Solution.value*a[dd]

        return (jdis_dd, jeq_dd, jcost_dd)   
        

    def perturb(self, cpself, dd):
        '''
        Perform the analysis after perturbing the dd-th design parameter
        '''

        print('computing %.3d derivative'%dd)
        
        #### set number of cores to be used
        cpself.VMOPTS.NumCores = ct.c_int( cpself.NumCoresVM[dd+1] )
        
        # perturb the design
        cpself.DESIGN.x[dd] = cpself.DESIGN.x[dd] + cpself.delta[dd]

        # change a few setting...
        cpself.SaveDict['OutputDir']=cpself.SaveDict['OutputDir']+'FD/'
        cpself.SaveDict['OutputFileRoot']=cpself.SaveDict['OutputFileRoot']+'FD%.3d'%dd
        cpself.SaveDict['SaveWake']=False
        #SaveDict['SaveProgress']=False
        
        cpself.save()
          
        # launch the wrapper...
        cpself.run_wrap( cpself.DESIGN.x )
        
        # save...
        cpself.save()
        
        # ... and compute the jacobian
        jdis_dd  = ( cpself.FUNC.gdis.val - self.FUNC.gdis.val )/cpself.delta[dd]
        jeq_dd   = ( cpself.FUNC.geq.val  - self.FUNC.geq.val  )/cpself.delta[dd]
        jcost_dd = ( cpself.FUNC.cost.val - self.FUNC.cost.val )/cpself.delta[dd]
        
        #del(cpself)      
        
        return (jdis_dd, jeq_dd, jcost_dd)          

        
        
        #---------------------------------------------------------------------- end compute Jacobian


    def eval_functionals(self):
        '''
        Evaluate cost/constraints
        '''
         
        # cost function (not list)
        if self._counter>-100:#0: #
            # by def. the cost will dep on the state
            args_values = self.get_args_values(self.FUNC.cost.args)
            self.FUNC.cost.val = self.FUNC.cost.fun( *args_values )
            self.FUNC.cost.fun_name = self.FUNC.cost.fun.__name__   # not sure why...
            
        # evaluate constraints
        for gset in [self.FUNC.geq, self.FUNC.gdis]:
            
            gsetval_list=[]
            gset.fun_name=[]
            Ng = len(gset.fun)
            
            for gg in range(Ng):
                
                # extract argument values
                ggargsval= self.get_args_values(gset.args[gg])
                
                # evaluate
                ggfun  = gset.fun[gg]
                ggval = ggfun( *ggargsval )
                gset.fun_name.append( ggfun.__name__ )
                
                # append values
                if np.isscalar(ggval): 
                    # output is a scalar
                    gsetval_list.append( ggval )
                else: 
                    # output is an array
                    Nelems = np.prod(ggval.shape)
                    gvec = np.reshape( ggval, ( Nelems ,) )
                    for gelem in gvec:
                        gsetval_list.append( gelem )
            gset.val = np.array(gsetval_list)               

        return


    def get_args_values(self,AttrList):
        '''
        Given a list of names of attribute attached to self, the method
        allows to retrieve their value.
        The method works even if the attributes are stored in subclasses.
        
        Eg: if self contains: self.A.B.c, the value of c is extracted
        by passing AttrList=['A.B.c']
        '''
        T=[]
        for attr in AttrList:
            sublist = attr.split('.')
            # reach attribute (even if contained in subclasses)
            subclass=self
            for sub in sublist:
                val = getattr(subclass,sub)
                subclass=val
            T.append(val)
        
        return tuple(T)   
        
                 
    def update_solver_input(self):
        ''' 
        Problem specific. Redefine this method in the optimisation 
        input file!
        '''
        pass
    

    def build_constr_dict(self):
        '''
        Builds list of dictionaries to define constraints for
        scipy.optimize.minimize
        
        constraints : dict or sequence of dict, optional
        Constraints definition (only for COBYLA and SLSQP). Each constraint is defined in a dictionary with fields:
        type : str
        Constraint type: ‘eq’ for equality, ‘ineq’ for inequality.
        fun : callable
        The function defining the constraint.
        jac : callable, optional
        The Jacobian of fun (only for SLSQP).
        args : sequence, optional
        Extra arguments to be passed to the function and Jacobian.  
        
        '''
        
        self.constr_dict=[]
        
        
        # for equalities
        base_dict = {'type': 'eq', 
                     'fun' : self.FUNC.geq.get_value, 
                     'jac' : self.get_geq_jac,
                     'args': ()    } 
        
        Ng = sum( self.FUNC.geq.len )
        list_eq_dict=[]
        
        for kk in range(Ng): 
            list_eq_dict.append(base_dict.copy())
            list_eq_dict[kk]['args']=(kk,)
            
        # for inequalities
        base_dict = {'type': 'ineq', 
                     'fun' : self.FUNC.gdis.get_value, 
                     'jac' : self.get_gdis_jac,
                     'args': ()    } 
        
        Ng = sum( self.FUNC.gdis.len )
        list_dis_dict=[]
        
        for kk in range(Ng): 
            list_dis_dict.append( base_dict.copy() )
            list_dis_dict[kk]['args']=(kk,)        
        
        
        '''
        # for equalities
        base_dict = {'type': 'eq', 
                     'fun' :  PyOpt.Utils.DerivedTypes.get_value_from_FuncProp_class , 
                     #'jac' : None,
                     'args': ()    } 
        
        Ng = sum( self.FUNC.geq.len )
        list_eq_dict=[]
        
        for kk in range(Ng): 
            list_eq_dict.append(base_dict.copy())
            list_eq_dict[kk]['args']=(self.FUNC.gdis,kk)
            
        # for inequalities
        base_dict = {'type': 'ineq', 
                     'fun' :  PyOpt.Utils.DerivedTypes.get_value_from_FuncProp_class , 
                    #'jac' : None,
                     'args': ()    } 
        
        Ng = sum( self.FUNC.gdis.len )
        list_dis_dict=[]
        
        for kk in range(Ng): 
            list_dis_dict.append( base_dict.copy() )
            list_dis_dict[kk]['args']=(self.FUNC.gdis,kk)                
        '''
        
        self.constr_dict=list_eq_dict + list_dis_dict
        
           
    def get_geq_jac(self,dummy,ii):
        '''dummy is to accommodate scipy.optimize.minimize standard 
        (1st argument of function returning the jacobian has to be
        the design)'''
        return self.FUNC.geq.jac[ii,:]
    def get_gdis_jac(self,dummy,ii):
        '''dummy is to accommodate scipy.optimize.minimize standard 
        (1st argument of function returning the jacobian has to be
        the design)'''
        return self.FUNC.gdis.jac[ii,:]           
    def get_cost_jac(self,dummy):
        '''dummy is to accommodate scipy.optimize.minimize standard 
        (1st argument of function returning the jacobian has to be
        the design)'''
        return self.FUNC.cost.jac


    def save(self):
        
        # save optimisation data
        PyLibs.io.save.h5file( self.SaveDict['OutputDir'], 
                               'opt_'+self.SaveDict['OutputFileRoot']+'.h5', 
                                self.DESIGN, self.PARAM, 
                                self.FUNC.cost, self.FUNC.geq, self.FUNC.gdis )
        


    def set_number_processes(self):
        
        '''
        method to redistribute the CPUs per process.
        method can be beneficial only if the amount of processors
        available is larger then the amount of design variables
        '''

        if self.parallelFDs is True:
            
            NumPrMin = self.DESIGN.Nx + 1
            PRratio = self.PROCESSORS//NumPrMin
            
            if PRratio == 0:
                # more then 1 loop required
                self.NumCoresVM = NumPrMin*[1]
                # check if during the last loop some processes will
                # not be used.
                Processes_out = self.PROCESSORS - NumPrMin%self.PROCESSORS
                Nloops = NumPrMin//self.PROCESSORS + 1
                ### method not optimum as
                ###    a. in initial loop more cpus then available are required
                ###    b. run time is reduced only if one of the loops run time
                ###        is reduced. Having a loop with some processes running
                ###        faster does not guarantee that overall time is cut down
                #for pp in range(Processes_out):
                #    self.NumCoresVM[pp] += 1 
                ###
                ###
            elif PRratio > 0:
                self.NumCoresVM = NumPrMin * [PRratio]
                # check processes out
                Processes_out = self.PROCESSORS - PRratio*NumPrMin 
                # If so, use more cpu for the VLM solution
                for pp in range(Processes_out):
                    self.NumCoresVM[pp] += 1 
            else:
                raise NameError("Ratio {PROCESSORS}/{Number Design Variable} not positive!")        
        
        else:
            self.NumCoresVM = NumPrMin*[self.VMOPTS.NumCores]
            