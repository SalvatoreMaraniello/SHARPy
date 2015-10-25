'''
@author:   sm6110
@date:     21 Oct 2015

@summary:  collects simple classes to define the input of a general
           optimisation problem
            
@todo:     In Design.__init__ there is the risk to overwrite exisisting
           methods/variables.
'''

import numpy as np



class Design():
    '''
    @summary: Class to store all the design variables for an optimisation 
    problem. The class can be built automatically by passing all the design
    variables as key words input.
    
    @param Vars: list containing all the design variables
    @param x: contains all the values of the design variables
               in a 1D array format
               
    @note: see example below. The class can also be built automatically
    
    @warning: can overwrite existing variables/methods!!!
    
    '''    
    
    def __init__(self, **kwargs):
          
        self.Vars=[]
        self.shape=[]
        self.len=[]
        
        for ww in kwargs:
            setattr(self, ww, kwargs[ww])
            self.Vars.append(ww)   
             
            # determine shape
            if np.isscalar(kwargs[ww]):
                self.shape.append( 1 )
                self.len.append( 1 )
            else:
                self.shape.append( kwargs[ww].shape )
                self.len.append( np.prod( kwargs[ww].shape ) )
               
        self.Nvars=len(self.Vars)
        
        # pack initial guess
        self.x=None
        self.pack()  
        self.Nx = len(self.x)   
        
        # Create Bound input       min | max
        self.bounds = self.Nx * [ (None,None), ]
        
    
    def set_bounds(self, **kwargs):
        '''
        quickly add bounds for arrays/scalars. If the target is 
        an array, the same bounds can be applied to each elements
        '''
         
        for ww in kwargs:
            # find position in Vars
            ii = self.Vars.index(ww)
            # find position in x, bounds
            kk = sum( self.len[:ii] )
            
            # allocate
            if type(kwargs[ww]) is tuple:
                # if an array, the same value is given to all
                self.bounds[kk:kk+self.len[ii]] = self.len[ii] * [ kwargs[ww] ]
                pass
            else:
                # if is a list, different bounds per array can be assigned
                if len(kwargs[ww])==self.len[ii]:
                    self.bounds[kk:kk+self.len[ii]] =  kwargs[ww]              
                else:
                    raise NameError('When the bounds of an array are specified through list,' +
                                    'this needs to have the same length as the array!!!'   )

        
    def pack(self):
        '''
        Goes from a list of design variables to a 1D array.
        '''
    
        x_list=[]
        #Nx = len(self.Vars)
        
        for ii in range(self.Nvars):
            
            val = getattr(self,self.Vars[ii])
            
            # append values
            if np.isscalar(val): 
                # output is a scalar
                x_list.append( val )
            else: 
                # output is an array
                Nelems = np.prod(self.shape[ii])
                valvec = np.reshape( val, ( Nelems ,) )
                for vv in valvec:
                    x_list.append( vv )
        
        self.x = np.array(x_list)  
          
    
    def unpack(self):
        '''
        updates the design variables from the 1D array x0
        '''
        
        kstart=0 # index in x
        for ii in range(self.Nvars):
            if self.shape[ii]==1:
                setattr(self, self.Vars[ii], self.x[kstart])
                kstart += 1
            else:
                Nelems = np.prod(self.shape[ii])
                setattr( self, self.Vars[ii], self.x[kstart:kstart+Nelems] )
                kstart += Nelems
                



class Functionals():
    '''
    Merely a container for cost/constraints data. To each of
    the set cost, geq, gdis is associated a FuncProp class 
    instance
    
    @param cost: cost function data.
        cost.fun: callable function
        cost.args: list of argument names in string format
        cost.val: scalar
    @param geq/gdis: equality and disequality constraints
        gxxx.fun: list of callable functions
        gxxx.args: lift of list of argument names in string format
        gxxx.val: list of numerical values (scalar/arrays)
    '''
    
    
    def __init__(self):
        
        self.cost=FuncProp('cost')
        self.geq=FuncProp('geq')
        self.gdis=FuncProp('gdis')
        
    
    def pack(self):
        pass


    def unpack(self):
        pass
        

  
        
class FuncProp():
    '''
    Container for functional properties. Used in Functional class
    initialisation.
    The class contains 3 attributes:
        fun: callable objective function or list of function required 
             to evaluate cost/constraints
        args: list, or list of lists, containing the input argument for
             the functions in fun
        val: numerical values
        shape: list containing the shape of the output of each fun
    '''
    
    def __init__(self,name=None):
    
        if name != None:
            self.name=name
    
        self.fun=None
        self.args=[]
        self.val=None
        #self.shape=[]
        self.len=[]
        
        
    def get_value(self,dummy,ii):
        '''
        Required to extract ii-th value from val array
        '''
        return self.val[ii]
    

    
    
class Parameters():
    '''
    Container for additional parameters required for the optimisation
    '''
    
    def __init__(self):
        pass
        
    


if __name__=='__main__':
    
    #  create Design class
    design=Design( scfv=np.zeros(5) ,
                   diad_angle=4.0  ,
                   a = 7, 
                   b = np.array([4,5,6]),
                   c = np.array([9, 10, 11]) )
    
    # variables have been attached as attributed
    print(design.scfv)
    print(design.diad_angle)
    # and also ordered according to
    print(design.Vars)
    # into
    print(design.x)
    print(design.shape)
    print(design.len)
    # the variables have so far no bounds
    print(design.bounds)
    # but these can be added as
    design.set_bounds(a=(-7,7),
                  c=(0,10),
                  scfv=(5,50) )
    print(design.bounds)
    # or to specify specific values for each element in the array
    design.set_bounds( b= [(-1,1), (-2,2), (-3,3)],
                       diad_angle=[(50,80)] )
    print(design.bounds) 
    
    
    
    
    
    
    