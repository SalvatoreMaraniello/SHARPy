'''
Created on 21 Sep 2015

@author: sm6110

Adapted from BeamLib/src/optimiser/wrapper/pysrc/lib.read

'''

import h5py
from warnings import warn


def h5file(filename):
    '''
    Read all entries of a HDF5 file and attaches them to class.
    Groups are saved as sub-classes, while dataset values are saved as class 
    attributes
    
    Important: though this method allows to read all the input/output classes 
    required to run the aeroelastic solutions, ctypes variable are not saved as
    such!
    
    '''
    
    class H: pass
    Hinst=H()
     
    # read and scan file
    hdfile=h5py.File(filename,'r')
    NamesList=[]                   # dataset names
    hdfile.visit(NamesList.append)
    
    for name in NamesList:
        #print('Found %s...'%name)
        if type(hdfile[name]) is h5py._hl.group.Group:  # @UndefinedVariable
            #print('    %s is a group!' %name)
            if hasattr(Hinst,'name') is False:
                setattr(Hinst,name,H())
        else:
            if '/' in name:
                subnames=name.split('/')
                if hasattr(Hinst,subnames[0]) is False:
                    setattr(Hinst,subnames[0],H())
                #print('    Extracting class %s'%subnames[0])
                subclass=getattr(Hinst,subnames[0])
                #print('    Allocating attribute %s'%subnames[1])
                setattr(subclass,subnames[1],hdfile[name].value)
                #print('    copying subclass %s back'%subnames[0])
                setattr(Hinst,subnames[0],subclass)
            else:
                setattr(Hinst,name,hdfile[name].value)

    '''
    #setattr(XBinst,'FollowerForce',hdfile['FollowerForce'].value)

    # parameters for constant beam span reconstruction 
    XBinst = conditional_reading(hdfile,XBinst,'cs_l2')

    setattr(XBinst,'',hdfile[''].value)
    XBinst = conditional_reading(hdfile,XBinst,'')
    '''
    
    # close and return
    hdfile.close()  
    
    return Hinst  
  


def conditional_reading(hdfile,obj,fieldname): 
    ''' 
    Given a hdf5 file 'hdfile' and the object obj, the routine:
        a. if the field fieldname is found and has a value, assigns it to the
           attribute obj.fieldname.  
        b. does nothing otherwise    
    '''
    
    try:
        val = hdfile[fieldname].value
        if val!='not found' and val!='no value':
            setattr(obj,fieldname,val)
    except:
        warn('Attribute "%s" not found!!!' %fieldname)
            
    return obj 

  
    
def h5series(rootname,attrlist):
    ''' 
    Given a list of datasets, creates a list of lists for all the solutions
    run for a DOE or optimisation. the output is in a list format
    '''
    
    outlist=[]
    
    cc=0
    go_on=True
    
    while go_on is True:
        cc_str =  '%.3d' % (cc)
        filename = rootname + cc_str + '.h5'
        
        try:
            print( 'Reading: %s' %(filename) )
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
            print( '%s not found. %s files read in total!' %(filename,cc_str) )
            go_on=False
    
    return outlist
     


def h5list(fileslist,attrlist):
    ''' 
    Equivalent to h5series but reads files from an user defined list (fileslist)
    
    All the attributes in 'attrlist' are read and stored in outlist.
    '''
    
    outlist=[]
    
    for filename in fileslist:
        

        print( 'Reading: %s' %(filename) )
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

    return outlist



if __name__=='__main__':
    
   filename='/home/sm6110/git/SHARPy_studies/aerocomp/pycoupled/Goland_NumNE3_NumEl12_M6_U20.0_sol312_dyn.h5'
   hd=h5file(filename) 

