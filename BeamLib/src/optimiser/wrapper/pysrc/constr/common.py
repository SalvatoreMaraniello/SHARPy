'''
Created on 12 Feb 2015

@author: sm6110


General Methods/Utility


'''


def return_args(C,AttrList):
    '''Given a class C and a list of attributes to the class, the function
    creates a tuple with the class attributes 
    
    Remark: duplicate of this function in module cost.py!!!
    '''
    
    T=[]
    for attr in AttrList:
        #print attr
        T.append(getattr(C,attr))
    
    return tuple(T)
