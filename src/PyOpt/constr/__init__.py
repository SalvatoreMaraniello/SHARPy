'''

@author: Salvatore Maraniello
@date: 20 Oct 2015

@summary: contains a collection of constraint methods for setting
          up optimisation problems
          
@note: by defaults, all disequality constrains functions return 
       a g value such that:
             g > 0
       if the constraint is verified as per scipy.optimize.minimize
       convention

'''