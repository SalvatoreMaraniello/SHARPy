
'''
SHARPy installation:

prerequisites: 
- Anaconda3-2.3.0-Linux-x86-64 or higher
- gfortran compiler
- modify UVLMLib/MakeFile
    EIGENDIR=PATHTO/eigen-eigen-3.2.5
    BOOSTDIR=PATHTO/boost_1_58_0 


'''

import os

CloneRepo = False

# clone repository
if CloneRepo: os.system('git clone https://github.com/SalvatoreMaraniello/SHARPy')
SHARPYDIR = os.getcwd() + '/SHARPy'

# add to bash fine
# export SHARPYDIR="/home/sm6110/git/SHARPy"


# compile Fortran libraries
os.chdir(SHARPYDIR+'/BeamLib/src/wrapper')
os.system('mkdir -p %s' %(SHARPYDIR+'/bin'))
os.system('make all')


# f2py run
os.chdir(SHARPYDIR+'/BeamLib/src/wrapper/f2py')
os.system('chmod u+x ./runf2py3x.sh')
os.system('./runf2py3x.sh')


# compile UVLM libraries
os.chdir(SHARPYDIR+'/UVLMLib')
os.system('make all')


# return to original directory
os.chdir(SHARPYDIR)

print('Installation almost completd...')
print('1. Add export SHARPYDIR="......./SHARPy" to your bash file')





