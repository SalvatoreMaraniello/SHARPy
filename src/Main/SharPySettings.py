'''@package PyBeam.Main.SharPySettings
@brief      System settings for SharPy
@author     Rob Simpson
@contact    r.simpson11@imperial.ac.uk
@version    0.0
@date       22/11/2012
@pre        None
@warning    None

-------------------------------------------------------------------------------
Modified: 25 Feb 2015 by S. Maraniello

Change list: 
    - new method to get home folder (use environmental variable)
    
Change in Progress:
    - use 
          " OutputDir, OutputfileRoot, 
          PlotTec, WriteOut, PlotOut, WriteUVLMdebug "    
      only as default values. these should be set in input/output files.  
        
'''
import getpass
import sys
import os


# Directories.
userid = getpass.getuser()

try: # if environmental variable is set-up
    SharPyProjectDir = os.environ["SHARPYDIR"]+'/' 
except:
    #SharPyProjectDir = '/home/' + userid + '/git/SHARPy/'
    SharPyProjectDir = os.path.expanduser('~') + '/git/SHARPy/'    

BeamLibDir = SharPyProjectDir + 'BeamLib/bin/'
BeamLibName = './BeamLib.so'

UVLMLibDir = SharPyProjectDir + 'UVLMLib/Debug/'
UVLMLibName = './UVLMLib.so'

OutputDir =  SharPyProjectDir + 'output/temp/'
OutputFileRoot = 'Foo'


# Python path.
sys.path.append(SharPyProjectDir + 'src')
sys.path.append(SharPyProjectDir + 'src/PyBeam')
sys.path.append(SharPyProjectDir + 'src/PyBeam/Utils')
sys.path.append(SharPyProjectDir + 'src/PyBeam/Solver')
sys.path.append(SharPyProjectDir + 'src/PyBeam/Main')
sys.path.append(SharPyProjectDir + 'src/PyFSI')
sys.path.append(SharPyProjectDir + 'src/PyFSI/Utils')
sys.path.append(SharPyProjectDir + 'src/PyAero')
sys.path.append(SharPyProjectDir + 'src/PyAero/UVLM')
sys.path.append(SharPyProjectDir + 'src/PyAero/UVLM/Utils')
sys.path.append(SharPyProjectDir + 'src/PyAero/UVLM/Solver')
sys.path.append(SharPyProjectDir + 'src/PyCoupled')
sys.path.append(SharPyProjectDir + 'src/PyCoupled/Utils')
sys.path.append(SharPyProjectDir + 'BeamLib/src/wrapper/f2py')

# Structural Code Constants.
MaxElNod = 3
DimMat = 24 # Memory allocated for sparse matrix storage in fortran solver.
RotEpsilon = 0.001 # Rotations below this are linearised.

# Tecplot.
PlotTec = True

# Output options.
WriteOut = True

# Live plotting options.
PlotOut = True

# Write UVLM info to text file for debugging.
WriteUVLMdebug = False


# Save Dictionary for Standard I/O:
# Imnport the dictionary in the input file pass it to the solver.
SaveDict={}
SaveDict['OutputDir'] = OutputDir
SaveDict['OutputFileRoot'] = OutputFileRoot
SaveDict['Format']='h5' # 'dat', 'all'
SaveDict['SaveProgress']=True
SaveDict['NumSavePoints']=10
SaveDict['SaveWake']=True
SaveDict['SaveWakeFreq']=50  # save every xx time-steps
SaveDict['SaveFD']=True      # FD 

if __name__ == '__main__':
    pass
