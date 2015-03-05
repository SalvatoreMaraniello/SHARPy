'''@package PyBeam.Main.SharPySettings
@brief      System settings for SharPy
@author     Rob Simpson
@contact    r.simpson11@imperial.ac.uk
@version    0.0
@date       22/11/2012
@pre        None
@warning    None

@Modified: 25 Feb 2015 by S. Maraniello
@Change list: change method to get home folder


'''
import getpass
import sys
import os


# Directories.
userid = getpass.getuser()

#sm
#SharPyProjectDir = '/home/' + userid + '/git/SHARPy/' # works on ubuntu
#SharPyProjectDir = os.path.expanduser('~') + '/git/SHARPy/' # works if installed in home
SharPyProjectDir = os.environ["SHARPYDIR"] + '/git/SHARPy/' # works if you setup an environmental variable during the install process

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

if __name__ == '__main__':
    pass
