
'''
SHARPy installation file

@author: Salvatore Maraniello
@contact: salvatore.maraniello10@imperial.ac.uk
@date: 1 Feb 2017
@brief: SHARPy installation file for Ubuntu. Read instructions below.


Prerequisites: 
--------------
- Anaconda3-2.3.0-Linux-x86-64 or higher (see below)
- gfortran compiler
- C++ libraries Eigen and Boost. If not availale, these can be downloaded using
this script (see Instruction before).


Instructions:
-------------
- Make sure you have installed Anaconda and check the name of the f2py 
executable shipped with it is 'f2py3'. If not, modify the 'f2py_executable'
variable accordingly. See also the troubleshooting section below.
- In the default setting, this script will allow you to install the latest
versions of the C++ libraries Eigen and Boost, required to compile the UVLM
solver. If you wish to use your own versions of these libraries, or prefer to
install them manually, set to False the 'DownloadEigen' and 'DownloadBoost' 
variables below and specify their location through the EIGENDIR and BOOSTDIR
variables.
- If you also want to download SHARPy, set CloneRepo=True
- The installer will try to modify the bash file ~/.bash_aliases. If this is
not your bash file, or you prefer to modify this manually, set ModifyBash=False
and follow the instructions below.
- Once you are happy with this setting, navigate into the SHARPy folder ( or 
to the path where you want to clone  the SHARPy folder) and run:
	python install.py
- The installer will try to add an environmental variable including the path to
SHARPY into your bash file ~/.bash_aliases'. If an error occurs at this state,
just add the following line 
    export SHARPYDIR="your/path/to/SHARPy"
to your bash file.


Uninstall:
----------
- Remove the SHARPy directory from your system. 


Installation Troubleshot: 
------------
- Finding the f2py executable name:
To determine the executable name, you can simply type 'f2py' in terminal, press
Tab and see at the suggested option. If required, use the 'whereis' command to 
find the location of the executable (on Unix systems, these are typically placed 
in /usr/bin or /usr/local/bin).
- Add path to SHARPy in bash file:
In Ubuntu systems you may want to add this to ~/.bash_aliases. To check that
the setting is correct, once you have modified your bash file, open a new 
terminal and type:
    echo $SHARPYDIR	


Installing Anaconda:
--------------------
Ref. https://docs.continuum.io/anaconda
From command line:
- Navigate to folder where you want to install Anaconda
- wget https://repo.continuum.io/archive/Anaconda3-2.3.0-Linux-x86_64.sh
- chmod u+x Anaconda3-2.3.0-Linux-x86_64.sh
- ./Anaconda3-2.3.0-Linux-x86_64.sh
- Follow installer instructions
--------------------------------------------------------------------------------
'''

import os


# Input (optional):
# -----------------
# Include here the paths to the C++ libraries Eigen and Boost. If these are not 
# installed on your system, the latest versions will be downloaded and installed
# under the SHARPy directory.
DownloadEigen=True
DownloadBoost=True
EIGENDIR=None
BOOSTDIR=None

# If you wish to clone SHARPy, set this to True 
CloneRepo = False

# Set the name of the f2py executable installed with Anaconda. 
f2py_executable='f2py3'

# Modify bash file
ModifyBash=True



# Installer: do not modify!
# ------------------------------------------------------------------------------


# Clone repository
# ----------------
if CloneRepo: 
    os.system('git clone https://github.com/SalvatoreMaraniello/SHARPy')
    SHARPYDIR = os.getcwd() + '/SHARPy'
else:
    SHARPYDIR = os.getcwd()
os.chdir(SHARPYDIR)


# Install C++ libraries
# ---------------------
# Determine paths
if EIGENDIR==None or BOOSTDIR==None: os.system('mkdir -p %s/ExtLib' %SHARPYDIR)
if EIGENDIR==None: EIGENDIR= SHARPYDIR + '/ExtLib/eigen'
if BOOSTDIR==None: BOOSTDIR= SHARPYDIR + '/ExtLib/boost_1_61_0'
# Eigen
if DownloadEigen:
    os.chdir(SHARPYDIR + '/ExtLib')
    os.system('hg clone https://bitbucket.org/eigen/eigen/')
    print('Eigen library successfully downloaded into %s' %EIGENDIR)
# Boost
if DownloadBoost:
    os.chdir(SHARPYDIR + '/ExtLib')
    os.system('wget -O boost_1_61_0.tar.gz '
              'http://sourceforge.net/projects/boost/files/boost/1.61.0/'
              'boost_1_61_0.tar.gz/download')
    os.system('tar zxvf boost_1_61_0.tar.gz')
    os.system('rm boost_1_61_0.tar.gz')
    print('Eigen library successfully downloaded into %s' %BOOSTDIR)	


# Compile UVLM solver
# -------------------
os.chdir(SHARPYDIR+'/UVLMLib')
# Read UVLM libraries Makefile
MakeFileName='./Makefile'
with open(MakeFileName, 'r') as file:
    Lines = file.readlines()
# Update paths to EIGEN/BOOST libraries
for ii in range(len(Lines)):
    if Lines[ii][:9]=="EIGENDIR=": 
    	print(Lines[ii])
    	Lines[ii]='EIGENDIR='+EIGENDIR+'\n'
    	print(Lines[ii])
    if Lines[ii][:9]=="BOOSTDIR=": 
    	Lines[ii]='BOOSTDIR='+BOOSTDIR+'\n'
# Overwrite Makefile
with open(MakeFileName, 'w') as file:
	for  line in Lines: file.write(line)
# compile UVLM libraries
os.system('make all')


# Compile Fortran libraries
# -------------------------
os.chdir(SHARPYDIR+'/BeamLib/src/wrapper')
os.system('mkdir -p %s' %(SHARPYDIR+'/BeamLib/bin'))
os.system('make all')


# Run f2py
# --------
os.chdir(SHARPYDIR+'/BeamLib/src/wrapper/f2py')
#os.system('chmod u+x ./runf2py3.sh')
#os.system('./runf2py3.sh')
os.system('%s -h lib_fem_f2py.pyf -m lib_fem lib_fem_f2py.f90'
	                                   ' --overwrite-signature'%f2py_executable)
os.system('%s -c -m lib_fem lib_fem_f2py.pyf lib_fem_f2py.f90' %f2py_executable)
os.system('%s -h lib_cbeam3_f2py.pyf -m lib_cbeam3 lib_cbeam3_f2py.f90'
	                                   ' --overwrite-signature'%f2py_executable)
os.system('%s -c -m lib_cbeam3 -I../ -L../../lib/src/ lib_cbeam3_f2py.pyf '
	      'lib_cbeam3_f2py.f90 ../../lib/src/lib_fem.o ../../lib/src/lib_rot.o'
	                             ' ../../lib/src/lib_rotvect.o'%f2py_executable)
os.system('%s -h lib_rotvect_f2py.pyf -m lib_rotvect lib_rotvect_f2py.f90 '
	                                    '--overwrite-signature'%f2py_executable)
os.system('%s -c -m lib_rotvect -I../ -L../../lib/src/ lib_rotvect_f2py.pyf'
	            ' lib_rotvect_f2py.f90 ../../lib/src/lib_rot.o'%f2py_executable)


# Export SHARPYDIR
# ----------------
# This step may require sudo permissions

# return to original directory
os.chdir(SHARPYDIR)

if ModifyBash:
    # Back up copy of bash file
    homedir = os.path.expanduser('~')
    try:
        os.system('cp %s/.bash_aliases %s/.bash_aliases-sharpy.backup'
                                                                   %(2*(homedir,)) )
        with open('%s/.bash_aliases' %homedir, 'a') as file:
    	    file.write('# SHARPY path:\n'
    			       'export SHARPYDIR="%s"' %SHARPYDIR)
    except:
        print('Warning: I could not export SHARPy path to your bash file !!!')
        print('Add the following line\n'
    		  'export SHARPYDIR="%s"\n'
    		  'to your bash file' %SHARPYDIR)

print('Installation completed!')


'''
Manually install packages
sudo apt-get install gfortran
sudo apt-get install python3-f2py
sudo apt-get install python3-dev

sudo apt install mercurial

sudo apt-get install python3-numpy
sudo apt-get install python3-scipy
sudo apt-get install python3-matplotlib
'''
