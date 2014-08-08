'''

Salvatore Maraniello
7 Aug 2014

Module containing setting paths for dependencies.

Modify the entry openmdao_abspath according to your installation!

'''

import os
import sys

# ----------------------------------------------------------------- Define Paths
# openMDAO main folder
openmdao_abspath='/home/sm6110/mylibs/openmdao-0.10.0'

#optimiser_abspath='/home/sm6110/git/SHARPy/BeamLib/src/optimiser'
wrapper_abspath=os.path.dirname(os.path.realpath(__file__))
optimiser_abspath=os.path.abspath(wrapper_abspath+'/..')

# -------------------------------------------------------------- Dynamic Library
wrapso_abspath=optimiser_abspath+'/bin/xbeamopt.so'

# ----------------------------------------------------------------- Append Paths
sys.path.append(optimiser_abspath+'/wrapper/pysrc')


''' 
OpenMDAO libraries... 
                                                   THIS MAY REQUIRE MANUAL CHANGE
'''
sys.path.append(openmdao_abspath+'/local/lib/python2.7/site-packages/setuptools-0.9.5-py2.7.egg') 
sys.path.append(openmdao_abspath+'/local/lib/python2.7/site-packages/Jinja2-2.4-py2.7.egg') 
sys.path.append(openmdao_abspath+'/local/lib/python2.7/site-packages/Pyevolve-0.6-py2.7.egg') 
sys.path.append(openmdao_abspath+'/local/lib/python2.7/site-packages/Pygments-1.3.1-py2.7.egg') 
sys.path.append(openmdao_abspath+'/local/lib/python2.7/site-packages/SetupDocs-1.0.5-py2.7.egg') 
sys.path.append(openmdao_abspath+'/local/lib/python2.7/site-packages/Sphinx-1.2.2-py2.7.egg') 
sys.path.append(openmdao_abspath+'/local/lib/python2.7/site-packages/bson-0.3.3-py2.7.egg') 
sys.path.append(openmdao_abspath+'/local/lib/python2.7/site-packages/cobyla-1.0.1-py2.7-linux-x86_64.egg') 
sys.path.append(openmdao_abspath+'/local/lib/python2.7/site-packages/conmin-1.0.1-py2.7-linux-x86_64.egg') 
sys.path.append(openmdao_abspath+'/local/lib/python2.7/site-packages/decorator-3.2.0-py2.7.egg') 
sys.path.append(openmdao_abspath+'/local/lib/python2.7/site-packages/docutils-0.10-py2.7.egg') 
sys.path.append(openmdao_abspath+'/local/lib/python2.7/site-packages/mock-1.0.1-py2.7.egg') 
sys.path.append(openmdao_abspath+'/local/lib/python2.7/site-packages/networkx-1.8.1-py2.7.egg') 
sys.path.append(openmdao_abspath+'/local/lib/python2.7/site-packages/newsumt-1.1.0-py2.7-linux-x86_64.egg') 
sys.path.append(openmdao_abspath+'/local/lib/python2.7/site-packages/nose-1.3.3-py2.7.egg') 
sys.path.append(openmdao_abspath+'/local/lib/python2.7/site-packages/openmdao.examples.bar3simulation-0.10.0-py2.7-linux-x86_64.egg') 
sys.path.append(openmdao_abspath+'/local/lib/python2.7/site-packages/openmdao.examples.enginedesign-0.10.0-py2.7-linux-x86_64.egg') 
sys.path.append(openmdao_abspath+'/local/lib/python2.7/site-packages/openmdao.examples.expected_improvement-0.10.0-py2.7.egg') 
sys.path.append(openmdao_abspath+'/local/lib/python2.7/site-packages/openmdao.examples.mdao-0.10.0-py2.7.egg') 
sys.path.append(openmdao_abspath+'/local/lib/python2.7/site-packages/openmdao.examples.nozzle_geometry_doe-0.10.0-py2.7.egg') 
sys.path.append(openmdao_abspath+'/local/lib/python2.7/site-packages/openmdao.examples.simple-0.10.0-py2.7.egg') 
sys.path.append(openmdao_abspath+'/local/lib/python2.7/site-packages/openmdao.lib-0.10.0-py2.7.egg') 
sys.path.append(openmdao_abspath+'/local/lib/python2.7/site-packages/openmdao.main-0.10.0-py2.7.egg') 
sys.path.append(openmdao_abspath+'/local/lib/python2.7/site-packages/openmdao.test-0.10.0-py2.7.egg') 
sys.path.append(openmdao_abspath+'/local/lib/python2.7/site-packages/openmdao.units-0.10.0-py2.7.egg')
sys.path.append(openmdao_abspath+'/local/lib/python2.7/site-packages/openmdao.util-0.10.0-py2.7.egg') 
sys.path.append(openmdao_abspath+'/local/lib/python2.7/site-packages/ordereddict-1.1-py2.7.egg') 
sys.path.append(openmdao_abspath+'/local/lib/python2.7/site-packages/pycrypto-2.3-py2.7-linux-x86_64.egg') 
sys.path.append(openmdao_abspath+'/local/lib/python2.7/site-packages/pyparsing-1.5.7-py2.7.egg') 
sys.path.append(openmdao_abspath+'/local/lib/python2.7/site-packages/slsqp-1.0.1-py2.7-linux-x86_64.egg') 
sys.path.append(openmdao_abspath+'/local/lib/python2.7/site-packages/traits-4.3.0-py2.7-linux-x86_64.egg') 
sys.path.append(openmdao_abspath+'/local/lib/python2.7/site-packages/PyYAML-3.10-py2.7-linux-x86_64.egg') 
sys.path.append(openmdao_abspath+'/local/lib/python2.7/site-packages/argh-0.15.1-py2.7.egg') 
sys.path.append(openmdao_abspath+'/local/lib/python2.7/site-packages/openmdao.gui-0.10.0-py2.7.egg')
sys.path.append(openmdao_abspath+'/local/lib/python2.7/site-packages/pathtools-0.1.2-py2.7.egg') 
sys.path.append(openmdao_abspath+'/local/lib/python2.7/site-packages/pyV3D-0.4.4-py2.7-linux-x86_64.egg') 
sys.path.append(openmdao_abspath+'/local/lib/python2.7/site-packages/pyzmq-13.1.0-py2.7-linux-x86_64.egg') 
sys.path.append(openmdao_abspath+'/local/lib/python2.7/site-packages/tornado-2.2.1-py2.7.egg') 
sys.path.append(openmdao_abspath+'/local/lib/python2.7/site-packages/watchdog-0.6.0-py2.7.egg') 
sys.path.append(openmdao_abspath+'/local/lib/python2.7/site-packages/EasyProcess-0.1.4-py2.7.egg') 
sys.path.append(openmdao_abspath+'/local/lib/python2.7/site-packages/PyVirtualDisplay-0.1.0-py2.7.egg')
sys.path.append(openmdao_abspath+'/local/lib/python2.7/site-packages/entrypoint2-0.0.5-py2.7.egg') 
sys.path.append(openmdao_abspath+'/local/lib/python2.7/site-packages/lazr.testing-0.1.2a-py2.7.egg') 
sys.path.append(openmdao_abspath+'/local/lib/python2.7/site-packages/mocker-1.1-py2.7.egg') 
sys.path.append(openmdao_abspath+'/local/lib/python2.7/site-packages/path.py-2.2.2-py2.7.egg') 
sys.path.append(openmdao_abspath+'/local/lib/python2.7/site-packages/selenium-2.35.0-py2.7.egg') 
sys.path.append(openmdao_abspath+'/local/lib/python2.7/site-packages/zope.exceptions-3.6.1-py2.7.egg')
sys.path.append(openmdao_abspath+'/local/lib/python2.7/site-packages/zope.testing-4.1.1-py2.7.egg')
sys.path.append(openmdao_abspath+'/local/lib/python2.7/site-packages/zope.testrunner-4.0.4-py2.7.egg')
sys.path.append(openmdao_abspath+'/lib/python2.7/site-packages/setuptools-0.9.5-py2.7.egg') 
sys.path.append(openmdao_abspath+'/lib/python2.7/site-packages/Jinja2-2.4-py2.7.egg')
sys.path.append(openmdao_abspath+'/lib/python2.7/site-packages/Pygments-1.3.1-py2.7.egg') 
sys.path.append(openmdao_abspath+'/lib/python2.7/site-packages/Sphinx-1.2.2-py2.7.egg') 
sys.path.append(openmdao_abspath+'/lib/python2.7/site-packages/decorator-3.2.0-py2.7.egg')
sys.path.append(openmdao_abspath+'/lib/python2.7/site-packages/docutils-0.10-py2.7.egg')
sys.path.append(openmdao_abspath+'/lib/python2.7/site-packages/pycrypto-2.3-py2.7-linux-x86_64.egg')
sys.path.append(openmdao_abspath+'/lib/python2.7/site-packages/pyparsing-1.5.7-py2.7.egg')
sys.path.append(openmdao_abspath+'/lib/python2.7') 
sys.path.append(openmdao_abspath+'/lib/python2.7/plat-linux2') 
sys.path.append(openmdao_abspath+'/lib/python2.7/lib-tk') 
sys.path.append(openmdao_abspath+'/lib/python2.7/lib-old') 
sys.path.append(openmdao_abspath+'/lib/python2.7/lib-dynload') 
sys.path.append(openmdao_abspath+'/local/lib/python2.7/site-packages') 
sys.path.append(openmdao_abspath+'/lib/python2.7/site-packages') 
sys.path.append(openmdao_abspath+'/local/lib/python2.7/site-packages/setuptools-0.9.5-py2.7.egg') 













