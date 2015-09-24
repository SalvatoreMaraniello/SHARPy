/usr/bin/f2py2.7 -h lib_fem_f2py.pyf -m lib_fem lib_fem_f2py.f90 --overwrite-signature
/usr/bin/f2py2.7 -c -m lib_fem lib_fem_f2py.pyf lib_fem_f2py.f90

/usr/bin/f2py2.7 -h lib_cbeam3_f2py.pyf -m lib_cbeam3 lib_cbeam3_f2py.f90 --overwrite-signature
/usr/bin/f2py2.7 -c -m lib_cbeam3 -I../ -L../../lib/src/ lib_cbeam3_f2py.pyf lib_cbeam3_f2py.f90 ../../lib/src/lib_fem.o ../../lib/src/lib_rot.o ../../lib/src/lib_rotvect.o

/usr/bin/f2py2.7 -h lib_rotvect_f2py.pyf -m lib_rotvect lib_rotvect_f2py.f90 --overwrite-signature
/usr/bin/f2py2.7 -c -m lib_rotvect -I../ -L../../lib/src/ lib_rotvect_f2py.pyf lib_rotvect_f2py.f90 ../../lib/src/lib_rot.o
