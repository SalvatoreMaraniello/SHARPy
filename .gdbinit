set stop-on-solib-events 0
set listsize 21

# file to load debug symbols
# /home/rjs10/SharPyProject/BeamLib/bin/
# add-symbol-file BeamLib.so [.. hex mem address ..] -readnow

# put python's source in gdb's source directories
dir /home/rjs10/python3.2-3.2/Python

# put beamlib source in gdb source dir
dir /home/rjs10/SharPyProject/BeamLib/src/fortran/install/src
dir /home/rjs10/SharPyProject/BeamLib/src/fortran/main
dir /home/rjs10/SharPyProject/BeamLib/src/lib/src
dir /home/rjs10/SharPyProject/BeamLib/src/wrapper

show directories

# from PYTHON_X.x_SRC_DIR/Misc/gdbinit
# If you use the GNU debugger gdb to debug the Python C runtime, you
# might find some of the following commands useful.  Copy this to your
# ~/.gdbinit file and it'll get loaded into gdb automatically when you
# start it up.  Then, at the gdb prompt you can do things like:
#
#    (gdb) pyo apyobjectptr
#    <module 'foobar' (built-in)>
#    refcounts: 1
#    address    : 84a7a2c
#    $1 = void
#    (gdb)
#
# NOTE: If you have gdb 7 or later, it supports debugging of Python directly
# with embedded macros that you may find superior to what is in here.
# See Tools/gdb/libpython.py and http://bugs.python.org/issue8032.

# Prints a representation of the object to stderr, along with the
# number of reference counts it current has and the hex address the
# object is allocated at.  The argument must be a PyObject*
define pyo
    # side effect of calling _PyObject_Dump is to dump the object's
    # info - assigning just prevents gdb from printing the
    # NULL return value
    set $_unused_void = _PyObject_Dump($arg0)
end

# Prints a representation of the object to stderr, along with the
# number of reference counts it current has and the hex address the
# object is allocated at.  The argument must be a PyGC_Head*
define pyg
    print _PyGC_Dump($arg0)
end

# print the local variables of the current frame
define pylocals
    set $_i = 0
    while $_i < f->f_code->co_nlocals
	if f->f_localsplus + $_i != 0
	    set $_names = co->co_varnames
	    set $_name = _PyUnicode_AsString(PyTuple_GetItem($_names, $_i))
	    printf "%s:\n", $_name
            pyo f->f_localsplus[$_i]
	end
        set $_i = $_i + 1
    end
end

# A rewrite of the Python interpreter's line number calculator in GDB's
# command language
define lineno
    set $__continue = 1
    set $__co = f->f_code
    set $__lasti = f->f_lasti
    set $__sz = ((PyVarObject *)$__co->co_lnotab)->ob_size/2
    set $__p = (unsigned char *)((PyBytesObject *)$__co->co_lnotab)->ob_sval
    set $__li = $__co->co_firstlineno
    set $__ad = 0
    while ($__sz-1 >= 0 && $__continue)
      set $__sz = $__sz - 1
      set $__ad = $__ad + *$__p
      set $__p = $__p + 1
      if ($__ad > $__lasti)
	set $__continue = 0
      else
        set $__li = $__li + *$__p
        set $__p = $__p + 1
      end
    end
    printf "%d", $__li
end

# print the current frame - verbose
define pyframev
    pyframe
    pylocals
end

define pyframe
    set $__fn = _PyUnicode_AsString(co->co_filename)
    set $__n = _PyUnicode_AsString(co->co_name)
    printf "%s (", $__fn
    lineno
    printf "): %s\n", $__n
### Uncomment these lines when using from within Emacs/XEmacs so it will
### automatically track/display the current Python source line
#    printf "%c%c%s:", 032, 032, $__fn
#    lineno
#    printf ":1\n"
end

### Use these at your own risk.  It appears that a bug in gdb causes it
### to crash in certain circumstances.

#define up
#    up-silently 1
#    printframe
#end

#define down
#    down-silently 1
#    printframe
#end

define printframe
    if $pc > PyEval_EvalFrameEx && $pc < PyEval_EvalCodeEx
	pyframe
    else
        frame
    end
end

# Here's a somewhat fragile way to print the entire Python stack from gdb.
# It's fragile because the tests for the value of $pc depend on the layout
# of specific functions in the C source code.

# Explanation of while and if tests: We want to pop up the stack until we
# land in Py_Main (this is probably an incorrect assumption in an embedded
# interpreter, but the test can be extended by an interested party).  If
# Py_Main <= $pc <= Py_GetArgcArv is true, $pc is in Py_Main(), so the while
# tests succeeds as long as it's not true.  In a similar fashion the if
# statement tests to see if we are in PyEval_EvalFrameEx().

# Note: The name of the main interpreter function and the function which
# follow it has changed over time.  This version of pystack works with this
# version of Python.  If you try using it with older or newer versions of
# the interpreter you may will have to change the functions you compare with
# $pc.

# print the entire Python call stack
define pystack
    while $pc < Py_Main || $pc > Py_GetArgcArgv
        if $pc > PyEval_EvalFrameEx && $pc < PyEval_EvalCodeEx
	    pyframe
        end
        up-silently 1
    end
    select-frame 0
end

# print the entire Python call stack - verbose mode
define pystackv
    while $pc < Py_Main || $pc > Py_GetArgcArgv
        if $pc > PyEval_EvalFrameEx && $pc < PyEval_EvalCodeEx
	    pyframev
        end
        up-silently 1
    end
    select-frame 0
end

# generally useful macro to print a Unicode string
def pu
  set $uni = $arg0 
  set $i = 0
  while (*$uni && $i++<100)
    if (*$uni < 0x80) 
      print *(char*)$uni++
    else
      print /x *(short*)$uni++
    end
  end
end
# end PYTHON_X.x_SRC_DIR/Misc/gdbinit


# used with pyattr to inspect an objects attribute list
# found at the bottom link on:
# http://wiki.python.org/moin/DebuggingWithGdb
define pyattrlist
 set $dict = *(PyDictObject**)(((int)$arg0)+$arg0.ob_type.tp_dictoffset)
 set $i = 0
 set $j = 0
 while $i < $dict.ma_mask
  if 0 != $dict.ma_table[$i].me_value
   echo \nattr#:
   p $j
   x/s ((PyStringObject*)$dict.ma_table[$i].me_key).ob_sval
   pyobjinfo $dict.ma_table[$i].me_value
   set $j = $j+1
  end
  set $i = $i+1
 end
end
document pyattrlist
show pythonic object attributes list
end

# used above
# found at the bottom link on:
# http://wiki.python.org/moin/DebuggingWithGdb
define pyattr
 set $dict = *(PyDictObject**)(((int)$arg0)+$arg0.ob_type.tp_dictoffset)
 set $i = 0
 set $j = 0
 while $i < $dict.ma_mask
  if 0 != $dict.ma_table[$i].me_value
   if $j == $arg1
    set $attr = $dict.ma_table[$i].me_value
    echo $attr:
    p $attr
   end
   set $j = $j+1
  end
  set $i = $i+1
 end
end
document pyattr
get pythonic object attribute
usage: pyattr <pyobj> <attr#>
atr# correlates to those shown in pyattrlist
end

#This function shows the pythonic backtrace.
define pbt
 set $i = 0
 set $j = 0
 while $i < 1000
  select $i
  if $eip >= &PyEval_EvalFrame
   if $eip < &PyEval_EvalCodeEx
    echo c frame #
    p $i
    echo py frame #
    p $j
    set $j = $j+1
    x/s ((PyStringObject*)f->f_code->co_filename)->ob_sval
    x/s ((PyStringObject*)f->f_code->co_name)->ob_sval
    echo line #
    p f->f_lineno
   end
  end
  set $i = $i+1
 end
end
document pbt
show python backtrace
end