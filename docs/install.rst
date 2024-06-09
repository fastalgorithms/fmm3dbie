.. role::matlab(code)
   :language: matlab

Get fmm3dBIE
=============

Installation from source
*************************

Installation directly from source is not recommended for Windows machines.
See source code with binaries below.

- Get the source code:

  * The source code can be downloaded from https://github.com/fastalgorithms/fmm3dbie

  * To get all dependencies along with the code run::

    git clone --recurse-submodules https://github.com/fastalgorithms/fmm3dbie.git

Install with Precompiled Windows Binaries
------------------------------------------

Windows compilation of fmm3dbie binaries and mex files tends to be a complicated
process. For intel systems, we have zip files with `pre-compiled Windows binaries
for matlab available <https://github.com/fastalgorithms/fmm3dbie/releases/tag/v1.0.0>`_.


Instructions:

- The x86 version should work on most intel machines. The avx2 specification is also commonly available and will be faster if it is compatible with your machine.

- Download the file and unzip. We recommend testing the fmm binaries by running the following in the matlab directory of fmm3dbie:

.. code:: matlab
   
   startup
   cd tests/
   runtests

Dependencies
************

This library is supported for unix/linux, Mac OSX, and Windows.

For the basic libraries

* Fortran compiler, such as ``gfortran`` packaged with GCC
* GNU make
* Blocked linear algebra routines and LAPACK (default used: Netlib BLAS
  on Linux machines, and framework accelerate on MacOS)

Optional:

* for building Python wrappers you will need ``python`` and ``pip``
* for building standard MATLAB wrappers: MATLAB
* for modifying MATLAB wrappers (experts only): ``mwrap``

Quick install instructions
*********************************************

Make sure you have dependencies downloaded, and `cd` into your fmm3dbie
directory. 

-  For linux, run ``make clean`` followed by ``make install``.
-  For Intel Mac OSX, run ``cp make.inc.macos.gnu make.inc`` followed by ``make clean`` and ``make install``.
-  For M1/M2/M3 Mac OSX, run ``cp make.inc.macos_arm.gnu make.inc`` followed by ``make clean`` and ``make install``.

This should compile the static library for FMM3D in ``FMM3D/lib-static/``, 
the static library for fmm3dbie in ``lib-static/``, 
the dynamic library for fmm3dbie in ``lib/`` and copy the dynamic 
library to ``$(HOME)/lib`` on Linux, and to ``/usr/local/lib`` on Mac OSX.
The location of the default installation directory can be changed by
running::

    make install PREFIX=(INSTALL_DIR)

If the FMM3D library was installed in a non-standard directory, you can
update the directory for where the library looks for ``libfmm3d.a`` by
running::
    
    make install PREFIX_FMM=(FMM_INSTALL_DIR)


In order to link against the dynamic library, you will have to update
the ``LD_LIBRARY_PATH`` environment
variable on Linux and ``DYLD_LIBRARY_PATH`` environment variable on Mac OSX
to the installation directory.
You may then link to the fmm3dbie library using the ``-L$(FMMBIE_INSTALL_DIR) -lfmm3dbie`` 
option.

.. note :: 
   On MacOSX, /usr/local/lib is included by default in the
   DYLD_LIBRARY_PATH.


To verify successful compilation of the program, run ``make test``
which compiles some fortran test drivers in ``test/`` linked against
the static library, after which it
runs the test programs. The last 9 lines of the terminal output should be::

   cat print_testres.txt
   Successfully completed 8 out of 8 tests in common testing suite
   Successfully completed 2 out of 2 tests in helm_wrappers testing suite
   Successfully completed 2 out of 2 tests in lap_wrappers testing suite
   Successfully completed 2 out of 2 tests in surface routs testing suite
   Successfully completed 27 out of 27 tests in tria_routs testing suite
   Successfully completed 24 out of 24 tests in quad_routs testing suite
   Successfully completed 1 out of 1 tests in quadratures testing suite
   rm print_testres.txt


To verify successful installation of the program, and the correct
setting for environment variables, run ``make test-dyn`` which compiles
some fortran test drivers in ``test/`` linked against the dynamic
library, after which it runs teh test prgram. The output ofshould be the
same as above.


.. note ::
   By default, ``make install`` creates the multithreaded version of the library. To
   compile the library in single-threaded mode, append
   ``OMP=OFF`` to the make task. For instance ``make install`` should be replaced by 
   ``make install OMP=OFF``. 
   

If ``make test`` fails, see more detailed instructions below. 

If ``make test-dyn`` fails with an error about not finding ``-lfmm3dbie`` 
make sure that the appropriate environment variables
have been set. If it fails with other issues, see more detailed
instructions below. 

Type ``make`` to see a list of other build options (language
interfaces, etc). Please see ``examples/`` for sample drivers.

If there is an error in testing on a standard set-up,
please file a bug report as a New Issue at https://github.com/fastalgorithms/fmm3dbie/issues

Examples
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``examples`` directory is a good place to see usage 
examples for Fortran. There are several sample drivers for Laplace, Helmholtz,
Stokes and Maxwell solvers.

For example, in the ``helmholtz`` directory, the sample drivers are
``helm_dir_fds_example.f``, ``helm_dir_iter_example.f``, and
``helm_dir_iter_example2.f``, and the corresponding makefiles
are ``helm_dir_fds_example.make``, ``helm_dir_iter_example.make``, and
``helm_dir_iter_example2.make``. These demonstrate how to link
to the dynamic library ``libfmm3dbie.so``. The first is an example for
using the matrix entry generators which are needed by fast direct
solvers. The last two are examples of iterative solvers - the geometry
specification is different in both examples. ``helm_dir_iter_example.f``
generates a triangulation of a surface which is known analytically,
using some basic routines from the ``xtri`` library, while
``helm_dir_iter_example2.f`` reads in a .go3 file of a triangulated
sphere. The repository comes with 2 .go3 triangulations of the sphere. 


Building Python wrappers
****************************

First make sure you have python3 and pip3 installed. 

You may then execute ``make python`` (after copying over the
operating system specific make.inc.* file to make.inc) which calls
pip for the install. 

See ``python/helm3d_dirsolver_demo.py`` to see
usage examples for the Python wrappers.


A few words about Python environments
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

There can be confusion and conflicts between various versions of Python and installed packages. It is therefore a very good idea to use virtual environments. Here's a simple way to do it (after installing python-virtualenv)::

  Open a terminal
  virtualenv -p /usr/bin/python3 env1
  . env1/bin/activate

Now you are in a virtual environment that starts from scratch. All pip installed packages will go inside the env1 directory. (You can get out of the environment by typing ``deactivate``)


Building the MATLAB wrappers
******************************

First make sure you have MATLAB installed. 

Then run ``make matlab`` (after copying over the operating
system specific make.inc.* file to make.inc) which links the .m files to
the .c file in the matlab folder.

To set the relevant paths, run ``startup`` in the ``matlab`` folder.

To run tests, you can run ``runtests`` in the ``matlab/tests`` 
directory and it should return::


   Totals:
   23 Passed, 0 Failed, 0 Incomplete.
   <time taken> seconds testing time.

Example codes for available in the ``matlab/demo`` folder.

Checkout our `MATLAB user guide <matlab_user_guide.html>`__  for info on how to use the MATLAB interface.

Tips for installing dependencies
**********************************

On Ubuntu linux
~~~~~~~~~~~~~~~~

On Ubuntu linux (assuming python3 as opposed to python)::

  sudo apt-get install make build-essential gfortran libopenblas-dev 


On Fedora/CentOS linux
~~~~~~~~~~~~~~~~~~~~~~~~

On a Fedora/CentOS linux system, these dependencies can be installed as 
follows::

  sudo yum install make gcc gcc-c++ gcc-gfortran libgomp openblas-devel 

.. _mac-inst:

On Mac OSX
~~~~~~~~~~~~~~~~~~~~~~~~

First setup Homebrew as follows. If you don't have Xcode, install
Command Line Tools by opening a terminal (from /Applications/Utilities/)
and typing::

  xcode-select --install

Then install Homebrew by pasting the installation command from
https://brew.sh

Then do::
  
  brew install gcc openblas 
  

Tips for installing optional dependencies
******************************************

Installing python3 and pip3
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

On Ubuntu linux
##################

::

  sudo apt-get install python3 python3-pip


On Mac OSX
############

Make sure you have homebrew installed. See `Tips for installing dependencies -> On Mac OSX <install.html#mac-inst>`__ 

::
  
  brew install python3

Then use `make python3` instead of `make python`. You will only need to
do this in case the default version of `python` and `pip` is not >=3.0 


Installing MWrap
~~~~~~~~~~~~~~~~~~

If you make any changes to the 
fortran code, you will need to regenerate the .c files
from the .mw files for which mwrap is required.
This is not needed for most users.
`MWrap <http://www.cs.cornell.edu/~bindel/sw/mwrap>`_
is a very useful MEX interface generator by Dave Bindel.

Make sure you have ``flex`` and ``bison`` installed.
Download version 0.33.5 or later from https://github.com/zgimbutas/mwrap, un-tar the package, cd into it, then::
  
  make
  sudo cp mwrap /usr/local/bin/


