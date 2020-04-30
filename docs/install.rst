Installation
============

Obtaining fmm3dbie
******************

The source code can be downloaded from https://gitlab.com/fastalgorithms/fmm3dbie 


Dependencies
************

This library is supported for unix/linux, and Mac OSX

For the basic libraries

* Fortran compiler, such as ``gfortran`` packaged with GCC
* GNU make
* `FMM3D <https://github.com/flatironinstitute/FMM3D>`_
* Blocked linear algebra routines and LAPACK (default used: Netlib BLAS
on Linux machines, and framework accelerate on MacOS)

Optional:

* for building Python wrappers you will need ``python`` and ``pip`` 

Quick install instructions
*********************************************

Make sure you have dependencies downloaded, and `cd` into your fmm3dbie
directory. 

-  For linux, run ``make test``.
-  For Mac OSX, run ``cp make.inc.macosx.gnu make.inc`` followed by ``make test``.

This should compile the static library
in ``lib-static/`` and some fortran test drivers in ``test/``, after which it
runs the test programs. The last 5 lines of the terminal output should be::

   cat print_testres.txt
   Successfully completed 3 out of 3 tests in common testing suite
   Successfully completed 2 out of 2 tests in helm_wrappers testing suite
   Successfully completed 27 out of 27 tests in tria_routs testing suite
   rm print_testres.txt

.. note ::
   By default, ``make test`` creates the multithreaded version of the library. To
   compile the library in single-threaded mode, append
   ``OMP=OFF`` to the make task. For instance ``make test`` should be replaced by 
   ``make test OMP=OFF``. 
   

If ``make test`` fails, see more detailed instructions below. If it succeeds, run
``make lib``, which creates the dynamic library (``libfmm3dbie.so``). You may then
link to the FMM library using the ``-lfmm3dbie`` option.

.. note ::
   On unix/Linux systems, you will need to include the location of
   libfmm3dbie.so in the environment variable LD_LIBRARY_PATH. 
   This can be set temporarily using::

      export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:<path-to-libfmm3dbie> 



.. note :: 
   On MacOSX, in order to link with the dynamic libraries, you will
   need to copy libfmm3d.so to ``usr/local/lib``. See any of the
   makefiles in the ``examples/`` directory for prototypes.

Type ``make`` to see a list of other build options (language
interfaces, etc). Please see ``examples/`` for sample drivers.

If there is an error in testing on a standard set-up,
please file a bug report as a New Issue at https://gitlab.com/fastalgorithms/fmm3dbie/issues

Examples
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``examples`` directory is a good place to see usage 
examples for Fortran.
There are three sample Fortran drivers for both iterative and direct
solvers for the Helmholtz Dirichlet problem. 

The sample drivers are
``helm_dir/helm_dir_fds_example.f``, ``helm_dir/helm_dir_iter_example.f``, and
``helm_dir/helm_dir_iter_example2.f``, and the corresponding makefiles
are ``helm_dir/helm_dir_fds_example.make``, ``helm_dir/helm_dir_iter_example.make``, and
``helm_dir/helm_dir_iter_example2.make``. These demonstrate how to link
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


