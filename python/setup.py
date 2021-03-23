import setuptools
import string
import os
from numpy.distutils.core import setup
from numpy.distutils.core import Extension
from sys import platform

pkg_name = "bie-solvers3dpy"

list_files=[]
list_files.append('../src/helm_wrappers/helm_comb_dir.f')
list_files.append('../src/common/sparse_reps.f')
list_files.append('../src/surface_routs/surf_routs.f90')
list_files.append('../src/surface_routs/vtk_routs.f90')
list_files.append('../src/stok_wrappers/stok_comb_vel.f')
list_files.append('../src/kernels/stok_kernels.f90')
list_files.append('../src/surface_routs/write_go3.f90')

FLIBS = os.getenv('FMMBIE_LIBS')
FLIBS = FLIBS.rstrip().split(' ')
FLIBS = list(filter(None,FLIBS))
FLIBS.append('../lib-static/libfmm3dbie.a')
if platform == "darwin":
    FLIBS.append('-L/usr/local/lib')
if platform == "linux" or platform == "linux2":
    FLIBS.append('-lopenblas')

helm = []
com = []
surf = []
stok = []
helm.append('helm_comb_dir_fds_csc_mem')
helm.append('helm_comb_dir_fds_csc_init')
helm.append('helm_comb_dir_fds_csc_matgen')
helm.append('helm_comb_dir_fds_block_mem')
helm.append('helm_comb_dir_fds_block_init')
helm.append('helm_comb_dir_fds_block_matgen')
helm.append('lpcomp_helm_comb_dir')
com.append('conv_to_csc')
surf.append('surf_vals_to_coefs')
surf.append('get_qwts')
surf.append('get_patch_id_uvs')
surf.append('surf_vtk_plot')
surf.append('surf_vtk_plot_scalar')
surf.append('surf_vtk_plot_vec')
stok.append('lpcomp_stok_comb_vel')
stok.append('stok_comb_vel_solver')
stok.append('stok_comb_vel_matgen')
stok.append('st3d_slp_vec')
stok.append('st3d_slp')
stok.append('st3d_dlp_vec')
stok.append('st3d_dlp')
stok.append('st3d_comb_vec')
stok.append('st3d_comb')
stok.append('st3d_strac_vec')
stok.append('st3d_strac')
surf.append('get_wtorus_geom')
surf.append('write_wtorus')
surf.append('get_sphere_geom')
surf.append('get_patch_distortion')


ext_helm = Extension(
    name='fmm3dbie',
    sources=list_files,
    f2py_options=['only:']+helm+com+surf+stok+[':'],
#    extra_f77_compile_args=FFLAGS,
#    extra_f90_compile_args=FFLAGS,
    extra_link_args=FLIBS
)

## TODO: fill in the info below
setup(
    name=pkg_name,
    version="0.1.0",
    author="Manas Rachh",
    author_email="mrachh@flatironinstitute.org",
    description="This pacakge contains basic routines Helmholtz dirichlet fast direct solver",
    url="",
    packages=setuptools.find_packages(),
    install_requires=[
        "pytest"
    ],
    ext_modules=[ext_helm],
    classifiers=(
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: Apache Software License",
        "Operating System :: OS Independent",
    )    
)
