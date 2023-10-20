# Python setup

from distutils.core import setup, Extension
from numpy import get_include

numpy_inc = get_include()		#  NumPy include path.
objs = "sht_init.o sht_kernels_a.o sht_kernels_s.o sht_fly.o sht_omp.o"
shtns_o = objs.split()			# transform to list of objects
libdir = "/home/wulff/local"
if len(libdir) == 0:
	libdir = []
else:
	libdir = [libdir+"/lib"]
cargs = "-std=c99 -fopenmp"
largs = " -L/opt/local/intel/ps_xe-2019.4/compilers_and_libraries_2019.4.243/linux/mkl/lib/intel64"
libs = "-Wl,--start-group -lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core -Wl,--end-group -lm "
libslist = libs.replace('-l','').split()	# transform to list of libraries

shtns_module = Extension('_shtns', sources=['shtns_numpy_wrap.c'],
	extra_objects=shtns_o, depends=shtns_o,
	extra_compile_args=cargs.split(),
	extra_link_args=largs.split(),
	library_dirs=libdir,
	libraries=libslist,
	include_dirs=[numpy_inc])

setup(name='SHTns',
	version='3.4.5',
	description='High performance Spherical Harmonic Transform',
	license='CeCILL',
	author='Nathanael Schaeffer',
	author_email='nathanael.schaeffer@univ-grenoble-alpes.fr',
	url='https://bitbucket.org/nschaeff/shtns',
	ext_modules=[shtns_module],
	py_modules=["shtns"],
	requires=["numpy"],
	)
