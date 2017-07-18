PATH_FORTRAN = '/usr/local/opt/gcc/lib/gcc/7'

cc = ['mex -v qld_interface.cpp qld.f ', ...
      '-L', PATH_FORTRAN, ' -lgfortran'];

eval(cc);
