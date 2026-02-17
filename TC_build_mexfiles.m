%% choose compitable fortran compiler for Matlab
% more information please refer to https://de.mathworks.com/support/compilers.html
mex -setup Fortran

%% build mex file for polyhedroal forward modelling coding in fortran
mex polyhedron_full.F
mex polyhedron_g.F
mex polyhedron_g_gx.F
mex polyhedron_gx.F
mex polyhedron_gxx.F
