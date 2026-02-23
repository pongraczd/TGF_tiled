%% choose compitable fortran compiler for Matlab
% more information please refer to https://de.mathworks.com/support/compilers.html
mex -setup Fortran

%% build mex file for polyhedroal forward modelling coding in fortran
mex src\polyhedron_full.F -outdir bin\
mex src\polyhedron_g.F -outdir bin\
mex src\polyhedron_g_gx.F -outdir bin\
mex src\polyhedron_gx.F -outdir bin\
mex src\polyhedron_gxx.F -outdir bin\
