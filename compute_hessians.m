% The following code is going to assemble Load Matrix using construct_A,
% solve using SPSD solvers and compute hessians using Quadratic fitting.
% This is all for now, refinement and interpolation will be added soon.
clear all;
close all;
file_name = input('Please input the root of mesh files', 'S');
mesh = readmesh2(file_name);
[A, b] = construct_A(mesh);
% L is lower triangular matrix from incomplete cholesky factorization
L = incomplete_cholesky(A);
% Preconditioned CGM will construct preconditioner M from L.
sol = preconditioned_CGM( A, b, L );

% Now I have to call Quadratic fit.

