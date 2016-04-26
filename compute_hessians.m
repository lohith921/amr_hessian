% The following code is going to assemble Load Matrix using construct_A,
% solve using SPSD solvers and compute hessians using Quadratic fitting.
% This is all for now, refinement and interpolation will be added soon.
% clear all;
close all;
file_name = input('Please input the root of mesh files', 'S');
mesh = readmesh2(file_name);
% construct_A is used to construct stiffness matrix
[A, b] = construct_A(mesh);
% L is lower triangular matrix from incomplete cholesky factorization
L = incomplete_cholesky(A);
% Preconditioned CGM will construct preconditioner M from L.
sol = preconditioned_CGM( A, b, L );
% M_tensor is the metric tensor at each node. We call Quadratic fit to compute it.
M_tensor = Quadratic_fit_modified_1(mesh.coords, mesh.tris, sol);

% nele is the # of triangles
% nele = length(mesh.tris);
% for i = 1:nele
%     [n1, n2, n3] = mesh.tris(:,i);
%     M1 = M_tensor(:,:,n1);
%     M2 = M_tensor(:,:,n2);
%     M3 = M_tensor(:,:,n3);
%     Mavg = (M1 + M2 + M3)/3;
%     if (det(Mavg) > limit)
%         mid1 = (n1 + n2)/2;
%         mid2 = (n2 + n3)/2;
%         mid3 = (n3 + n1)/2;
%         mesh.tris(i,:) = [n1 mid1 mid3];
%         mesh.tris(nele + 1, :) = [mid1 mid2 mid3];
%         mesh.tris(nele + 2, :) = [mid1 n2 mid2];
%         mesh.tris(nele + 3, :) = [mid2 n3 mid3];
%         nele = nele +3;
%     end
% end



