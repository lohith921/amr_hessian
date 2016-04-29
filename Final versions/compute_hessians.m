% The following code is going to assemble stiffness Matrix using construct_A, solve using SPSD solvers and compute hessians using 
% Quadratic fitting. This is all for now, refinement and interpolation will be added soon.
% clear all;
close all;
%file_name = input('Please input the root of mesh files', 'S');
mesh = readmesh2('simple');
% construct_A is used to construct stiffness matrix
[A, b] = construct_A(mesh);
% L is lower triangular matrix from incomplete cholesky factorization
L = incomplete_cholesky(A);
% Preconditioned CGM will construct preconditioner M from L.
sol = preconditioned_CGM(A, b, L);
% M_tensor is the metric tensor at each node. We call Quadratic fit to compute it.
M_tensor = Quadratic_fit_modified_1(mesh.coords, mesh.tris, sol);
% nele is the # of triangles
nele = length(mesh.tris);
nnodes = length(mesh.coords);
limit = 10^10;
num_e = nele;
for i = 1:num_e
    %[n1, n2, n3] = mesh.tris(i,:);
    n1 = mesh.tris(i,1); n2 = mesh.tris(i,2); n3 = mesh.tris(i,3);
%     n1 
%     n2 
%     n3
    M1 = M_tensor(:,:,n1);
    M2 = M_tensor(:,:,n2);
    M3 = M_tensor(:,:,n3);
    Mavg = (M1 + M2 + M3)/3;
    if (det(Mavg) > limit)
        m1 = nnodes + 1; m2 = nnodes + 2; m3 = nnodes + 3;
        mid1 = (mesh.coords(n1,:) + mesh.coords(n2,:))/2;
        mid2 = (mesh.coords(n2,:) + mesh.coords(n3,:))/2;
        mid3 = (mesh.coords(n3,:) + mesh.coords(n1,:))/2;
                
        mesh.coords(m1,:) = mid1;
        mesh.coords(m2,:) = mid2;
        mesh.coords(m3,:) = mid3;
        
        M_tensor(:,:,m1) = (M1 + M2)/2;
        M_tensor(:,:,m2) = (M2 + M3)/2;
        M_tensor(:,:,m3) = (M3 + M1)/2;
        
        nnodes = nnodes + 3;
        
        [nb1, ot1] = find_neighbor_with(n2, n1, mesh.tris);
        [nb2, ot2] = find_neighbor_with(n3, n2, mesh.tris);
        [nb3, ot3] = find_neighbor_with(n1, n3, mesh.tris);
        
        if ~isempty(nb1)
            nele = nele + 1;
            nb13 = mesh.tris(nb1, ot1);
            mesh.tris(nb1, :) = [m1 n1 nb13];
            mesh.tris(nele, :) = [n2 m1 nb13];
        end
        
        if ~isempty(nb2)
            nele = nele + 1;
            nb21 = mesh.tris(nb2, ot2);
            mesh.tris(nb2, :) = [nb21 m2 n2];
            mesh.tris(nele, :) = [nb21 n3 m2];
        end
        
        if ~isempty(nb3)
            nele = nele + 1;
            nb32 = mesh.tris(nb3, ot3);
            mesh.tris(nb3, :) = [n3 nb32 m3];
            mesh.tris(nele, :) = [n3 nb32 n1];
        end
        
        mesh.tris(i,:) = [n1 m1 m3]; 
        mesh.tris(nele + 1, :) = [m1 m2 m3];
        mesh.tris(nele + 2, :) = [m1 n2 m2];
        mesh.tris(nele + 3, :) = [m2 n3 m3];
        nele = nele + 3;
    end
end

% removed code
% node_file = sprintf('%s.M.%s',file_name,'node');
% ele_file = sprintf('%s.M.%s',file_name,'ele');
% n_file = fopen(node_file,'w');
% e_file = fopen(ele_file,'w');
% fwrite(n_file, [nnodes 
% fwrite(n_file,mesh.coords(n1)
%         mesh.coords(m1,:) = (mesh.coords(n1,:) + mesh.coords(n2,:))/2;
%         mesh.coords(m2,:) = (mesh.coords(n2,:) + mesh.coords(n3,:))/2;
%         mesh.coords(m3,:) = (mesh.coords(n3,:) + mesh.coords(n1,:))/2;
        
%         %d1 = find(mesh.coords(:,:) == mid1(1,:))'%  && mesh.coords(:,2) == mid1(:,2));

%         mesh.tris(i,:) = [n1 m1 m3];
%         
%         nnodes = nnodes + 3;
%         for 1:
