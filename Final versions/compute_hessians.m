function [mesh, old_mesh, sol] = compute_hessians(file_name)
% The following code is going to assemble stiffness Matrix using construct_A, solve using SPSD solvers and compute hessians using 
% Quadratic fitting. This is all for now, refinement and interpolation will be added soon.
% clear all;
    close all;
    %file_name = input('Please input the root of mesh files', 'S');
    mesh = readmesh2(file_name);
    [q_vector] = compute_Quality(mesh);
    figure(1);
    histogram(q_vector);
    title('Histogram of Quality before refinement');
    % construct_A is used to construct stiffness matrix
    % [A,boundary_edges, b] = construct_A(mesh);
    [A, b] = construct_A(mesh);
    % L is lower triangular matrix from incomplete cholesky factorization
    L = incomplete_cholesky(A);
    % Preconditioned CGM will construct preconditioner M from L.
    sol = preconditioned_CGM(A, b, L);
    figure(2);
    title('Solution Plot');
    trisurf((mesh.tris),mesh.coords(:,1),mesh.coords(:,2),0*mesh.coords(:,1),sol,'edgecolor','k','facecolor','interp');
    view(2), axis equal,colorbar;
    
    % M_tensor is the metric tensor at each node. We call Quadratic fit to compute it.
    M_tensor = Quadratic_fit_modified_1(mesh.coords, mesh.tris, sol);
    det_Hess = zeros(length(M_tensor),1);
%     for i = 1:length(M_tensor)
%         det_Hess(i) = det(M_tensor(:,:,i));
%     end
    % nele is the # of triangles
    figure(3);
    title('Mesh before refining');
    drawmesh2_simple(mesh);
    old_mesh = mesh;
    nele = length(mesh.tris);
    nnodes = length(mesh.coords);
    % build a boundary_edge datastructure.
    
    % num_e = nele;
    for i = 1:nele
        n1 = mesh.tris(i,1);        n2 = mesh.tris(i,2);        n3 = mesh.tris(i,3);
        M1 = M_tensor(:,:,n1);      M2 = M_tensor(:,:,n2);      M3 = M_tensor(:,:,n3);
        Mavg = (M1 + M2 + M3)/3;
        if (norm(Mavg,2) > 1)
            m1 = nnodes + 1;        m2 = nnodes + 2;            m3 = nnodes + 3;
        
            mid1 = (mesh.coords(n1,:) + mesh.coords(n2,:))/2;
            mid2 = (mesh.coords(n2,:) + mesh.coords(n3,:))/2;
            mid3 = (mesh.coords(n3,:) + mesh.coords(n1,:))/2;
                
            mesh.coords(m1,:) = mid1;
            mesh.coords(m2,:) = mid2;
            mesh.coords(m3,:) = mid3;
            
            mesh.border([m1; m2; m3]) = 0;
            M_tensor(:,:,m1) = (M1 + M2)/2;     M_tensor(:,:,m2) = (M2 + M3)/2;     M_tensor(:,:,m3) = (M3 + M1)/2;
        
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
                mesh.tris(nb2, :) = [n2 nb21 m2];
                mesh.tris(nele, :) = [n3 m2 nb21];
            end
        
            if ~isempty(nb3)
                nele = nele + 1;
                nb32 = mesh.tris(nb3, ot3);
                mesh.tris(nb3, :) = [n3 nb32 m3];
                mesh.tris(nele, :) = [m3 nb32 n1];
            end
        
            mesh.tris(i,:) = [n1 m1 m3]; 
            mesh.tris(nele + 1, :) = [m1 m2 m3];
            mesh.tris(nele + 2, :) = [m1 n2 m2];
            mesh.tris(nele + 3, :) = [m2 n3 m3];
            nele = nele + 3;
        end
    end
    figure(4)
    title('Mesh after refinement');
    drawmesh2_simple(mesh);
%     [q_new] = compute_Quality(mesh);
%     figure(5);
%     title('Histogram of Quality after refinement');
%     histogram(q_new);
    %
  %  write_file(mesh, file_name);
    % New mesh has been created, now lets solve for u on the new mesh using old mesh.
    %sol_new = interpolate_p(old_mesh, sol, mesh);
    % figure(6);
    % trisurf((mesh.tris),mesh.coords(:,1),mesh.coords(:,2),0*mesh.coords(:,1),u,'edgecolor','k','facecolor','interp');
    % view(2), axis equal,colorbar;
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
%         mesh.border(m1) = mesh.border(n1,1) & mesh.border(n2,1);
%         mesh.border(m2) = mesh.border(n2,1) & mesh.border(n3,1);
%         mesh.border(m3) = mesh.border(n3,1) & mesh.border(n1,1);
%             flag = find(boundary_edges(:,1:2) == [n1 n2]);
%             if (flag > 0)
%                 mesh.border(m1) = mesh.border(n1);
%             end
%             if (find(boundary_edges(:,1:2) == [n2 n3]))
%                 mesh.border(m2) = mesh.border(n2);
%             end
%             if (find(boundary_edges(:,1:2) == [n3 n1]))
%                 mesh.border(m3) = mesh.border(n3);
%             end
%     %% compute boundaries in new mesh
%     edge_structure = zeros(nele, 6);
%     for i = 1:nele
%         g1 = mesh.tris(i,1);        g2 = mesh.tris(i,2);        g3 = mesh.tris(i,3);
%         edge_structure(i, 1:6) = [g1 g2 g2 g3 g3 g1];
%     end
%     mesh.border = zeros(length(mesh.coords),1);
%     for i = 1:nele
%         g1 = mesh.tris(i,1);        g2 = mesh.tris(i,2);        g3 = mesh.tris(i,3);
%         for j= 1:3 % for each edge e1 e2 e3. e1: g1-g2, e2: g2-g3, e3: g3-g1
%             if ~(find(((edge_structure(:,mod(j+1,3)) == g1) & (edge_structure(:,mod(i+2,3)) == g2))))
%                 mesh.border(g1) = 1;            mesh.border(g2) = 1;
%             end
%             if ~(find(((edge_structure(:,mod(i+1,3)) == g2) & (edge_structure(:,mod(i+2,3)) == g3))))
%                 mesh.border(g2) = 1;            mesh.border(g3) = 1;
%             end
%             if ~(find(((edge_structure(:,mod(i+1,3)) == g3) & (edge_structure(:,mod(i+2,3)) == g1))))
%                 mesh.border(g3) = 1;            mesh.border(g1) = 1;
%             end
%         end
%     end
