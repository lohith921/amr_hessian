function [u_new] = interpolate_p(o_mesh, o_sol, n_mesh)
% This function interpolates mesh solution from coarse mesh to the
% refined mesh. It takes 2 inputs, 1: coarse mesh, 2: refined mesh.

n_ele = length(n_mesh.tris);                % # of elements in new mesh.
n_nodes = length(n_mesh.coords);            % # of nodes in new mesh.
N = zeros(n_nodes,n_nodes);                     % Nodal basis functions.
u_new = zeros(n_nodes,1);                	% Solution vector for new mesh.
u_old = -1*ones(n_nodes,1);                 % Solution vector for old mesh
W = zeros(n_nodes, n_nodes);                % Weight matrix
% Compute local stiffness matrices and then assemble then into global stiffness matrix.
for k = 1:n_ele
  tri_nodes = n_mesh.tris(k,:);               % contains node #'s of this triangle.
  g1 = tri_nodes(1);    g2 = tri_nodes(2);      g3 = tri_nodes(3);
  
  X = n_mesh.coords(tri_nodes,1);             % x coordinates for vertices.
  Y = n_mesh.coords(tri_nodes,2);             % y coordinates for vertices.
  x1 = X(1);  y1 = Y(1);
  x2 = X(2);  y2 = Y(2);
  x3 = X(3);  y3 = Y(3);
  % Compute xg , yg coordinates for Gauss points.
  xg1 = (x1 + x2)/2; yg1 = (y1 + y2)/2;
  xg2 = (x2 + x3)/2; yg2 = (y2 + y3)/2;
  xg3 = (x3 + x1)/2; yg3 = (y3 + y1)/2;
  
  dx21 = (x2 - x1); dy21 = (y2 - y1);
  dx31 = (x3 - x1); dy31 = (y3 - y1);
  dx32 = (x3 - x2); dy32 = (y3 - y2);
  
  A = (dx21*dy31 - dx31*dy21)/2;
  denom = 4*A;
  % Find shape function values for Gauss points.
  N1Xg1 = 1 - (dx21*dy32 - dy21*(dx31 + dx32))/denom;
  N1Xg2 = 1 - ((dx21 + dx31) + (dy21 + dy31)*(x1 + x2))/denom;
  N1Xg3 = 1 + (dx31*dy21 + (dy31*dx32))/denom;
  
  N2Xg1 = (dx21*dy31 - dx31*dy21)/denom;
  N2Xg2 = ((dx21 + dx31)*dy31 - (dy21+dy31)*dx31)/denom;
  N2Xg3 = 0;
  
  N3Xg1 = -(dy31*dx32 + dx31*dy21)/denom;
  N3Xg2 = -((dx21 + dx31)*dy21 + (dy31 + dy21)*dx32)/denom;
  N3Xg3 = -(dy21*dx31 + dy31*dx32)/denom;
  
  % computing N matrix
  N(g1,g1) = N(g1,g1) + N1Xg1;      N(g1,g2) = N(g1,g2) + N1Xg2;       N(g1,g3) = N(g1,g3) + N1Xg3;
  N(g2,g1) = N(g2,g1) + N2Xg1;      N(g2,g2) = N(g2,g2) + N2Xg2;       N(g2,g3) = N(g2,g3) + N2Xg3;
  N(g3,g1) = N(g3,g1) + N3Xg1;      N(g3,g2) = N(g3,g2) + N3Xg2;       N(g3,g3) = N(g3,g3) + N3Xg3;
  
  W(g1,g1) = W(g1,g1) + 1/3;        W(g1,g2) = W(g1,g2) + 1/3;         W(g1,g3) = W(g1,g3) + 1/3;
  W(g2,g1) = W(g2,g1) + 1/3;        W(g2,g2) = W(g1,g2) + 1/3;         W(g2,g3) = W(g2,g3) + 1/3;
  W(g3,g1) = W(g3,g1) + 1/3;        W(g3,g2) = W(g3,g2) + 1/3;         W(g3,g3) = W(g3,g3) + 1/3;
  
  if( u_old(g1) == -1)
      u_old(g1) = compute_uold(o_mesh, o_sol,  [xg1; yg1]);
  end
  if(u_old(g2) == -1)
      u_old(g2) = compute_uold(o_mesh, o_sol, [xg2; yg2]);
  end
  if(u_old(g3) == -1)
      u_old(g3) = compute_uold(o_mesh, o_sol, [xg3; yg3]);
  end
end

M1 = (N.')*W*N;
Q = (N.')*W*u_old;
u_new = M1\Q;
end

% % Determine interior nodes 
% interior = find(n_mesh.border < 1);
% counter = length(interior);
% i = counter;
% % Determine boundary nodes 
% boundary = find(n_mesh.border > 0);
% b = length(boundary);
% Nodal basis functions.
% N = zeros(i+b, i+b);
% u_new = zeros(i+b, 1);
% u_old = -1*ones(i+b,1);
% % Weight matrix
% W = zeros(i+b, i+b);
% % Global ordering for nodes.
% glbal = [interior' boundary'];

%   % Get global ordering for this tri.
%   g1 = find(glbal == tri_nodes(1));
%   g2 = find(glbal == tri_nodes(2));
%   g3 = find(glbal == tri_nodes(3));
