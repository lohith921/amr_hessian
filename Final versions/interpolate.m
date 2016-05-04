function [u_new] = interpolate(o_mesh, n_mesh)
% This function interpolates mesh solution from coarse mesh to the
% refined mesh. It takes 2 inputs, 1: coarse mesh, 2: refined mesh.

% Determine interior nodes 
interior = find(n_mesh.border < 1);
counter = length(interior);
i = counter;
% Determine boundary nodes 
boundary = find(n_mesh.border > 0);
b = length(boundary);
% Nodal basis functions.
N = zeros(i+b, i+b);
u_new = zeros(i+b, 1);
% Weight matrix
W = zeros(i+b, i+b);
% Global ordering for nodes.
glbal = [interior' boundary'];
% Compute local stiffness matrices and then assemble then into 
% global stiffness matrix.
for k = 1:length(n_mesh.tris)
  % Get this tri nodes. mytris~~ my triangle.
  tri_nodes = mesh.tris(k,:);
  % Get global ordering for this tri.
  g1 = find(glbal == tri_nodes(1));
  g2 = find(glbal == tri_nodes(2));
  g3 = find(glbal == tri_nodes(3));
  % Get x, y coordinates for vertices.
  X = mesh.coords(tri_nodes,1);
  Y = mesh.coords(tri_nodes,2);
  x1 = X(1);  y1 = Y(1);
  x2 = X(2);  y2 = Y(2);
  x3 = X(3);  y3 = Y(3);
  % Compute xg , yg coordinates for Gauss points.
%   xg1 = (x1 + x2)/2; yg1 = (y1 + y2)/2;
%   xg2 = (x2 + x3)/2; yg2 = (y2 + y3)/2;
%   xg3 = (x3 + x1)/2; yg3 = (y3 + y1)/2;
  
  dx21 = (x2 - x1); dy21 = (y2 - y1);
  dx31 = (x3 - x1); dy31 = (y3 - y1);
  dx32 = (x3 - x2); dy32 = (y3 - y2);
  
  m32 = dy32/dx32;          m31 = dy31/dx31;
  m21 = dy21/dx21;          
  
  A = (dx21*dy31 - dx31*dy21)/2;
  denom = 4*A;
  % Find shape function values for Gauss points.
  N1Xg1 = 1 + ((m32 - m21)*dx21*dx32)/denom;
  N1Xg2 = 1 - dx32*(dx21*(m32 - m21) + dx31*(m32 - m31))/denom;
  N1Xg3 = 1 - (dx31*dx32*(m32 - m31))/denom;
  
  N2Xg1 = dx21*dx31*(m31 - m21)/denom;
  N2Xg2 = N2Xg1;
  N2Xg3 = 0;
  
  N3Xg1 = -dx31*dy21/denom;
  N3Xg2 = -(2*dx31*dy21 + dy31*dx32)/denom;
  N3Xg3 = -(dy21*dx31 + dy31*dx32)/denom;
  
  % computing N matrix
  N(g1,g1) = N(g1,g1) + N1Xg1;      N(g1,g2) = N(g1,g2) + N1Xg2;       N(g1,g3) = N(g1,g3) + N1Xg3;
  N(g2,g1) = N(g2,g1) + N2Xg1;      N(g2,g2) = N(g2,g2) + N2Xg2;       N(g2,g3) = N(g2,g3) + N2Xg3;
  N(g3,g1) = N(g3,g1) + N3Xg1;      N(g3,g2) = N(g3,g2) + N3Xg2;       N(g3,g3) = N(g3,g3) + N3Xg3;
  
  W(g1,g1) = W(g1,g1) + 1/3;        W(g1,g2) = W(g1,g2) + 1/3;         W(g1,g3) = W(g1,g3) + 1/3;
  W(g2,g1) = W(g2,g1) + 1/3;        W(g2,g2) = W(g1,g2) + 1/3;         W(g2,g3) = W(g2,g3) + 1/3;
  W(g3,g1) = W(g3,g1) + 1/3;        W(g3,g2) = W(g3,g2) + 1/3;         W(g3,g3) = W(g3,g3) + 1/3;
  % Compute shape functions for this triangle.
  % Compute phi1, phi2, phi3.
  denom = x3*y1 - x3*y2 - x2*y1 - x1*y3 + x1*y2 + x2*y3;
  m1 = (-y3+y2)/denom;
  n1 = -(x2-x3)/denom;
  m2 = (-y1+y3)/denom;
  n2 = (-x3+x1)/denom;
  m3 = -(-y1+y2)/denom;
  n3 = -(x1-x2)/denom;
  
  alpha_aniso = 1;
  beta_aniso = 1;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    
  v1 = [alpha_aniso.*m1; beta_aniso.*n1];
  v2 = [alpha_aniso.*m2; beta_aniso.*n2];
  v3 = [alpha_aniso.*m3; beta_aniso.*n3];
  % Compute area.
  delta = abs(det([x1 y1 1; x2 y2 1; x3 y3 1]));
  % Compute dot products of gradients.
  A(g1,g1) = A(g1,g1)+dot(v1,v1)*0.5*delta;
  A(g1,g2) = A(g1,g2)+dot(v1,v2)*0.5*delta;
  A(g1,g3) = A(g1,g3)+dot(v1,v3)*0.5*delta;
  A(g2,g2) = A(g2,g2)+dot(v2,v2)*0.5*delta;
  A(g2,g3) = A(g2,g3)+dot(v2,v3)*0.5*delta;
  A(g3,g3) = A(g3,g3)+dot(v3,v3)*0.5*delta;
  A(g2,g1) = A(g1,g2);
  A(g3,g1) = A(g1,g3);
  A(g3,g2) = A(g2,g3);

    

end

