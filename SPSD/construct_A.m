% Script to compute the weights using discretization of Laplace's equation.
clear all;
close all;
warning off;
% Get input
%filebase=input('Enter the name of the mesh:  ');
%filebase=input('shape');
% Read-in initial mesh
% mesh=readmesh2('bar_i');
%mesh=readmesh2('a_shape.1');
mesh = readmesh2('annulus_initial');
% mesh=readmesh2('four_circles');
%mesh=readmesh2('cylinder_boundary_smoothed_initial');
%figure;
%drawmesh2(mesh);
%inverted = evalquality2(mesh,filebase,1,1);
% Number of nodes (200)
n = length(mesh.coords);
% Determine interior nodes (135)
interior = find(mesh.border < 1);
counter = length(interior);
i = counter;
% Determine boundary nodes (65)
boundary = find(mesh.border > 0);
b = length(boundary);
% Create space for A.
A = zeros(i+b,i+b);
B = zeros(i+b,1);
% Global ordering for nodes.
glbal = [interior' boundary'];
% Compute local stiffness matrices and then assemble then into 
% global stiffness matrix.
for k=1:length(mesh.tris)
  % Get this tri.
  mytris = mesh.tris(k,:);
  % Get global ordering for this tri.
  g1 = find(glbal == mytris(1));
  g2 = find(glbal == mytris(2));
  g3 = find(glbal == mytris(3));
  % Get x, y coordinates for vertices.
  x = mesh.coords(mytris,1);
  y = mesh.coords(mytris,2);
  x1 = x(1);
  y1 = y(1);
  x2 = x(2);
  y2 = y(2);
  x3 = x(3);
  y3 = y(3);
  % Compute shape functions for this triangle.
  % Compute phi1, phi2, phi3.
  denom = x3*y1-x3*y2-x2*y1-x1*y3+x1*y2+x2*y3;
  m1 = (-y3+y2)/denom;
  n1 = -(x2-x3)/denom;
  m2 = (-y1+y3)/denom;
  n2 = (-x3+x1)/denom;
  m3 = -(-y1+y2)/denom;
  n3 = -(x1-x2)/denom;
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % jacobian matrix
%   ele_7=[x2-x1 x3-x1;y2-y1 y3-y1];
%   
%   if ( det(ele_7) > 0)
%   [a7,v7, q7,delta7]=tmp(ele_7);
%   asp_ratio(k)=delta7(1,1)/delta7(2,2); 
% 
%   else
%     ele_7=[x3-x1 x2-x1;y3-y1 y2-y1];
%     asp_ratio(k)=delta7(1,1)/delta7(2,2); 
% 
%         
%   end
%   % case 1
%   if (asp_ratio(k) <=1)
%       alpha_aniso=asp_ratio(k);
%       beta_aniso=1;
%   else
%       alpha_aniso=1;
%       beta_aniso=asp_ratio(k);
 % end

  % case 2 (bad)
%    if (asp_ratio(k) <=1)
%       beta_aniso=asp_ratio(k);
%       alpha_aniso=1;
%   else
%       beta_aniso=1;
%       alpha_aniso=asp_ratio(k);
%   end
  
      
    % alpha_aniso=delta7(2,2);
    % beta_aniso=delta7(1,1);
  
  
  alpha_aniso=1;
  beta_aniso=1;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    
  v1 = [alpha_aniso.*m1; beta_aniso.*n1];
  v2 = [alpha_aniso.*m2; beta_aniso.*n2];
  v3 = [alpha_aniso.*m3; beta_aniso.*n3];
  % Compute delta.
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
  
  coeff_1=0.01;
%   
  B(g1) = 2*new_func(x1,y1, coeff_1)/12 + new_func(x2,y2, coeff_1)/12 + new_func(x3,y3, coeff_1)/12;
  B(g2) = new_func(x1,y1, coeff_1)/12 + 2*new_func(x2,y2, coeff_1)/12 + new_func(x3,y3, coeff_1)/12;
  B(g3) = new_func(x1,y1, coeff_1)/12 + new_func(x2,y2, coeff_1)/12 + 2*new_func(x3,y3, coeff_1)/12;
%   B(g1)=0;
%   B(g2)=0;
%   B(g3)=0;

end;
% i
%drawmesh2(mesh);