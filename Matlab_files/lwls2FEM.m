% Script to compute the weights using discretization of Laplace's equation.
%close all, clear all, clc;
%warning off;
close all;
% Get input
%filebase=input('Enter the name of the mesh:  ');
%filebase=input('shape');

% Read-in initial mesh
mesh=readmesh2('circ_2_4514.1');
%figure;
%drawmesh2(mesh);
%inverted = evalquality2(mesh,filebase,1,1);

% Number of nodes (200)
n = length(mesh.coords);

% Determine interior nodes (135)
interior = find(mesh.border < 1);
counter = length(interior);
i = counter;
%interior
% Determine boundary nodes (65)
boundary = find(mesh.border > 0);
b = length(boundary);
%boundary
% Create space for A.
%fprintf('interior + boundary: %d\n',i+b);
A = zeros(i+b,i+b);

% Global ordering for nodes.
glbal = [interior' boundary'];
%glbal
% Compute local stiffness matrices and then assemble then into 
% global stiffness matrix.
for k=1:length(mesh.tris)

  % Get this tri.
  mytris = mesh.tris(k,:);
  mytris
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

  alpha_aniso=15;
  beta_aniso=2;
  
  
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
  
  coeff_1=-1e7;
  
  B(g1)=2*new_func(x1,y1, coeff_1)/12 + new_func(x2,y2, coeff_1)/12 + new_func(x3,y3, coeff_1)/12;
  B(g2)=new_func(x1,y1, coeff_1)/12 + 2*new_func(x2,y2, coeff_1)/12 + new_func(x3,y3, coeff_1)/12;
  B(g3)=new_func(x1,y1, coeff_1)/12 + new_func(x2,y2, coeff_1)/12 + 2*new_func(x3,y3, coeff_1)/12;
  
end;


% Extract A_I and A_B.
A_I = A(1:i,1:i);
A_B = A(1:i,i+1:end);
%clear A;

% Move mesh.
newmesh = mesh;
ele=length(newmesh.tris);
node=length(newmesh.coords);

% Move boundary of mesh.

% Sample boundary motion.
boundary = find(newmesh.border > 0);

%newmesh.coords(boundary,1) = newmesh.coords(boundary,1)+3;
%newmesh.coords(boundary,2) = newmesh.coords(boundary,2)-2;

%filebase_1=input('Enter the name of the mesh:  ');
%filebase=input('shape');

% Read-in initial mesh
%newmesh=readmesh2('a_shape.1');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% max_x_1=max(mesh.coords(:,1));
% max_y_1=max(mesh.coords(:,2));
% 
% 
% qwqw=1.0;
% 
% newmesh.coords(boundary,1) = -1 .* ( 1+ mesh.coords(boundary,2)./max_y_1) .*cos(pi.*mesh.coords(boundary,1)./max_x_1);
% newmesh.coords(boundary,2) = (1 + (2.*qwqw-1).* mesh.coords(boundary,2)./max_y_1).*sin(pi.*mesh.coords(boundary,1)./max_x_1);
% 

newmesh.coords(boundary,1) =   mesh.coords(boundary,1).^2.* mesh.coords(boundary,2).^2 ;
newmesh.coords(boundary,2) =    mesh.coords(boundary,2) ;



%node_new_xi]=-1*(1+node_y[i])*cos(PI*node_x[i]);
%node_new_y[i]=(1+(2*aaaa-1)*node_y[i])*sin(PI*node_x[i]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
newmesh.tris=mesh.tris;


% Solve linear system
A_I = sparse(A_I);
R = chol(A_I);
A_B = sparse(A_B);
b = -A_B*[newmesh.coords(boundary,1) newmesh.coords(boundary,2)];
solution = R\(R'\(b));
newmesh.coords(interior,1) = solution(:,1);
newmesh.coords(interior,2) = solution(:,2);

% Check for inverted elements.
%inverted = evalquality2(newmesh,filebase,3,1);

% draw initial mesh

figure2=figure;

for i=1:ele
   xxx_1 = mesh.coords( mesh.tris(i,1), 1);
   yyy_1 = mesh.coords( mesh.tris(i,1), 2);
    
   xxx_2 = mesh.coords( mesh.tris(i,2), 1);
   yyy_2 = mesh.coords( mesh.tris(i,2), 2);
   
   xxx_3 = mesh.coords( mesh.tris(i,3), 1);
   yyy_3 = mesh.coords( mesh.tris(i,3), 2);
   
   ppp=[xxx_1 xxx_2 xxx_3; yyy_1 yyy_2 yyy_3];
   figure2=fill(ppp(1,:), ppp(2,:), 'w');
   hold on;
end

saveas(figure2, 'b_femwarp.jpg');

figure5=figure;

for i=1:ele
   new_xxx_1 = newmesh.coords( mesh.tris(i,1), 1);
   new_yyy_1 = newmesh.coords( mesh.tris(i,1), 2);
    
   new_xxx_2 = newmesh.coords( mesh.tris(i,2), 1);
   new_yyy_2 = newmesh.coords( mesh.tris(i,2), 2);
   
   new_xxx_3 = newmesh.coords( mesh.tris(i,3), 1);
   new_yyy_3 = newmesh.coords( mesh.tris(i,3), 2);
   
   new_ppp=[new_xxx_1 new_xxx_2 new_xxx_3; new_yyy_1 new_yyy_2 new_yyy_3];
   
   new_area_1 = [ new_xxx_2 - new_xxx_1 new_xxx_3 - new_xxx_1; new_yyy_2-new_yyy_1   new_yyy_3 - new_yyy_1];
   new_area_2(i) = new_area_1(1,1)*new_area_1(2,2) - new_area_1(1,2)*new_area_1(2,1);
   
   if (new_area_2(i) > 0)
   figure5=fill(new_ppp(1,:), new_ppp(2,:), 'w');
   new_index(i)=0;
   else
   figure5=fill(new_ppp(1,:), new_ppp(2,:), 'r');
   new_index(i)=1;
   end
   hold on;
end

saveas(figure5, 'a_femwarp.jpg');
sum(new_index)



