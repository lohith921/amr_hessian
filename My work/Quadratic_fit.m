% mesh = struct('coords',nodes,'tris',ele,'border',boundary);
% Read-in initial mesh
mesh=readmesh2_modified('trial.1');
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
%A = zeros(i+b,i+b);

% Global ordering for nodes.
glbal = [interior' boundary'];

% Compute local stiffness matrices and then assemble then into 
% global stiffness matrix.

%for k=1:length(mesh.coords)
  % Get this node.
  mynode = mesh.coords(1,:);

  % Get global ordering for this tri.
  %g1 = find(glbal == mytris(1));
  %g2 = find(glbal == mytris(2));
  %g3 = find(glbal == mytris(3));

  % Get x, y coordinates for vertices.
  %x = mesh.coords(mytris,1);
  %y = mesh.coords(mytris,2);
 
  x=mynode(1);
  y=mynode(2);
  neighbs=mesh.neighbour(1,:);
  %neighbs
  n1=neighbs(1);
  n2=neighbs(2);
  n3=neighbs(3);
  %[x1 y1] = mesh.coords(n1,:);
  %[x2 y2] = mesh.coords(n2,:);
  if(n1 ~= -1)
      X=mesh.coords(n1,:);
      x1=X(1); y1=X(2);
  else
      x1=0; y1=0;
  end
  if(n2 ~= -1)
      Y=mesh.coords(n2,:);
      x2=Y(1); y2=Y(2);
  else
      x2=0; y2=0;
  end
  if(n3 ~= -1)
      Z = mesh.coords(n3,:);
      x3=Z(1); y3=Z(2);
  else
      x3=0; y3=0;
  end
  A=zeros(4,4);
  %A(1,1)=x^2;
  A(:,1)=[x^2 x1^2 x2^2 x3^2];
  A(:,2)=[y^2 y1^2 y2^2 y3^2];
  A(:,3)=[x*y x1*y1 x2*y2 x3*y3];
  A(:,4)=[1 1 1 1];
  A
%end
  
  
  

  