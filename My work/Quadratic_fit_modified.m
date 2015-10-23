function [Hess1 ]=Quadratic_fit_modified(p, e, t, sol)
% p-points, e-edges, t-triangles
nnodes = length(p);
Hess1 = zeros(nnodes,1);
Hess2 = zeros(nnodes,1);
for k=1:nnodes
 % checking which triangles have node k as first vertex
  t1 = find(t(1,:) == k);
 % triangle t1 has node k as first vertex, lets look for other vertices
  nb1 = find(t(2,t1));
  nb2 = find(t(3,t1));
  
 % checking which triangles have node k as second vertex
  t2 = find(t(2,:) == k);
 % t2 has node k as 2nd vertex, lets look for other vertices
  nb3 = find(t(1,t2));
  nb4 = find(t(3,t2));
  
 % checking which triangles have node k as third vertex
 t3 = find(t(3,:) == k);
 % t3 has node k as 23rd vertex, lets look for other vertices
  nb5 = find(t(1,t3));
  nb6 = find(t(3,t3));
 % neighbs- neighbour vector 
  neighbs = [nb1 nb2 nb3 nb4 nb5 nb6];
  neighbs1 = unique(neighbs);
 % Nodes- node vector to hold upto 6 nodes including the current one to fit
 % the curve.
  Nodes=zeros(6,1);
 % X-vector to hold x coordinates of all 6 nodes selected.
  X=zeros(6,1);
 % Y-vector to hold y coordinates of all 6 nodes selected.
  Y=zeros(6,1);
 % U-solution vector for the nodes selected.
  U=zeros(6,1);
  Nodes(1,1)=k;
  i=2;
  % code for constructing node vector 6x1
  for j=1:length(neighbs1)
      if( neighbs1(j)~=0 && i<=6)
          Nodes(i,1)=neighbs1(j);
          i=i+1;
      end
  end
  if i<=6
      for j=i:6
          Nodes(j,1)=0;
      end
  end
  % node vector is done.
  
  % now lets extract X, Y values for nodes and construct the matrix system
  % Ac=u
  A=zeros(6,6);
  A(:,6)=1;
  for i=1:6
      if(Nodes(i,1)~=0)
          n=Nodes(i,1);
         % X(i,1)=p(1,n);
          %Y(i,1)=p(2,n);
          x = p(1,n);
          y = p(2,n);
          U(i,1)=sol(n,1);
          A(i,1) = x*x;
          A(i,2) = y*y;
          A(i,3) = x*y;
          A(i,4) = x;
          A(i,5) = y;
      else
          x = 0;
          y = 0;
          U(i,1)=0;
          A(i,1) = x*x;
          A(i,2) = y*y;
          A(i,3) = x*y;
          A(i,4) = x;
          A(i,5) = y;
      end
  end
  % there is some rank issue sort it out.
  c = U\A;
  hess = [2*c(1) c(3); c(3) 2*c(2)];
  hess
  % Hess1(k)=norm(hess);
  %Hess2(k)=norm(hess,2);
end
%end
  
  
  

  