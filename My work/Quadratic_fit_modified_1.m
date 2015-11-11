function [Hess1 ]=Quadratic_fit_modified_1(p, e, t, sol)
% This function incorporates methods for under and over determined systems.
% p-points, e-edges, t-triangles
nnodes = length(p);
Hess1 = zeros(nnodes,1);
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
    neighbs = [nb1  nb2  nb3  nb4  nb5  nb6];
    neighbs = neighbs';
    neighbs1 = unique(neighbs);
    % s is the # of neighbours
    [s, z] = size(neighbs1);
    % Nodes- node vector to hold nodes including the current one to fit  the curve.
    Nodes = zeros(s,1);
    %  X-vector to hold x coordinates of nodes selected.
    %  X = zeros(s,1);
    %  Y-vector to hold y coordinates of all 6 nodes selected.
    %  Y = zeros(s,1);
    %  U-solution vector for the nodes selected.
    U = zeros(s,1);
    Nodes(1,1) = k;
    i = 2;
    % code for constructing node vector sx1
    for j = 1:s
        if( neighbs1(j)~=0 && i<=6)
            Nodes(i,1) = neighbs1(j);
            i=i+1;
        end
    end
    %   node vector is done.
    %   now lets extract X, Y values for nodes and construct the matrix system
    %   Ac=u
    A=zeros(s,6);
    A(:,6)=1;
    for i=1:s
        if(Nodes(i,1)~=0)
            n=Nodes(i,1);
            % X(i,1)=p(1,n);
            % Y(i,1)=p(2,n);
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
    % format long
    % there is some rank issue sort it out.
    if(s<6) % under determined system x = M'inv(MM')Y
        M = A*A';
        if M == M' % checking if M is symmetric
            [T,p] = chol(M); % trying to get the cholesky factorization.
            if p==0 % means M is positive definite, T is valid
                w = T\U;
                z = T'\w;
                c = A'*z; % solving for coefficient matrix c
            end
        end
    elseif(s>6) % over determined system x=inv(M'M)M'y
        rk = rank(A);
        if(rk<min(s,6))
%       [Q,R]=qr(A);
            R = QR_HOUSE(A);
            Q = Q_HOUSE(R);
            y = U\Q;
            c = y\R;
        end
        else
            c = U\A;
        end
        hess = [2*c(1) c(3); c(3) 2*c(2)];
        % Frobenius norm
        %Hess1(k) = norm(hess);
        Hess1(k) = norm(hess,Inf);
        % Hess1(k)=norm(hess);
        %******************************************************************5
        % lets begin the refinement part %
%        write_file(p,e,t,Hess,'Elliptic');
    end
    %end
    %A1 = inv(A);
        %c = A1*U;
        %lets try c=(A'A)inverse*A'*U
        %A1 = A'*A;
        %A2=inverse(A1);
        %c=A2*A'*U;
  
  
  

  