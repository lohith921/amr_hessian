function [Hess1]=Quadratic_fit_centers(p, e, t, sol)
% This function fits the polynomial at triangle centers and computes the
% hessian there.Under and over determined systems are also solved.
% p-points, e-edges, t-triangles
ntris = length(t);
Hess1 = zeros(ntris,1);
X = zeros(3,2);
x = zeros(3,2);
b = zeros(3,1);
ctrd = zeros(2,1);
for k=1:ntris
    %     p1 = p(t(1,k)); % p1 holds x,y of 1st vertex
    %     p2 = p(t(2,k)); % p2 holds x,y of 2nd vertex
    %     p3 = p(t(3,k)); % p3 holds x,y of 3rd vertex
    for i = 1:3
        pt = p(t(i,k));
        X(i,:) = pt;
        b(i) = sol(t(i,k));
    end
    % computing the midpoints
    %     for i = 1:3
    %         x(i,:) = (X(i,:) + X(mod(i+1,3),:))/2;
    %     end
    % computing the centroid of the triangle
    ctrd(1) = sum(X(:,1))/3;
    ctrd(2) = sum(X(:,2))/3;
    % a1x^2 + a2y^2 + a3x + a4y + a5xy + a6 is our polynomial
    % Ac = b is my matrix system.
    A = zeros(3,6);
    c = zeros(6,1);
    A(:,6)=1;
    for i=1:3
        A(i,1) = X(i,1)*X(i,1); % xi^2
        A(i,2) = X(i,2)*X(i,2); % yi62
        A(i,3) = X(i,1); % xi
        A(i,4) = X(i,2); % yi
        A(i,5) = X(i,1)*X(i,2); % xi*yi
    end
    % Ac = b is an under determined system.
%     M = A*A';
At = A';
M = At*A;
    b = At*b;
    if M == M' % checking if M is symmetric
       [T,p1] = chol(M); % trying to get the cholesky factorization.
%        [L,U] = lu(M);
        if p1 == 0 % means M is positive definite, T is valid
            w = L\b;
            z = U\w;
            c = z;
            %c = A'*z; % solving for coefficient matrix c
        end
    end
    hess = [2*c(1) c(5); c(5) 2*c(2)];
    Hess1(k) = norm(hess,Inf);
end
%     % code for constructing node vector sx1
%     for j = 1:s
%         if( neighbs1(j)~=0 && i<=6)
%             Nodes(i,1) = neighbs1(j);
%             i=i+1;
%         end
%     end
%     %   node vector is done.
%     %   now lets extract X, Y values for nodes and construct the matrix system
%     %   Ac=u
% if(Nodes(i,1)~=0)
%             n = Nodes(i,1);
%             x = p(1,n);
%             y = p(2,n);
%             b(i,1)=sol(n,1);
%             A(i,1) = x*x;
%             A(i,2) = y*y;
%             A(i,3) = x*y;
%             A(i,4) = x;
%             A(i,5) = y;
%         else
%             x = 0;
%             y = 0;
%             b(i,1)=0;
%             A(i,1) = x*x;
%             A(i,2) = y*y;
%             A(i,3) = x*y;
%             A(i,4) = x;
%             A(i,5) = y;
%         end
% elseif(s>6) % over determined system x=inv(M'M)M'y
%         rk = rank(A);
%         if(rk >= min(s,6)) % A has full rank.
%             % [Q,R]=qr(A); QR factorization.
%             R = QR_HOUSE(A);
%             Q = Q_HOUSE(R);
%             u = Q'*b;
%             % solving for c by backward substitution
%             c(6) = b(6)/R(6,6);
%             for i = n:1
%                 for j= i+1:n
%                     su = su + R(i,j)*c(i);
%                 end
%                 c(i) = (b(i)-su)/R(i,i);
%             end
%         end
% Hess1(k)=norm(hess);
    %******************************************************************5
    % lets begin the refinement part %
    %        write_file(p,e,t,Hess,'Elliptic');

