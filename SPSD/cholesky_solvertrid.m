function [ x ] = cholesky_solvertrid( A,b )
% This program solves the system Ax=b where A is a tridiagonal SPD, and b
% is a vector.
[n,~] = size(A(1,:));
% computing the cholesky factorization.
for k=1:n
    A(k,k) = sqrt(A(k,k));
    A(k+1:n,k) = A(k+1:n,k)/A(k,k);
    for j = k+1:n
        A(j:n,j) = A(j:n,j) - A(j:n,k)*A(j,k);
    end
end
L = A; % lower triangular matrix
Lt = A.'; % upper triangular matrix
% LLtx = b is the new system, lets refactor it to Ly = b & Ltx = y
% first we will solve for y using forward substitution and then for x using
% back ward substitution
[n,~] = size(L);
y = zeros(n,1);
y(1) = b(1)/L(1,1);
for i=2:n
    j = i-1;
    y(i) = (b(i)-L(i,j)*y(j))/L(i,i);
end

% back ward solving for x
x = zeros(n,1);
x(n) = y(n) / Lt(n,n);
for i = n-1:-1:1
    j = i+1;
    x(i) = (y(i)-L(i,j)*y(i))/L(i,i);
end


end

