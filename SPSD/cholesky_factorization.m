function [ A ] = cholesky_factorization(A)
% This program computes the cholesky factorization 
% A = GG' of given input matrix
% and returns the lower triangular matrix. 
n = size(A(1,:));
for k=1:n
    A(k,k) = sqrt(A(k,k));
    A(k+1:n,k) = A(k+1:n,k)/A(k,k);
    for j = k+1:n
        A(j:n,j) = A(j:n,j) - A(j:n,k)*A(j,k);
    end
end
end

