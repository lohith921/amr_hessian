function [ A ] = incomplete_cholesky(A)
% This program computes the incomplete cholesky factorization of input
% matrix and returns the lower triangular L matrix
n = size(A(1,:));
% B = A; 
for k = 1:n
    A(k,k) = sqrt(A(k,k));
    for i = k+1:n
        if A(i,k) ~= 0
            A(i,k) = A(i,k)/A(k,k);
        end
    end
    for j=k+1:n
        for i = j:n
            if A(i,j) ~= 0
                A(i,j) = A(i,j) - A(i,k)*A(j,k);
            end
        end
    end
end

