function [ y ] = forward_substitution( L,b )
% This program computes the y in Ly=b using forward substitution method.
[n,~] = size(L);
y = zeros(n,1);
y(1) = b(1)/L(1,1);
for i=2:n
    j = i-1;
    y(i) = (b(i)-L(i,j)*y(j))/L(i,i);
end
end

