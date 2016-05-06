function [ x ] = back_substitution( Lt,y )
% This program implements the back substitution algorithm for cholesky
% factorization method . L'x = y. Lt is the upper triangular matrix and y
% is vector.
[n,~] = size(Lt);
x = zeros(n,1);
% x(n) = y(n) / Lt(n,n);
% for i = n-1:-1:1
%     j = i+1;
%     x(i) = (y(i)-Lt(i,j)*y(i))/Lt(i,i);
% end
for i = n:-1:1
   x(i) = ( y(i) - Lt(i, :)*x )/Lt(i, i);
end
end

