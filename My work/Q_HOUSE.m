function [Q] = Q_HOUSE(A);

% input:   A - (n x m) matrix, QR - decomposition containing the 
%              necessary information of the Householder vectors v
%              in the lower triangle and R in the upper triangle
% output   Q - orthonormal matrix Q = Q_{n-1}* .... *Q_1
%              with Q_k = I-2/(v'*v)*(v*v')

% Author     : Stefan Hüeber
% Date       : May, 12, 2003
% Institution: University of Stuttgart,
%              Institut for Applied Analysis and Numerical Mathematics,
%              High Performance Scientific Computing
% Version    : 1.0

n = size(A,1);
Q = eye(n);
for k = 1:n-1
	v = ones(n+1-k,1);
	v(2:n+1-k) = A(k+1:n,k);
	Qk = eye(n);
	Qk(k:n,k:n) = eye(n+1-k) - (2/(v'*v))*(v*v');
	Q = Qk*Q;
end