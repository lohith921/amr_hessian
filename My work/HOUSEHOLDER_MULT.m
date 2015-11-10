function [A] = HOUSEHOLDER_MULT(A,v);

% input:   A - matrix
%          v - Householder vector 
% output:  A - transformed matrix 
%              A = A - \(\frac{\text{2}}{\text{v'v}}(\text{vv'})\)A 

% Author     : Stefan Hüeber
% Date       : May, 8, 2003
% Institution: University of Stuttgart,
%              Institut for Applied Analysis and Numerical Mathematics,
%              High Performance Scientific Computing
% Version    : 1.0

beta = -2/(v'*v);
w = v'*A;          % w is a line vector
A = A + beta*v*w;
