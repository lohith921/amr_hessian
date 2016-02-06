function [ x ] = deflated_PCGM( A, b, M, P )
% This program computes the solution of the deflated system MPAx=MPb using the PCGM
% ( M as preconditioner, P is deflation vector). It returns vector x. Here A is
% non-singular
n = size(b); % dimension of system.
x_hat = zeros(n,1); % x0 =0
r = b; % r0 = b;
r_hat = P*r;
s = cholesky_solvertrid(M,r_hat); % z0 = inv(M)r0
p = s; % p0 = z0
while r_hat ~= 0
    alpha = (s.')*r_hat/(p.')*P*A*p;
    x1_hat = x_hat + alpha*p; % Xj+1 = Xj + alphaPj
    r1_hat = r_hat - alpha*P*A*p; % rj+1 = rj + alphaPj
    s1 = cholesky_solvertrid(M,r1_hat);
    beta = (s1.')*r1_hat/(s.')*r_hat; % beta = (rj+1,zj+1)/(rj,zj)
    p = s1 + beta*p;
    r_hat = r1_hat;
    s = s1;
end
x = x1;
end



