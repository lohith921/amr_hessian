function [ x ] = preconditioned_CGM(A, b, L)
% This program computes the solution of the system Ax=b using the CGM
% method with M as preconditioner. It returns vector x. Here A is
% non-singular. M is constructed from my L as M = L(L.')
M = L*transpose(L);
[n,~] = size(b); % dimension of system.
x = zeros(n,1); % x0 =0
x1 = x;
r = b; % r0 = b;
%z = cholesky_solvertrid(M,rk); % z0 = inv(M)r0
Z = cholesky_factorization(M);
y = forward_substitution(Z,r);
z = back_substitution(transpose(Z),y);
p = z; % p0 = z0
Z1 = Z;
error = 0;
while (error <= (10^-3) )
    alpha = (transpose(r)*z)/(transpose(p)*(A*p));
    x1 = x + alpha*p; % Xj+1 = Xj + alphaPj
    r1 = r - alpha*A*p; % rj+1 = rj + alphaPj
    %Z1 = cholesky_factorization(M);
    y1 = forward_substitution(Z1,r1);
    z1 = back_substitution(transpose(Z1),y1);
%z1 = cholesky_solvertrid(M,r1);
    beta = (transpose(r1)*z1)/(transpose(r)*z); % beta = (rj+1,zj+1)/(rj,zj)
    p1 = z1 + beta*p;
    error = norm(r1-r,2);
    error
    r = r1;
    z = z1;
    x = x1;
    p = p1;
    %norm(r,2)
end
x = x1;
end



% rk = r0;
% % rk1 = rk;
% rk2 = zeros(n,1);
% zk2 = zeros(n,1);
% while(rk ~= 0)
%     rk1 = rk;
%     % solve for Mzk = rk
%     zk = cholesky_solvertrid(M,rk);
%     k = k+1;
%     if k==1
%         pk = zk;
%         rk2 = rk;
%         zk1 = zk; % for next round of iteration.
%         xk1 = x0;
%     else
%         betak = rk1.'*zk1/rk.'*zk1;
%         pk = zk1 + betak*pk1;
%     end
%     alphak = rk1.'*zk1/pk.'*A*pk;
%     xk = xk1 + alphak*pk;
% %     rk = rk1 - alphak*A*pk;
% %     zk1 = zk;
% % end
% r = zeros(n,1) % represents rk in algorithm
% r = b; % before the start of algorithm
% r1 = zeros(n,1); % represents rk-1 in algorithm
% r2 = zeros(n,1); % represents rk-2 in algorithm
% z1 = zeros(n,1); % represents zk-1 in algorithm
% z2 = zeros(n,1); % represents zk-2 in algorithm
% z = zeros(n,1); % represents zk in algorithm
% p = zeros(n,1); % represents pk in algorithm
% x = zeros(n,1); % represents xk in the algorithm
% k = 0;
% while(r ~= 0)
%     % solve for Mzk = rk
%     z = cholesky_solvertrid(M,r);
%     k = k+1;
%     if k==1
%         p = z;
%          r1 = r; % for next round.
%          z1 = z; % for next round of iteration.
% %         xk1 = x0;
%     else
%         betak = r1.'*z1/r2.'*z2;
%         pk = zk1 + betak*pk1;
%     end
%     alphak = rk1.'*zk1/pk.'*A*pk;
%     xk = xk1 + alphak*pk;
%     rk = rk1 - alphak*A*pk;
%     zk1 = zk;
% end
% x = xk;
% end

