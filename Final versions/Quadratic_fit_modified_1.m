function metric_t = Quadratic_fit_modified_1(p, t, solution)
% This function incorporates methods for under and over determined systems.
% p-points, e-edges, t-triangles
n_nodes = length(p);
% 3D array to hold metric tensors, 2x2 matrix for each node.
metric_t = zeros(2,2,n_nodes);
for k = 1:n_nodes
    % checking which triangles have node k as first vertex
    t1 = find(t(1,:) == k);
    % triangle t1 has node k as first vertex, lets look for other vertices
    nb1 = t(2,t1);
    nb2 = t(3,t1);
    % checking which triangles have node k as second vertex
    t2 = find(t(2,:) == k);
    % t2 has node k as 2nd vertex, lets look for other vertices
    nb3 = t(1,t2);
    nb4 = t(3,t2);
    % checking which triangles have node k as third vertex
    t3 = find(t(3,:) == k);
    % t3 has node k as 23rd vertex, lets look for other vertices
    nb5 = t(1,t3);
    nb6 = t(3,t3);
    % neighbs- neighbour vector
    neighbs = [nb1 nb2 nb3 nb4 nb5 nb6];
    neighbs = neighbs';
    neighbs1 = unique(neighbs);
    % s is the # of neighbours
    [s, ~] = size(neighbs1);
    % Nodes- node vector to hold nodes including the current one to fit  the curve.
    Nodes = zeros(s,1);
       %  b-solution vector for the nodes selected. choosing for easy notation
    %  purpose while using LU
    rsv = zeros(s,1);
    Nodes(1,1) = k;
    i = 2;
    % code for constructing node vector sx1
    for j = 1:s
        if( (neighbs1(j) ~= 0) && (i <= 6))
            Nodes(i,1) = neighbs1(j);
            i=i+1;
        end
    end
    %   node vector is done.
    %   now lets extract X, Y values for nodes and construct the matrix system
    %   Ac=u
    A = zeros(s,6);
    c = zeros(6,1);
    A(:,6) = 1;
    for i = 1:s
        if (Nodes(i,1) ~= 0)
            n = Nodes(i,1);
            x = p(n,1);
            y = p(n,2);
            rsv(i,1) = solution(n,1);
            A(i,1) = x*x;
            A(i,2) = y*y;
            A(i,3) = x*y;
            A(i,4) = x;
            A(i,5) = y;
        else
            x = 0;
            y = 0;
            rsv(i,1)=0;
            A(i,1) = x*x;
            A(i,2) = y*y;
            A(i,3) = x*y;
            A(i,4) = x;
            A(i,5) = y;
        end
    end
    % format long
    % there is some rank issue sort it out.
    if(s < 6) % under determined system x = M'inv(MM')Y
        M = A*(A.');
        if M == (M.') % checking if M is symmetric
            [T,p1] = chol(M); % trying to get the cholesky factorization.
            if p1==0 % means M is positive definite, T is valid
                w = double(T\rsv);
                z = double(T.'\w);
                c = double(A.'*z); % solving for coefficient matrix c
            end
        end
    elseif(s >= 6) % over determined system x=inv(M'M)M'y
        rk = rank(A);
        if(rk >= min(s,6)) % A has full rank.
            % [Q,R]=qr(A); QR factorization.
            R = QR_HOUSE(A);
            Q = Q_HOUSE(R);
            u = (Q.')*rsv;
            % solving for c by backward substitution
            c(6) = rsv(6)/R(6,6);
            for i = n:1
                for j= i+1:n
                    su = su + R(i,j)*c(i);
                end
                c(i) = (u(i)-su)/R(i,i);
            end
        end
    else
        c = rsv\A;
    end
    hess = [2*c(1) c(3); c(3) 2*c(2)];
    % Eigen value decomposition of hess
    [V,D] = eig(hess);
    D = abs(D);
    % might have to change it to VDinv(V)
    % met = V*D*inv(V);
    met = V*D*(V.');
    metric_t(:,:,k) = met;
end
end
%         else % A doesnot have full rank.
%             [ Q R E]=qr(A);
%         end
%     Hess1(k) = norm(hess,Inf);
%     Hess1(k)=norm(hess);   

% lets consider swapping edges
% ntris = length(t);
% nedge = length(e);
% code for swapping edges
%% 
% for i = 1:ntris-1
%     x = t(1,i);
%     y = t(2,i);
%     z = t(3,i);
%     X = p(:,x);
%     Y = p(:,y);
%     Z = p(:,z);
%     w = 0;
%     Nbs = [x y z];
%     j = 1;
%     while(j~=ntris)
%         if j==i
%             j=j+1;
%             continue;
%         else
% we will find nodes which form a quadrilateral with each edge. 
%         if( any(Nbs(:)==t(3,j)) && any(Nbs(:)==t(1,j)))
%            flag = 1;
%            w = t(2,j);
%            break;
%         elseif( any(Nbs(:)==t(2,j)) && any(Nbs(:)==t(3,j)))
%             flag = 2;
%             w = t(1,j);
%             break;
%         elseif( any(Nbs(:)==t(1,j)) && any(Nbs(:)==t(2,j)))
%             flag = 3;
%             w = t(3,j);
%             break;
%         else
%             j = j+1;
%         end
%         end
%     end
%     W = p(:,w);
%     T = find((t(1,:)==x) & (t(3,:)==y));
%     Mavg = (metric_t(:,:,x) + metric_t(:,:,y) + metric_t(:,:,z) + metric_t(:,:,w))/4;
%     L1 = Y-X;
%     L2 = Z-X;
%     L3 = Z-W;
%     L4 = Y-W;
%     t1 = L1(1,:)*L2(2,:)-L1(2,:)*L2(1,:);
%     t2 = (Z-W).' * Mavg * (Y-W);
%     t3 = (Y-X).' * Mavg * (Z-X);
%     t4 = L3(1,:)*L4(2,:)-L3(2,:)*L4(1,:);
%     if (t1*t2 + t3*t4)<0
%         ed = 1;
%         while(ed ~= nedge)
%             k = ed;
%          if (e(1,k)==y && e(2,k)==z) 
%              e(1,k) = w;
%              e(2,k) = x;
%              break;
%          elseif (e(2,k)==y && e(1,k)==z)
%              e(2,k) = w;
%              e(1,k) = x;
%              break;
%          else
%              ed = ed + 1;
%          end
%         end
%         t(3,i) = w;
%         t(1,j) = w;
%         t(2,j) = y;
%         t(3,j) = z;
%             t(2,k) = y;            
%             t(3,i) = w;
%             t(1,i) = z;
%             t(1,j) = w;
%             t(2,j) = z;
%             t(3,j) = x;        
%         
%             t(1,i) = w;
%             
%     end
% end
%%
% swapping edges ends here.
%  % lets begin the refinement part, working on edges %
%     nedges = length(e);
%     for i=1:nedges
%         n1 = e(1,i);
%         n2 = e(2,i);
%         v1 = p(:,n1); % vertex 1
%         v2 = p(:,n2); % vertex 2
%         M1 = metric_t(:,:,n1);
%         M2 = metric_t(:,:,n2);
%         Mavg = (M1+M2)/2;
%         v = v1 - v2;
%         vt = v.';
%         temp = Mavg*v;
%         L = sqrt((vt)*temp);       
%         fprintf('\nlength of edge %d and %d is',n1,n2);
%         L 
%     end

 %  X-vector to hold x coordinates of nodes selected.
    %  X = zeros(s,1);
    %  Y-vector to hold y coordinates of all 6 nodes selected.
    %  Y = zeros(s,1);



%  l1 = sqrt(v1'*metric_t(:,:,n1)*v1);
%         l2 = sqrt(v2'*metric_t(:,:,n2)*v2);
%         if (l1+l2)~=0
%           L = 2*(l1^2+l1*l2+l2^2)/(l1+l2)*3;
%         else 
%           L = 0;
%         end




