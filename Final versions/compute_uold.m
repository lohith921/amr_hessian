% This function computes the solution at the gauss point using linear
% interpolation. Inputs are old mesh, solution vector, gauss point(vector)
function [val] = compute_uold(o_mesh, sol, Xg)
    val = 0;
    for k = 1:length(o_mesh.tris)
        tri_nodes = o_mesh.tris(k,:);
        X1 = o_mesh.coords(tri_nodes(1),:);
        X2 = o_mesh.coords(tri_nodes(2),:);
        X3 = o_mesh.coords(tri_nodes(3),:);
        left = [X1 1; X2 1; X3 1];
        right = [Xg ; 1];
        weights = (left')\right;
        if( sum((weights >= 0)) == 3)
        % if (bool == 1) % This point lies on or inside current triangle.
                 A = abs(det([(X2' - X1') (X3' - X1')])/2);			% area of the whole triangle
                A1 = abs(det([(Xg - X2') (Xg - X3')])/2);			%  area of triangle formed by Xg, X2, X3
                A2 = abs(det([(Xg - X1') (Xg - X3')])/2);			%  area of triangle formed by Xg, X1, X3
                A3 = abs(det([(Xg - X1') (Xg - X2')])/2);			%  area of triangle formed by Xg, X2, X1
                %if ((A1 + A2 + A3) == A)
                    %phi1 = A1/A;      phi2 = A2/A;      phi3 = A3/A;
                phi = [A1 A2 A3]/A;								% row vector of interpolation function values
                val = phi*sol(tri_nodes');
                break;
        end
   end
end
%% The following function is used to check if a point lies inside or on the triangle.
% function b = checkinside(A, B, C, P)
%     %Compute vectors        
%     v0 = C - A;
%     v1 = B - A;
%     v2 = P - A;
% 
%     % Compute dot products
%     dot00 = dot(v0, v0);
%     dot01 = dot(v0, v1);
%     dot02 = dot(v0, v2);
%     dot11 = dot(v1, v1);
%     dot12 = dot(v1, v2);
% 
%     % Compute barycentric coordinates
%     invDenom = 1 / (dot00 * dot11 - dot01 * dot01);
%     u = (dot11 * dot02 - dot01 * dot12) * invDenom;
%     v = (dot00 * dot12 - dot01 * dot02) * invDenom;
% 
%     % Check if point is in triangle
%     if (u >= 0) && (v >= 0) && (u + v <= 1)
%         b = 1;
%     else
%         b = 0;
%     end 
% end
%         xmin = min([X1(1) X2(1) X3(1)]);
%         xmax = max([X1(1) X2(1) X3(1)]);
%         ymin = min([X1(2) X2(2) X3(2)]);
%         ymax = max([X1(2) X2(2) X3(2)]);
       % check = ((Xg(1) < xmin) || (Xg(1) > xmax)) && ((Xg(2) < ymin) || (Xg(2) > ymax));
        %if ~check
           % bool = checkinside((X1'), (X2'), (X3'), Xg);

