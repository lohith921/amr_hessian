% This function is used to find the neighbor triangle which shares edge
% with current triangle. p1, p2 are vertices of the edge. tri is the array
% of triangles. The triangle number is returned in n and the index of third
% vertex is returned in ot
function [n, ot] = find_neighbor_with(p1, p2, tri)
    n = find((tri(:,1) == p1) & (tri(:,2) == p2));
    ot = 3;
    if isempty(n)
        n = find((tri(:,2) == p1) & (tri(:,3) == p2));
        ot = 1;
        if isempty(n)
            n = find((tri(:,3) == p1) & (tri(:,1) == p2));
            ot = 2;
        end
    end
end