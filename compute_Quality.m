% The function computes the quality of mesh.
function [qual] = compute_Quality(mesh)
nele = length(mesh.tris);
quality = zeros(nele,1);
c = ones(3,1);
for i = 1:nele
    a = zeros(3,3);
    P1 = mesh.coords(mesh.tris(i,1),:);
    P2 = mesh.coords(mesh.tris(i,2),:);
    P3 = mesh.coords(mesh.tris(i,3),:);
    a = [P1 1; P2 1; P3 1];
    A = det(a);
    h1 = norm(P1-P2,2);
    h2 = norm(P2-P3,2);
    h3 = norm(P3-P1,2);
    quality(i) = 4*A*sqrt(3)/(h1^2 + h2^2 + h3^3);
end
qual = quality;
% end