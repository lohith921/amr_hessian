function [ ] = write_file( mesh, filebase )
% write_file writes the mesh information into .node .ele files etc
%   Detailed explanation goes here
    type1 = 'r.ele';
    type2 = 'r.node';
    file1 = sprintf('%s.%s', filebase, type1);
    file2 = sprintf('%s.%s', filebase, type2);
    nele = length(mesh.tris);
    nnodes = length(mesh.coords);
    % beginning writing to ele file
    fele = fopen(file1,'wt');
    fprintf(fele, '%d\t %d\t %d\t %d\n', nele, 3, 0, 0);
    for i = 1:nele
        fprintf(fele, ' %d\t %d\t %d\t %d\n', i, mesh.tris(i,1), mesh.tris(i,2), mesh.tris(i,3));
    end
    fclose(fele);
    % end of writing to ele file
    
    % beginning writing to node fil
    fnode = fopen(file2, 'wt');
    fprintf(fnode, '%d\t %d\t %d\t %d\n', nnodes, 2, 0, 0);
    for i = 1:nnodes
        fprintf(fnode, ' %d\t %f\t %f\t %d\t %d\n', i, mesh.coords(i,1), mesh.coords(i,2), 0, 0);
    end
    fclose(fnode);
    % end of writing to node file
end

