function mesh = readmesh2_modified(filebase)
% Filenames
type1 = 'node';
type2 = 'ele';
type3 = 'neigh';
filebase
% Read in nodes
file = sprintf('%s.%s',filebase,type1);
[num_nodes,attr] = textread(file,'%d %*d %d %*d',1);
nodes = zeros(num_nodes,2);
if (attr==0)
  [nodes(:,1), nodes(:,2) boundary] = textread(file,'%*d %f %f %d',...
     num_nodes,'headerlines',1);
elseif (attr==1)
  [nodes(:,1), nodes(:,2) boundary] = textread(file,'%*d %f %f %*f %d',...
     num_nodes,'headerlines',1);
end;

%Read in element connectivity
file = sprintf('%s.%s',filebase,type2);
[num_ele] = textread(file,'%d %*d %*d',1);
ele = zeros(num_ele,3);
[ele(:,1), ele(:,2), ele(:,3)] = textread(file,'%*d %d %d %d',...
   num_ele,'headerlines',1);
%%%%%%%%%%%%%
% comment it for suzanne
% for hydrocephalus..need to comment it
ele=ele+1;

% lets read neighbour file
file = sprintf('%s.%s',filebase,type3);
[num_neighs] = textread(file,'%d %*d %*d',1);
neighs = zeros(num_neighs,3);
[neighs(:,1), neighs(:,2), neighs(:,3)] = textread(file,'%*d %d %d %d',...
   num_neighs,'headerlines',1);

%%%%%%%%%%%%%
mesh = struct('coords',nodes,'tris',ele,'border',boundary,'neighbour',neighs);
