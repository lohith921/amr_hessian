function mesh = readmesh2(filebase)
% Filenames

type1 = 'node';
type2 = 'ele';
% Read in nodes
file = sprintf('%s.%s',filebase,type1);
[num_nodes,attr] = textread(file,'%d %*d %d %*d',1);
nodes = zeros(num_nodes,2);
if (attr == 0)
  [nodes(:,1), nodes(:,2) boundary] = textread(file,'%*d %f %f %d',...
     num_nodes,'headerlines',1);
elseif (attr == 1)
  [nodes(:,1), nodes(:,2) boundary] = textread(file,'%*d %f %f %*d %d',...
     num_nodes,'headerlines',1);
end;

% Read in element connectivity
file = sprintf('%s.%s',filebase,type2);
[num_ele] = textread(file,'%d %*d %*d',1);
ele = zeros(num_ele,3);
[ele(:,1), ele(:,2), ele(:,3)] = textread(file,'%*d %d %d %d',...
   num_ele,'headerlines',1);

%%%%%%%%%%%%%
% comment it for suzanne
% for hydrocephalus..need to comment it
%ele=ele+1;

%%%%%%%%%%%%%
mesh = struct('coords',nodes,'tris',ele,'border',boundary);
