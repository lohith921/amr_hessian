% This script is written and read by pdetool and should NOT be edited.
% There are two recommended alternatives:
 % 1) Export the required variables from pdetool and create a MATLAB script
 %    to perform operations on these.
 % 2) Define the problem completely using a MATLAB script. See
 %    http://www.mathworks.com/help/pde/examples/index.html for examples
 %    of this approach.
function pdemodel
[pde_fig,ax]=pdeinit;
pdetool('appl_cb',1);
set(ax,'DataAspectRatio',[1 1.2 1]);
set(ax,'PlotBoxAspectRatio',[2.5 1.6666666666666667 1]);
set(ax,'XLim',[-2 3]);
set(ax,'YLim',[-2 2]);
set(ax,'XTickMode','auto');
set(ax,'YTickMode','auto');

% Geometry description:
pdeellip(1.309806835066865,-0.035661218424962726,1.4933135215453199,1.1054977711738494,...
0,'E1');
pdeellip(1.2875185735512638,0.044576523031203852,0.38632986627043087,0.72213967310549787,...
0,'E2');
set(findobj(get(pde_fig,'Children'),'Tag','PDEEval'),'String','E1-E2')

% Boundary conditions:
pdetool('changemode',0)
pdesetbd(8,...
'dir',...
1,...
'1',...
'0')
pdesetbd(7,...
'dir',...
1,...
'1',...
'0')
pdesetbd(6,...
'dir',...
1,...
'1',...
'0')
pdesetbd(5,...
'dir',...
1,...
'1',...
'0')
pdesetbd(4,...
'dir',...
1,...
'1',...
'0')
pdesetbd(3,...
'dir',...
1,...
'1',...
'0')
pdesetbd(2,...
'dir',...
1,...
'1',...
'0')
pdesetbd(1,...
'dir',...
1,...
'1',...
'0')

% Mesh generation:
setappdata(pde_fig,'Hgrad',1.3);
setappdata(pde_fig,'refinemethod','regular');
setappdata(pde_fig,'jiggle',char('on','mean',''));
setappdata(pde_fig,'MesherVersion','preR2013a');
pdetool('initmesh')
pdetool('jiggle')
pdetool('jiggle')
pdetool('jiggle')

% PDE coefficients:
pdeseteq(1,...
'1.0',...
'0.0',...
'10',...
'1.0',...
'0:10',...
'0.0',...
'0.0',...
'[0 100]')
setappdata(pde_fig,'currparam',...
['1.0';...
'0.0';...
'10 ';...
'1.0'])

% Solve parameters:
setappdata(pde_fig,'solveparam',...
char('0','1000','10','pdeadworst',...
'0.5','longest','0','1E-4','','fixed','Inf'))

% Plotflags and user data strings:
setappdata(pde_fig,'plotflags',[1 1 1 1 1 1 1 1 0 1 0 1 1 0 0 0 0 1]);
setappdata(pde_fig,'colstring','');
setappdata(pde_fig,'arrowstring','');
setappdata(pde_fig,'deformstring','');
setappdata(pde_fig,'heightstring','');

% Solve PDE:
pdetool('solve')
