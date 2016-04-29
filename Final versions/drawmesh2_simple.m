% This function is used to draw the mesh without any border highliting
function []=drawmesh2_simple(mesh)

n=length(mesh.tris);


for i=1:n

  tri1=mesh.tris(i,1);
  tri2=mesh.tris(i,2);
  tri3=mesh.tris(i,3);
  xnode1 = mesh.coords(tri1,1);
  xnode2 = mesh.coords(tri2,1);
  xnode3 = mesh.coords(tri3,1);

  ynode1=mesh.coords(tri1,2);
  ynode2=mesh.coords(tri2,2);
  ynode3=mesh.coords(tri3,2);
  hold on;
  plot([xnode1,xnode2],[ynode1,ynode2],'k');
  plot([xnode2,xnode3],[ynode2,ynode3],'k');
  plot([xnode1,xnode3],[ynode1,ynode3],'k');
  
  end;

axis('square');
%   b1 = mesh.border(tri1);
%   b2 = mesh.border(tri2);
%   b3 = mesh.border(tri3);
% fprintf('triangle %d, node numbers are %d %d %d\n',i,tri1,tri2,tri3);


%   if (b1 & b2)
%     c1=1;
%   else 
%     c1=0;
%   end;
% 
%   if (b2 & b3)
%     c2=1;
%   else
%     c2=0;
%   end;
% 
%   if (b1 & b3)
%     c3=1;
%   else
%     c3=0;
%   end;

%   % WILL WANT TO SET THIS TO 1 LATER.
%   w = 3;
% 
%   if (~c1)
%     %line([xnode1,xnode2],[ynode1,ynode2]);
%     plot([xnode1,xnode2],[ynode1,ynode2],'k');
%   else
%     plot([xnode1,xnode2],[ynode1,ynode2],'r','LineWidth',w);
% %     plot([xnode1,xnode2],[ynode1,ynode2],'k');
%   end;
% 
%   if (i==1)
%     hold on;
%   end;
%   if (~c2)
%     %line([xnode2,xnode3],[ynode2,ynode3]);
%     plot([xnode2,xnode3],[ynode2,ynode3],'k');
%   else
%     plot([xnode2,xnode3],[ynode2,ynode3],'r','LineWidth',w);
% %     plot([xnode2,xnode3],[ynode2,ynode3],'k');
%   end;
%   if (~c3)
%     %line([xnode1,xnode3],[ynode1,ynode3]);
%     plot([xnode1,xnode3],[ynode1,ynode3],'k');
%   else
%     plot([xnode1,xnode3],[ynode1,ynode3],'r','LineWidth',w);
% %     plot([xnode1,xnode3],[ynode1,ynode3],'k');
%   end;


