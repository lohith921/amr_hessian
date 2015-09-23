x=linspace(0,2*pi,360);
y=cos(x);
fileid=fopen('circ.poly','w');
fprintf(fileid,'%d %d %d %d\n',360,2,1,0);
for i=1:360
    fprintf(fileid,'%d %f %f\n',i,x(i),y(i));
end