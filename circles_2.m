angle=linspace(0,2*pi,45);
x=cos(angle);
y=sin(angle);
x1=x*0.5;
y1=y*0.5;
fileid=fopen('circ_2_4514.node','w');
fprintf(fileid,'%d %d %d %d\n',90,2,1,1);
for i=1:45
    fprintf(fileid,'%d %f %f %d %d\n',i,x(i),y(i),0,1);
end
j=46;
for i=1:45
    fprintf(fileid,'%d %f %f %d %d\n',j,x1(i),y1(i),0,0);
    j=j+1;
end
fprintf(fileid,'\0');
fclose(fileid);