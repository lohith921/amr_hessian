
angle=linspace(0,2*pi,45);
x=cos(angle);
y=sin(angle);
fileid=fopen('./Files/circ1_45.poly','w');
fprintf(fileid,'%d %d %d %d\n',45,2,1,0);
for i=1:45
    fprintf(fileid,'%d %f %f %f\n',i,x(i),y(i),0);
end
fprintf(fileid,'%d %d\n',45,0);
fprintf(fileid,'%d %d %d\n',1,45,1);
for i=2:45
    fprintf(fileid,'%d %d %d\n',i,i-1,i);
end
fprintf(fileid,'\0');
fclose(fileid);