function write_xuap(xuap,frame)
% WRITE_XUAP(FRAME)
% writes an element of the XUAP structure (xuap.x, xuap.y, z,u,v,w,ax,ay,az,id)
% into a XUAP.FRAME
% 

np = length(xuap.x);
x = xuap.x;
y = xuap.y;
z = xuap.z;
u = xuap.u;
v = xuap.v;
w = xuap.w;
ax = xuap.ax;
ay = xuap.ay;
az = xuap.az;
id = xuap.id;


filename = sprintf('xuap.%d',frame);
fid = fopen(filename,'w');
fprintf(fid,'%d\n',np);
fprintf(fid,'%9.4f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f %5d\n',[x,y,z,u,v,w,ax,ay,az,id]');
fclose(fid);
