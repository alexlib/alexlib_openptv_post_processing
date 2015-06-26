function write_rt_is(filename)



xx=x;
yy=y;
zz=z;
c1=c1;
c2=c2;
totalpx=total_pix;
np = length(total_pix);
pnr = 1:np;
xpx=xpix;
ypx=ypix;
grval=grv;

fid = fopen(filename,'w');
fprintf(fid,'%d\n',np);
fprintf(fid,'%4d %9.4f %9.4f %5d %5d %5d %5d %5d\n',[pnr',xx,yy,zz,c1,c2,totalpx,xpx,ypx,grval]');
fclose(fid);
