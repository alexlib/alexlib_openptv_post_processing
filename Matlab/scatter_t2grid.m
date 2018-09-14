
%This function is used to create a Eulieran data of the 3D position and 
%velocities %of the particles from the scatter data obtained by traj2scatter t.m 
%Its required inputs are the scatter data and the grid paramters  for  the  positions .  
%The  optional  variable ,  q,  can  be  
%used  to  specify  the  number  of  data  files  used .  
%This function  will  output  
%the  3D position  and  velocity  grid  (xg ,  yg ,  zg ,  ug ,  vg ,  
%wg)  and  the  average  
%velocities  (u avg ,  v avg ,  w avg) .  

%Input:[scatter t,x i,x interval,x f,y i,y interval,y f ,z i,z interval ,z f,number of files(optional)] 
%Output:[xg,yg,zg,ug,vg,wg,u avg,v avg,w avg] 
%Jin Tae Kim 2014 University of Illinois at Urbana Champaign 

function [varargout] = scatter_t2grid( scatter_t ) 
 
x1 = min(cat(1,scatter_t.xf));
x2 = max(cat(1,scatter_t.xf));

y1 = min(cat(1,scatter_t.yf));
y2 = max(cat(1,scatter_t.yf));


z1 = min(cat(1,scatter_t.zf));
z2 = max(cat(1,scatter_t.zf));

[nx,ny,nz] = deal(5);

[xg,yg,zg] = ndgrid(linspace(x1,x2,nx),linspace(y1,y2,ny),linspace(z1,z2,nz));

[ug,vg,wg] = deal(0*xg);

for i=1:length(scatter_t) % uses every data file in the specified directory and tracks how fast the program is running
    x = scatter_t(i).x; %%to avoid the interpolation around zero
    y = scatter_t(i).y;
    z = scatter_t(i).z;
    u = scatter_t(i).u;
    v = scatter_t(i).v;
    w = scatter_t(i).w;
    
    % ind = (find(abs(u) > 10e5 & abs(v) > 10e5 & abs(w) > 10e5)) ; % for outliers
    % if length(ind) > 100
    ug(:,:,:,i) = griddata(x,y,z,u,xg,yg,zg,'natural') ; % interpolates scattered data and creates grid for the velocities
    vg(:,:,:,i) = griddata(x,y,z,v,xg,yg,zg,'natural') ;
    wg(:,:,:,i) = griddata(x,y,z,w,xg,yg,zg,'natural') ;
    % end
end

figure
quiver3(xg,yg,zg,mean(ug,4),mean(vg,4),mean(wg,4))