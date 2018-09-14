
function grid = scatter_to_grid( scatter_t, n)
% 
%

x1 = min(cat(1,scatter_t.x));
x2 = max(cat(1,scatter_t.x));

y1 = min(cat(1,scatter_t.y));
y2 = max(cat(1,scatter_t.y));


z1 = min(cat(1,scatter_t.z));
z2 = max(cat(1,scatter_t.z));

[nx,ny,nz] = deal(n);

grid = repmat(struct( 'x' ,[] , 'y' ,[] , 'z' ,[] , 'u' ,[] , 'v' ,[] , 'w' ,[]), size(scatter_t,4),1);


[grid.x,grid.y,grid.z] = ndgrid(linspace(x1,x2,nx),linspace(y1,y2,ny),linspace(z1,z2,nz));

[grid.u,grid.v,grid.w] = deal(repmat(grid.x,1,1,1,length(scatter_t)));

for i=1:length(scatter_t)
    F =  scatteredInterpolant(scatter_t(i).x,scatter_t(i).y,scatter_t(i).z,...
        scatter_t(i).u) ;
    grid.u(:,:,:,i) = F(grid.x,grid.y,grid.z);
    F =  scatteredInterpolant(scatter_t(i).x,scatter_t(i).y,scatter_t(i).z,...
        scatter_t(i).v) ;
    grid.v(:,:,:,i) = F(grid.x,grid.y,grid.z);
    F =  scatteredInterpolant(scatter_t(i).x,scatter_t(i).y,scatter_t(i).z,...
        scatter_t(i).w) ;
    grid.w(:,:,:,i) = F(grid.x,grid.y,grid.z);
end
