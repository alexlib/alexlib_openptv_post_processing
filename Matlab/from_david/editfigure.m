%EDIT a figure
%update the file name and the file path

%uiopen('D:\3DPTV\Reutcavity-29dec08\matlab-process\matlab-res-scene25Re50050fps\Aug_present\3DmapvelocityNoedit.fig',1)
x = get(get(gca,'children'),'xdata');
y = get(get(gca,'children'),'ydata');
z = get(get(gca,'children'),'zdata');
u = get(get(gca,'children'),'udata');
v = get(get(gca,'children'),'vdata');
w = get(get(gca,'children'),'wdata');
 
xmin = min(x(:)); 
ymin = min(y(:)); 
zmin = min(z(:));
 
xmax = max(x(:)); 
ymax= max(y(:)); 
zmin = max(z(:));
 
[sx,sy,sz] = meshgrid(0,-0.05:0.01:0.05,-0.05:0.01:0.05);
hlines = streamline(x,y,z,u,v,w,sx,sy,sz);

%set(hlines, 'LineWidth',2, 'Color', 'r','LineStyle','-->')
set(hlines,'LineWidth',2,'Color','r')
 
view(3)
daspect([2,2,2])
axis tight
h=coneplot(x, y, z, u, v, w, sx, sy, sz,y, 1,'quiver');
set(h, 'Linewidth',2, 'color','p')
lighting gouraud

% and
%x-y view
figure
quiver(squeeze(x(:,:,8)),squeeze(y(:,:,8)),squeeze(u(:,:,8)),squeeze(v(:,:,8)),'LineWidth',2)
hold on
[sx,sy,sz] = meshgrid(-0.05:0.01:0.05,-0.05:0.01:0.05, 0);
hlines = streamline(x,y,z,u,v,w,sx,sy,sz);
set(hlines,'LineWidth',1,'Color','r')
xlabel('x(m)'); ylabel('y(m)');
%x-z view
figure
quiver(squeeze(x(8,:,:)),squeeze(z(8,:,:)),squeeze(u(8,:,:)),squeeze(w(8,:,:)),'LineWidth',2)
%hold on
%[sx,sy,sz] = meshgrid(-0.05:0.01:0.05,0,-0.05:0.01:0.05);
%hlines = streamline(x,y,z,u,v,w,sx,sy,sz);
%set(hlines,'LineWidth',1,'Color','r')
xlabel('x(m)'); ylabel('z(m)');
%y-z view
figure
quiver(squeeze(y(:,8,:)),squeeze(z(:,8,:)),squeeze(v(:,8,:)),squeeze(w(:,8,:)),'LineWidth',2)
%hold on
%[sx,sy,sz] = meshgrid(0,-0.05:0.01:0.05,-0.05:0.01:0.05);
%hlines = quiver(x,y,z,u,v,w,sx,sy,sz);%or plot-check
%set(hlines,'LineWidth',1,'Color','r')
xlabel('y(m)'); ylabel('z(m)');


