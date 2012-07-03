function [newtraj] = link_trajectories_cubicspline(traj)
% should work for scattered data
x = cat(1,traj.xf);
y = cat(1,traj.yf);
z = cat(1,traj.zf);
t = cat(1,traj.t);

x = x(:)';
y = y(:)';
z = z(:)';
t = t(:)';

dt = min(diff(t));

ti = linspace(t(1),t(end),length(t));



[xi,p] = csaps(t,x,[],ti);
[yi,p] = csaps(t,y,p,ti);
[zi,p] = csaps(t,z,p,ti);

newtraj.xf = xi;
newtraj.yf = yi;
newtraj.zf = zi;


newtraj.uf = gradient5(newtraj.xf,ti*dt);
newtraj.vf = gradient5(newtraj.yf,ti*dt);
newtraj.wf = gradient5(newtraj.zf,ti*dt);

newtraj.axf = gradient5(newtraj.uf,ti*dt);
newtraj.ayf = gradient5(newtraj.vf,ti*dt);
newtraj.azf = gradient5(newtraj.wf,ti*dt);

newtraj.t = ti;

newtraj.trajid = traj(1).trajid(1)*ones(length(xi),1);
newtraj = ensure_order_traj(newtraj);

