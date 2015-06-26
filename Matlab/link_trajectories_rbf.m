function [newtraj] = link_trajectories_rbf(traj,smoothness)
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


% figure
% plot(t,x,'o',t,y,'s',t,z,'^');
% hold on

if nargin < 2
    smoothness = 0.1;
end


coeffx = rbfcreate(t,x,'RBFsmooth',smoothness);
coeffy = rbfcreate(t,y,'RBFsmooth',smoothness);
coeffz = rbfcreate(t,z,'RBFsmooth',smoothness);



xi = rbfinterp(ti,coeffx);
yi = rbfinterp(ti,coeffy);
zi = rbfinterp(ti,coeffz);


% plot(ti,xi,'r--',ti,yi,'g--',ti,zi,'b--')


newtraj.xf = xi;
newtraj.yf = yi;
newtraj.zf = zi;
%{
newtraj.uf = gradient(newtraj.xf,ti/50);
newtraj.vf = gradient(newtraj.yf,ti/50);
newtraj.wf = gradient(newtraj.zf,ti/50);

newtraj.axf = gradient(newtraj.uf,ti/50);
newtraj.ayf = gradient(newtraj.vf,ti/50);
newtraj.azf = gradient(newtraj.wf,ti/50);
%}

newtraj.uf = gradient5(newtraj.xf,ti*dt);
newtraj.vf = gradient5(newtraj.yf,ti*dt);
newtraj.wf = gradient5(newtraj.zf,ti*dt);

newtraj.axf = gradient5(newtraj.uf,ti*dt);
newtraj.ayf = gradient5(newtraj.vf,ti*dt);
newtraj.azf = gradient5(newtraj.wf,ti*dt);

newtraj.t = ti;

newtraj.trajid = traj(1).trajid(1)*ones(length(xi),1);
newtraj = ensure_order_traj(newtraj);

