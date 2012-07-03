function [newtraj] = link_2D_trajectories_rbf(traj,smoothness)
% should work for scattered data
x = cat(1,traj.xf);
y = cat(1,traj.yf);
t = cat(1,traj.t);

x = x(:)';
y = y(:)';
t = t(:)';
tinit = t(1);
t = t - tinit;

dt = min(diff(t));

ti = linspace(t(1),t(end),length(t));


% figure
% plot(t,x,'o',t,y,'s',t,z,'^');
% hold on

if nargin < 2
    smoothness = 0.1;
end


coeffx = rbfcreate(t,x,'RBFsmooth',smoothness,'RBFFunction','cubic');
coeffy = rbfcreate(t,y,'RBFsmooth',smoothness,'RBFFunction','cubic');



xi = rbfinterp(ti,coeffx);
yi = rbfinterp(ti,coeffy);


% plot(ti,xi,'r--',ti,yi,'g--',ti,zi,'b--')


newtraj.xf = xi;
newtraj.yf = yi;
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

newtraj.axf = gradient5(newtraj.uf,ti*dt);
newtraj.ayf = gradient5(newtraj.vf,ti*dt);

newtraj.t = ti+tinit;

newtraj.trajid = traj(1).trajid(1)*ones(length(xi),1);
newtraj = ensure_order_traj(newtraj);

