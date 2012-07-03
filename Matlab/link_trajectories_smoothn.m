function [traj] = link_trajectories_smoothn(traj,s)
% TRAJ = LINK_TRAJECTORIES_SMOOTHN(TRAJ,S)

if nargin < 2
    s = [];
end

x = cat(1,traj.xf);
y = cat(1,traj.yf);
z = cat(1,traj.zf);
t = cat(1,traj.t);

x = x(:)';
y = y(:)';
z = z(:)';
t = t(:)';

dt = t(2)-t(1);

% ti = linspace(t(1),t(end),length(t));

traj.xf = smoothn(x,s,'robust');
traj.yf = smoothn(y,s,'robust');
traj.zf = smoothn(z,s,'robust');

% plot3(x,y,z,'r.',xs,ys,zs,'b-',xss,yss,zss,'k--')

traj.uf = gradient5(traj.xf,t);
traj.vf = gradient5(traj.yf,t);
traj.wf = gradient5(traj.zf,t);

traj.axf = gradient5(traj.uf,t);
traj.ayf = gradient5(traj.vf,t);
traj.azf = gradient5(traj.wf,t);

traj = ensure_order_traj(traj);

