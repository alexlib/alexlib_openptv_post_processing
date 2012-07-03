function [traj] = link_2d_trajectories_smoothn(traj,s)
% TRAJ = LINK_TRAJECTORIES_SMOOTHN(TRAJ,S)

if nargin < 2
    s = 1000;
end

x = cat(1,traj.xf);
y = cat(1,traj.yf);
t = cat(1,traj.t);

x = x(:)';
y = y(:)';
t = t(:)';

dt = min(diff(t));

% ti = linspace(t(1),t(end),length(t));

traj.xf = smoothn(x,s,'robust');
traj.yf = smoothn(y,s,'robust');

% plot3(x,y,z,'r.',xs,ys,zs,'b-',xss,yss,zss,'k--')

% traj.uf = gradient5(traj.xf,t);
% traj.vf = gradient5(traj.yf,t);
% 
% traj.axf = gradient5(traj.uf,t);
% traj.ayf = gradient5(traj.vf,t);

traj = ensure_order_traj(traj);

