function plot_long_trajectories(data,minLength)
% Plots long trajectories
% Usage:
% PLOT_LONG_TRAJECTORIES(SCENE5,<MINLENGTH>)
% where SCENE5 is a structure that includes 
% field 'traj' and optional input MINLENGTH is the 
% shortest length of the trajectory to plot.
%
% SEE ALSO: building_trajectories.m

% Author: Alex Liberzon
% Copyright (c) 2007, TAU

trajID = cat(2,data.traj);
trajLen = trajID*0;

figure, hold on
minLength = 20;
xf = cat(1,data.xf);
yf = cat(1,data.yf);
zf = cat(1,data.zf);

for i = 1:max(trajID)
    traj = find(trajID == i);
    trajLen(i) = length(traj);
    if length(traj) >= minLength
        plot3(xf(traj(1)),yf(traj(1)),zf(traj(1)),'ro');
        plot3(xf(traj),yf(traj),zf(traj),'.');
    end
end
