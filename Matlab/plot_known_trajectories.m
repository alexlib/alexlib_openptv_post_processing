function varargout = plot_known_trajectories(data,trajID2Plot)
% Plots known trajectories
% Usage:
% PLOT_known_trajectories(data,trajID2Plot)
% where DATA is a structure that includes 
% field 'traj' and  input trajID2Plot is the 
% indexes of the needed to plot trajectories.
%
% SEE ALSO: building_trajectories.m

% Author: Alex Liberzon
% Copyright (c) 2007, TAU

trajID = cat(2,data.traj);
% trajLen = trajID*0;

figure, hold on
xlabel('x'),ylabel('y'),zlabel('z')
box on
view(3)
xf = cat(1,data.xf);
yf = cat(1,data.yf);
zf = cat(1,data.zf);

for i = 1:length(trajID2Plot)
    traj = find(trajID == trajID2Plot(i));
    plot3(xf(traj(2)),yf(traj(2)),zf(traj(2)),'ro');
    plot3(xf(traj(3:end-1)),yf(traj(3:end-1)),zf(traj(3:end-1)),'.','MarkerSize',4);
  end
% varargout{1} = trajLen;