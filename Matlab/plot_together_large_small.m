function varargout = plot_together_large_small(large,small)
% Plots trajectories of two datasets together
% Usage:
% PLOT_LONG_TRAJECTORIES(TRAJ,<MINLENGTH>)
% where SCENE5 is a structure that includes
% field 'traj' and optional input MINLENGTH is the
% shortest length of the trajectory to plot.
%
% SEE ALSO: building_trajectories.m

% Author: Alex Liberzon
% Copyright (c) 2011, TAU

small_t_span = cat(1,small.t);
% small_trajid = cat(1,small.trajid);
small_traj_ind = [];
for i = 1:length(small)
    small_traj_ind = cat(1,small_traj_ind,i*ones(length(small(i).xf),1));
end




large_trajlen = zeros(length(large),1);
going_up = [];
for i = 1:length(large)
    large_trajlen(i) = length(large(i).xf);
    if any(large(i).yf > -20) && any(large(i).yf(2:end) > large(i).yf(1))
        going_up = cat(1,going_up,i);
    end
end
% longest = find(large_trajlen >= prctile(large_trajlen,0));
% large = large(intersect(longest,going_up));
large = large(going_up);


% hf = figure;
% hold on;

for i = 1:length(large) % for all trajectories of the large particles
    tspan = large(i).t;
    figure;
    hold on;
    plot3d_traj(large(i),1);
    
    relevant = [];
    for j = 1:length(tspan) % last 10 frames are probably not interesting
        ind = small_t_span == tspan(j);
        relevant = cat(1,relevant,unique(small_traj_ind(ind)));
    end
    relevant = unique(relevant);
    
    for k = 1:length(relevant)
        plot3d_traj(small(relevant(k)),0);
    end
    
end
% drawnow;
% pause
end


