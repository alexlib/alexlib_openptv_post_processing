function varargout = plot_long_trajectories(traj,minLength,figureNum)
% Plots long trajectories
% Usage:
%   PLOT_LONG_TRAJECTORIES(TRAJ,<MINLENGTH>,<FIGUREHANDLE>)
% where TRAJ is a structure that includes
% field 'traj' and optional input MINLENGTH is the
% shortest length of the trajectory to plot.
% If FigureHandles is an open figure, adds the trajectories to that plot
% with some random color and without markers.
%
% SEE ALSO: building_trajectories.m

% Author: Alex Liberzon
% Copyright (c) 2007,2015 TAU

% default minimum length to plot, saves memory and space
if nargin < 2, minLength = 10; end

trajLen = zeros(length(traj),1);
only_long = isfield(traj,'trajid');
if only_long % we can select according to length
    for i = 1:length(traj)
        trajLen(i) = length(traj(i).trajid);
        if traj(i).trajid == -999
            trajLen(i) = 0;
        end
    end
    id = find(trajLen >= minLength);
else % we cannot, take everything?
    id = 1:length(traj);
end

if nargin < 3 % new figure, use markers
    figure1 = figure;
    axes('Parent',figure1);
    box('on');
    grid('on');
    hold('all');
    
    % Create labels
    xlabel('$x$','Interpreter','Latex');
    ylabel('$z$','Interpreter','Latex');
    zlabel('$y$','Interpreter','Latex');
    
    for i = 1:length(id)
        %     plot3(traj(id(i)).xf/1000,traj(id(i)).yf/1000,traj(id(i)).zf/1000,...
        plot3(traj(id(i)).xf,traj(id(i)).zf,traj(id(i)).yf,...
            'MarkerFaceColor','w','Marker','o',...
            'MarkerSize',6,...
            'UserData',traj(id(i)).trajid(1));
        plot3(traj(id(i)).xf(1),traj(id(i)).zf(1),traj(id(i)).yf(1),...
            'MarkerFaceColor','w','Marker','>',...
            'MarkerSize',8);
    end
else % old figure, adding trajectories single color, single marker
    figure1 = figure(figureNum);
    hold on
    color = rand(3,1);
    for i = 1:length(id)
        plot3(traj(id(i)).xf,traj(id(i)).zf,traj(id(i)).yf,...
            'MarkerEdgeColor',color,...
            'MarkerFaceColor',color,...
            'Marker','.',...
            'MarkerSize',6,...
            'UserData',traj(id(i)).trajid(1));
    end
end

% campos([-0.2208   -0.2416    0.6716]);
% camup([-0.1478    0.9178    0.4739]);
% camtarget([-0.0119    0.0014    0.0014]);

if nargout == 1
    varargout{1} = figure1;
elseif nargout == 2
    varargout{1} = figure1;
    varargout{2} = TrajLen;
end
view(3)
