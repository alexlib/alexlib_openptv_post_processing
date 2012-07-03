function varargout = plot_2d_trajectories(traj,minLength,figureNum)
% plot_2d_trajectories
%
% SEE ALSO: building_trajectories.m

% Author: Alex Liberzon
% Copyright (c) 2007, TAU

if nargin > 2
    figure1 = figureNum;
    color = 'r';
else
    % Create figure
    figure1 = figure;
    color = 'b';
    % Create axes
    % axes('Parent',figure1,...
    %     'PlotBoxAspectRatio',[1 1 1],...
    %     'FontSize',12,...
    %     'DataAspectRatio',[1 1.6 1.6],...
    %     'Color',[0.9412 0.9412 0.9412],...
    %     'CameraViewAngle',10.34,...
    %     'CameraUpVector',[0 1 0],...
    %     'CameraPosition',[0.2148 -0.2713 0.4102]);
    axes('Parent',figure1);
    box('on');
    grid('on');
    hold('all');
end



% Create xlabel
xlabel('$x$','Interpreter','Latex');

% Create ylabel
ylabel('$z$','Interpreter','Latex');


if nargin < 2, minLength = 10; end

trajLen = zeros(length(traj),1);
only_long = isfield(traj,'trajid');
idx = 1:length(traj);
switch only_long
    case 1
        for i = 1:length(traj)
            trajLen(i) = length(traj(i).trajid);
            if traj(i).trajid == -999
                trajLen(i) = 0;
            end
        end
        id = find(trajLen >= minLength);
    case 0
        id = 1:length(traj);
end

if nargin < 3
    for i = 1:length(id)
        %     plot3(traj(id(i)).xf/1000,traj(id(i)).yf/1000,traj(id(i)).zf/1000,...
        plot(traj(id(i)).xf,traj(id(i)).yf,...
            'MarkerFaceColor','w','Marker','o',...
            'MarkerSize',3,...
            'UserData',id(i));
    end
else
    
    for i = 1:length(id)
        plot(traj(id(i)).xf,traj(id(i)).yf,...
            'MarkerFaceColor','w','Marker','*',...
            'MarkerSize',3,...
            'UserData',id(i));
    end
end

if nargout
    varargout{1} = trajLen;
end
set(gca,'ydir','reverse')
