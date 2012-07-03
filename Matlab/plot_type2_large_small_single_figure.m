function [hf1] = plot_type2_large_small_single_figure(filename,large_quantity,small_quantity, R)
% PLOT_TYPE2_LARGE_SMALL(MATFILENAME,LARGE_QUANTITY,SMALL_QUANTITY, R)
% plots the data of LARGE_QUANTITY (e.g. 'uf') and of SMALL-QUANTITY (e.g. 'wf') for the event
% given by the MATFILENAME (.e.g 'large1934980')
%
% Following the mark_type2_events script
% we can load largeXXXXXX.mat and eventsXXXXXX.mat files
% and plot it in a variety of ways
% R is the distance in which we want to search the small particles.
% Inputs:
%       FILENAME - 'large1934980'
%       LARGE_QUANTITY,SMALL_QUANTITY - e.g. 'yf','uf','axf'
%       R = radius of the sphere in millimeters to include the small
%       particles
%
% Example:
%   plot_type2_large_small('large1934980','yf','yf')
%

% Copyright (c) 2011, Turbulence Structure Laboratory

dist = @(x,y) sqrt((x(1,:)-y(1)).^2 + (x(2,:)-y(2)).^2 + (x(3,:)-y(3)).^2);

if nargin == 2
    small_quantity = large_quantity;
    R = 15; % mm, maximum distance of the particles to plot
elseif nargin == 3
    R = 15; % mm
end


if ismac
    matdirectory = '/Users/alex/Dropbox/resuspension/2011/trajectories(186-194)';
else
    matdirectory ='C:\Users\hadar\Dropbox\resuspension\2011\trajectories(186-194)';
end
% large_ones = dir(fullfile(matdirectory,'large*'));

[path,filename,ext] = fileparts(filename);

% Load the largeXXXX.mat
large = load(fullfile(matdirectory,filename));
% and rename it to 'large'
large = large.(filename);

% Load the eventsXXXXX.mat
events = load(fullfile(matdirectory,strrep(filename,'large','events')));

% Load the small particle data and rename to 'small'
small = load(fullfile(matdirectory,strrep(filename,'large','small')));
small = small.(strrep(filename,'large','small'));

numEvents = size(events.events,1);
numTraj = length(small);


trajIds = unique(cat(1,large.trajid));


for i = 1:numEvents
    % Select the trajectory and the time of event
    time_of_event = events.events(i,2);
    trajind = trajIds == events.events(i,1);
    
    data = large(trajind);
    ind = find(data.t == time_of_event);
    
    hf1 = figure; hold on
    title(sprintf('Data: %s',filename))
    ylabel(sprintf('large quantity: %s',large_quantity));
    xlabel('Time [frames]')
    %     subplot(211), hold on
    %     subplot(212); hold on
    
    
    % Plot the large particles
    %     subplot(211)
    plot(data.t(3:end-2),data.(large_quantity)(3:end-2),'LineWidth',1,'DisplayName',sprintf('%d',data.trajid(1)));
    plot(data.t(ind),data.(large_quantity)(ind),'kx');%marks the time of event.
    
    %     subplot(212)
    for j = 1:numTraj
        [tf,loc] = ismember(time_of_event,small(j).t(3:end-2));
        if tf % i.e. if this small trajectory is also appearing at that time
            % position of the small particle at the time_of_event
            xf = small(j).xf(loc);
            yf = small(j).yf(loc);
            zf = small(j).zf(loc);
            
            % position of the large particle at the time of event:
            xp = data.xf(ind);
            yp = data.yf(ind);
            zp = data.zf(ind);
            
            distance = dist([xf,yf,zf]',[xp,yp,zp]);
            if distance <= R
                plot(small(j).t(3:end-2),small(j).(small_quantity)(3:end-2),'LineWidth',1,'Color','r','DisplayName',sprintf('%d',small(j).trajid(1)));
                plot(small(j).t(loc),small(j).(small_quantity)(loc),'kx');
            end
        end
    end
end


