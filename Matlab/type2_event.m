function [events] = type2_event(large,verbose)
% type2_event  Marks the events in terms of vertical position
%
% [EVENTS] = TYPE2(LARGE,VERBOSE)
%
% Inputs:
%  large - data structure of large particles
%  verbose  - 0 for a silent mode (default if empty), 1 for graphics
%
% See also: MARK_TYPE2_EVENT
%
%

% Copyright (c) 2011, Turbulence Structure Laboratory, Tel Aviv University
%
% Author: Hadar Traugott and Alex Liberzon
%


events = [];
numTraj = length(large);

if nargin < 2
    verbose = 0; % 1 if you want to see the plots
end


% We define automatically the position of the ground
% we take the first locations of the particles
% Warning !!!
% this approach works only if there are more particles that really start on
% the ground and resuspend as compared to all kind of others. If there is
% one particle on the ground and others somewhere, it might break down.

ids = cat(1,large.trajid);
yfs = cat(1,large.yf);
[junk,j] = unique(ids,'first');
lower_bound = mode(round(yfs(j)));
upper_bound = mode(ceil(yfs(j)));

large_particles_on_the_floor = [];
k = 1;
for i = 1:numTraj
    if large(i).yf(1) >= lower_bound  && large(i).yf(1) <= upper_bound % the first point of the trajectory is on the ground, empirically found
        large_particles_on_the_floor(k) = i;
        k = k + 1;
    end
end

if verbose,
    hf1 = figure; hold on
    for i = 1:numTraj,
        if ismember(i,large_particles_on_the_floor)
            plot(large(i).t,large(i).yf,'LineWidth',2,'DisplayName',sprintf('%d',large(i).trajid(1))); 
xlabel(sprintf('t (sec) %s'));
ylabel(sprintf('y (mm) %s'));
        else
            plot(large(i).t,large(i).yf,'--');
        end
    end
    hold off
end


large = large(large_particles_on_the_floor);
numTraj = length(large);

% figure, hold on
for i = 1:numTraj,
    % event is when the y changes above the 3rd digit
    smoothy = round(large(i).yf*100);
    event = find(diff(smoothy) >=2,1,'first');
    
    % event = find(diff(round(large(i).yf*100),1)==0,1,'last');
    % if the movement is downwards, it's not an event
    if diff(large(i).yf(event:event+1)) < 0, event = []; end
    % if we loose the particle after less than 10 frames, it's not an event
    if large(i).t(end) <= large(i).t(event) + 10, event = []; end;
    
    % new condition for the jumping case (trajid=1, large1896330)
    % if the particle is back below the upper_bound, it's out
    if all(large(i).yf(event:end) < upper_bound), event = []; end
    
    
    if ~isempty(event)
        
        % so we have an event here
        events = cat(1,events,[large(i).trajid(1),large(i).t(event)]); %the first row in the events is the index, the second row is the time of the event.
        
        
        if verbose,
            figure(hf1);hold on;
            plot(large(i).t,large(i).yf,'LineWidth',1,'Color','g','DisplayName',sprintf('%d',large(i).trajid(1)));
            plot(large(i).t(event),large(i).yf(event),'ro');%marks the event 

            hold off
        end
        
        
    end
    
end






