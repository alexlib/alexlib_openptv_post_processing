function [events,large] = type2_event_v2(large,verbose)
% type2_event_v2  Marks the events in terms of vertical position
%
% [EVENTS] = TYPE2_EVENT_V2(LARGE,VERBOSE)
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
% Last modified:
% 14-12-2012 - new event detection


events = [];
numTraj = length(large);



% We define automatically the position of the ground
% we take the first locations of the particles
% Warning !!!
% this approach works only if there are more particles that really start on
% the ground and resuspend as compared to all kind of others. If there is
% one particle on the ground and others somewhere, it might break down.

yfs = zeros(numTraj,1);
for i = 1:numTraj
    yfs(i) = large(i).yf(1);
end
lower_bound = mode(round(yfs));
upper_bound = mode(ceil(yfs)) + .5 ; % half a millimeter above the most common (mode) value

large_particles_on_the_floor = [];
k = 1;
for i = 1:numTraj
    if large(i).yf(1) > lower_bound  && large(i).yf(1) <= upper_bound  && max(large(i).yf) > upper_bound + 2 % the first point of the trajectory is on the ground, empirically found
        large_particles_on_the_floor(k) = i;
        k = k + 1;
    end
end

if verbose,
    hf1 = figure; hold on
    for i = 1:length(large),
        if ismember(i,large_particles_on_the_floor)
            plot(large(i).t,large(i).yf,'Color','r','LineWidth',1,'DisplayName',sprintf('%d',large(i).trajid(1)));
            xlabel(sprintf('t (sec) %s'));
            ylabel(sprintf('y (mm) %s'));

        else
            plot(large(i).t,large(i).yf,'--','DisplayName',sprintf('%d',large(i).trajid(1)));
        end
    
    end
    title(sprintf('Data: %d',large(1).t(1)))
    hold off
end


large = large(large_particles_on_the_floor);
numTraj = length(large);

% figure, hold on
for i = 1:numTraj,
    
    % let's go from the edge of the upper_bound to the moment of the event,
    % backwards
    
    % when it crossed first time the upper bound:
    first_time_above_upper_bound = find(large(i).yf > upper_bound,1,'first');
    
    % now going backwards we find the last time the particle crossed
    % "jumped" above it's initial position which is varying
    
    relevant_part = 1:first_time_above_upper_bound;
    
    
    
    
    % event is when the y changes above the 3rd digit
    smoothy = smooth(large(i).yf(relevant_part));
    % event = find(diff(smoothy) >= 0 ,1, 'first');
    event = find((smoothy - min(smoothy)) > .1,1,'first');
    
    % event = find(diff(round(large(i).yf*100),1)==0,1,'last');
    % if the movement is downwards, it's not an event
    % HERE WE CAN IMPROVE BY ADDING AVERAGE VELOCITY OVER SOME PERIOD OF
    % TIME
    % if mean(diff(large(i).yf(event:event+5))) <= 0, event = []; end
    % if we loose the particle after less than 5 frames, it's not an event
    if large(i).t(end) <= large(i).t(event) + 5, event = []; end;
    
    
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
% keyboard






