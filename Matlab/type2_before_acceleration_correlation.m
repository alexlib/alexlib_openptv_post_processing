function [data_large,data_small] = type2_before_acceleration_correlation(large,small,time_interval,R, verbose,figureNum)
% type2_before  is similar to type2  and it  defines  the velocity of the  small particles  before the movement of the large one. 
%TYPE2 should solve the problem of complex selections of events and the
% respective flow fields.
% [DATA_LARGE,DATA_SMALL] = ...
%   TYPE2(LARGE,SMALL,TIME_INTERVAL,CLOSE_ENOUGH,VERBOSE)
% Inputs:
%  large - data structure of large particles
%  small - -//- of tracers
%  time_interval - how many time frames before the event we include in the
%  analysis
%  close_enough - how big is the sphere of the "influence" in mm
%  R - radius of sphera 
%  verbose  - 0 for a silent mode (default if empty), 1 for graphics
% Example:
%
%   >> load large_10680
%   >> load small_10680
%   >> selection = type2_before(large_10680,small_10680,10,30,30,1)
%

% lets' define type2 event
% plot the vertical position of the date to see if we can see the events
% there:

data_large = [];
data_small = [];

dist = @(x,y) sqrt((x(1,:)-y(1)).^2 + (x(2,:)-y(2)).^2 + (x(3,:)-y(3)).^2);

numTraj = length(large);
numSmallTraj = length(small);

% close_enough = 15; % 5 mm from the particle
% time_interval = 10;
if nargin < 5
    verbose = 0; % 1 if you want to see the plots
end

% hf2 = figure; hold on;


%figure, hold on
%for i = 1:numTraj,
% plot(large(i).t,large(i).yf,'DisplayName',sprintf('%d',large(i).trajid(1)));
%ylabel('y');
%end; hold off
% clearly we see there particles that were too high
% we see those that were at the level of the wall and then raised up

% velocity is less clear in this case, one has to know the position 


% manually fond for large_10680
% lower_bound = -22; % mm
% upper_bound = -21.5; % mm

% let's try to define automatically the position of the ground
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
 else
   plot(large(i).t,large(i).yf,'--');
    end
 end
hold off
end


% Now, let's say we have the rolling/sliding particles on the ground
% identified
% we should find now the type2 event locations for each
% lets' define Type2 event as
% detachment

% take only the rolling ones
% and using the two more significant digits (x 100) and rounding it
% try to find the place where the location really changes, more than 0.5 mm
% or so
% It works apparently for the date from 10680

%{
clear all
load large_10680
load small_10680
type2(large_10680,small_10680)
%}


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
    % if the particle is back below the upper_bound within  10 frames, it's
    % out
    if all(large(i).yf(event:end) < upper_bound), event = []; end
    
    
    if ~isempty(event)  && event > time_interval,
        %plot(large(i).t,large(i).yf,large(i).t(event),large(i).yf(event),'o');
        
        % so we have an event here
        time_of_event = large(i).t(event); % we can decide if we take two steps back or so
        x_of_event = large(i).xf(event-time_interval);
        y_of_event = large(i).yf(event-time_interval);
        z_of_event = large(i).zf(event-time_interval);
        v_of_event = large(i).vf(event-time_interval);
        u_of_event = large(i).uf(event-time_interval);
        w_of_event = large(i).wf(event-time_interval);
        % we take first the u component of the particle at the event
        % as is, with its direction
        % WE SHALL CHANGE IT FOR OTHER CONDITIONS
        data_large = cat(1,data_large,max(large(i).ayf(event-time_interval:event)));
        
        
        if verbose,
            figure(hf1);hold on;
            plot(large(i).t,large(i).yf,'LineWidth',2,'Color','r','DisplayName',sprintf('%d',large(i).trajid(1)));
            hold off
        end

        
        % if we need statistics of events, we can collect it here as well
        % e.g. u,v,w,ax,ay,az, ....
        
        
        % A bit before the event, let's take 20 time steps
        % 160 fps => 0.125 sec = 1/8 sec
        before_event = time_of_event-time_interval:time_of_event;
        
        % figure, plot_together_large_small(large(i),small);
        if verbose
%             hf = figure;
           [junk,hf] =  plot_long_trajectories(large(i));
        end
        
        [xs,zs]= range_before_event (x_of_event, z_of_event,u_of_event,w_of_event,R) ;
        
        data = [];
        % let's find the small ones
        for ii = 1:numSmallTraj % for all small ones
            relevant = find(small(ii).t >= before_event(1) & small(ii).t <= time_of_event);
            if ~isempty(relevant)
                
                tmp = small(ii);
                tmp.xf = tmp.xf(relevant);
                tmp.yf = tmp.yf(relevant);
                tmp.zf = tmp.zf(relevant);
                tmp.t = tmp.t(relevant);
                tmp.trajid = tmp.trajid(relevant);
                
                % Plot all the relevant ones:
%                 if verbose,  plot_long_trajectories(tmp,1,hf); end
                
                % let's find those that are close_enough
                
                %                 for j = 1:length(relevant)
                %                     if sqrt((small(ii).xf(relevant(j))-x_of_event)^2 +...
                %                             (small(ii).yf(relevant(j))-y_of_event)^2+...
                %                             (small(ii).zf(relevant(j))-z_of_event)^2) <= close_enough
                %                         data(k) = small(ii).uf(relevant(j));
                %                         k = k + 1;
                %                     end
                %                 end
                distances = dist([tmp.xf,tmp.yf,tmp.zf]',[xs,y_of_event,zs]);
                if any(distances <= R)
                    figure(figureNum);
                    plot(distances(distances<=R),tmp.ayf(distances<=R)*max(large(i).ayf(event-time_interval:event)),'*');
                    drawnow
                    % Plot only the closest relevant ones
                    if verbose,plot_long_trajectories(tmp,1,hf); end
                    % data = cat(1,data,tmp.uf(distances <= R));
                        data = cat(1,data,tmp.ayf(distances <= R));
                end
            end
            
        end
        % data_small = cat(1,data_small,mean(data));
        data_small = cat(1,data_small,data);
    end
    
end

disp(sprintf('No. of events = %d\n',length(data_small)));








