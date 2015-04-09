clc, clear all

dist = @(x,y) sqrt((x(1,:)-y(1)).^2 + (x(2,:)-y(2)).^2 + (x(3,:)-y(3)).^2);

if ismac
    matdirectory = '/Users/alex/Dropbox/resuspension/2011/trajectories(186-194)';
else
    matdirectory ='C:\Users\hadar\Dropbox\resuspension\2011\trajectories(186-194)\events_110';
end

large_ones = dir(fullfile(matdirectory,'large*'));

time_of_event = {};
%data_large = {};
%data_small = {};
data_force = {};

R = 25; %[mm]

density_l = 0.998 ; %density of liquid [gr/cm^3]
volume_liquid = 4* pi* (R/10)^3 / 3; % [cm^3]
mass_liquid = density_l * volume_liquid/1000; % [Kg]

%for all the large particles
for kFiles = 1:length(large_ones)  
    
    % get the name of the data
    filename = large_ones(kFiles).name;
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
    numSmallTraj = length(small);
    
    trajIds = unique(cat(1,large.trajid));
    
    
    
    for iEvents = 1:numEvents
        % Select the trajectory and the time of event
        time_of_event{kFiles,iEvents} = events.events(iEvents,2);
        trajind = trajIds == events.events(iEvents,1);
        
        % Take the single trajectory in which we had an event aside into
        % the temprorary 'data'
        
        data = large(trajind);
        
        
        % Along that trajectory, at each time step
        
        for ind = 1:length(data.t)
            
            xp = data.xf(ind);
            yp = data.yf(ind);
            zp = data.zf(ind);
           
            counter = 0;
            force = 0;
            
            for j = 1:numSmallTraj % for all small ones
                
                distance = dist([small(j).xf,small(j).yf,small(j).zf]',[xp,yp,zp])';
                
                relevant = find(small(j).t == data.t(ind) &  distance <= R);
                
                if ~isempty(relevant)
                    
                    tmp = small(j);
                    tmp.axf = tmp.axf(relevant);
                    tmp.ayf = tmp.ayf(relevant);
                    tmp.azf = tmp.azf(relevant);
                                      
                    force = force + mass_liquid * tmp.ayf;
                    
                    counter = counter + 1;
                    
                 
                end % relevant
            end % smalltraj
            
            %data_small{kFiles,iEvents,ind} =  tmp_velocity;
            data_force{kFiles,iEvents,ind} =  norm (force) / counter; 
            time_data_small{kFiles, iEvents,ind} = data.t(ind);
            
          
            
        end
    end
end

save force_R25_liquid_80_y.mat


%%
figure, hold on,
for kFiles = 1:size(data_force,1)
    for iEvents = 1:size(data_force,2)
        if ~isempty([time_of_event{kFiles,iEvents}])
            plot([time_data_small{kFiles,iEvents,:}]-[time_data_small{kFiles,iEvents,1}],[data_force{kFiles,iEvents,:}]); 
            tmp = [data_force{kFiles,iEvents,:}];
            t1 = max([time_of_event{kFiles,iEvents}]-[time_data_small{kFiles,iEvents,1}]-1,1);
            plot(t1,tmp(t1+1),'ro','MarkerSize',10); % time of the event
           
        end
    end
end
%xlim = get(gca,'xlim');
%plot(xlim,[f_drag_critical,f_drag_critical], 'r'); % make a line
%plot(xlim,[f_drag_critical_1,f_drag_critical_1], 'g'); % make a line
%hold off
