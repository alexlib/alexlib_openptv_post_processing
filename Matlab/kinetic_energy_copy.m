clear all

dist = @(x,y) sqrt((x(1,:)-y(1)).^2 + (x(2,:)-y(2)).^2 + (x(3,:)-y(3)).^2);

if ismac
    matdirectory = '/Users/alex/Dropbox/resuspension/2011/trajectories(186-194)';
else
    matdirectory ='C:\Users\hadar\Dropbox\resuspension\2011\trajectories(186-194)\events_90';
end

large_ones = dir(fullfile(matdirectory,'large*'));

time_of_event = {};
data_large = {};
data_small = {};


R = 25;

 %for all the large particles
for kFiles = 1: length(large_ones)   
    
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
            up = data.uf(ind);
            vp = data.vf(ind);
            wp = data.wf(ind);
            kinetic_energy_p = 1/2* (up^2 + vp^2 + wp^2 );
            
            % storage
            data_large{kFiles,iEvents,ind} = kinetic_energy_p;
            
            counter = 0;
            tmp_kinetic_energy = 0;
            
            for j = 1:numSmallTraj % for all small ones
                
                distance = dist([small(j).xf,small(j).yf,small(j).zf]',[xp,yp,zp])';
                
                relevant = find(small(j).t == data.t(ind) &  distance <= R);
                
                if ~isempty(relevant)
                    
                    tmp = small(j);
                    tmp.uf = tmp.uf(relevant);
                    tmp.vf = tmp.vf(relevant);
                    tmp.wf = tmp.wf(relevant);
                    
                    %tmp_kinetic_energy = tmp_kinetic_energy +1/2 *(tmp.uf^2+ tmp.vf^2 + tmp.wf^2); % sum of all the relevant small particles
                    %tmp_kinetic_energy = tmp_kinetic_energy +1/2 *(tmp.uf^2+ tmp.wf^2); % sum of all the relevant small particles
                    tmp_kinetic_energy = tmp_kinetic_energy +1/2 *( tmp.vf^2 ); % sum of all the relevant small particles
                    counter = counter + 1;
                    
                    % if verbose,plot_long_trajectories(tmp(relevant),1,hf);   end
                    
                end % relevant
            end % smalltraj
            
            data_small{kFiles,iEvents,ind} =  tmp_kinetic_energy / counter; % average the small particles kinetic energy
            time_data_small{kFiles, iEvents,ind} = data.t(ind);
            
            % i= number of events
            %ind= time of event
            
        end
    end
end

save kinetic_energy_R25_90.mat


%%
figure, hold on,
for kFiles = 1:size(data_small,1)
    for iEvents = 1:size(data_small,2)
        if ~isempty([time_of_event{kFiles,iEvents}])
            plot([time_data_small{kFiles,iEvents,:}]-[time_data_small{kFiles,iEvents,1}],[data_small{kFiles,iEvents,:}]); % x = time, y= kinetic energy small
            tmp = [data_small{kFiles,iEvents,:}];
            t1 = max([time_of_event{kFiles,iEvents}]-[time_data_small{kFiles,iEvents,1}]-1,1);
            plot(t1,tmp(t1+1),'ro','MarkerSize',10); % time of the event
        end
    end
end
hold off










