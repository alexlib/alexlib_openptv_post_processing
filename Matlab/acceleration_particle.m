clc, clear all

dist = @(x,y) sqrt((x(1,:)-y(1)).^2 + (x(2,:)-y(2)).^2 + (x(3,:)-y(3)).^2);

if ismac
    matdirectory = '/Users/alex/Dropbox/resuspension/2011/trajectories_2012(186-194)';
else
    matdirectory ='C:\Users\hadar\Dropbox\resuspension\Matlab\trajectories_alex_2012(186-194)\events_110';
end

large_ones = dir(fullfile(matdirectory,'large*'));


 [data_acceleration_f1, data_acceleration_p_f1,time_data_large...
   time_of_event,data_acceleration_p] = deal({});

R = 25; %mm


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
        data = large(trajind);
        
        
        % Along that trajectory, at each time step
        
        for ind = 1:length(data.t)
            
            xp = data.xf(ind);
            yp = data.yf(ind);
            zp = data.zf(ind);

            axp = data.xf(ind)/1000;
            ayp = data.yf(ind)/1000;
            azp = data.zf(ind)/1000;
            
            acceleration_p = [axp,ayp,azp];% [m/s^2]
%             
%             data_large_y{kFiles,iEvents,ind} = ayp ;
%             data_large_x{kFiles,iEvents,ind} = axp ; 
%             data_large_z{kFiles,iEvents,ind} = azp ;
%             
            time_data_large{kFiles, iEvents,ind} = data.t(ind);
     
            acceleration_f1 = zeros(1,3);
            acceleration_p_f1= zeros(1,3);
            
            counter = 0;
            
            for j = 1:numSmallTraj % for all small ones
                counter = counter + 1; 
                distance = dist([small(j).xf,small(j).yf,small(j).zf]',[xp,yp,zp])';
                
                relevant = find(small(j).t == data.t(ind) &  distance <= R);
                
                if isempty(relevant), continue, end
                
                tmp = small(j);
                
                if tmp.yf(relevant) <= -21.5 && abs(tmp.yf(relevant)-yp )>= 10 ,continue, end
               
                
                axf = tmp.axf(relevant)/1000;
                ayf= tmp.ayf(relevant)/1000;
                azf= tmp.azf(relevant)/1000;
 
                acceleration_f = [axf,ayf,azf]; %[m/s^2]
                acceleration_f1= acceleration_f1+ acceleration_f;
                
                acceleration_p_f = [axp-axf,ayp-ayf,azp-azf];
                acceleration_p_f1 = acceleration_p_f1 + acceleration_p_f;
                
               

            end % smalltraj
            
            if counter > 0
                % Averaging in the sphere of radius R

                data_acceleration_f1{kFiles,iEvents,ind} = acceleration_f1/counter;
                data_acceleration_p_f1{kFiles,iEvents,ind} = acceleration_p_f1/counter;
            else

                data_acceleration_f1{kFiles,iEvents,ind} = zeros(1,3);
                data_acceleration_p_f1{kFiles,iEvents,ind} = zeros(1,3);
            end
            
            data_acceleration_p{kFiles,iEvents,ind} = acceleration_p;

        end %ind
    end %iEvents
end %kFiles

save acceleration_110.mat


%%
figure, hold on,
for kFiles = 1:size(data_acceleration_p,1)
    for iEvents = 1:size(data_acceleration_p,2)
        if ~isempty([time_of_event{kFiles,iEvents}]);
            tmp = cat(1,data_acceleration_p{kFiles,iEvents,:});
            vertical_tmp = tmp(:,2);
            plot([time_data_large{kFiles,iEvents,:}]-[time_data_large{kFiles,iEvents,1}],vertical_tmp);
            t1 = max([time_of_event{kFiles,iEvents}]-[time_data_large{kFiles,iEvents,1}]-1,1);
            plot(t1,vertical_tmp(t1+1),'r*','MarkerSize',10); % time of the event
            
        end
    end
end
xlabel('time from the first event [frame]')
ylabel('vertical acceleration particle [m/s^2]');

figure, hold on,
for kFiles = 1:size(data_acceleration_p,1)
    for iEvents = 1:size(data_acceleration_p,2)
        if ~isempty([time_of_event{kFiles,iEvents}]);
            tmp = cat(1,data_acceleration_p{kFiles,iEvents,:});
            horizontal_tmp = tmp(:,1);
            plot([time_data_large{kFiles,iEvents,:}]-[time_data_large{kFiles,iEvents,1}],horizontal_tmp);
            t1 = max([time_of_event{kFiles,iEvents}]-[time_data_large{kFiles,iEvents,1}]-1,1);
            plot(t1,horizontal_tmp(t1+1),'r*','MarkerSize',10); % time of the event
            
        end
    end
end
xlabel('time from the first event [frame]')
ylabel('horizontal acceleration particle_x [m/s^2]');

                
figure, hold on,
for kFiles = 1:size(data_acceleration_p,1)
    for iEvents = 1:size(data_acceleration_p,2)
        if ~isempty([time_of_event{kFiles,iEvents}]);
            tmp = cat(1,data_acceleration_p{kFiles,iEvents,:});
            horizontal_tmp = tmp(:,3);
            plot([time_data_large{kFiles,iEvents,:}]-[time_data_large{kFiles,iEvents,1}],horizontal_tmp);
            t1 = max([time_of_event{kFiles,iEvents}]-[time_data_large{kFiles,iEvents,1}]-1,1);
            plot(t1,horizontal_tmp(t1+1),'r*','MarkerSize',10); % time of the event
            
        end
    end
end
xlabel('time from the first event [frame]')
ylabel('horizontal acceleration fluid_z [m/s^2]');

%%%%%%%%%%%%

figure, hold on,
for kFiles = 1:size(data_acceleration_p,1)
    for iEvents = 1:size(data_acceleration_p,2)
        if ~isempty([time_of_event{kFiles,iEvents}]);
            tmp = cat(1,data_acceleration_f1{kFiles,iEvents,:});
            vertical_tmp = tmp(:,2);
            plot([time_data_large{kFiles,iEvents,:}]-[time_data_large{kFiles,iEvents,1}],vertical_tmp);
            t1 = max([time_of_event{kFiles,iEvents}]-[time_data_large{kFiles,iEvents,1}]-1,1);
            plot(t1,vertical_tmp(t1+1),'r*','MarkerSize',10); % time of the event
            
        end
    end
end
xlabel('time from the first event [frame]')
ylabel('vertical acceleration fluid [m/s^2]');

figure, hold on,
for kFiles = 1:size(data_acceleration_p,1)
    for iEvents = 1:size(data_acceleration_p,2)
        if ~isempty([time_of_event{kFiles,iEvents}]);
            tmp = cat(1,data_acceleration_f1{kFiles,iEvents,:});
            horizontal_tmp = tmp(:,1);
            plot([time_data_large{kFiles,iEvents,:}]-[time_data_large{kFiles,iEvents,1}],horizontal_tmp);
            t1 = max([time_of_event{kFiles,iEvents}]-[time_data_large{kFiles,iEvents,1}]-1,1);
            plot(t1,horizontal_tmp(t1+1),'r*','MarkerSize',10); % time of the event
            
        end
    end
end
xlabel('time from the first event [frame]')
ylabel('horizontal acceleration fluid_x [m/s^2]');

figure, hold on,
for kFiles = 1:size(data_acceleration_p,1)
    for iEvents = 1:size(data_acceleration_p,2)
        if ~isempty([time_of_event{kFiles,iEvents}]);
            tmp = cat(1,data_acceleration_f1{kFiles,iEvents,:});
            horizontal_tmp = tmp(:,3);
            plot([time_data_large{kFiles,iEvents,:}]-[time_data_large{kFiles,iEvents,1}],horizontal_tmp);
            t1 = max([time_of_event{kFiles,iEvents}]-[time_data_large{kFiles,iEvents,1}]-1,1);
            plot(t1,horizontal_tmp(t1+1),'r*','MarkerSize',10); % time of the event
            
        end
    end
end
xlabel('time from the first event [frame]')
ylabel('horizontal acceleration fluid_z [m/s^2]');



