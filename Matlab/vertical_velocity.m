clc, clear all

dist = @(x,y) sqrt((x(1,:)-y(1)).^2 + (x(2,:)-y(2)).^2 + (x(3,:)-y(3)).^2);

if ismac
    matdirectory = '/Users/alex/Dropbox/resuspension/2011/trajectories_2012(186-194)';
else
    matdirectory ='C:\Users\dell\Dropbox\resuspension\Matlab\trajectories_alex_2012(186-194)\events_100';
end

large_ones = dir(fullfile(matdirectory,'large*'));

time_of_event = {};

data_small_V = {};
data_small_U_W = {};

time_data_small = {};

Ws=10.8;
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
           
           
            
            counter = 0;
            V_f=0;
            Uf_Wf1=0;
            
            
            for j = 1:numSmallTraj % for all small ones
                
                distance = dist([small(j).xf,small(j).yf,small(j).zf]',[xp,yp,zp])';
                
                relevant = find(small(j).t == data.t(ind)&  distance <= R); %  & small(j).yf <= -5 );
                                 
                               
                if ~isempty(relevant)  
                    
                  
                    tmp = small(j);
                    tmp.yf = tmp.yf(relevant);
                    
                    if tmp.yf <= -21 ,continue, end
                   
                    tmp.uf = tmp.uf(relevant); %[mm/s]
                    tmp.vf = tmp.vf(relevant);
                    tmp.wf = tmp.wf(relevant);
                    
                    uf_wf = [tmp.uf,tmp.wf];
                    uf_wf = norm(uf_wf);
                    
                    Uf_Wf1= Uf_Wf1+ uf_wf;
                    
                    V_f= V_f+ tmp.vf;
                  
                   counter = counter + 1;
                    
                    % if verbose,plot_long_trajectories(tmp(relevant),1,hf);   end
               
                end % relevant
            end % smalltraj
            
            %storge
            data_small_V{kFiles,iEvents,ind} = V_f / counter/Ws; 
            data_small_U_W{kFiles,iEvents,ind} = Uf_Wf1/ counter/Ws; 
            
            time_data_small{kFiles, iEvents,ind} = data.t(ind);
           
           
        end
        
        % d_tmp_kinetic_energy{kFiles,iEvents,ind} = diff(data_small{kFiles,iEvents});
    end
end

save average_fluid_velocity_100.mat


%%
figure, hold on,
for kFiles = 1:size(data_small_V,1)
    for iEvents = 1:size(data_small_V,2)
        if ~isempty([time_of_event{kFiles,iEvents}])
            %plot([time_data_small{kFiles,iEvents,:}]-[time_data_small{kFiles,iEvents,1}],[data_small{kFiles,iEvents,:}]); % x = time
            tmp = cat(1,data_small_V{kFiles,iEvents,:});
            t1 = max([time_of_event{kFiles,iEvents}]-[time_data_small{kFiles,iEvents,1}]-1,1);
            plot(t1,tmp(t1+1),'ro','MarkerSize',10); % time of the event
        end
    end
end
xlabel('time[frame]')
ylabel('V_f/W_s');



figure, hold on,
for kFiles = 1:size(data_small_V,1)
    for iEvents = 1:size(data_small_V,2)
        if ~isempty([time_of_event{kFiles,iEvents}])
            %plot([time_data_small{kFiles,iEvents,:}]-[time_data_small{kFiles,iEvents,1}],[data_small{kFiles,iEvents,:}]); % x = time
            tmp = cat(1,data_small_U_W{kFiles,iEvents,:});
            t1 = max([time_of_event{kFiles,iEvents}]-[time_data_small{kFiles,iEvents,1}]-1,1);
            plot(t1,tmp(t1+1),'ro','MarkerSize',10); % time of the event
        end
    end
end
xlabel('time[frame]')
ylabel('(U_f^2+W_f^2)^0^.^5)');






