clc, clear all

dist = @(x,y) sqrt((x(1,:)-y(1)).^2 + (x(2,:)-y(2)).^2 + (x(3,:)-y(3)).^2);

if ismac
    matdirectory = '/Users/alex/Dropbox/resuspension/2011/trajectories_2012(186-194)';
else
    matdirectory ='C:\Users\dell\Dropbox\resuspension\Matlab\trajectories_alex_2012(186-194)\events_110';
end

large_ones = dir(fullfile(matdirectory,'large*'));

time_of_event = {};
U_f = {};
V_f = {};
W_f = {};
t_kinetic_energy={};



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
           
            
            counter = 0;
           
            u_f=0;
            v_f=0;
            w_f=0;
            
%             u_t1{kFiles,iEvents,ind}=0;
%             v_t1{kFiles,iEvents,ind}=0;
%             w_t1{kFiles,iEvents,ind}=0;
%             
            t_kinetic_energy{kFiles,iEvents,ind} = 0;
            
         
            
            
            for j = 1:numSmallTraj % for all small ones
                
                distance = dist([small(j).xf,small(j).yf,small(j).zf]',[xp,yp,zp])';
                
                relevant = find(small(j).t == data.t(ind)&  distance <= R); 
                                 
                               
                if ~isempty(relevant)  
                    
                   
                    tmp = small(j);
                    tmp.yf = tmp.yf(relevant);
                    
                    if tmp.yf <= -21 && abs(tmp.yf-yp )>= 10,continue, end
                   
                    
                    tmp.uf = tmp.uf(relevant); %[mm/s]
                    tmp.vf = tmp.vf(relevant);
                    tmp.wf = tmp.wf(relevant);
     
                   
                   u_f= u_f+ tmp.uf;
                   v_f= v_f+ tmp.vf;
                   w_f= w_f+ tmp.wf;
                  
                    counter = counter + 1;
                    
                   
                end % relevant
            end % smalltraj
            
                   U_f{kFiles,iEvents,ind} = u_f / counter;
                   V_f{kFiles,iEvents,ind} = v_f / counter;
                   W_f{kFiles,iEvents,ind} = w_f / counter;
                   avg_kinetic_energy{kFiles,iEvents,ind} = 1/2 *(U_f{kFiles,iEvents,ind}^2+ V_f{kFiles,iEvents,ind}^2 + W_f{kFiles,iEvents,ind}^2);
                     
            counter = 0;
                   for j = 1:numSmallTraj % for all small ones
                
                distance = dist([small(j).xf,small(j).yf,small(j).zf]',[xp,yp,zp])';
                
                relevant = find(small(j).t == data.t(ind)&  distance <= R); 
                                 
                               
                if ~isempty(relevant)  
                    
                   
                    tmp = small(j);
                    tmp.yf = tmp.yf(relevant);
                    
                    if tmp.yf <= -21 ,continue, end
                   
                    
                    tmp.uf = tmp.uf(relevant); %[mm/s]
                    tmp.vf = tmp.vf(relevant);
                    tmp.wf = tmp.wf(relevant);
                    
                   u_t = tmp.uf - U_f{kFiles,iEvents,ind};
                   v_t = tmp.vf - V_f{kFiles,iEvents,ind}; 
                   w_t = tmp.wf - W_f{kFiles,iEvents,ind};
                   
                   t_kinetic_energy{kFiles,iEvents,ind} = t_kinetic_energy{kFiles,iEvents,ind} + 1/2 *(u_t^2+ v_t^2 + w_t^2);
                   
%                    u_t1{kFiles,iEvents,ind}= u_t1{kFiles,iEvents,ind}+ u_t; 
%                    v_t1{kFiles,iEvents,ind}= v_t1{kFiles,iEvents,ind}+ v_t; 
%                    w_t1{kFiles,iEvents,ind}= w_t1{kFiles,iEvents,ind}+ w_t; 
%                    
                   
                  counter = counter + 1;
                    
                   
               end % relevant
                  end % smalltraj
                  
            %storge
          
             
            data_small{kFiles,iEvents,ind} = t_kinetic_energy{kFiles,iEvents,ind} / counter /1000^2; % average the small particles kinetic energy [m/s]
            data_small_1{kFiles,iEvents,ind} = avg_kinetic_energy{kFiles,iEvents,ind} /1000^2; % 
            %u_t2{kFiles,iEvents,ind}= u_t1{kFiles,iEvents,ind}/counter/1000;
            %v_t2{kFiles,iEvents,ind}= v_t1{kFiles,iEvents,ind}/counter/1000;
            
            time_data_small{kFiles, iEvents,ind} = data.t(ind);
           
           
            % i= number of events
            %ind= time of event
            
        end %ind
        
        
    end %iEvents
end %kFiles

save turbulent_kinetic_energy_110_25_10.mat


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

figure, hold on,
for kFiles = 1:size(data_small_1,1)
    for iEvents = 1:size(data_small_1,2)
        if ~isempty([time_of_event{kFiles,iEvents}])
            plot([time_data_small{kFiles,iEvents,:}]-[time_data_small{kFiles,iEvents,1}],[data_small_1{kFiles,iEvents,:}]); % x = time, y= kinetic energy small
            tmp = [data_small_1{kFiles,iEvents,:}];
            t1 = max([time_of_event{kFiles,iEvents}]-[time_data_small{kFiles,iEvents,1}]-1,1);
            plot(t1,tmp(t1+1),'ro','MarkerSize',10); % time of the event
        end
    end
end
hold off


% figure, hold on,
% for kFiles = 1:size(data_small_1,1)
%     for iEvents = 1:size(data_small_1,2)
%         if ~isempty([time_of_event{kFiles,iEvents}])  
%             tmp = [u_t2{kFiles,iEvents,:}];
%             tmp1 = [v_t2{kFiles,iEvents,:}];
%             t1 = max([time_of_event{kFiles,iEvents}]-[time_data_small{kFiles,iEvents,1}]-1,1);
%             plot(tmp(t1+1),tmp1(t1+1),'ro','MarkerSize',10); % time of the event
%             %scatter(tmp1,tmp2);
%         end
%     end
% end
% hold off


% 
% 
% figure, hold on,
% for kFiles = 1:size(data_small_1,1)
%     for iEvents = 1:size(data_small_1,2)
%         if ~isempty([time_of_event{kFiles,iEvents}])  
%             tmp = [uv{kFiles,iEvents,:}];
%             t1 = max([time_of_event{kFiles,iEvents}]-[time_data_small{kFiles,iEvents,1}]-1,1);
%             plot(t1,tmp(t1+1),'ro','MarkerSize',10); % time of the event
%               
%         end
%     end
% end
% hold off
% 
% 
% 




