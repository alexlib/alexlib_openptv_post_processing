clc, clear all

dist = @(x,y) sqrt((x(1,:)-y(1)).^2 + (x(2,:)-y(2)).^2 + (x(3,:)-y(3)).^2);

if ismac
    matdirectory = '/Users/alex/Dropbox/resuspension/2011/trajectories_2012(186-194)';
else
    matdirectory ='C:\Users\dell\Dropbox\resuspension\Matlab\trajectories_alex_2012(186-194)\events_100';
end

large_ones = dir(fullfile(matdirectory,'large*'));

time_of_event = {};
U_f = {};
V_f = {};
W_f = {};
U_W={};

R = 25;
Ws=10.8;

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
            uf_wf1=0;
            
             u_t1{kFiles,iEvents,ind}=0;
             v_t1{kFiles,iEvents,ind}=0;
             w_t1{kFiles,iEvents,ind}=0;
             ut_wt1{kFiles,iEvents,ind}=0; 
             v_t2{kFiles,iEvents,ind}=0;
            
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
                    uf_wf = [tmp.uf,tmp.wf];
                    uf_wf = norm(uf_wf);
                    
                    
                   u_f= u_f+ tmp.uf;
                   v_f= v_f+ tmp.vf;
                   w_f= w_f+ tmp.wf;
                   
                   uf_wf1= uf_wf1+ uf_wf;
                 
                   counter = counter + 1;
                    
                   
                end % relevant
            end % smalltraj
            
                   U_f{kFiles,iEvents,ind} = u_f / counter;
                   V_f{kFiles,iEvents,ind} = v_f / counter;
                   W_f{kFiles,iEvents,ind} = w_f / counter;
                   U_W{kFiles,iEvents,ind} = uf_wf1/ counter; 
                     
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
                   ut_wt = [u_t,w_t];
                   ut_wt = norm(ut_wt);               
                   
                   u_t1{kFiles,iEvents,ind}= u_t1{kFiles,iEvents,ind}+ u_t; 
                   v_t1{kFiles,iEvents,ind}= v_t1{kFiles,iEvents,ind}+ v_t; 
                   w_t1{kFiles,iEvents,ind}= w_t1{kFiles,iEvents,ind}+ w_t; 
                   ut_wt1{kFiles,iEvents,ind}= ut_wt1{kFiles,iEvents,ind}+ ut_wt;
                     
                   v_t2{kFiles,iEvents,ind}= sqrt(v_t1{kFiles,iEvents,ind}^2);
                  counter = counter + 1;
                    
                   
               end % relevant
                  end % smalltraj
                  
            %storge
          
             
            data_small{kFiles,iEvents,ind} = v_t1{kFiles,iEvents,ind} / counter /Ws; % average the small particles kinetic energy [mm/s]
            data_small_1{kFiles,iEvents,ind} = V_f{kFiles,iEvents,ind}/Ws ;
            data_small_2{kFiles,iEvents,ind} = v_t2{kFiles,iEvents,ind} / counter /Ws;
            
            data_small_U_W{kFiles,iEvents,ind} = U_W{kFiles,iEvents,ind}/Ws; 
            u_t2{kFiles,iEvents,ind}= u_t1{kFiles,iEvents,ind}/counter/Ws;
            v_t2{kFiles,iEvents,ind}= w_t1{kFiles,iEvents,ind}/counter/Ws;
            Ut_Wt{kFiles,iEvents,ind} = ut_wt1{kFiles,iEvents,ind}/ counter/Ws; 
            
            time_data_small{kFiles, iEvents,ind} = data.t(ind);
           
           
            % i= number of events
            %ind= time of event
            
        end %ind
        
        
    end %iEvents
end %kFiles

save vertical_velocity_100_25_10.mat


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


figure, hold on,
for kFiles = 1:size(data_small_U_W,1)
    for iEvents = 1:size(data_small_U_W,2)
        if ~isempty([time_of_event{kFiles,iEvents}])
            plot([time_data_small{kFiles,iEvents,:}]-[time_data_small{kFiles,iEvents,1}],[data_small_U_W{kFiles,iEvents,:}]); % x = time, y= kinetic energy small
            tmp = [data_small_U_W{kFiles,iEvents,:}];
            t1 = max([time_of_event{kFiles,iEvents}]-[time_data_small{kFiles,iEvents,1}]-1,1);
            plot(t1,tmp(t1+1),'ro','MarkerSize',10); % time of the event
        end
    end
end
hold off

figure, hold on,
for kFiles = 1:size(Ut_Wt,1)
    for iEvents = 1:size(Ut_Wt,2)
        if ~isempty([time_of_event{kFiles,iEvents}])
            plot([time_data_small{kFiles,iEvents,:}]-[time_data_small{kFiles,iEvents,1}],[Ut_Wt{kFiles,iEvents,:}]); % x = time, y= kinetic energy small
            tmp = [Ut_Wt{kFiles,iEvents,:}];
            t1 = max([time_of_event{kFiles,iEvents}]-[time_data_small{kFiles,iEvents,1}]-1,1);
            plot(t1,tmp(t1+1),'ro','MarkerSize',10); % time of the event
        end
    end
end
hold off

figure, hold on,
for kFiles = 1:size(data_small_2,1)
    for iEvents = 1:size(data_small_2,2)
        if ~isempty([time_of_event{kFiles,iEvents}])
            plot([time_data_small{kFiles,iEvents,:}]-[time_data_small{kFiles,iEvents,1}],[data_small_2{kFiles,iEvents,:}]); % x = time, y= kinetic energy small
            tmp = [data_small_2{kFiles,iEvents,:}];
            t1 = max([time_of_event{kFiles,iEvents}]-[time_data_small{kFiles,iEvents,1}]-1,1);
            plot(t1,tmp(t1+1),'ro','MarkerSize',10); % time of the event
        end
    end
end
hold off

