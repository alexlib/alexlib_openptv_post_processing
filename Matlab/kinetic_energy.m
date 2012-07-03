

if ismac
    matdirectory = '/Users/alex/Dropbox/resuspension/2011/trajectories(186-194)';
else
    matdirectory ='C:\Users\hadar\Dropbox\resuspension\2011\trajectories(186-194)';
end
large_ones = dir(fullfile(matdirectory,'large*'));

for k = 1:3 %length(large_ones)
data_large = 0;
data_small={};
filename = large_ones(k).name;
[path,filename,ext] = fileparts(filename);

time_interval = 4; 
R = 50;


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
dist = @(x,y) sqrt((x(1,:)-y(1)).^2 + (x(2,:)-y(2)).^2 + (x(3,:)-y(3)).^2);





for i = 1:5%numEvents
    % Select the trajectory and the time of event
    time_of_event = events.events(i,2);
    trajind = trajIds == events.events(i,1);
    before_event = (time_of_event-time_interval):time_of_event;
    data = large(trajind);
    
    %ind = find(data.t == time_of_event);
    %ind = find(data.t == (time_of_event- time_interval));
    for ind = 1:5 % length(data.t)
       
        xp = data.xf(ind);
        yp = data.yf(ind);
        zp = data.zf(ind);
        up = data.uf(ind);
        vp = data.vf(ind);
        wp = data.wf(ind);
        kinetic_energy_p = 1/3* (up^2 + vp^2 + wp^2 );
        data_large(i) = kinetic_energy_p;
          
        counter=0;
        tmp_kinetic_energy = 0;
        for j = 1:numSmallTraj % for all small ones
            
            distance = dist([small(j).xf,small(j).yf,small(j).zf]',[xp,yp,zp])';
           
            relevant = find(small(j).t == data.t(ind) &  distance <= R);
            
            if relevant
                
                tmp = small(j);
                tmp.uf = tmp.uf(relevant);
                tmp.vf = tmp.vf(relevant);
                tmp.wf = tmp.wf(relevant);
              
                tmp_kinetic_energy = tmp_kinetic_energy +1/3 *(tmp.uf^2+ tmp.vf^2 + tmp.wf^2);
                
              counter=counter+1;
              
               % if verbose,plot_long_trajectories(tmp(relevant),1,hf);   end
                
            end % relevent
        end % smalltraj
        
           data_small{ind,i,k} =  tmp_kinetic_energy / counter;
            % i= number of events
            %ind= time of event
            
        end
    end
end


figure, hold on, for i = 1:size(data_small,3), for j = 1:size(data_small,2), plot([data_small{:,j,i}]); %plot(events.events(j,2),1e-3,'*');
    end, 
end










