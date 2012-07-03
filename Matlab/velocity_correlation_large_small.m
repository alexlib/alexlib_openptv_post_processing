function collected_data = velocity_correlation_large_small(filename)
% Usage:
% temp = velocity_correlation_large_small('large1933608');
% plot (temp(:,1), temp(:,2)./temp(1,2),'o')

if ismac
    matdirectory = '/Users/alex/Dropbox/resuspension/2011/trajectories(186-194)';
else
    matdirectory ='C:\Users\hadar\Dropbox\resuspension\2011\trajectories(186-194)';
end

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
numTraj = length(small);

trajIds = unique(cat(1,large.trajid));

%time_interval = 0; 

collected_data = zeros(10000,4);
counter = 0;


for i = 1:numEvents
    % Select the trajectory and the time of event
    time_of_event = events.events(i,2);
    trajind = trajIds == events.events(i,1);
    
    data = large(trajind);
    
  ind = find(data.t == time_of_event);
  
   %ind = find(data.t == (time_of_event- time_interval));
   
  if ind
    
for j = 1:numTraj
    [tf,loc] = ismember(time_of_event,small(j).t );
   
    if tf
        % small particle at the time of event
        xf = small(j).xf(loc);
        yf = small(j).yf(loc);
        zf = small(j).zf(loc);
        uf = small(j).uf(loc);
        vf = small(j).vf(loc);
        wf = small(j).wf(loc);
        if yf <= -21.5 ,continue, end
        
        velocity_f = [uf,vf,wf]; %  full velocity of the small particles in the time of the event
        
        % large particle at the time of event:
        xp = data.xf(ind);
        yp = data.yf(ind);
        zp = data.zf(ind);
        up = data.uf(ind);
        vp = data.vf(ind);
        wp = data.wf(ind);
        velocity_p = [up,vp,wp]; % full velocity of the large particle in the time of the event.
        %distance vector
        distance_r = [xf-xp,yf-yp,zf-zp];
        
        distance_r1 = distance_r/norm(distance_r);
        % velocity in direction of distance vector
        v_f11 = velocity_f * distance_r1'; % projection of small particle velocity on distance vector
        v_p11 = velocity_p * distance_r1'; % projection of large particle velocity on distance vector
        correlation_11 = v_f11 * v_p11;
        
        % correlation in y direction 
        correlation_22 = wf * wp ;
          
        % correlation in x direction 
        correlation_33 = up * uf ;
        
        
        counter = counter + 1;
        collected_data(counter,1) = norm(distance_r);
        collected_data(counter,2) = correlation_11;
        collected_data(counter,3) = correlation_22;
        collected_data(counter,4) = correlation_33;
    end   
end
  end
end

collected_data = collected_data(1:counter,:);

