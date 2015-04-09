function [collected_data] = collect_velocity_large(fullname,dt)
% collect_velocity_large collects the horizontal and vertical velocities
% of the large particles at the event moment

collected_data = [];

[filepath,filename,ext] = fileparts(fullname);

% Load the largeXXXX.mat
large = load(fullfile(filepath,filename));
% and rename it to 'large'
large = large.(filename);

% Load the eventsXXXXX.mat
events = load(fullfile(filepath,strrep(filename,'large','events')));

% Load the small particle data and rename to 'small'
% small = load(fullfile(matdirectory,strrep(filename,'large','small')));
% small = small.(strrep(filename,'large','small'));

numEvents = size(events.events,1);
% numTraj = length(small);


trajIds = unique(cat(1,large.trajid));


for i = 1:numEvents
    % Select the trajectory and the time of event
    time_of_event = events.events(i,2);
    trajind = trajIds == events.events(i,1);
    
    data = large(trajind);
    ind = find(data.t == time_of_event+dt,1);
    
    if ~isempty(ind)
        
        dax = abs(data.axf(ind)) * cosine(data.axf(ind),data.uf(ind));
        day = abs(data.ayf(ind)) * cosine(data.ayf(ind),data.vf(ind));
        
        
        % collected_data = cat(1,collected_data,[data.axf(ind),data.ayf(ind),data.azf(ind)]);
        collected_data = cat(1,collected_data,[dax,day,data.ayf(ind)]);
        
    end
    
end


