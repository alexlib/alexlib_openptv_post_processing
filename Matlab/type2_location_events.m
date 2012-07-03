function [collected_data] = type2_location_events(filename,quantity1,quantity2)

% Copyright (c) 2011, Turbulence Structure Laboratory

collected_data = [];

if ismac
    matdirectory = '/Users/alex/Dropbox/resuspension/2011/trajectories(186-194)';
else
    matdirectory ='C:\Users\hadar\Dropbox\resuspension\2011\trajectories(186-194)';
end
% large_ones = dir(fullfile(matdirectory,'large*'));

[path,filename,ext] = fileparts(filename);

% Load the largeXXXX.mat
large = load(fullfile(matdirectory,filename));
% and rename it to 'large'
large = large.(filename);

% Load the eventsXXXXX.mat
events = load(fullfile(matdirectory,strrep(filename,'large','events')));

numEvents = size(events.events,1);

trajIds = unique(cat(1,large.trajid));


for i = 1:numEvents
    % Select the trajectory and the time of event
    time_of_event = events.events(i,2);
    trajind = trajIds == events.events(i,1);
    
    data = large(trajind);
    ind = find(data.t == time_of_event);
         
            % position of the large particle at the time of event:
            xp = data.xf(ind);
            zp = data.zf(ind);
            
   collected_data = cat(1,collected_data,[data.(quantity1)(ind),data.(quantity2)(ind)]);
end
       
