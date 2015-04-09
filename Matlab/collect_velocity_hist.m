function [collected_data] = collect_velocity_hist(filename)

collected_data = [];


if ismac
    matdirectory = '/Users/alex/Dropbox/resuspension/2011/trajectories(186-194)';
else
    matdirectory ='C:\Users\hadar\Dropbox\resuspension\Matlab\trajectories_alex_2012(186-194)';
end
% large_ones = dir(fullfile(matdirectory,'large*'));

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

%numEvents = size(events.events,1);
numTraj = length(small);


%trajIds = unique(cat(1,large.trajid));


% for i = 1:numEvents
%     % Select the trajectory and the time of event
%     time_of_event = events.events(i,2);
%     trajind = trajIds == events.events(i,1);
%     
%      data = large(trajind);
%     ind = find(data.t == time_of_event);
%     
%       if ind
    for j = 1:numTraj
      % [tf,loc] = ismember(time_of_event,small(j).t);
        
        
        %if tf % i.e. if this small trajectory is also appearing at that time
            
            yf = small(j).yf;
            zf = small(j).zf;
            location = find( yf >= -21 & yf <= -15 );
            
            if ~isempty(location) 
                
            uf= small(j).uf(location); % [mm/s]
            vf= small(j).vf(location);
            wf= small(j).wf(location);
            
%             uf_wf = [uf,wf];
%             uf_wf = norm(uf_wf);
          
            collected_data = cat(1,collected_data,[uf,vf,wf]);
            %collected_data = cat(1,collected_data,[uf_wf,vf]);
            end
            
           % end
        end  
%     end
%       end
% end


