function run_plot_particle_trajectory
% run_.... is the outer loop for any function ...
% that performs the same operation on all the data
% collected from 3D-PTV in the resuspension experiment
%
if ismac
    matdirectory = '/Users/alex/Dropbox/resuspension/2011/trajectories_2012(186-194)';    
else
    matdirectory ='C:\Users\hadar\Dropbox\resuspension\2011\trajectories_2012(186-194)';
end

directories = {'events_80','events_90','events_100','events_110'};

styles = {'ro','bs','gd','m>'};

hf = figure;
hold on;


for k = 1:4
    
    large_ones = dir(fullfile(matdirectory,directories{k},'large*'));
    
    
    for i = 1:length(large_ones)
        disp(large_ones(i).name) % to visualize the run
        plot_particle_trajectory(fullfile(matdirectory,large_ones(i).name),hf,styles{k});
    end
end
hold off
end




function plot_particle_trajectory(fullname,hf,styles)
% plots large particles trajectories, shifting the event point to
% the 0,0,0 position
% the function should visualize the "coherent" behaviour of the large
% particles
% Inputs:
%       filename =  full (incl. the file path) name of the large
% particles dataset
%       hf = figure handle to which we add the plot


[filepath,filename,ext] = fileparts(fullname);

% Load the largeXXXX.mat
large = load(fullfile(filepath,filename));
% and rename it to 'large'
large = large.(filename);

% Load the eventsXXXXX.mat
events = load(fullfile(filepath,strrep(filename,'large','events')));
numEvents = size(events.events,1);

trajIds = unique(cat(1,large.trajid));


for i = 1:numEvents
    % Select the trajectory and the time of event
    time_of_event = events.events(i,2);
    trajind = trajIds == events.events(i,1);
    
    data = large(trajind);
    ind = find(data.t == time_of_event,1); % ind = index of the event
    
    if ~isempty(ind)
        
        % Once when all "detach" at the same point:
                x = data.xf - data.xf(ind);
                y = data.yf - data.yf(ind);
                z = data.zf - data.zf(ind);
        
        % Once without shifts
        
%                 x = data.xf;
%                 y = data.yf;
%                 z = data.zf;
        
        
        % Once when all "start" at the same point:
%         x = data.xf - data.xf(1);
%         y = data.yf - data.yf(1);
%         z = data.zf - data.zf(1);
        
        figure(hf);
        
        plot3(x,z,y,styles);
        %         ...
        %             'MarkerFaceColor','k','Marker','o',...
        %             'MarkerEdgeColor','r',...
        %             'MarkerSize',4,...
        %             'LineStyle',':');
        plot3(x(1),z(1),y(1),...
            'MarkerFaceColor','c','Marker','s',...
            'MarkerEdgeColor','k',...
            'MarkerSize',4);
        
    end
    
end
end



