function [hf1] = plot_type2_event(filename,quantity)
% PLOT_TYPE2_EVENT(MATFILENAME,QUANTITY)
% plot_type2_event
% Following the mark_type2_events script
% we can load largeXXXXXX.mat and eventsXXXXXX.mat files
% and plot it in a variety of ways
% note: eventsxxxx only marks the moment of event according to type2
% definitions. here we can load the large particle trajectories, 
%mark the moment of event (with eventsxxxx  files)and plot the desired quantity vs the absolute time of the event.   
% Inputs:
%       FILENAME - 'large1934980' 
%       QUANTITY - e.g. 'yf','uf','axf'
%
% Example:
%   plot_type2_event('large1934980','yf')

if ismac
    matdirectory = '/Users/alex/Dropbox/resuspension/2011/trajectories(186-194)';
else
    matdirectory ='C:\Users\hadar\Dropbox\resuspension\2011\trajectories(186-194)';
end
% large_ones = dir(fullfile(matdirectory,'large*'));

% Load the largeXXXX.mat
large = load(fullfile(matdirectory,filename));
% use only the proper data, in some files there are too many variables
large = large.(filename);
% Load the eventsXXXXX.mat
load(fullfile(matdirectory,strrep(filename,'large','events'))); % replace 'large' with 'events' in the name of the function

numTraj = length(large);

hf1 = figure; hold on
for i = 1:numTraj,
    % If in events:
    if ismember(i,events(:,1))
        ind = events(:,1) == i;
        time_of_event = events(ind,2);% time of event is the second row of ind 
        ind = find(large(i).t == time_of_event);% 
        plot(large(i).t,large(i).(quantity),'LineWidth',1,'DisplayName',sprintf('%d',large(i).trajid(1)));
        plot(large(i).t(ind),large(i).(quantity)(ind),'ro'); %marks the event.
        % plot(large(i).t,large(i).vf,'LineWidth',1,'DisplayName',sprintf('%d',large(i).trajid(1)));
        % plot(large(i).t(ind),large(i).vf(ind),'ro');
    else
        % plot(large(i).t,large(i).vf,'k--');
plot(large(i).t,large(i).(quantity),'k--');
    end
end
hold off
