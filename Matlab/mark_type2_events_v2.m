% mark type2 events_v2
% Last modified:
% 14-12-2012
%
% we are running the event identification every time we need to get some
% data
% it's just a waste of time. we need to run this script once for the
% directory and mark the events for every data in largeXXXXXX.mat file
%

% load the file
% find the events
% save the largeXXXXXX_events.mat file with the simple information
% in the following table:
% trajId, time_of_event
% trajId, time_of_event

% next every time we need the data, we simply use this markup.
% the markup can change with time, but then events could be saved
% together with the respective type2_event definition function


if ismac
    % matdirectory = '/Users/alex/Dropbox/resuspension/Matlab/trajectories_2012(186-194)';
    matdirectory = '/Users/alex/Desktop/trajectories_2012(186-194)';
else
    matdirectory ='C:\Users\hadar\Desktop\events_110';
end
large_ones = dir(fullfile(matdirectory,'large*'));

% if Reynolds stresses are required, one can use:
name_of_the_function ='type2_event_v2'; % Only to mark the events

% Plot the results?
verbose = 1;


% hf2 = figure; hold on;
for i = 1:length(large_ones)
    % Get the number of the run, it resembles the frames in Streams 5
    num = large_ones(i).name(6:end-4);
    % load the largeXXXXX.mat file
    large = load(fullfile(matdirectory,large_ones(i).name));
    
    tmp = large.(['large',num]);
    
    removelist = 1;
    while ~isempty(removelist)
        [tmp,removelist] = haitao_linking_criteria_v3(tmp,3,5);
        % disp(removelist)
    end
    

       
        % CHANGE THE NAME OF THE FILE THAT IS RUNNING
        events = feval(name_of_the_function,...
            tmp,verbose);
        
        save(fullfile(matdirectory,['events',num2str(num)]),'events');
        
        large.(['large',num]) = tmp;
        
        save(fullfile(matdirectory,large_ones(i).name),'-struct','large');

end
