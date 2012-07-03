% mark type2 events
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
    matdirectory = '/Users/alex/Dropbox/resuspension/2011/trajectories(186-194)';
else
    matdirectory ='C:\Users\hadar\Dropbox\resuspension\2011\trajectories(186-194)';
end
large_ones = dir(fullfile(matdirectory,'large*'));

% if Reynolds stresses are required, one can use:
name_of_the_function ='type2_event'; % Only to mark the events

% Plot the results? 
verbose = 0;


% hf2 = figure; hold on;
for i = 1:length(large_ones)
    % Get the number of the run, it resembles the frames in Streams 5
    num = large_ones(i).name(6:end-4);
    % load the largeXXXXX.mat file
    load(fullfile(matdirectory,large_ones(i).name));
    try
        disp(['large',num])
        % CHANGE THE NAME OF THE FILE THAT IS RUNNING
        events = feval(name_of_the_function,...
            eval(['large',num]),verbose);
       
        save(fullfile(matdirectory,['events',num2str(num)]),'events');
    catch
        continue
    end
    %keyboard
    
end
