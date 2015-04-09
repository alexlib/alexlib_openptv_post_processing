clear all


if ismac
    matdirectory = '/Users/alex/Dropbox/resuspension/2011/trajectories(186-194)';
else
    matdirectory ='C:\Users\hadar\Dropbox\resuspension\2011\trajectories(186-194)/events_110';
end

large_ones = dir(fullfile(matdirectory,'large*'));

%time_of_event = ();
%data_large1 = {};
%data_large2 = {};
%time_data= {};

counter1 = 0;
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
    
    numEvents = size(events.events,1);    
    trajIds = unique(cat(1,large.trajid));
    
    for iEvents = 1:numEvents
        % Select the trajectory and the time of event
        time_of_event{kFiles,iEvents} = events.events(iEvents,2);
        trajind = trajIds == events.events(iEvents,1);
        % Take the single trajectory in which we had an event aside into
        % the temprorary 'data'
        
        data = large(trajind);
        event = find(data.t == time_of_event{kFiles,iEvents}); 
        first = find(data.t(1));
        
        figure1 = figure; 
        color = 'b'; axes('Parent',figure1); box('on');grid('on'); hold('all');
      
        xlabel('$x$','Interpreter','Latex');ylabel('$z$','Interpreter','Latex');

   for i = 1:length(data)
        plot(data.xf/295,data.zf/295,...
            'MarkerFaceColor','w','Marker','o',...
            'MarkerSize',3); 
         plot(data.xf(first)/295,data.zf(first)/295,...
            'MarkerFaceColor','r','Marker','*',...
            'MarkerSize',10);
         plot(data.xf(event)/295,data.zf(event)/295,...
            'MarkerFaceColor','r','Marker','*',...
            'MarkerSize',10);
        
   end
     
         counter1= counter1+1;
    end   
end


