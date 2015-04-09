clc, clear all

dist = @(x,y) sqrt((x(1,:)-y(1)).^2 + (x(2,:)-y(2)).^2 + (x(3,:)-y(3)).^2);

if ismac
    matdirectory = '/Users/alex/Dropbox/resuspension/2011/trajectories(186-194)';
else
    matdirectory ='C:\Users\hadar\Desktop\trajectories_alex_2012(186-194)';
end

large_ones = dir(fullfile(matdirectory,'large*'));

 data = {} ;
 
%for all the large particles
for kFiles = 1:length(large_ones)
    
    % get the name of the data
    filename = large_ones(kFiles).name;
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
    numSmallTraj = length(small);
    
    trajIds = unique(cat(1,large.trajid));
    
    
    
    for iEvents = 1:numEvents
        % Select the trajectory and the time of event
       trajind = trajIds == events.events(iEvents,1);
       data = large(trajind);
       
       %time_of_event{kFiles,iEvents} = events.events(iEvents,2); 
        
       %ind = find(data.t == time_of_event);
        
       for j = 1:numSmallTraj
           
       %loc = find (small(j).t==time_of_event);

            %yf= small(j). yf(loc);
            %uf = small(j).uf(loc);
            %vf = small(j).vf(loc);
            %wf = small(j).wf(loc);
            
            yf= small(j). yf;
            uf = small(j).uf;
            vf = small(j).vf;
            wf = small(j).wf;
            
            
           if (yf >= -21.5 && yf <= -6.5)
            collected_data = cat(1,collected_data,[uf,vf,wf]);
           end 
          
       end 
        
    end
     
     data = cat(1,data,collected_data);
end

% hf1=figure;
% subplot(311), hold on
% nhist(data(:,1),41),hold on
% title(sprintf('small quantity: %s',quantity1))
% 
% hold on
% subplot(312),
% nhist(data(:,2),41)
% title(sprintf('small quantity: %s',quantity2))
% 
% hold on
% subplot(313),
% nhist(data(:,3),41)
% title(sprintf('small quantity: %s',quantity3))
% 
