dist = @(x,y) sqrt((x(1,:)-y(1)).^2 + (x(2,:)-y(2)).^2 + (x(3,:)-y(3)).^2);

if ismac
    matdirectory = '/Users/alex/Dropbox/resuspension/2011/trajectories_2012(186-194)';
else
    matdirectory ='C:\Users\hadar\Dropbox\resuspension\Matlab\trajectories_alex_2012(186-194)\events_80';
end

large_ones = dir(fullfile(matdirectory,'large*'));

time_of_event = {};
data_large = {};
data_small = {};


R = 25;

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
    
    % Load the small particle data and rename to 'small'
    small = load(fullfile(matdirectory,strrep(filename,'large','small')));
    small = small.(strrep(filename,'large','small'));
    
    numEvents = size(events.events,1);
    numSmallTraj = length(small);
    
    trajIds = unique(cat(1,large.trajid));
    
    
    
    for iEvents = 1:numEvents
        % Select the trajectory and the time of event
        time_of_event{kFiles,iEvents} = events.events(iEvents,2);
        trajind = trajIds == events.events(iEvents,1);
        
       
        data = large(trajind);
      
        for ind = 1:length(data.t==time_of_event{kFiles,iEvents})
            
            xp = data.xf(ind);
            yp = data.yf(ind);
            zp = data.zf(ind);
            
            axp = data.axf(ind)/1000;
            ayp = data.ayf(ind)/1000;
            azp = data.azf(ind)/1000;
            
            acceleration_p = [axp,ayp,azp];% [m/s^2]
            
            acceleration_f1 = zeros(1,3);
            counter = 0;

            for j = 1:numSmallTraj % for all small ones
                
                distance = dist([small(j).xf,small(j).yf,small(j).zf]',[xp,yp,zp])';
                
                relevant = find(small(j).t == data.t(ind)&  distance <= R); %  & small(j).yf <= -5 );
                                 
                               
                if ~isempty(relevant)  % && ~isempty(small(j).t(ind-1))
                    tmp = small(j);
                    tmp.yf = tmp.yf(relevant);
                    
                    if tmp.yf <= -21.5 && tmp.yf <= -11,continue, end
               
                    
                    axf = tmp.axf(relevant)/1000; 
                    ayf = tmp.ayf(relevant)/1000;
                    azf = tmp.azf(relevant)/1000;
                    
                    acceleration_f = [axf,ayf,azf]; %[m/s^2]
                    acceleration_f1= acceleration_f1+ acceleration_f;
                   
                   counter = counter + 1;
    
                end % relevant
            end % smalltraj
            
            data_acceleration_p{kFiles,iEvents,ind} = acceleration_p;
            data_acceleration_f1{kFiles,iEvents,ind} = acceleration_f1/counter;
  
            time_data_small{kFiles, iEvents,ind} = data.t(ind);
  
            
        end
    end
end


            tmp1 = [data_acceleration_p{:}];
            tmp2 = [data_acceleration_f1{:}];
            
            tmp1a = tmp1(2:3:end);
    tmp2a = tmp2(2:3:end);

figure, hold on, nhist(tmp1a(:),101,'r-'); nhist(tmp2a(:),101,'b--');


figure, hold on, for i = -20:20, c = corrcoef(circshift(tmp1a(:),[-i 0]),tmp2a(:)); plot(i,c(1,2),'o'); end


figure, hold on, plot(tmp2a(:),'b'); plot(tmp1a(:),'k'), plot(circshift(tmp1a(:),[9,0]),'r')

figure, hold on, scatter(circshift(tmp1a(:),[0,0]),tmp2a(:),'r');  scatter(circshift(tmp1a(:),[-4,0]),tmp2a(:),'b');

            
figure
scatter(tmp1,tmp2);
xlabel('particle')
ylabel('tracers')
