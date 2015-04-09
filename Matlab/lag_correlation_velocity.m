clc, clear all

dist = @(x,y) sqrt((x(1,:)-y(1,:)).^2 + (x(2,:)-y(2,:)).^2 + (x(3,:)-y(3,:)).^2);

if ismac
    matdirectory = '/Users/alex/Dropbox/resuspension/2011/trajectories(186-194)';
else
    matdirectory ='C:\Users\hadar\Dropbox\resuspension\2011\trajectories(186-194)\events_80';
end

large_ones = dir(fullfile(matdirectory,'large*'));

collected_data = {};

R = 25;

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
        
        % Take the single trajectory in which we had an event aside into
        % the temprorary 'data'
        
        data = large(trajind);
        
       time_of_event = events.events(iEvents,2);
       event = find(data.t == time_of_event); 
       first = find(data.t(1));
       
            
            ind_large= first:event;
            for j = 1:numSmallTraj % for all small ones
             % for  delta_t = 1:1:100 ;
                  
               % if (event + delta_t) > length(data.t), continue,end % needs to be fixed
                 
                %ind_small = find (small(j).t <= data.t (event+delta_t) & ...
                %small(j).t >= data.t (first+delta_t));
                ind_small =  small(j).t (ind_large);
                
                if isempty(ind_small), continue, end  % if there are no small particles in the correct time
                if length(ind_small) ~= length(ind_large),continue, end % needs to be fixed
               
              
                %ind_large_shift = find (data.t <= data.t (event+delta_t) & ...
                %data.t >= data.t (first+delta_t)); % for the calculation of distance, shift the large particle to the same time as the small one
           
               distance = dist([small(j).xf(ind_small),small(j).yf(ind_small),small(j).zf(ind_small)]'...
               ,[data.xf(ind_large),data.yf(ind_large),data.zf(ind_large)]');
           
        
               relevant =  distance <= R & (small(j).yf(ind_small))' >= -21.5; % time and location
               
          
               time_rev_large = ind_large;% (relevant);
               time_rev_small = ind_small (relevant); 
               
               up = data.uf(time_rev_large);
               vp = data.vf(time_rev_large);
               wp = data.wf(time_rev_large);
               
               tmp.uf = small(j).uf(time_rev_small);
               tmp.vf = small(j).vf(time_rev_small);
               tmp.wf = small(j).wf(time_rev_small); 
               
                        acorrelation_uu = xcorr(up,tmp.uf,50); 
                        acorrelation_vv = xcorr(vp,tmp.vf,50);
                        acorrelation_ww = xcorr(wp,tmp.wf,50);
                        %acorrelation_vp = xcorr(velocity_p,velocity_f,50);
                   
             collected_data{kFiles,iEvents,j,delta_t, 1}= cat(1,acorrelation_uu);
             collected_data{kFiles,iEvents,j,delta_t, 2}= cat(1,acorrelation_vv);
             collected_data{kFiles,iEvents,j,delta_t, 3}= cat(1,acorrelation_ww);
             %collected_data{kFiles,iEvents,j,delta_t, 4}= cat(1,acorrelation_vp);
                
                     
             % end %delta_t        
            end %numsmalltraj
                        
           
    end % iEvents
    
end % kfiles


%temp  = sortrows(collected_data{kFiles,iEvents,j,delta_t, 1}, 1); % organizes tha data from small to large

save lag_correlations.mat

%[bin_centers,average_in_bins,n] = plot_average_in_bins(temp,0:5:100);

 %figure
 %plot(bin_centers,average_in_bins(1:end-1)/max(average_in_bins(1:end-1)),'o')

%figure
%bar(bin_centers,n(1:end-1))

