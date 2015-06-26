% lag_correlation_velocity_org.m is the original file created by Hadar
% Traugott
% see the next version in lag_correlation_velocity.m

clc, clear all

dist = @(x,y) sqrt((x(1,:)-y(1)).^2 + (x(2,:)-y(2)).^2 + (x(3,:)-y(3)).^2);

if ismac
    matdirectory = '/Users/alex/Dropbox/resuspension/2011/trajectories(186-194)';
else
    matdirectory ='C:\Users\hadar\Dropbox\resuspension\2011\trajectories(186-194)\events_80';
end

large_ones = dir(fullfile(matdirectory,'large*'));



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
        data_small = small(trajind); 
        
        % Along that trajectory, at each time step
        
        for ind = 1:length(data.t)
            
            xp = data.xf(ind);
            yp = data.yf(ind);
            zp = data.zf(ind);
            up = data.uf(ind);
            vp = data.vf(ind);
            wp = data.wf(ind);
            
            velocity_p = [up,vp,wp];
            
            
            counter = 0;
            
            for j = 1:numSmallTraj % for all small ones
              
            counter_1 = 0;
                
            for k = 1:length(data_small.t) % for the small trajectory, at each time step
                
                distance = dist([small(k).xf,small(k).yf,small(k).zf]',[xp,yp,zp])';
                
                relevant = find(distance <= R)  ;
                
                if ~isempty(relevant)
                    
                    tmp = small(k);
                    tmp.x = tmp.xf(relevant);
                    tmp.y = tmp.yf(relevant);
                    tmp.z = tmp.zf(relevant);
                    
                   if tmp.y <= -21.5 ,continue, end
                      
                    tmp.uf = tmp.uf(relevant);
                    tmp.vf = tmp.vf(relevant);
                    tmp.wf = tmp.wf(relevant);
                    
                    velocity_f= [tmp.uf,tmp.vf,tmp.wf];
                    
                     
                    %distance_r = [tmp.x-xp,tmp.y-yp,tmp.z-zp];
                    %distance_r1 = distance_r/norm(distance_r);
                    
                    
                end % relevant
                
                correlation_uu  = up * tmp.uf;
                correlation_vv  = vp * tmp.vf;
                correlation_ww  = wp * tmp.wf; 
                    
                counter_1 = counter_1 + 1;
                    
           end %data_small.t
              
               counter = counter + 1;
               
               collected_data{kFiles,iEvents,ind,counter,counter_1, 1}= correlation_uu;
               collected_data{kFiles,iEvents,ind,counter,counter_1, 2}= correlation_vv;
               collected_data{kFiles,iEvents,ind,counter,counter_1, 3}= correlation_ww;
               
            end % smalltraj
     
            
        end %
        
        
    end % iEvents
    
    
end % kfiles

%collected_data = collected_data{kFiles,iEvents,ind,1:counter,1:counter_1, :};

save lag_correlations.mat

  %data_correlation = cat(1,data_correlation,collected_data);
  %temp  = sortrows(data_correlation , 1); % organizes tha data from small to large 
 
 %[bin_centers,average_in_bins,n] = plot_average_in_bins(temp,0:10:150);
 
 % figure
 % plot(bin_centers,average_in_bins(1:end-1)/max(average_in_bins(1:end-1)),'o')
  
  %figure
 % bar(bin_centers,n(1:end-1))

