clc, clear all

dist = @(x,y) sqrt((x(1,:)-y(1)).^2 + (x(2,:)-y(2)).^2 + (x(3,:)-y(3)).^2);

if ismac
    matdirectory = '/Users/alex/Dropbox/resuspension/2011/trajectories(186-194)';
else
    matdirectory ='C:\Users\hadar\Dropbox\resuspension\Matlab\trajectories_alex_2012(186-194)\events_110';
end

large_ones = dir(fullfile(matdirectory,'large*'));

time_of_event = {};
data_velocity_p_f = {};
basset_force1= {};
time_data_small = {};

R = 25;


D = 550 * 10^(-6); %diameter of the solid particle [m]
r = D/2;
A= pi * r^2 ;% the projected grain area perpendicular to flow direction
V = 4/3 * pi * r^3; %[m^3]
k_v=1E-6;

density_l = 0.998*1000 ; % density of water ,temperature= 22 [gr/cm^3]
density_s = 1.065*1000 ; %density of solid patricles [gr/cm^3]
g = 9.81;


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
        time_of_event{kFiles,iEvents} = events.events(iEvents,2);
        trajind = trajIds == events.events(iEvents,1);
        
        % Take the single trajectory in which we had an event aside into
        % the temprorary 'data'
        
        data = large(trajind);
        
        for i = 1:length(data.t)
            basset_force1{kFiles,iEvents,i} = zeros(1,3);
            %velocity_p_f1{kFiles,iEvents,i} = zeros(1,3);
            data_velocity_p_f{kFiles,iEvents,i} = zeros(1,3);
        end
        
        % Along that trajectory, at each time step
    
        
        for ind = 1:length(data.t)
    
            xp = data.xf(ind);
            yp = data.yf(ind);
            zp = data.zf(ind);
            
            up = data.uf(ind)/1000;
            vp = data.vf(ind)/1000;
            wp = data.wf(ind)/1000;
            
            velocity_p = [up,vp,wp];
            velocity_p_f1=zeros(1,3);
            counter = 0;
            
            for j = 1:numSmallTraj % for all small ones
                
                distance = dist([small(j).xf,small(j).yf,small(j).zf]',[xp,yp,zp])';
                
                relevant = find(small(j).t == data.t(ind) &  distance <= R);
                
                if ~isempty(relevant)
                    
                    tmp = small(j);
                    
                    if tmp.yf <= -21 ,continue, end
                    
                    tmp.uf = tmp.uf(relevant)/1000;
                    tmp.vf = tmp.vf(relevant)/1000;
                    tmp.wf = tmp.wf(relevant)/1000;
                    
                    velocity_f=[tmp.uf,tmp.vf,tmp.wf];
                    
                    
                    velocity_p_f =  [up-tmp.uf,vp-tmp.vf,wp-tmp.wf];
                    velocity_p_f1= velocity_p_f1 + velocity_p_f;
                    
                    counter = counter + 1;
                    
                    
                end % relevant
            end % smalltraj
            
           if length(velocity_p_f1) ~= 3 || any(isnan(velocity_p_f1)), length(velocity_p_f1), end
           data_velocity_p_f{kFiles,iEvents,ind} = velocity_p_f1/counter; 
           
           time_data_small{kFiles, iEvents,ind} = data.t(ind); 
  
        end %ind
        
        
        
        dt = 1/160;
        tau = 6:-1:1;
        
        
        % basset force loop
        for ind = 6:length(data.t)
            
            urel = [data_velocity_p_f{kFiles,iEvents,ind-5:ind}]; % 6 x [U,V,W]
            
            for k = 1:3
                % dt in seconds, urel in m/sec
                du = gradient(urel(k:3:end),dt);
                
                y = du./ sqrt(pi * k_v * tau * dt);
                
                tmp1 =  6*pi * r^2* k_v * density_l * dt * trapz(y);
                
                basset_force1{kFiles,iEvents,ind}(k) = tmp1;
                
                
            end
              
        end % basset
        
    end %iEvents
end %Kfiles


save basset_force_110
% 
figure, hold on,
for kFiles = 1:size(basset_force1,1)
    for iEvents = 1:size(basset_force1,2)
        if ~isempty([time_of_event{kFiles,iEvents}]);
            tmp2 = cat(1,basset_force1{kFiles,iEvents,:});
            norm_basset_force1 = sum(tmp2.^2,2).^(0.5);
            plot([time_data_small{kFiles,iEvents,1:end}]-[time_data_small{kFiles,iEvents,1}],norm_basset_force1);
            t1 = max([time_of_event{kFiles,iEvents}]-[time_data_small{kFiles,iEvents,1}]-1,1);
            plot(t1,norm_basset_force1(t1+1),'r*','MarkerSize',10); % time of the event
            
        end
    end
end


