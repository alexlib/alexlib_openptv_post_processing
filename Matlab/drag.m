clc, clear all

dist = @(x,y) sqrt((x(1,:)-y(1)).^2 + (x(2,:)-y(2)).^2 + (x(3,:)-y(3)).^2);

if ismac
    matdirectory = '/Users/alex/Dropbox/resuspension/2011/trajectories(186-194)';
else
    matdirectory ='C:\Users\hadar\Dropbox\resuspension\2011\trajectories(186-194)\events_110';
end

large_ones = dir(fullfile(matdirectory,'large*'));

time_of_event = {};
%data_large = {};
%data_small = {};
data_drag = {};

R = 25;


D = 550 * 10^(-6); %diameter of the solid particle [m]
r = D/2;
A= pi * r^2 ;% the projected grain area perpendicular to flow direction

density_l = 0.998 ; % density of water ,temperature= 22 [gr/cm^3]
density_s = 1.065 ; %density of solid patricles [gr/cm^3]

CD = 2; % drag coefficient

%u_critical = 0.0655;
%f_drag_critical= 0.5* density_l*1000 * CD * A * (u_critical)^2;

% according to the article: impulse an particle dislodgement undr turbulent flow conditions
ws = (density_s -density_l)* 1000 * pi  * 9.81 *4 / 3 * r^3 ;%submerged particle weight 
fv = (1+0.5*(density_l/(density_s-density_l))); % 
f_drag_critical =  fv * ws  ; 
%u_square_critical = 2*fv*ws/(density_l* CD *A);

%according to shields parameter. article: grid stirred turbulence/Medina,Redondo
%s = density_s/density_l;
%theta_critical = 0.06; 
%u_critical= sqrt(theta_critical*(s-1)* 9.81 * D); %critical bed shear velocity
%f_drag_critical= 0.5* density_l*1000 * CD * A * (u_critical)^2;

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
        
        
        % Along that trajectory, at each time step
        
        for ind = 1:length(data.t)
            
            xp = data.xf(ind);
            yp = data.yf(ind);
            zp = data.zf(ind);
            up = data.uf(ind);
            vp = data.vf(ind);
            wp = data.wf(ind);
            %velocity_p = [up,vp,wp];
            
            % storage
            %data_large{kFiles,iEvents,ind} = velocity_p;
            
            counter = 0;
            drag_force = 0;
            
            for j = 1:numSmallTraj % for all small ones
                
                distance = dist([small(j).xf,small(j).yf,small(j).zf]',[xp,yp,zp])';
                
                relevant = find(small(j).t == data.t(ind) &  distance <= R);
                
                if ~isempty(relevant)
                    
                    tmp = small(j);
                    tmp.uf = tmp.uf(relevant);
                    tmp.vf = tmp.vf(relevant);
                    tmp.wf = tmp.wf(relevant);
                    %tmp_velocity = [tmp.uf, tmp.vf, tmp.wf];
                    
                    velocity_p_f = [tmp.uf-up,tmp.vf-vp,tmp.wf-wp];
                    %velocity_p_f = [up-tmp.uf,vp-tmp.vf,wp-tmp.wf];
                    drag_force = drag_force + 0.5* density_l*1000* CD * A * norm(velocity_p_f) * velocity_p_f;% sum of all the relevant small particles
                    
                    counter = counter + 1;
                    
                    % if verbose,plot_long_trajectories(tmp(relevant),1,hf);   end
                    
                end % relevant
            end % smalltraj
            
            %data_small{kFiles,iEvents,ind} =  tmp_velocity;
            data_drag{kFiles,iEvents,ind} =  norm (drag_force) / counter; % average  of the   absolute drag force
            time_data_small{kFiles, iEvents,ind} = data.t(ind);
            
          
            
        end
    end
end

save drag_force_R25_110_article.mat


%%
figure, hold on,
for kFiles = 1:size(data_drag,1)
    for iEvents = 1:size(data_drag,2)
        if ~isempty([time_of_event{kFiles,iEvents}])
            plot([time_data_small{kFiles,iEvents,:}]-[time_data_small{kFiles,iEvents,1}],[data_drag{kFiles,iEvents,:}]); 
            tmp = [data_drag{kFiles,iEvents,:}];
            t1 = max([time_of_event{kFiles,iEvents}]-[time_data_small{kFiles,iEvents,1}]-1,1);
            plot(t1,tmp(t1+1),'ro','MarkerSize',10); % time of the event
           
        end
    end
end
xlim = get(gca,'xlim');
plot(xlim,[f_drag_critical,f_drag_critical], 'r'); % make a line
%plot(xlim,[f_drag_critical_1,f_drag_critical_1], 'g'); % make a line
hold off
