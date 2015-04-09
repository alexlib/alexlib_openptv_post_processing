clc, clear all

dist = @(x,y) sqrt((x(1,:)-y(1)).^2 + (x(2,:)-y(2)).^2 + (x(3,:)-y(3)).^2);

if ismac
    matdirectory = '/Users/alex/Dropbox/resuspension/2011/trajectories(186-194)/events_80';
else
    matdirectory ='C:\Users\hadar\Dropbox\resuspension\2011\trajectories(186-194)\events_80';
end

large_ones = dir(fullfile(matdirectory,'large*'));

time_of_event = {};

%forces:
data_drag = {};
data_inertia = {};
data_lift = {};
data_pressure = {};

% R = 15;


D = 550E-6; %diameter of the solid particle [m]
r = D/2;
A= pi * r^2 ;% the projected grain area perpendicular to flow direction
V = pi*D^3/6;

density_l = 0.998 ; % density of water ,temperature= 22 [gr/cm^3]
density_s = 1.065 ; %density of solid patricles [gr/cm^3]

mass_particle = density_s *1000 * V; % [Kg]
 

CD = 3; % drag coefficient
%CI = 1.5; % inertia coefficient

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
            
            axp = data.axf(ind);
            ayp = data.ayf(ind);
            azp = data.azf(ind);
            
            velocity_p = [up,vp,wp];
            acceleration_p = [axp,ayp,azp];
            
            %acceleration_p_1 = acceleration_p/norm(acceleration_p);
            
            force_r = acceleration_p * mass_particle;
            
            % storage
            %data_large{kFiles,iEvents,ind} = velocity_p;
            
            counter = 0;
            
            drag_force = 0;
            
            closest_distance = Inf;
            closest_point_velocity = zeros(1,3);
            closest_point_acceleration = zeros(1,3);
            
            for j = 1:numSmallTraj % for all small ones
                
                relevant = find(small(j).t == data.t(ind)); % &  distance <= R);
                
                distance = dist([small(j).xf(relevant),small(j).yf(relevant),small(j).zf(relevant)]',[xp,yp,zp])';
                
                if distance < closest_distance & small(j).yf(relevant) > -21.5 % not on the ground
                    closest_distance = distance;
                    closest_point_velocity  = [small(j).uf(relevant),...
                        small(j).vf(relevant),small(j).wf(relevant)];
                    closest_point_acceleration = [small(j).axf(relevant), ...
                        small(j).ayf(relevant), small(j).azf(relevant)];
                    
                end
            end
            
            tmp.uf = closest_point_velocity(1);
            tmp.vf = closest_point_velocity(2);
            tmp.wf = closest_point_velocity(3);
            
            tmp.axf = closest_point_acceleration(1);
            tmp.ayf = closest_point_acceleration(2);
            tmp.azf = closest_point_acceleration(3);
            
            acceleration_f = [tmp.axf,tmp.ayf,tmp.azf];
            
            velocity_p_f = [tmp.uf-up,tmp.vf-vp,tmp.wf-wp]; % relative velocity
            velocity_p_f_1 = velocity_p_f/ norm(velocity_p_f);
            
            
           % drag_force = 0.5* density_l*1000* CD * A * norm(velocity_p_f) * velocity_p_f;
            
            tmp.axf = closest_point_acceleration(1);
            tmp.ayf = closest_point_acceleration(2);
            tmp.azf = closest_point_acceleration(3);
            
             acceleration_p_f = [axp-tmp.axf,ayp-tmp.ayf,azp-tmp.azf];
             
           added_mass = 0.5; % according to "drag and lift on microscopic bubbles entrained by a vortex" (1995)
            
           inertia_force =  -4/3 * pi * r^3 *(density_s * 1000 * acceleration_p + density_l * 1000 * added_mass  * acceleration_p_f);            
           
           pressure_force = 4/3 * density_l * pi * r^3 * acceleration_f;
           buoyancy_force = 4/3 *(density_s - density_l)* 1000 * 9.81 * pi * r^3;
           buoyancy_force_1 = [0,0,buoyancy_force];
           
           total_force = inertia_force + pressure_force + buoyancy_force_1;
           drag_force = total_force * velocity_p_f_1'; % total_force_horizontal
           lift_force = total_force - drag_force; % total_force_verticle
           
           %lift_force = force_r - drag_force;
            data_force_r {kFiles,iEvents,ind} =  force_r;
            data_total_force{kFiles,iEvents,ind} =  total_force;
            data_drag{kFiles,iEvents,ind} =  drag_force; % norm (drag_force) / counter; % average  of the   absolute drag force
            data_inertia{kFiles,iEvents,ind} = inertia_force;
            data_lift{kFiles,iEvents,ind} = lift_force;
            %data_buoyancy{kFiles,iEvents,ind} = buoyancy_force;
            data_pressure{kFiles,iEvents,ind} = pressure_force;
            
            time_data_small{kFiles, iEvents,ind} = data.t(ind);

        end
    end
end

save forces_closest_particle_80.mat

%%force_r
figure, hold on,
for kFiles = 1:size(data_drag,1)
    for iEvents = 1:size(data_drag,2)
        if ~isempty([time_of_event{kFiles,iEvents}]);
            tmp = cat(1,data_force_r{kFiles,iEvents,:});
            norm_force_r = sum(tmp.^2,2).^(0.5);
            plot([time_data_small{kFiles,iEvents,:}]-[time_data_small{kFiles,iEvents,1}],norm_force_r);
            % tmp = [data_drag{kFiles,iEvents,:}];
            t1 = max([time_of_event{kFiles,iEvents}]-[time_data_small{kFiles,iEvents,1}]-1,1);
            plot(t1,norm_force_r(t1+1),'ro','MarkerSize',10); % time of the event
            
        end
    end
end

hold off

%%total_force
figure, hold on,
for kFiles = 1:size(data_drag,1)
    for iEvents = 1:size(data_drag,2)
        if ~isempty([time_of_event{kFiles,iEvents}]);
            tmp = cat(1,data_total_force{kFiles,iEvents,:});
            norm_total_force = sum(tmp.^2,2).^(0.5);
            plot([time_data_small{kFiles,iEvents,:}]-[time_data_small{kFiles,iEvents,1}],norm_total_force);
            % tmp = [data_drag{kFiles,iEvents,:}];
            t1 = max([time_of_event{kFiles,iEvents}]-[time_data_small{kFiles,iEvents,1}]-1,1);
            plot(t1,norm_total_force(t1+1),'ro','MarkerSize',10); % time of the event
            
        end
    end
end

hold off

%%drag
figure, hold on,
for kFiles = 1:size(data_drag,1)
    for iEvents = 1:size(data_drag,2)
        if ~isempty([time_of_event{kFiles,iEvents}]);
            tmp = cat(1,data_drag{kFiles,iEvents,:});
            norm_drag = sum(tmp.^2,2).^(0.5);
            plot([time_data_small{kFiles,iEvents,:}]-[time_data_small{kFiles,iEvents,1}],norm_drag);
            % tmp = [data_drag{kFiles,iEvents,:}];
            t1 = max([time_of_event{kFiles,iEvents}]-[time_data_small{kFiles,iEvents,1}]-1,1);
            plot(t1,norm_drag(t1+1),'ro','MarkerSize',10); % time of the event
            
        end
    end
end
hold off

%lift
figure, hold on,
for kFiles = 1:size(data_drag,1)
    for iEvents = 1:size(data_drag,2)
        if ~isempty([time_of_event{kFiles,iEvents}]);
            tmp = cat(1,data_lift{kFiles,iEvents,:});
            norm_lift = sum(tmp.^2,2).^(0.5);
            plot([time_data_small{kFiles,iEvents,:}]-[time_data_small{kFiles,iEvents,1}],norm_lift);
            % tmp = [data_drag{kFiles,iEvents,:}];
            t1 = max([time_of_event{kFiles,iEvents}]-[time_data_small{kFiles,iEvents,1}]-1,1);
            plot(t1,norm_lift(t1+1),'ro','MarkerSize',10); % time of the event
            
        end
    end
end



%% inertia
figure, hold on,
for kFiles = 1:size(data_drag,1)
    for iEvents = 1:size(data_drag,2)
        if ~isempty([time_of_event{kFiles,iEvents}]);
            tmp = cat(1,data_inertia{kFiles,iEvents,:});
            norm_inertia = sum(tmp.^2,2).^(0.5);
            plot([time_data_small{kFiles,iEvents,:}]-[time_data_small{kFiles,iEvents,1}],norm_inertia);
            % tmp = [data_drag{kFiles,iEvents,:}];
            t1 = max([time_of_event{kFiles,iEvents}]-[time_data_small{kFiles,iEvents,1}]-1,1);
            plot(t1,norm_inertia(t1+1),'ro','MarkerSize',10); % time of the event
            
        end
    end
end

hold off



%pressure
figure, hold on,
for kFiles = 1:size(data_drag,1)
    for iEvents = 1:size(data_drag,2)
        if ~isempty([time_of_event{kFiles,iEvents}]);
            tmp = cat(1,data_pressure{kFiles,iEvents,:});
            norm_pressure = sum(tmp.^2,2).^(0.5);
            plot([time_data_small{kFiles,iEvents,:}]-[time_data_small{kFiles,iEvents,1}],norm_pressure);
            % tmp = [data_drag{kFiles,iEvents,:}];
            t1 = max([time_of_event{kFiles,iEvents}]-[time_data_small{kFiles,iEvents,1}]-1,1);
            plot(t1,norm_pressure(t1+1),'ro','MarkerSize',10); % time of the event
            
        end
    end
end

hold off







%% Sum of both
%figure, hold on,
%for kFiles = 1:size(data_drag,1)
 %   for iEvents = 1:size(data_drag,2)
  %      if ~isempty([time_of_event{kFiles,iEvents}]);
   %         tmpI = cat(1,data_inertia{kFiles,iEvents,:});
    %        tmpD = cat(1,data_drag{kFiles,iEvents,:});
     %       tmp = tmpI + tmpD; % vector sum
      %      norm_total = sum(tmp.^2,2).^(0.5);
       %     plot([time_data_small{kFiles,iEvents,:}]-[time_data_small{kFiles,iEvents,1}],norm_total);
            % tmp = [data_drag{kFiles,iEvents,:}];
        %    t1 = max([time_of_event{kFiles,iEvents}]-[time_data_small{kFiles,iEvents,1}]-1,1);
         %   plot(t1,norm_total(t1+1),'ro','MarkerSize',10); % time of the event
            
        %end
   % end
%end
%xlim = get(gca,'xlim');
%plot(xlim,[f_drag_critical,f_drag_critical], 'r'); % make a line
%plot(xlim,[f_drag_critical_1,f_drag_critical_1], 'g'); % make a line
%hold off