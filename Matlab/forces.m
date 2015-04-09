clc, clear all

dist = @(x,y) sqrt((x(1,:)-y(1)).^2 + (x(2,:)-y(2)).^2 + (x(3,:)-y(3)).^2);

if ismac
    matdirectory = '/Users/alex/Dropbox/resuspension/2011/trajectories(186-194)';
else
    matdirectory ='C:\Users\hadar\Dropbox\resuspension\Matlab\trajectories_alex_2012(186-194)\events_110';
end

large_ones = dir(fullfile(matdirectory,'large*'));


time_data_small = {};
data_drag = {};
data_inertia = {};
data_lift = {};
data_pressure = {};
%data_force_r= {};
data_total_force = {};

time_of_event = {};

R = 25;


D = 550E-6; %diameter of the solid particle [m]
r = D/2;
A= pi * r^2 ;% the projected grain area perpendicular to flow direction
V = 4/3 * pi * r^3; %[m^3]

density_l = 0.998 * 1000; % kg/m^3 ; % density of water ,temperature= 22 
density_s = 1.062 * 1000; % kg/m^3 ; %density of solid patricles [gr/cm^3]

mass_particle = density_s * V; % [Kg]
%Cd = 4; % drag coefficient

added_mass = 0.5; % according to "drag and lift on microscopic bubbles entrained by a vortex" (1995)
g = 9.81; % m/s^2 



% according to the article: impulse an particle dislodgement undr turbulent flow conditions
%ws = (density_s -density_l)* 1000 * pi  * 9.81 *4 / 3 * r^3 ;%submerged particle weight
%fv = (1+0.5*(density_l/(density_s-density_l))); %
%f_drag_critical =  fv * ws  ;
%u_square_critical = 2*fv*ws/(density_l* CD *A);


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
            up = data.uf(ind)/1000; %[m/s]
            vp = data.vf(ind)/1000;
            wp = data.wf(ind)/1000;
            
            axp = data.axf(ind)/1000;
            ayp = data.ayf(ind)/1000;
            azp = data.azf(ind)/1000;
            
            velocity_p = [up,vp,wp]; %[m/s]
            acceleration_p = [axp,ayp,azp];% [m/s^2]
            %acceleration_p_1 = acceleration_p/norm(acceleration_p);
            
            %force_r = acceleration_p * mass_particle; % f_r = m*a_particle
            
            
            % storage
            %data_large{kFiles,iEvents,ind} = velocity_p;
            
            counter = 0;

            
                    [inertia_force, pressure_force,...
                        total_force, drag_force, lift_force,drag_Urev ] = deal(zeros(1,3));
                    
            
            
            for j = 1:numSmallTraj % for all small ones
                
                distance = dist([small(j).xf,small(j).yf,small(j).zf]',[xp,yp,zp])';
                
                relevant = find(small(j).t == data.t(ind) &  distance <= R);
                
                if isempty(relevant), continue, end
                    
                    tmp = small(j);
                    
                    
                    if tmp.yf(relevant) <= -21 ,continue, end
                    
                    uf = tmp.uf(relevant)/1000; %[m/s]
                    vf = tmp.vf(relevant)/1000;
                    wf = tmp.wf(relevant)/1000;
                    
                    axf = tmp.axf(relevant)/1000; %[m/s^2]
                    ayf = tmp.ayf(relevant)/1000;
                    azf = tmp.azf(relevant)/1000;
                    
                    acceleration_f = [axf,ayf,azf]; %[m/s^2]
                    
                    acceleration_p_f = [axp-axf,ayp-ayf,azp-azf]; 
                    
                    
                    velocity_f_p = [uf-up,vf-vp,wf-wp];
                    velocity_f_p_1 = velocity_f_p / norm(velocity_f_p);
                    
                    drag_Urev = drag_Urev + 0.5* density_l* A * norm(velocity_f_p) * velocity_f_p;
                    
                    
                    buoyancy_force = -(density_s - density_l) * V * g;
                    buoyancy_force_1 = [0,buoyancy_force,0]; % i.e. parallel to g
                    
                    inertia_force =  inertia_force + (-V *(density_s  * acceleration_p + density_l * added_mass  * acceleration_p_f));
                    pressure_force = pressure_force + (V * density_l * acceleration_f);
                    
                    total_force = total_force + (-(inertia_force + pressure_force + buoyancy_force_1));
                    
                    
                    %drag_force = drag_force + (force_r * velocity_f_p_1'); % total_force_horizontal
                    %lift_force = lift_force + (force_r - drag_force); % total_force_verticle
                    
                    
                    drag_force = drag_force + (total_force * velocity_f_p_1'); % total_force_horizontal
                    lift_force = lift_force + (total_force - drag_force); % total_force_verticle
                    
                    
                    counter = counter + 1;
                    
                    
                    
%                 end % relevant
            end % smalltraj
            
            
            %data_force_r {kFiles,iEvents,ind} =  force_r;
            data_total_force{kFiles,iEvents,ind} = total_force/counter;
            data_drag{kFiles,iEvents,ind} =  drag_force/counter;
            data_inertia{kFiles,iEvents,ind} = inertia_force/counter;
            data_lift{kFiles,iEvents,ind} = lift_force/counter;
            data_pressure{kFiles,iEvents,ind} = pressure_force/counter;      
            data_drag_Urev{kFiles,iEvents,ind} =  norm (drag_Urev) / counter; % average  of the   absolute drag force
            
            time_data_small{kFiles, iEvents,ind} = data.t(ind);
            
            
            
            
        end
    end
end

save forces_R25_110.mat


% %%force_r
% figure, hold on,
% for kFiles = 1:size(data_drag,1)
%     for iEvents = 1:size(data_drag,2)
%         if ~isempty([time_of_event{kFiles,iEvents}]);
%             tmp = cat(1,data_force_r{kFiles,iEvents,:});
%             norm_force_r = sum(tmp.^2,2).^(0.5);
%             plot([time_data_small{kFiles,iEvents,:}]-[time_data_small{kFiles,iEvents,1}],norm_force_r);
%             % tmp = [data_drag{kFiles,iEvents,:}];
%             t1 = max([time_of_event{kFiles,iEvents}]-[time_data_small{kFiles,iEvents,1}]-1,1);
%             plot(t1,norm_force_r(t1+1),'ro','MarkerSize',10); % time of the event
%             
%         end
%     end
% end


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



% %%force_r vs total force
% figure, hold on,
% for kFiles = 1:size(data_drag,1)
%     for iEvents = 1:size(data_drag,2)
%         if ~isempty([time_of_event{kFiles,iEvents}]);
%             tmp1 = cat(1,data_force_r{kFiles,iEvents,:});
%             norm_force_r = sum(tmp1.^2,2).^(0.5);
%             tmp2 = cat(1,data_total_force{kFiles,iEvents,:});
%             norm_total_force = sum(tmp2.^2,2).^(0.5);
%              t1 = max([time_of_event{kFiles,iEvents}]-[time_data_small{kFiles,iEvents,1}]-1,1);
% 
% scatter(norm_force_r(t1+1),norm_total_force(t1+1),100);
% % scatter(tmp1,tmp2);
%             
%         end
%     end
% end
% xlabel('F = ma')
% ylabel('sum of forces');


%%pressure vs total force
figure, hold on,
for kFiles = 1:size(data_drag,1)
    for iEvents = 1:size(data_drag,2)
        if ~isempty([time_of_event{kFiles,iEvents}]);
            tmp1 = cat(1,data_pressure{kFiles,iEvents,:});
            norm_pressure = sum(tmp1.^2,2).^(0.5);
            tmp2 = cat(1,data_total_force{kFiles,iEvents,:});
            norm_total_force = sum(tmp2.^2,2).^(0.5);
             t1 = max([time_of_event{kFiles,iEvents}]-[time_data_small{kFiles,iEvents,1}]-1,1);

scatter(norm_pressure(t1+1),norm_total_force(t1+1),100);
% scatter(tmp1,tmp2);
            
        end
    end
end
xlabel('Pressure force')
ylabel('sum of forces');


%%drag vs lift force
figure, hold on,
for kFiles = 1:size(data_drag,1)
    for iEvents = 1:size(data_drag,2)
        if ~isempty([time_of_event{kFiles,iEvents}]);
            tmp1 = cat(1,data_drag{kFiles,iEvents,:});
            norm_drag = sum(tmp1.^2,2).^(0.5);
            tmp2 = cat(1,data_lift{kFiles,iEvents,:});
            norm_lift = sum(tmp2.^2,2).^(0.5);
             t1 = max([time_of_event{kFiles,iEvents}]-[time_data_small{kFiles,iEvents,1}]-1,1);

scatter(norm_drag(t1+1),norm_lift(t1+1),100);
% scatter(tmp1,tmp2);
            
        end
    end
end
xlabel('Drag force')
ylabel('Lift force');


%drag vs urev^2 force
figure, hold on,
for kFiles = 1:size(data_drag,1)
  for iEvents = 1:size(data_drag,2)
       if ~isempty([time_of_event{kFiles,iEvents}]);
          tmp1 = cat(1,data_drag_Urev{kFiles,iEvents,:});
           norm_drag_Urev = sum(tmp1.^2,2).^(0.5);
           tmp2 = cat(1,data_drag{kFiles,iEvents,:});
           norm_drag = sum(tmp2.^2,2).^(0.5);
            t1 = max([time_of_event{kFiles,iEvents}]-[time_data_small{kFiles,iEvents,1}]-1,1);

scatter(norm_drag_Urev(t1+1),norm_drag(t1+1),100);
 %scatter(tmp1,tmp2);
            
       end
   end
end
xlabel('urev^2')
ylabel('drag force');



