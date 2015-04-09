clc, clear all

dist = @(x,y) sqrt((x(1,:)-y(1)).^2 + (x(2,:)-y(2)).^2 + (x(3,:)-y(3)).^2);

if ismac
    matdirectory = '\Users\alex\Dropbox\resuspension\Matlab\trajectories_alex_2012(186-194)\events_80';
else
    matdirectory ='C:\Users\dell\Dropbox\resuspension\Matlab\trajectories_alex_2012(186-194)\events_110';
end

large_ones = dir(fullfile(matdirectory,'large*'));

[   time_data, time_of_event, drag_Urev, drag_force,  data_acceleration_p...
    inertia_force, inertia_force_part1, lift_force, pressure_force, ...
    total_force, buoyancy_force_1, inertia_force_part2,norm_velocity_f_p1,  ...
    data_velocity_f_p1, data_velocity_p_f1,data_velocity_f1, norm_velocity_f_p, ...
    data_acceleration_f1, data_acceleration_p_f1,...
    Xp,Yp,Zp,Up,Vp,Wp,Re] = deal({});



R = 25; %mm


D = 550E-6; %diameter of the solid particle [m]
r = D/2;
A= pi * r^2 ;% the projected grain area perpendicular to flow direction
V = 4/3 * pi * r^3; %[m^3]
k_v=1E-6;

density_l = 0.998 * 1000; % kg/m^3 ; % density of water ,temperature= 22
density_s = 1.062 * 1000; % kg/m^3 ; %density of solid patricles [gr/cm^3]

mass_particle = density_s * V; % [Kg]

added_mass = 0.5; % according to "drag and lift on microscopic bubbles entrained by a vortex" (1995)
g = 9.81; % m/s^2

dt = 1/160;
%tau = 6:-1:1;


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
    
     disp(kFiles)
    
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
            
            Xp{kFiles,iEvents,ind} = xp;
            Yp{kFiles,iEvents,ind} = yp;
            Zp{kFiles,iEvents,ind} = zp;
            
            Up{kFiles,iEvents,ind} = up;
            Vp{kFiles,iEvents,ind} = vp;
            Wp{kFiles,iEvents,ind} = wp;
            
            
            axp = data.axf(ind)/1000; %[m/s^2]
            ayp = data.ayf(ind)/1000;
            azp = data.azf(ind)/1000;
            
            velocity_p = [up,vp,wp]; %[m/s]
            acceleration_p = [axp,ayp,azp];% [m/s^2]
            
            velocity_p_f1 = zeros(1,3);
            velocity_f_p1 = zeros(1,3);
            acceleration_f1 = zeros(1,3);
            acceleration_p_f1= zeros(1,3);
            velocity_f1=zeros(1,3);
            
            counter = 0;
            
            for j = 1:numSmallTraj % for all small ones
                
                distance = dist([small(j).xf,small(j).yf,small(j).zf]',[xp,yp,zp])';
                
                relevant = find(small(j).t == data.t(ind) &  distance <= R);
                
                if isempty(relevant), continue, end
                
                tmp = small(j);
                
                if tmp.yf(relevant) <= -21.5 && abs(tmp.yf(relevant)-yp )>= 10  ,continue, end
                counter = counter + 1;
                
                uf = tmp.uf(relevant)/1000; %[m/s]
                vf = tmp.vf(relevant)/1000;
                wf = tmp.wf(relevant)/1000;
                
                axf = tmp.axf(relevant)/1000;
                ayf= tmp.ayf(relevant)/1000;
                azf= tmp.azf(relevant)/1000;
                
                velocity_f = [uf,vf,wf];
                velocty_f1= velocity_f1+velocity_f;
                
                velocity_p_f = [up-uf,vp-vf,wp-wf];
                velocity_p_f1= velocity_p_f1 + velocity_p_f;
                
                velocity_f_p = [uf-up,vf-vp,wf-wp];
                velocity_f_p1= velocity_f_p1 + velocity_f_p;
                
                
                acceleration_f = [axf,ayf,azf]; %[m/s^2]
                acceleration_f1= acceleration_f1+ acceleration_f;
                
                acceleration_p_f = [axp-axf,ayp-ayf,azp-azf];
                acceleration_p_f1 = acceleration_p_f1 + acceleration_p_f;
  
                
            end % smalltraj
            
            if counter > 0
                % Averaging in the sphere of radius R
                data_velocity_f1{kFiles,iEvents,ind} = velocity_f1/counter;
                
                data_velocity_p_f1{kFiles,iEvents,ind} = velocity_p_f1/counter;
                data_velocity_f_p1{kFiles,iEvents,ind} = velocity_f_p1/counter;
                
                norm_velocity_f_p1 {kFiles,iEvents,ind} = data_velocity_f_p1{kFiles,iEvents,ind} / ...
                    norm(data_velocity_f_p1{kFiles,iEvents,ind});
                
                data_acceleration_f1{kFiles,iEvents,ind} = acceleration_f1/counter;
                data_acceleration_p_f1{kFiles,iEvents,ind} = acceleration_p_f1/counter;
            else
                data_velocity_f1{kFiles,iEvents,ind} = zeros(1,3);
                data_velocity_p_f1{kFiles,iEvents,ind} = zeros(1,3);
                data_velocity_f_p1{kFiles,iEvents,ind} = zeros(1,3);
                
                norm_velocity_f_p1 {kFiles,iEvents,ind} = zeros(1,3);
                
                data_acceleration_f1{kFiles,iEvents,ind} = zeros(1,3);
                data_acceleration_p_f1{kFiles,iEvents,ind} = zeros(1,3);
            end
            
            
            data_velocity_p_f1{kFiles,iEvents,ind} = velocity_p_f1/counter;
            data_velocity_f_p1{kFiles,iEvents,ind} = velocity_f_p1/counter;
            
            norm_velocity_f_p1 {kFiles,iEvents,ind}= data_velocity_f_p1{kFiles,iEvents,ind} / ...
            norm(data_velocity_f_p1{kFiles,iEvents,ind});
            
            
            data_acceleration_f1{kFiles,iEvents,ind} = acceleration_f1/counter;
            data_acceleration_p_f1{kFiles,iEvents,ind} = acceleration_p_f1/counter;
            
            data_acceleration_p{kFiles,iEvents,ind} = acceleration_p;
            time_data{kFiles, iEvents,ind} = data.t(ind);
           
            
 
        
%         %basset_force1 = cell(kFiles,iEvents,length(data.t));
%         for i = 1:length(data.t)
%             basset_force1{kFiles,iEvents,i} = zeros(1,3);
%         end
%         
%         for ind = 6:length(data.t)
%             
%             urel = [data_velocity_p_f1{kFiles,iEvents,ind-5:ind}]; % 6 x [U,V,W]
%             
%             
%             for k = 1:3
%                 % dt in seconds, urel in m/sec
%                 du = gradient(urel(k:3:end),dt);
%                 
%                 y = du./ sqrt(pi * k_v * tau * dt);
%                 
%                 tmp =  6*pi * r^2* k_v * density_l * dt * trapz(y);
%                 
%                 basset_force1{kFiles,iEvents,ind}(k) = tmp;
%                 
%             end % k
            
            
            drag_Urev{kFiles,iEvents,ind} = 0.5* density_l* A * ...
                norm(data_velocity_f_p1{kFiles,iEvents,ind})* data_velocity_f_p1{kFiles,iEvents,ind};
            
            buoyancy_force_1 {kFiles,iEvents,ind}= [0,-(density_s - density_l) * V * g,0]; 
            
            inertia_force {kFiles,iEvents,ind}=  -V * (density_s  * data_acceleration_p{kFiles,iEvents,ind} +....
                density_l * added_mass  * data_acceleration_p_f1{kFiles,iEvents,ind});
            
            inertia_force_part1 {kFiles,iEvents,ind}=  -V * density_s  * data_acceleration_p{kFiles,iEvents,ind};
               
            inertia_force_part2 {kFiles,iEvents,ind}=  -V * density_l * ...
            added_mass* data_acceleration_p_f1{kFiles,iEvents,ind};
            
            pressure_force {kFiles,iEvents,ind}= V * density_l *...
                data_acceleration_f1{kFiles,iEvents,ind};
            
            total_force{kFiles,iEvents,ind} = -(inertia_force{kFiles,iEvents,ind} +....
                pressure_force{kFiles,iEvents,ind} + buoyancy_force_1{kFiles,iEvents,ind});
           
  
            
            drag_force{kFiles,iEvents,ind} = dot(total_force{kFiles,iEvents,ind}, ...
                data_velocity_f_p1{kFiles,iEvents,ind}).* norm_velocity_f_p1{kFiles,iEvents,ind};
             
            
            lift_force {kFiles,iEvents,ind}= total_force{kFiles,iEvents,ind} - ...
                drag_force{kFiles,iEvents,ind}; % total_force_verticle
            
            Re{kFiles,iEvents,ind} = D*(data_velocity_f_p1{kFiles,iEvents,ind})/k_v;
             
            drag_Vs_Urev{kFiles,iEvents,ind}= drag_force{kFiles,iEvents,ind}./ drag_Urev{kFiles,iEvents,ind};
            
            lift_Vs_Urev{kFiles,iEvents,ind}= lift_force{kFiles,iEvents,ind}./ drag_Urev{kFiles,iEvents,ind};
           
            
             end % ind 
    end % iEvents
end % Kfiles

save forces_R25_nobasset_10_110.mat

%%
% all forces togther
for kFiles= 1:size(total_force,1);
   for iEvents=1:size(total_force,2);  
    if ~isempty([time_of_event{kFiles,iEvents}]);
       figure, hold on, 
        tmp = cat(1,total_force{kFiles,iEvents,:});
        norm_total_force = sum(tmp.^2,2).^(0.5);
%         tmp1 = cat(1,basset_force1{kFiles,iEvents,:});
%         norm_basset_force1 = sum(tmp1.^2,2).^(0.5);
        tmp3 = cat(1,inertia_force{kFiles,iEvents,:});
        norm_inertia = sum(tmp3.^2,2).^(0.5);
        tmp4 = cat(1,inertia_force_part1{kFiles,iEvents,:});
        norm_inertia_part1 = sum(tmp4.^2,2).^(0.5);
        tmp5 = cat(1,inertia_force_part2{kFiles,iEvents,:});
        norm_inertia_part2 = sum(tmp5.^2,2).^(0.5);
        tmp6 = cat(1,pressure_force{kFiles,iEvents,:});
        norm_pressure = sum(tmp6.^2,2).^(0.5);
        tmp7 = cat(1,buoyancy_force_1 {kFiles,iEvents,:});
        norm_buoyancy = sum(tmp7.^2,2).^(0.5);
        tmp8 = cat(1,drag_force {kFiles,iEvents,:});
        norm_drag = sum(tmp8.^2,2).^(0.5);
        tmp9 = cat(1,lift_force {kFiles,iEvents,:});
        norm_lift = sum(tmp9.^2,2).^(0.5);
      
        plot([time_data{kFiles,iEvents,1:end}]-[time_data{kFiles,iEvents,1}],norm_total_force,'b');
        %plot([time_data{kFiles,iEvents,1:end}]-[time_data{kFiles,iEvents,6}],norm_basset_force1,'r');
        plot([time_data{kFiles,iEvents,1:end}]-[time_data{kFiles,iEvents,1}],norm_inertia,'g');

        plot([time_data{kFiles,iEvents,1:end}]-[time_data{kFiles,iEvents,1}],norm_pressure,'m');
        plot([time_data{kFiles,iEvents,1:end}]-[time_data{kFiles,iEvents,1}], norm_buoyancy,'c');
        
        plot([time_data{kFiles,iEvents,1:end}]-[time_data{kFiles,iEvents,1}], norm_drag,'y');
        plot([time_data{kFiles,iEvents,1:end}]-[time_data{kFiles,iEvents,1}], norm_lift,'c');
         
        t1 = max([time_of_event{kFiles,iEvents}]-[time_data{kFiles,iEvents,1}]-1,1);
        plot(t1,norm_total_force(t1+1),'r*','MarkerSize',10); % time of the event 
        %plot(t1,norm_basset_force1(t1+6),'r*','MarkerSize',10); 
        plot(t1,norm_inertia(t1+1),'r*','MarkerSize',10); 
        plot(t1,norm_pressure(t1+1),'r*','MarkerSize',10); 
        plot(t1,norm_buoyancy(t1+1),'r*','MarkerSize',10);
        plot(t1,norm_drag(t1+1),'r*','MarkerSize',10); 
        plot(t1,norm_lift(t1+1),'r*','MarkerSize',10);
        
        
        xlabel('time from the first event [frame]')
        ylabel('forces [N]');
    end
  end
end
%%%%%%%%%%%%%%%%%
% all forces togther_vertical
for kFiles= 1:size(total_force,1);
   for iEvents=1:size(total_force,2);  
    if ~isempty([time_of_event{kFiles,iEvents}]);
       figure, hold on, 
        tmp = cat(1,total_force{kFiles,iEvents,:});
        vertical_total_force = tmp(:,2);
%         tmp1 = cat(1,basset_force1{kFiles,iEvents,:});
%         vertical_basset_force1 = tmp1(:,2);
        tmp3 = cat(1,inertia_force{kFiles,iEvents,:});
        vertical_inertia = tmp3(:,2);
        tmp4 = cat(1,inertia_force_part1{kFiles,iEvents,:});
        vertical_inertia_part1 = tmp4(:,2);
        tmp5 = cat(1,inertia_force_part2{kFiles,iEvents,:});
        vertical_inertia_part2 = tmp5(:,2);
        tmp6 = cat(1,pressure_force{kFiles,iEvents,:});
        vertical_pressure = tmp6(:,2);
        tmp7 = cat(1,buoyancy_force_1 {kFiles,iEvents,:});
        buoyancy = tmp7(:,2);
        tmp8 = cat(1,drag_force {kFiles,iEvents,:});
        vertical_drag = tmp8(:,2);
        tmp9= cat(1,lift_force {kFiles,iEvents,:});
        vertical_lift = tmp9(:,2);
%         
        plot([time_data{kFiles,iEvents,1:end}]-[time_data{kFiles,iEvents,1}],vertical_total_force,'b');
       % plot([time_data{kFiles,iEvents,1:end}]-[time_data{kFiles,iEvents,6}],vertical_basset_force1,'r');
        plot([time_data{kFiles,iEvents,1:end}]-[time_data{kFiles,iEvents,1}],vertical_inertia,'g');
        plot([time_data{kFiles,iEvents,1:end}]-[time_data{kFiles,iEvents,1}], vertical_drag,'y');
        plot([time_data{kFiles,iEvents,1:end}]-[time_data{kFiles,iEvents,1}], vertical_lift,'c');
%         plot([time_data{kFiles,iEvents,6:end}]-[time_data{kFiles,iEvents,6}],vertical_inertia_part1,'y');
%         plot([time_data{kFiles,iEvents,6:end}]-[time_data{kFiles,iEvents,6}],vertical_inertia_part2,'c');
        plot([time_data{kFiles,iEvents,1:end}]-[time_data{kFiles,iEvents,1}], vertical_pressure,'m');
         plot([time_data{kFiles,iEvents,1:end}]-[time_data{kFiles,iEvents,1}], buoyancy,'c');
        
        t1 = max([time_of_event{kFiles,iEvents}]-[time_data{kFiles,iEvents,1}]-1,1);
        plot(t1,vertical_total_force(t1+1),'r*','MarkerSize',10); % time of the event 
%        plot(t1,vertical_basset_force1(t1+6),'r*','MarkerSize',10); 
        plot(t1,vertical_inertia(t1+1),'r*','MarkerSize',10); 
%         plot(t1,vertical_inertia_part1(t1+1),'r*','MarkerSize',10);
%         plot(t1,vertical_inertia_part2(t1+1),'r*','MarkerSize',10); 
        plot(t1,vertical_pressure(t1+1),'r*','MarkerSize',10); 
        plot(t1,buoyancy(t1+1),'r*','MarkerSize',10);
        xlabel('time from the first event [frame]')
        ylabel('vertical forces [N]');
    end
  end
end

%%
% %all forces togther_horizontal
% for kFiles= 1:size(total_force,1);
%    for iEvents=1:size(total_force,2);  
%     if ~isempty([time_of_event{kFiles,iEvents}]);
%        figure, hold on, 
%         tmp = cat(1,total_force{kFiles,iEvents,:});
%         horizontal_total_force = sum(tmp(:,[1 3]).^2,2).^(0.5);
% %         tmp1 = cat(1,basset_force1{kFiles,iEvents,:});
% %         horizontal_basset_force1 = sum(tmp1(:,[1 3]).^2,2).^(0.5);
%         tmp3 = cat(1,inertia_force{kFiles,iEvents,:});
%         horizontal_inertia = sum(tmp3(:,[1 3]).^2,2).^(0.5);
%         tmp4 = cat(1,inertia_force_part1{kFiles,iEvents,:});
%         %horizontal_inertia_part1 = sum(tmp4(:,[1 3]).^2,2).^(0.5);
%         %tmp5 = cat(1,inertia_force_part2{kFiles,iEvents,:});
%         %horizontal_inertia_part2 = sum(tmp5(:,[1 3]).^2,2).^(0.5);
%         tmp6 = cat(1,pressure_force{kFiles,iEvents,:});
%         horizontal_pressure = sum(tmp6(:,[1 3]).^2,2).^(0.5);
%         tmp8 = cat(1,drag_force {kFiles,iEvents,:});
%         horizontal_drag = sum(tmp8(:,[1 3]).^2,2).^(0.5);
%         tmp9= cat(1,lift_force {kFiles,iEvents,:});
%         horizontal_lift = sum(tmp9(:,[1 3]),2).^(0.5);
% 
%         
%         plot([time_data{kFiles,iEvents,1:end}]-[time_data{kFiles,iEvents,1}],horizontal_total_force,'b');
%        % plot([time_data{kFiles,iEvents,1:end}]-[time_data{kFiles,iEvents,6}],horizontal_basset_force1,'r');
%         plot([time_data{kFiles,iEvents,1:end}]-[time_data{kFiles,iEvents,1}],horizontal_inertia,'g');
% %         plot([time_data{kFiles,iEvents,6:end}]-[time_data{kFiles,iEvents,6}],horizontal_inertia_part1,'y');
% %         plot([time_data{kFiles,iEvents,6:end}]-[time_data{kFiles,iEvents,6}],horizontal_inertia_part2,'c');
%         plot([time_data{kFiles,iEvents,1:end}]-[time_data{kFiles,iEvents,1}], horizontal_pressure,'m');
%         plot([time_data{kFiles,iEvents,1:end}]-[time_data{kFiles,iEvents,1}], horizontal_drag,'y');
%         plot([time_data{kFiles,iEvents,1:end}]-[time_data{kFiles,iEvents,1}], horizontal_lift,'c');
% 
%         
%         t1 = max([time_of_event{kFiles,iEvents}]-[time_data{kFiles,iEvents,1}]-1,1);
%         plot(t1,horizontal_total_force(t1+1),'r*','MarkerSize',10); % time of the event 
% %        plot(t1,horizontal_basset_force1(t1+6),'r*','MarkerSize',10); 
%         plot(t1,horizontal_inertia(t1+1),'r*','MarkerSize',10); 
% %         plot(t1,horizontal_inertia_part1(t1+1),'r*','MarkerSize',10);
% %         plot(t1,horizontal_inertia_part2(t1+1),'r*','MarkerSize',10); 
%         plot(t1,horizontal_pressure(t1+1),'r*','MarkerSize',10); 
%         xlabel('time from the first event [frame]')
%         ylabel(' horizontal forces [N]');
% 
%     end
%   end
% end

%all forces togther_horizontal_fix!
for kFiles= 1:size(total_force,1);
   for iEvents=1:size(total_force,2);  
    if ~isempty([time_of_event{kFiles,iEvents}]);
       figure, hold on, 
        tmp = cat(1,total_force{kFiles,iEvents,:});
        horizontal_total_force = sum(tmp(:,[1 3]),2);
%         tmp1 = cat(1,basset_force1{kFiles,iEvents,:});
%         horizontal_basset_force1 = sum(tmp1(:,[1 3]),2);
        tmp3 = cat(1,inertia_force{kFiles,iEvents,:});
        horizontal_inertia = sum(tmp3(:,[1 3]),2);
        tmp4 = cat(1,inertia_force_part1{kFiles,iEvents,:});
        tmp6 = cat(1,pressure_force{kFiles,iEvents,:});
        horizontal_pressure = sum(tmp6(:,[1 3]),2);
        tmp8 = cat(1,drag_force {kFiles,iEvents,:});
        horizontal_drag = sum(tmp8(:,[1 3]),2);
        tmp9= cat(1,lift_force {kFiles,iEvents,:});
        horizontal_lift = sum(tmp9(:,[1 3]),2);

        
        plot([time_data{kFiles,iEvents,1:end}]-[time_data{kFiles,iEvents,1}],horizontal_total_force,'b');
%        plot([time_data{kFiles,iEvents,1:end}]-[time_data{kFiles,iEvents,6}],horizontal_basset_force1,'r');
        plot([time_data{kFiles,iEvents,1:end}]-[time_data{kFiles,iEvents,1}],horizontal_inertia,'g');
%         plot([time_data{kFiles,iEvents,6:end}]-[time_data{kFiles,iEvents,6}],horizontal_inertia_part1,'y');
%         plot([time_data{kFiles,iEvents,6:end}]-[time_data{kFiles,iEvents,6}],horizontal_inertia_part2,'c');
        plot([time_data{kFiles,iEvents,1:end}]-[time_data{kFiles,iEvents,1}], horizontal_pressure,'m');
        plot([time_data{kFiles,iEvents,1:end}]-[time_data{kFiles,iEvents,1}], horizontal_drag,'y');
        plot([time_data{kFiles,iEvents,1:end}]-[time_data{kFiles,iEvents,1}], horizontal_lift,'c');

        
        t1 = max([time_of_event{kFiles,iEvents}]-[time_data{kFiles,iEvents,1}]-1,1);
        plot(t1,horizontal_total_force(t1+1),'r*','MarkerSize',10); % time of the event 
%        plot(t1,horizontal_basset_force1(t1+6),'r*','MarkerSize',10); 
        plot(t1,horizontal_inertia(t1+1),'r*','MarkerSize',10); 
%         plot(t1,horizontal_inertia_part1(t1+1),'r*','MarkerSize',10);
%         plot(t1,horizontal_inertia_part2(t1+1),'r*','MarkerSize',10); 
        plot(t1,horizontal_pressure(t1+1),'r*','MarkerSize',10); 
        xlabel('time from the first event [frame]')
        ylabel(' horizontal forces [N]');

    end
  end
end




%% total_force
figure, hold on,
for kFiles = 1:size(total_force,1)
    for iEvents = 1:size(total_force,2)
        if ~isempty([time_of_event{kFiles,iEvents}]);
            tmp = cat(1,total_force{kFiles,iEvents,:});
            norm_total_force = sum(tmp.^2,2).^(0.5);
            plot([time_data{kFiles,iEvents,1:end}]-[time_data{kFiles,iEvents,1}],norm_total_force);
            t1 = max([time_of_event{kFiles,iEvents}]-[time_data{kFiles,iEvents,1}]-1,1);
            plot(t1,norm_total_force(t1+1),'r*','MarkerSize',10); % time of the event
            
        end
    end
end
xlabel('time from the first event [frame]')
ylabel('total force [N]');


%%total_force Vs Urev
figure, hold on,
for kFiles = 1:size(total_force,1)
    for iEvents = 1:size(total_force,2)
        if ~isempty([time_of_event{kFiles,iEvents}]);
            tmp = cat(1,total_force{kFiles,iEvents,:});
            norm_total_force = sum(tmp.^2,2).^(0.5);
            tmp1 = cat(1,drag_Urev{kFiles,iEvents,:});
            norm_Urev = sum(tmp1.^2,2).^(0.5); 
            t1 = max([time_of_event{kFiles,iEvents}]-[time_data{kFiles,iEvents,1}]-1,1);
            scatter(norm_Urev(t1+1),abs(tmp(t1+1)),100);
            
        end
    end
end
xlabel('0.5*p_l*A*u_r^2')
ylabel('total force [N]');

%%total_force Vs Re_r
figure, hold on,
for kFiles = 1:size(total_force,1)
    for iEvents = 1:size(total_force,2)
        if ~isempty([time_of_event{kFiles,iEvents}]);
            tmp = cat(1,total_force{kFiles,iEvents,:});
            norm_total_force = sum(tmp.^2,2).^(0.5);
            tmp1 = cat(1,Re{kFiles,iEvents,:});
            norm_Re = sum(tmp1.^2,2).^(0.5);
            t1 = max([time_of_event{kFiles,iEvents}]-[time_data{kFiles,iEvents,1}]-1,1);
            scatter(norm_Re(t1+1),norm_total_force(t1+1),100);
            
        end
    end
end
xlabel('Re_r')
ylabel('total force [N]');


%total force vertical_time
figure, hold on,
for kFiles = 1:size(total_force,1)
    for iEvents = 1:size(total_force,2)
        if ~isempty([time_of_event{kFiles,iEvents}]);
            tmp = cat(1,total_force{kFiles,iEvents,:});
            vertical_tmp = tmp(:,2);
            plot([time_data{kFiles,iEvents,1:end}]-[time_data{kFiles,iEvents,1}],vertical_tmp);
            t1 = max([time_of_event{kFiles,iEvents}]-[time_data{kFiles,iEvents,1}]-1,1);
            plot(t1,vertical_tmp(t1+1),'r*','MarkerSize',10); % time of the event
            
        end
    end
end
xlabel('time from the first event [frame]')
ylabel('vertical total force [N]');
% 
% 
% %total force horizontal_time
% figure, hold on,
% for kFiles = 1:size(total_force,1)
%     for iEvents = 1:size(total_force,2)
%         if ~isempty([time_of_event{kFiles,iEvents}]);
%             tmp = cat(1,total_force{kFiles,iEvents,:});
%             norm_horizontal_tmp = sum(tmp(:,[1 3]).^2,2).^(0.5);
%             plot([time_data{kFiles,iEvents,1:end}]-[time_data{kFiles,iEvents,1}],norm_horizontal_tmp);
%             t1 = max([time_of_event{kFiles,iEvents}]-[time_data{kFiles,iEvents,1}]-1,1);
%             plot(t1,norm_horizontal_tmp(t1+1),'r*','MarkerSize',10); % time of the event
%             
%         end
%     end
% end
% xlabel('time from the first event [frame]')
% ylabel('horizontal total force [N]');



%total force vertical
figure, hold on,
for kFiles = 1:size(total_force,1)
    for iEvents = 1:size(total_force,2)
        if ~isempty([time_of_event{kFiles,iEvents}]);
            tmp = cat(1,total_force{kFiles,iEvents,:});
            tmp1 = cat(1,drag_Urev{kFiles,iEvents,:});
            norm_Urev = abs(tmp1(:,2)); 
            t1 = max([time_of_event{kFiles,iEvents}]-[time_data{kFiles,iEvents,1}]-1,1);
            scatter(norm_Urev(t1+1),abs(tmp(t1+1,2)),100);
            
        end
    end
end
xlabel('0.5*p_l*A*u_r^2')
ylabel(' vertical total force [N]');



%total force horizontal
figure, hold on,
for kFiles = 1:size(total_force,1)
    for iEvents = 1:size(total_force,2)
        if ~isempty([time_of_event{kFiles,iEvents}]);
            tmp = cat(1,total_force{kFiles,iEvents,:});
            norm_horizontal_tmp = sum(tmp(:,[1 3]).^2,2).^(0.5);
            tmp1 = cat(1,drag_Urev{kFiles,iEvents,:});
            norm_Urev = abs(tmp1(:,[1 3])); 
            t1 = max([time_of_event{kFiles,iEvents}]-[time_data{kFiles,iEvents,1}]-1,1);
            scatter(norm_Urev(t1+1),norm_horizontal_tmp(t1+1),100);
            
        end
    end
end
xlabel('0.5*p_l*A*u_r^2')
ylabel('total force horizontal [N]');


% %basset_force
% figure, hold on,
% for kFiles = 1:size(basset_force1,1)
%     for iEvents = 1:size(basset_force1,2)
%         if ~isempty([time_of_event{kFiles,iEvents}]);
%             tmp2 = cat(1,basset_force1{kFiles,iEvents,:});
%             norm_basset_force1 = sum(tmp2.^2,2).^(0.5);
%             plot([time_data{kFiles,iEvents,1:end}]-[time_data{kFiles,iEvents,1}],norm_basset_force1);
%             t1 = max([time_of_event{kFiles,iEvents}]-[time_data{kFiles,iEvents,1}]-1,1);
%             plot(t1,norm_basset_force1(t1+1),'r*','MarkerSize',10); % time of the event
%             
%         end
%     end
% end
% xlabel('time from the first event [frame]')
% ylabel('basset force [N]');
% 

% %basset_force_vertical
% figure, hold on,
% for kFiles = 1:size(basset_force1,1)
%     for iEvents = 1:size(basset_force1,2)
%         if ~isempty([time_of_event{kFiles,iEvents}]);
%             tmp2 = cat(1,basset_force1{kFiles,iEvents,:});
%             vertical_basset_force1 = tmp2(:,2);
%             plot([time_data{kFiles,iEvents,1:end}]-[time_data{kFiles,iEvents,1}],vertical_basset_force1);
%             t1 = max([time_of_event{kFiles,iEvents}]-[time_data{kFiles,iEvents,1}]-1,1);
%             plot(t1,vertical_basset_force1(t1+1),'r*','MarkerSize',10); % time of the event
%             
%         end
%     end
% end
% xlabel('time from the first event [frame]')
% ylabel('  vertical basset force [N]');
% 
% %basset_force_horizontal
% figure, hold on,
% for kFiles = 1:size(basset_force1,1)
%     for iEvents = 1:size(basset_force1,2)
%         if ~isempty([time_of_event{kFiles,iEvents}]);
%             tmp2 = cat(1,basset_force1{kFiles,iEvents,:});
%             horizontal_basset_force1 = sum(tmp2(:,[1 3]).^2,2).^(0.5);
%             plot([time_data{kFiles,iEvents,1:end}]-[time_data{kFiles,iEvents,1}],horizontal_basset_force1);
%             t1 = max([time_of_event{kFiles,iEvents}]-[time_data{kFiles,iEvents,1}]-1,1);
%             plot(t1,horizontal_basset_force1(t1+1),'r*','MarkerSize',10); % time of the event
%             
%         end
%     end
% end
% xlabel('time from the first event [frame]')
% ylabel('horizontal basset force [N]');


%% inertia
figure, hold on,
for kFiles = 1:size(inertia_force,1)
    for iEvents = 1:size(inertia_force,2)
        if ~isempty([time_of_event{kFiles,iEvents}]);
            tmp3 = cat(1,inertia_force{kFiles,iEvents,:});
            norm_inertia = sum(tmp3.^2,2).^(0.5);
            plot([time_data{kFiles,iEvents,1:end}]-[time_data{kFiles,iEvents,1}],norm_inertia);
            % tmp = [data_drag{kFiles,iEvents,:}];
            t1 = max([time_of_event{kFiles,iEvents}]-[time_data{kFiles,iEvents,1}]-1,1);
            plot(t1,norm_inertia(t1+1),'r*','MarkerSize',10); % time of the event
            
        end
    end
end
xlabel('time from the first event [frame]')
ylabel('inertia force [N]');

% %% inertia_vertical
% figure, hold on,
% for kFiles = 1:size(inertia_force,1)
%     for iEvents = 1:size(inertia_force,2)
%         if ~isempty([time_of_event{kFiles,iEvents}]);
%             tmp3 = cat(1,inertia_force{kFiles,iEvents,:});
%             vertical_inertia = tmp3(:,2);
%             plot([time_data{kFiles,iEvents,1:end}]-[time_data{kFiles,iEvents,1}],vertical_inertia);
%             % tmp = [data_drag{kFiles,iEvents,:}];
%             t1 = max([time_of_event{kFiles,iEvents}]-[time_data{kFiles,iEvents,1}]-1,1);
%             plot(t1,vertical_inertia(t1+1),'r*','MarkerSize',10); % time of the event
%             
%         end
%     end
% end
% xlabel('time from the first event [frame]')
% ylabel('vertical inertia force [N]');
% 
% %% inertia _horizontal
% figure, hold on,
% for kFiles = 1:size(inertia_force,1)
%     for iEvents = 1:size(inertia_force,2)
%         if ~isempty([time_of_event{kFiles,iEvents}]);
%             tmp3 = cat(1,inertia_force{kFiles,iEvents,:});
%             horizontal_inertia = sum(tmp3(:,[1 3]).^2,2).^(0.5);
%             plot([time_data{kFiles,iEvents,1:end}]-[time_data{kFiles,iEvents,1}],horizontal_inertia);
%             % tmp = [data_drag{kFiles,iEvents,:}];
%             t1 = max([time_of_event{kFiles,iEvents}]-[time_data{kFiles,iEvents,1}]-1,1);
%             plot(t1,horizontal_inertia(t1+1),'r*','MarkerSize',10); % time of the event
%             
%         end
%     end
% end
% xlabel('time from the first event [frame]')
% ylabel('horizontal inertia force [N]');

%% inertia_part1
figure, hold on,
for kFiles = 1:size(inertia_force_part1,1)
    for iEvents = 1:size(inertia_force_part1,2)
        if ~isempty([time_of_event{kFiles,iEvents}]);
            tmp4 = cat(1,inertia_force_part1{kFiles,iEvents,:});
            norm_inertia_part1 = sum(tmp4.^2,2).^(0.5);
            plot([time_data{kFiles,iEvents,1:end}]-[time_data{kFiles,iEvents,1}],norm_inertia_part1);
            % tmp = [data_drag{kFiles,iEvents,:}];
            t1 = max([time_of_event{kFiles,iEvents}]-[time_data{kFiles,iEvents,1}]-1,1);
            plot(t1,norm_inertia_part1(t1+1),'r*','MarkerSize',10); % time of the event
            
        end
    end
end
xlabel('time from the first event [frame]')
ylabel('ma_p [N]');

% %% inertia_part1_vertical
% figure, hold on,
% for kFiles = 1:size(inertia_force_part1,1)
%     for iEvents = 1:size(inertia_force_part1,2)
%         if ~isempty([time_of_event{kFiles,iEvents}]);
%             tmp4 = cat(1,inertia_force_part1{kFiles,iEvents,:});
%             vertical_inertia_part1 = tmp4(:,2);
%             plot([time_data{kFiles,iEvents,1:end}]-[time_data{kFiles,iEvents,1}],vertical_inertia_part1);
%             % tmp = [data_drag{kFiles,iEvents,:}];
%             t1 = max([time_of_event{kFiles,iEvents}]-[time_data{kFiles,iEvents,1}]-1,1);
%             plot(t1,vertical_inertia_part1(t1+1),'r*','MarkerSize',10); % time of the event
%             
%         end
%     end
% end
% xlabel('time from the first event [frame]')
% ylabel('vertical ma_p [N]');
% %% inertia_part1_horizontal
% figure, hold on,
% for kFiles = 1:size(inertia_force_part1,1)
%     for iEvents = 1:size(inertia_force_part1,2)
%         if ~isempty([time_of_event{kFiles,iEvents}]);
%             tmp4 = cat(1,inertia_force_part1{kFiles,iEvents,:});
%             horizontal_inertia_part1 = sum(tmp4(:,[1 3]).^2,2).^(0.5);
%             plot([time_data{kFiles,iEvents,1:end}]-[time_data{kFiles,iEvents,1}],horizontal_inertia_part1);
%             % tmp = [data_drag{kFiles,iEvents,:}];
%             t1 = max([time_of_event{kFiles,iEvents}]-[time_data{kFiles,iEvents,1}]-1,1);
%             plot(t1,horizontal_inertia_part1(t1+1),'r*','MarkerSize',10); % time of the event
%             
%         end
%     end
% end
% xlabel('time from the first event [frame]')
% ylabel(' horizontal ma_p [N]');



%% inertia_part2
figure, hold on,
for kFiles = 1:size(inertia_force_part2,1)
    for iEvents = 1:size(inertia_force_part2,2)
        if ~isempty([time_of_event{kFiles,iEvents}]);
            tmp5 = cat(1,inertia_force_part2{kFiles,iEvents,:});
            norm_inertia_part2 = sum(tmp5.^2,2).^(0.5);
            plot([time_data{kFiles,iEvents,1:end}]-[time_data{kFiles,iEvents,1}],norm_inertia_part2);
            % tmp = [data_drag{kFiles,iEvents,:}];
            t1 = max([time_of_event{kFiles,iEvents}]-[time_data{kFiles,iEvents,1}]-1,1);
            plot(t1,norm_inertia_part2(t1+1),'r*','MarkerSize',10); % time of the event
            
        end
    end
end
xlabel('time from the first event [frame]')
ylabel('inertia (a_p-a_f) [N]');

% %% inertia_part2_vertical
% figure, hold on,
% for kFiles = 1:size(inertia_force_part2,1)
%     for iEvents = 1:size(inertia_force_part2,2)
%         if ~isempty([time_of_event{kFiles,iEvents}]);
%             tmp5 = cat(1,inertia_force_part2{kFiles,iEvents,:});
%             vertical_inertia_part2 = tmp5(:,2);
%             plot([time_data{kFiles,iEvents,1:end}]-[time_data{kFiles,iEvents,1}],vertical_inertia_part2);
%             % tmp = [data_drag{kFiles,iEvents,:}];
%             t1 = max([time_of_event{kFiles,iEvents}]-[time_data{kFiles,iEvents,1}]-1,1);
%             plot(t1,vertical_inertia_part2(t1+1),'r*','MarkerSize',10); % time of the event
%             
%         end
%     end
% end
% xlabel('time from the first event [frame]')
% ylabel('vertical inertia (a_p-a_f) [N]');
% 
% %% inertia_part2_horizontal
% figure, hold on,
% for kFiles = 1:size(inertia_force_part2,1)
%     for iEvents = 1:size(inertia_force_part2,2)
%         if ~isempty([time_of_event{kFiles,iEvents}]);
%             tmp5 = cat(1,inertia_force_part2{kFiles,iEvents,:});
%             horizontal_inertia_part2 = sum(tmp5(:,[1 3]).^2,2).^(0.5);
%             plot([time_data{kFiles,iEvents,1:end}]-[time_data{kFiles,iEvents,1}],horizontal_inertia_part2);
%             % tmp = [data_drag{kFiles,iEvents,:}];
%             t1 = max([time_of_event{kFiles,iEvents}]-[time_data{kFiles,iEvents,1}]-1,1);
%             plot(t1,horizontal_inertia_part2(t1+1),'r*','MarkerSize',10); % time of the event
%             
%         end
%     end
% end
% xlabel('time from the first event [frame]')
% ylabel('horizontal inertia (a_p-a_f) [N]');

%pressure
figure, hold on,
for kFiles = 1:size(pressure_force,1)
    for iEvents = 1:size(pressure_force,2)
        if ~isempty([time_of_event{kFiles,iEvents}]);
            tmp6 = cat(1,pressure_force{kFiles,iEvents,:});
            norm_pressure = sum(tmp6.^2,2).^(0.5);
            plot([time_data{kFiles,iEvents,1:end}]-[time_data{kFiles,iEvents,1}],norm_pressure);
            % tmp = [data_drag{kFiles,iEvents,:}];
            t1 = max([time_of_event{kFiles,iEvents}]-[time_data{kFiles,iEvents,1}]-1,1);
            plot(t1,norm_pressure(t1+1),'r*','MarkerSize',10); % time of the event
            
        end
    end
end
xlabel('time from the first event [frame]')
ylabel('pressure force [N]');

% %pressure_vertical
% figure, hold on,
% for kFiles = 1:size(pressure_force,1)
%     for iEvents = 1:size(pressure_force,2)
%         if ~isempty([time_of_event{kFiles,iEvents}]);
%             tmp6 = cat(1,pressure_force{kFiles,iEvents,:});
%             vertical_pressure = tmp6(:,2);
%             plot([time_data{kFiles,iEvents,1:end}]-[time_data{kFiles,iEvents,1}],vertical_pressure);
%             % tmp = [data_drag{kFiles,iEvents,:}];
%             t1 = max([time_of_event{kFiles,iEvents}]-[time_data{kFiles,iEvents,1}]-1,1);
%             plot(t1,vertical_pressure(t1+1),'r*','MarkerSize',10); % time of the event
%             
%         end
%     end
% end
% xlabel('time from the first event [frame]')
% ylabel('vertical pressure force [N]');
% 
% %pressure_horizontal
% figure, hold on,
% for kFiles = 1:size(pressure_force,1)
%     for iEvents = 1:size(pressure_force,2)
%         if ~isempty([time_of_event{kFiles,iEvents}]);
%             tmp6 = cat(1,pressure_force{kFiles,iEvents,:});
%            horizontal_pressure = sum(tmp6.^2,2).^(0.5);
%             plot([time_data{kFiles,iEvents,1:end}]-[time_data{kFiles,iEvents,1}],horizontal_pressure);
%             % tmp = [data_drag{kFiles,iEvents,:}];
%             t1 = max([time_of_event{kFiles,iEvents}]-[time_data{kFiles,iEvents,1}]-1,1);
%             plot(t1,horizontal_pressure(t1+1),'r*','MarkerSize',10); % time of the event
%             
%         end
%     end
% end
% xlabel('time from the first event [frame]')
% ylabel('horizontal pressure force [N]');

%%drag
figure, hold on,
for kFiles = 1:size(drag_force,1)
    for iEvents = 1:size(drag_force,2)
        if ~isempty([time_of_event{kFiles,iEvents}]);
            tmp7 = cat(1,drag_force{kFiles,iEvents,:});
            norm_drag = sum(tmp7.^2,2).^(0.5);
            plot([time_data{kFiles,iEvents,1:end}]-[time_data{kFiles,iEvents,1}],norm_drag);
            % tmp = [data_drag{kFiles,iEvents,:}];
            t1 = max([time_of_event{kFiles,iEvents}]-[time_data{kFiles,iEvents,1}]-1,1);
            plot(t1,norm_drag(t1+1),'r*','MarkerSize',10); % time of the event
            
        end
    end
end
xlabel('time from the first event [frame]')
ylabel('drag force [N]');

%lift
figure, hold on,
for kFiles = 1:size(lift_force,1)
    for iEvents = 1:size(lift_force,2)
        if ~isempty([time_of_event{kFiles,iEvents}]);
            tmp8 = cat(1,lift_force{kFiles,iEvents,:});
            norm_lift = sum(tmp8.^2,2).^(0.5);
            plot([time_data{kFiles,iEvents,1:end}]-[time_data{kFiles,iEvents,1}],norm_lift);
            % tmp = [data_drag{kFiles,iEvents,:}];
            t1 = max([time_of_event{kFiles,iEvents}]-[time_data{kFiles,iEvents,1}]-1,1);
            plot(t1,norm_lift(t1+1),'r*','MarkerSize',10); % time of the event
            
        end
    end
end
xlabel('time from the first event [frame]')
ylabel('lift force [N]');

%%drag vs lift force
figure, hold on,
for kFiles = 1:size(drag_force,1)
    for iEvents = 1:size(drag_force,2)
        if ~isempty([time_of_event{kFiles,iEvents}]);
            tmp1 = cat(1,drag_force{kFiles,iEvents,:});
            norm_drag = sum(tmp1.^2,2).^(0.5);
            tmp2 = cat(1,lift_force{kFiles,iEvents,:});
            norm_lift = sum(tmp2.^2,2).^(0.5);
            t1 = max([time_of_event{kFiles,iEvents}]-[time_data{kFiles,iEvents,1}]-1,1);
            
            scatter(norm_drag(t1+1),norm_lift(t1+1),100);
            % scatter(tmp1,tmp2);
            
        end
    end
end
xlabel('Drag force')
ylabel('Lift force');


%drag vs urev^2 
figure, hold on,
for kFiles = 1:size(drag_force,1)
    for iEvents = 1:size(drag_force,2)
        if ~isempty([time_of_event{kFiles,iEvents}]);
            tmp1 = cat(1,drag_Urev{kFiles,iEvents,:});
            norm_drag_Urev = sum(tmp1.^2,2).^(0.5);
            tmp2 = cat(1,drag_force{kFiles,iEvents,:});
            norm_drag = sum(tmp2.^2,2).^(0.5);
            t1 = max([time_of_event{kFiles,iEvents}]-[time_data{kFiles,iEvents,1}]-1,1);
            
            scatter(norm_drag_Urev(t1+1),norm_drag(t1+1),100);
            %scatter(tmp1,tmp2);
            
        end
    end
end
xlabel('0.5*p_l*A*u_r^2')
ylabel('drag force');


% lift vs urev^2 
figure, hold on,
for kFiles = 1:size(lift_force,1)
    for iEvents = 1:size(lift_force,2)
        if ~isempty([time_of_event{kFiles,iEvents}]);
            tmp1 = cat(1,drag_Urev{kFiles,iEvents,:});
            norm_drag_Urev = sum(tmp1.^2,2).^(0.5);
            tmp2 = cat(1,lift_force{kFiles,iEvents,:});
            norm_drag = sum(tmp2.^2,2).^(0.5);
            t1 = max([time_of_event{kFiles,iEvents}]-[time_data{kFiles,iEvents,1}]-1,1);
            
            scatter(norm_drag_Urev(t1+1),norm_drag(t1+1),100);
            %scatter(tmp1,tmp2);
            
        end
    end
end
xlabel('0.5*p_l*A*u_r^2')
ylabel('lift force');


%drag/urev^2 Vs Re_relative
figure, hold on,
for kFiles = 1:size(drag_force,1)
    for iEvents = 1:size(drag_force,2)
        if ~isempty([time_of_event{kFiles,iEvents}]);
            tmp1 = cat(1,Re{kFiles,iEvents,:});
            norm_Re = sum(tmp1.^2,2).^(0.5);
            tmp2 = cat(1,drag_Vs_Urev{kFiles,iEvents,:});
            norm_drag_Vs_Urev = sum(tmp2.^2,2).^(0.5);
            t1 = max([time_of_event{kFiles,iEvents}]-[time_data{kFiles,iEvents,1}]-1,1);
            
            scatter( norm_Re(t1+1),norm_drag_Vs_Urev(t1+1),100);
            %scatter(tmp1,tmp2);
            
        end
    end
end
xlabel('Re_r');
ylabel('drag force/0.5*p_l*A*u_r^2')


%lift/urev^2 Vs Re_relative
figure, hold on,
for kFiles = 1:size(lift_force,1)
    for iEvents = 1:size(lift_force,2)
        if ~isempty([time_of_event{kFiles,iEvents}]);
            tmp1 = cat(1,Re{kFiles,iEvents,:});
            norm_Re = sum(tmp1.^2,2).^(0.5);
            tmp2 = cat(1,lift_Vs_Urev{kFiles,iEvents,:});
            norm_lift_Vs_Urev = sum(tmp2.^2,2).^(0.5);
            t1 = max([time_of_event{kFiles,iEvents}]-[time_data{kFiles,iEvents,1}]-1,1);
            
            scatter( norm_Re(t1+1),norm_lift_Vs_Urev(t1+1),100);
            %scatter(tmp1,tmp2);
            
        end
    end
end
xlabel('Re_r');
ylabel('lift force/0.5*p_l*A*u_r^2')


