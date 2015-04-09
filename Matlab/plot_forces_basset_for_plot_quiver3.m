%% quiver of the large particle trajectories with
% their own velocity in red
% total force in blue
% relative velocity in green

% Data is produced by forces_basset_v2

load forces_R25_10_90.mat

for kFiles = 1:size(total_force,1)
    for iEvents = 1:size(total_force,2)
        
        x = [Xp{kFiles,iEvents,:}]; if isempty(x), continue, end
        
        figure, hold on; view(3); grid on; box on;
        
        z = [Zp{kFiles,iEvents,:}];
        y = [Yp{kFiles,iEvents,:}];
        u = [Up{kFiles,iEvents,:}];
        v = [Vp{kFiles,iEvents,:}];
        w = [Wp{kFiles,iEvents,:}];
        
        q = mean((u.^2 + v.^2 + w.^2).^0.5);
        % red is the particle velocity
        quiver3(x,z,y,u/q,w/q,v/q,1.2,'r');
        title('red = Up, blue = total force, magenta = Urel, green=Uf')
        
        f = cat(1,total_force{kFiles,iEvents,:});
        nf = mean(sum(f.^2,2).^(0.5));
        
        
        % if the data is not created by forces_basset_v2
        % then the length of x,y,z is different from f
        % we shift it by 6 due to basset_force, check the other file
        if length(x) ~= length(f)
            x1 = x(6:end);
            y1 = y(6:end);
            z1 = z(6:end);
        else
            x1 = x;
            y1 = y;
            z1 = z;
        end
        
        
        quiver3(x1',z1',y1',f(:,1)/nf,f(:,3)/nf,f(:,2)/nf,1.2,'b');
        
        f = cat(1,data_velocity_f_p1{kFiles,iEvents,:}); 
        %nf = mean(sum(f.^2,2).^(0.5));
        quiver3(x',z',y',f(:,1)/nf,f(:,3)/nf,f(:,2)/nf,1.2,'m');
        
        f = cat(1,data_velocity_f1{kFiles,iEvents,:}); 
        %nf = mean(sum(f.^2,2).^(0.5));
        quiver3(x',z',y',f(:,1)/nf,f(:,3)/nf,f(:,2)/nf,1.2,'g');
        
        event = find([time_data{kFiles,iEvents,:}] == time_of_event{kFiles,iEvents});
        
        plot3(x(event),z(event),y(event),'k*','MarkerSize',10);
    end
end


%%
for kFiles = 1:size(total_force,1)
    for iEvents = 1:size(total_force,2)
        
        x = [Xp{kFiles,iEvents,:}]; if isempty(x), continue, end
        
        figure, hold on; view(3); grid on; box on;
        title('red = total, green = drag, blue = lift')
        
        z = [Zp{kFiles,iEvents,:}];
        y = [Yp{kFiles,iEvents,:}];
        u = [Up{kFiles,iEvents,:}];
        v = [Vp{kFiles,iEvents,:}];
        w = [Wp{kFiles,iEvents,:}];
        
        q = mean((u.^2 + v.^2 + w.^2).^0.5);
        
        % quiver3(x,z,y,u/q,w/q,v/q,1.2,'r');
        
       if length(x) ~= length(f)
            x1 = x(6:end);
            y1 = y(6:end);
            z1 = z(6:end);
        else
            x1 = x;
            y1 = y;
            z1 = z;
        end
        
        
        
        
        f = cat(1,total_force{kFiles,iEvents,:});
        nf = mean(sum(f.^2,2).^(0.5));
        
        quiver3(x1',z1',y1',f(:,1)/nf,f(:,3)/nf,f(:,2)/nf,1.2,'r');
        
        
        
        f = cat(1,drag_force{kFiles,iEvents,:});
        % nf = mean(sum(f.^2,2).^(0.5));
        
        quiver3(x1',z1',y1',f(:,1)/nf,f(:,3)/nf,f(:,2)/nf,1.2,'g');
        
        f = cat(1,lift_force{kFiles,iEvents,:});
        % nf = mean(sum(f.^2,2).^(0.5));
        
        quiver3(x1',z1',y1',f(:,1)/nf,f(:,3)/nf,f(:,2)/nf,1.2,'b');
        
        
        
        event = find([time_data{kFiles,iEvents,:}] == time_of_event{kFiles,iEvents});
        
        plot3(x(event),z(event),y(event),'k*','MarkerSize',10);
    end
end

%%
%%
for kFiles = 1:size(total_force,1)
    for iEvents = 1:size(total_force,2)
        
        x = [Xp{kFiles,iEvents,:}]; if isempty(x), continue, end
        
        figure, hold on; view(3); grid on; box on;
        title('red = total, green = pressure, blue = basset, magenta = inertia')
        
        z = [Zp{kFiles,iEvents,:}];
        y = [Yp{kFiles,iEvents,:}];
        u = [Up{kFiles,iEvents,:}];
        v = [Vp{kFiles,iEvents,:}];
        w = [Wp{kFiles,iEvents,:}];
        
        
        
        
        q = mean((u.^2 + v.^2 + w.^2).^0.5);
        
        % quiver3(x,z,y,u/q,w/q,v/q,1.2,'r');
        
        
        f = cat(1,total_force{kFiles,iEvents,:});
        nf = mean(sum(f.^2,2).^(0.5));
        
        
        
        if length(x) ~= length(f)
            x1 = x(6:end);
            y1 = y(6:end);
            z1 = z(6:end);
        else
            x1 = x;
            y1 = y;
            z1 = z;
        end
        
        quiver3(x1',z1',y1',f(:,1)/nf,f(:,3)/nf,f(:,2)/nf,1.2,'r');
        
        
        
        f = cat(1,pressure_force{kFiles,iEvents,:});
        % nf = mean(sum(f.^2,2).^(0.5));
        quiver3(x1',z1',y1',f(:,1)/nf,f(:,3)/nf,f(:,2)/nf,1.2,'g');
        
        try
            f = cat(1,basset_force1{kFiles,iEvents,6:end});
            % nf = mean(sum(f.^2,2).^(0.5));
            
            quiver3(x1',z1',y1',f(:,1)/nf,f(:,3)/nf,f(:,2)/nf,1.2,'b');
        catch
            ;
        end
        
        f = cat(1,inertia_force{kFiles,iEvents,:});
        quiver3(x1',z1',y1',f(:,1)/nf,f(:,3)/nf,f(:,2)/nf,1.2,'m');
        
        
        event = find([time_data{kFiles,iEvents,:}] == time_of_event{kFiles,iEvents});
        
        plot3(x(event),z(event),y(event),'k*','MarkerSize',10,'LineWidth',2);
    end
end



% 
% %% Plot other forces vs time
% plot_results = false;
% 
% if plot_results
%     
%     
%     
%     %%total_force
%     figure, hold on,
%     for kFiles = 1:size(total_force,1)
%         for iEvents = 1:size(total_force,2)
%             if ~isempty([time_of_event{kFiles,iEvents}]);
%                 tmp = cat(1,total_force{kFiles,iEvents,:});
%                 norm_total_force = sum(tmp.^2,2).^(0.5);
%                 plot([time_data{kFiles,iEvents,:}]-[time_data{kFiles,iEvents,1}],norm_total_force);
%                 % tmp = [data_drag{kFiles,iEvents,:}];
%                 t1 = max([time_of_event{kFiles,iEvents}]-[time_data{kFiles,iEvents,1}]-1,1);
%                 plot(t1,norm_total_force(t1+1),'r*','MarkerSize',10); % time of the event
%                 
%             end
%         end
%     end
%     xlabel('time from the first event [frame]')
%     ylabel('total force [N]');
%     
%     
%     %basset_force
%     figure, hold on,
%     for kFiles = 1:size(basset_force1,1)
%         for iEvents = 1:size(basset_force1,2)
%             if ~isempty([time_of_event{kFiles,iEvents}]);
%                 tmp2 = cat(1,basset_force1{kFiles,iEvents,:});
%                 if isempty(tmp2), continue, end
%                 norm_basset_force1 = sum(tmp2.^2,2).^(0.5);
%                 plot([time_data{kFiles,iEvents,:}]-[time_data{kFiles,iEvents,1}],norm_basset_force1);
%                 t1 = max([time_of_event{kFiles,iEvents}]-[time_data{kFiles,iEvents,1}]-1,1);
%                 plot(t1,norm_basset_force1(t1+1),'r*','MarkerSize',10); % time of the event
%                 
%             end
%         end
%     end
%     xlabel('time from the first event [frame]')
%     ylabel('basset force [N]');
%     
%     
%     
%     %%drag
%     figure, hold on,
%     for kFiles = 1:size(drag_force,1)
%         for iEvents = 1:size(drag_force,2)
%             if ~isempty([time_of_event{kFiles,iEvents}]);
%                 tmp = cat(1,drag_force{kFiles,iEvents,:});
%                 norm_drag = sum(tmp.^2,2).^(0.5);
%                 plot([time_data{kFiles,iEvents,:}]-[time_data{kFiles,iEvents,1}],norm_drag);
%                 % tmp = [data_drag{kFiles,iEvents,:}];
%                 t1 = max([time_of_event{kFiles,iEvents}]-[time_data{kFiles,iEvents,1}]-1,1);
%                 plot(t1,norm_drag(t1+1),'r*','MarkerSize',10); % time of the event
%                 
%             end
%         end
%     end
%     xlabel('time from the first event [frame]')
%     ylabel('drag force [N]');
%     
%     %lift
%     figure, hold on,
%     for kFiles = 1:size(lift_force,1)
%         for iEvents = 1:size(lift_force,2)
%             if ~isempty([time_of_event{kFiles,iEvents}]);
%                 tmp = cat(1,lift_force{kFiles,iEvents,:});
%                 norm_lift = sum(tmp.^2,2).^(0.5);
%                 plot([time_data{kFiles,iEvents,:}]-[time_data{kFiles,iEvents,1}],norm_lift);
%                 % tmp = [data_drag{kFiles,iEvents,:}];
%                 t1 = max([time_of_event{kFiles,iEvents}]-[time_data{kFiles,iEvents,1}]-1,1);
%                 plot(t1,norm_lift(t1+1),'r*','MarkerSize',10); % time of the event
%                 
%             end
%         end
%     end
%     xlabel('time from the first event [frame]')
%     ylabel('lift force [N]');
%     
%     %% inertia
%     figure, hold on,
%     for kFiles = 1:size(inertia_force,1)
%         for iEvents = 1:size(inertia_force,2)
%             if ~isempty([time_of_event{kFiles,iEvents}]);
%                 tmp = cat(1,inertia_force{kFiles,iEvents,:});
%                 norm_inertia = sum(tmp.^2,2).^(0.5);
%                 plot([time_data{kFiles,iEvents,:}]-[time_data{kFiles,iEvents,1}],norm_inertia);
%                 % tmp = [data_drag{kFiles,iEvents,:}];
%                 t1 = max([time_of_event{kFiles,iEvents}]-[time_data{kFiles,iEvents,1}]-1,1);
%                 plot(t1,norm_inertia(t1+1),'r*','MarkerSize',10); % time of the event
%                 
%             end
%         end
%     end
%     xlabel('time from the first event [frame]')
%     ylabel('inertia force [N]');
%     
%     %% inertia_part1
%     figure, hold on,
%     for kFiles = 1:size(inertia_force_part1,1)
%         for iEvents = 1:size(inertia_force_part1,2)
%             if ~isempty([time_of_event{kFiles,iEvents}]);
%                 tmp = cat(1,inertia_force_part1{kFiles,iEvents,:});
%                 norm_inertia_part1 = sum(tmp.^2,2).^(0.5);
%                 plot([time_data{kFiles,iEvents,:}]-[time_data{kFiles,iEvents,1}],norm_inertia_part1);
%                 % tmp = [data_drag{kFiles,iEvents,:}];
%                 t1 = max([time_of_event{kFiles,iEvents}]-[time_data{kFiles,iEvents,1}]-1,1);
%                 plot(t1,norm_inertia_part1(t1+1),'r*','MarkerSize',10); % time of the event
%                 
%             end
%         end
%     end
%     xlabel('time from the first event [frame]')
%     ylabel('ma_p [N]');
%     
%     %% inertia_part2
%     figure, hold on,
%     for kFiles = 1:size(inertia_force_part2,1)
%         for iEvents = 1:size(inertia_force_part2,2)
%             if ~isempty([time_of_event{kFiles,iEvents}]);
%                 tmp = cat(1,inertia_force_part2{kFiles,iEvents,:});
%                 norm_inertia_part2 = sum(tmp.^2,2).^(0.5);
%                 plot([time_data{kFiles,iEvents,:}]-[time_data{kFiles,iEvents,1}],norm_inertia_part2);
%                 % tmp = [data_drag{kFiles,iEvents,:}];
%                 t1 = max([time_of_event{kFiles,iEvents}]-[time_data{kFiles,iEvents,1}]-1,1);
%                 plot(t1,norm_inertia_part2(t1+1),'r*','MarkerSize',10); % time of the event
%                 
%             end
%         end
%     end
%     xlabel('time from the first event [frame]')
%     ylabel('inertia (a_p-a_f) [N]');
%     
%     %pressure
%     figure, hold on,
%     for kFiles = 1:size(pressure_force,1)
%         for iEvents = 1:size(pressure_force,2)
%             if ~isempty([time_of_event{kFiles,iEvents}]);
%                 tmp = cat(1,pressure_force{kFiles,iEvents,:});
%                 norm_pressure = sum(tmp.^2,2).^(0.5);
%                 plot([time_data{kFiles,iEvents,:}]-[time_data{kFiles,iEvents,1}],norm_pressure);
%                 % tmp = [data_drag{kFiles,iEvents,:}];
%                 t1 = max([time_of_event{kFiles,iEvents}]-[time_data{kFiles,iEvents,1}]-1,1);
%                 plot(t1,norm_pressure(t1+1),'r*','MarkerSize',10); % time of the event
%                 
%             end
%         end
%     end
%     xlabel('time from the first event [frame]')
%     ylabel('pressure force [N]');
%     
%     %%drag vs lift force
%     figure, hold on,
%     for kFiles = 1:size(drag_force,1)
%         for iEvents = 1:size(drag_force,2)
%             if ~isempty([time_of_event{kFiles,iEvents}]);
%                 tmp1 = cat(1,drag_force{kFiles,iEvents,:});
%                 norm_drag = sum(tmp1.^2,2).^(0.5);
%                 tmp2 = cat(1,lift_force{kFiles,iEvents,:});
%                 norm_lift = sum(tmp2.^2,2).^(0.5);
%                 t1 = max([time_of_event{kFiles,iEvents}]-[time_data{kFiles,iEvents,1}]-1,1);
%                 
%                 scatter(norm_drag(t1+1),norm_lift(t1+1),100);
%                 % scatter(tmp1,tmp2);
%                 
%             end
%         end
%     end
%     xlabel('Drag force')
%     ylabel('Lift force');
%     
%     
%     %drag vs urev^2
%     figure, hold on,
%     for kFiles = 1:size(drag_force,1)
%         for iEvents = 1:size(drag_force,2)
%             if ~isempty([time_of_event{kFiles,iEvents}]);
%                 tmp1 = cat(1,drag_Urev{kFiles,iEvents,:});
%                 norm_drag_Urev = sum(tmp1.^2,2).^(0.5);
%                 tmp2 = cat(1,drag_force{kFiles,iEvents,:});
%                 norm_drag = sum(tmp2.^2,2).^(0.5);
%                 t1 = max([time_of_event{kFiles,iEvents}]-[time_data{kFiles,iEvents,1}]-1,1);
%                 
%                 scatter(norm_drag_Urev(t1+1),norm_drag(t1+1),100);
%                 %scatter(tmp1,tmp2);
%                 
%             end
%         end
%     end
%     xlabel('0.5*p_l*A*u_r^2')
%     ylabel('drag force');
%     
%     
%     %drag/urev^2 Vs Re_relative
%     figure, hold on,
%     for kFiles = 1:size(drag_force,1)
%         for iEvents = 1:size(drag_force,2)
%             if ~isempty([time_of_event{kFiles,iEvents}]);
%                 tmp1 = cat(1,Re{kFiles,iEvents,:});
%                 norm_Re = sum(tmp1.^2,2).^(0.5);
%                 tmp2 = cat(1,drag_Vs_Urev{kFiles,iEvents,:});
%                 norm_drag_Vs_Urev = sum(tmp2.^2,2).^(0.5);
%                 t1 = max([time_of_event{kFiles,iEvents}]-[time_data{kFiles,iEvents,1}]-1,1);
%                 
%                 scatter( norm_Re(t1+1),norm_drag_Vs_Urev(t1+1),100);
%                 %scatter(tmp1,tmp2);
%                 
%             end
%         end
%     end
%     xlabel('Re_r');
%     ylabel('drag force/0.5*p_l*A*u_r^2')
%     
%     
%     %lift/urev^2 Vs Re_relative
%     figure, hold on,
%     for kFiles = 1:size(lift_force,1)
%         for iEvents = 1:size(lift_force,2)
%             if ~isempty([time_of_event{kFiles,iEvents}]);
%                 tmp1 = cat(1,Re{kFiles,iEvents,:});
%                 norm_Re = sum(tmp1.^2,2).^(0.5);
%                 tmp2 = cat(1,lift_Vs_Urev{kFiles,iEvents,:});
%                 norm_lift_Vs_Urev = sum(tmp2.^2,2).^(0.5);
%                 t1 = max([time_of_event{kFiles,iEvents}]-[time_data{kFiles,iEvents,1}]-1,1);
%                 
%                 scatter( norm_Re(t1+1),norm_lift_Vs_Urev(t1+1),100);
%                 %scatter(tmp1,tmp2);
%                 
%             end
%         end
%     end
%     xlabel('Re_r');
%     ylabel('lift force/0.5*p_l*A*u_r^2')
    
    
    
%end