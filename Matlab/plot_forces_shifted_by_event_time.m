 load forces_R25_10_80.mat

%%
% all forces togther_vertical

figure, hold on

k = 0;

for kFiles= 1:size(total_force,1);
   for iEvents=1:size(total_force,2);  
    if ~isempty([time_of_event{kFiles,iEvents}]);
       % figure, hold on, 
      tmp = cat(1,total_force{kFiles,iEvents,:});
%     vertical_total_force = tmp(:,2);
%      tmp = cat(1,basset_force1{kFiles,iEvents,6:end});
%      tmp = cat(1,inertia_force{kFiles,iEvents,:});
%      tmp = cat(1,pressure_force{kFiles,iEvents,:});
%      vertical_inertia = tmp(:,2);
%      tmp = cat(1,inertia_force_part1{kFiles,iEvents,:});
%      vertical_inertia_part1 = tmp(:,2);
%     tmp = cat(1,inertia_force_part2{kFiles,iEvents,:});
%      vertical_inertia_part2 = tmp(:,2);
%     tmp = cat(1,data_acceleration_p{kFiles,iEvents,:});
%      vertical_acceleration_p = tmp(:,2);
%       tmp = cat(1,drag_force{kFiles,iEvents,:});
%     tmp = cat(1,lift_force{kFiles,iEvents,:});
%      tmp = cat(1,data_acceleration_f1{kFiles,iEvents,:});
%      tmp = cat(1,data_acceleration_p_f1{kFiles,iEvents,:});

%        plot([time_data{kFiles,iEvents,1:end}]-[time_data{kFiles,iEvents,6}],vertical_basset_force1,'r');
%        plot([time_data{kFiles,iEvents,6:end}]-[time_data{kFiles,iEvents,6}],vertical_inertia,'g');
%         plot([time_data{kFiles,iEvents,6:end}]-[time_data{kFiles,iEvents,6}],vertical_inertia_part1,'y');
%         plot([time_data{kFiles,iEvents,6:end}]-[time_data{kFiles,iEvents,6}],vertical_inertia_part2,'c');
%        plot([time_data{kFiles,iEvents,6:end}]-[time_data{kFiles,iEvents,6}], vertical_pressure,'m');
%         plot([time_data{kFiles,iEvents,6:end}]-[time_data{kFiles,iEvents,6}], buoyancy,'c');

       % Vertical is tmp(:,2); or horizontal_x (tmp(:,1)),horizontal_z (tmp(:,3))
        %tmp = tmp(:,2).';
        %tmp = tmp(:,1).'; 
        %tmp = tmp(:,3).';
        %tmp = tmp(:,[1 3]).';
        tmp = sum(tmp(:,[1 3]),2).';

        t1 = max([time_of_event{kFiles,iEvents}]-[time_data{kFiles,iEvents,1}]-1,1);
        k  = k + 1;
        xd{k} = [time_data{kFiles,iEvents,6:end}]-[time_data{kFiles,iEvents,6}]-t1;
        %xd{k} =[time_data{kFiles,iEvents,1:end}]-[time_data{kFiles,iEvents,1}]-t1;%
        % use for acceleration particle and acceleration fluid
        
        yd{k} = tmp;
        
        plot(xd{k},yd{k},'Color',[.5 .5 .5]);
        plot(t1-t1,tmp(t1+1),'ro','MarkerSize',10,'MarkerFaceColor','r'); % time of the event 
%        plot(t1,vertical_basset_force1(t1+6),'r*','MarkerSize',10); 
%        plot(t1,vertical_inertia(t1+1),'r*','MarkerSize',10); 
%         plot(t1,vertical_inertia_part1(t1+1),'r*','MarkerSize',10);
%         plot(t1,vertical_inertia_part2(t1+1),'r*','MarkerSize',10); 
%        plot(t1,vertical_pressure(t1+1),'r*','MarkerSize',10); 
%        plot(t1,buoyancy(t1+1),'r*','MarkerSize',10);
        xlabel('time [frame]')
        ylabel(' horizontal force [N]');
    end
  end
end





% ch = findobj(gca,'type','line','color','b');
% yd = zeros(1,401);
% counter = ones(1,401); 
% for i = 1:length(ch)
%     tmp = get(ch(i),'ydata');
%     x = get(ch(i),'xdata');
%     ind = find(x == 0); 
%     ii = 200-ind+1:200+length(tmp)-ind;
%     yd(ii) = yd(ii) + tmp;
%     counter(ii) = counter(ii) + 1;
% end

% figure, hold on;
xx = -200:200;
yy = xx*0;
counter = yy + 1;

for k = 1:length(xd)
    tmp = interp1(xd{k},yd{k},xx,'nearest',0); 
    if ~any(isnan(tmp)),
         yy = yy + tmp;
         counter(find(tmp)) = counter(find(tmp)) + 1;
    end
end

plot(xx,yy./counter,'LineWidth',4,'Color','b');

return

