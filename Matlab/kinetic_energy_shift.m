load vertical_velocity_110_25_10.mat


figure, hold on

k = 0;

for kFiles= 1:size(data_small_2,1);
   for iEvents=1:size(data_small_2,2);  
    if ~isempty([time_of_event{kFiles,iEvents}]);
        
       tmp = cat(1,data_small_2{kFiles,iEvents,:});
      
        t1 = max([time_of_event{kFiles,iEvents}]-[time_data_small{kFiles,iEvents,1}]-1,1);
        k  = k + 1;
        
        xd{k} =[time_data_small{kFiles,iEvents,1:end}]-[time_data_small{kFiles,iEvents,1}]-t1;
       
        
        yd{k} = tmp;
        
        plot(xd{k},yd{k},'Color',[.5 .5 .5]);
        plot(t1-t1,tmp(t1+1),'ro','MarkerSize',10,'MarkerFaceColor','r'); % time of the event 
        xlabel('time [frame]')
        ylabel('(<v_t>^2)^0^.^5/w_s');
    end
  end
end






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

