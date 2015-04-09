clc, clear all

dist = @(x,y) sqrt((x(1,:)-y(1)).^2 + (x(2,:)-y(2)).^2 + (x(3,:)-y(3)).^2);

if ismac
    matdirectory = '\Users\alex\Dropbox\resuspension\Matlab\trajectories_alex_2012(186-194)\events_80';
else
    matdirectory ='C:\Users\dell\Dropbox\resuspension\Matlab\trajectories_alex_2012(186-194)\events_80';
end

large_ones = dir(fullfile(matdirectory,'large*'));

[   time_data, time_of_event,data_acceleration_p,error_velocity_2...
    Xp,Yp,Zp, Up,Vp,Wp,data_velocity_f1,std_velocity_2,error_acceleration_2...
    data_acceleration_f1,data_velocity_p,std_acceleration_2,Axp,Ayp,Azp, acceleration_p_f...
    data_acceleration_p_f_1,acceleration_pf2,std_acceleration_pf2,error_acceleration_pf2...
     basset_force2,data_velocity_p_f1,basset_force3,data_basset_force4,basset_force4, basset, basset1...
     std_basset1,error_basset1,data_velocity_p_f,data_velocity_pf1, data_velocity_p_f2] = deal({});



R = 25; %mm

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
            
            Axp{kFiles,iEvents,ind} = axp;
            Ayp{kFiles,iEvents,ind} = ayp;
            Azp{kFiles,iEvents,ind} = azp;
            
            velocity_p = [up,vp,wp]; %[m/s]
            acceleration_p = [axp,ayp,azp];% [m/s^2]
            
           
            acceleration_f1 = zeros(1,3);
            velocity_f1 = zeros(1,3);
            acceleration_p_f_1 = zeros(1,3);
            velocity_p_f1 = zeros(1,3);
           
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
                velocty_f1= velocity_f1 + velocity_f;
                
                velocity_p_f = [up-uf,vp-vf,wp-wf];
                velocity_p_f1= velocity_p_f1 + velocity_p_f;
               
                acceleration_f = [axf,ayf,azf]; %[m/s^2]
                acceleration_f1= acceleration_f1+ acceleration_f;
              
                acceleration_p_f = [axp-axf,ayp-ayf,azp-azf]; %[m/s^2]
                acceleration_p_f_1 = acceleration_p_f_1 + acceleration_p_f; %[m/s^2]
                %data_velocity_p_f{kFiles,iEvents,ind,j} = velocity_p_f;

            end % smalltraj
             
            
          
            if counter > 0
                % Averaging in the sphere of radius R
                data_velocity_f1{kFiles,iEvents,ind} = velocity_f1/counter;
                data_velocity_p_f1{kFiles,iEvents,ind} = velocity_p_f1/counter;
                data_acceleration_f1{kFiles,iEvents,ind} = acceleration_f1/counter;
                
                data_velocity_p_f2{kFiles,iEvents,ind} = velocity_p-data_velocity_f1{kFiles,iEvents,ind};
                
                data_acceleration_p_f_1{kFiles,iEvents,ind} = acceleration_p_f_1/counter;
                
                data_velocity_p_f1{kFiles,iEvents,ind} = velocity_p_f1/counter;
                     
            else
                data_velocity_f1{kFiles,iEvents,ind} = zeros(1,3);
               data_velocity_p_f1{kFiles,iEvents,ind} = zeros(1,3);
                data_acceleration_f1{kFiles,iEvents,ind} = zeros(1,3);
                
                data_acceleration_p_f_1{kFiles,iEvents,ind} = zeros(1,3);
                
                data_velocity_p_f2{kFiles,iEvents,ind}=zeros(1,3);
                
                data_velocity_p_f1{kFiles,iEvents,ind} = zeros(1,3);
            end
            
    
           
            
          %%%%%%%%%%%%%%%%%%%%%%%%
         %standard deviation%
            acceleration_2 = zeros(1,3);
            velocity_2 = zeros(1,3);
            acceleration_pf2 = zeros(1,3);
            velocity_pf1=zeros(1,3);
            %counter = 0;
            
            for j = 1:numSmallTraj % for all small ones
                
                distance = dist([small(j).xf,small(j).yf,small(j).zf]',[xp,yp,zp])';
                
                relevant = find(small(j).t == data.t(ind) &  distance <= R);
                
                if isempty(relevant), continue, end
                
                tmp = small(j);
                
                if tmp.yf(relevant) <= -21.5 && abs(tmp.yf(relevant)-yp )>= 10  ,continue, end
               % counter = counter + 1;
                
                uf = tmp.uf(relevant)/1000; %[m/s]
                vf = tmp.vf(relevant)/1000;
                wf = tmp.wf(relevant)/1000;
                
                velocity_f = [uf,vf,wf];
                 
                velocity_1= (velocity_f - data_velocity_f1{kFiles,iEvents,ind}).^2;
                velocity_2= velocity_2+velocity_1;
                
                velocity_p_f = [up-uf,vp-vf,wp-wf];
                velocity_pf= (velocity_p_f - data_velocity_p_f1{kFiles,iEvents,ind});
                velocity_pf1= velocity_pf + velocity_pf1;
                
                axf = tmp.axf(relevant)/1000;
                ayf= tmp.ayf(relevant)/1000;
                azf= tmp.azf(relevant)/1000;
               
                acceleration_f = [axf,ayf,azf]; %[m/s^2]
                acceleration_p_f = [axp-axf,ayp-ayf,azp-azf]; %[m/s^2]
                
                acceleration_1 = (acceleration_f - data_acceleration_f1{kFiles,iEvents,ind}).^2;
                acceleration_2 = acceleration_2 + acceleration_1;
                
               acceleration_pf1= (acceleration_p_f - data_acceleration_p_f_1{kFiles,iEvents,ind}).^2;
               acceleration_pf2= acceleration_pf1 + acceleration_pf2;
               
              
               
            end % smalltraj 
            
             if counter > 0
                % Averaging in the sphere of radius R
                std_velocity_2{kFiles,iEvents,ind} = sqrt(velocity_2/(counter-1));
                error_velocity_2{kFiles,iEvents,ind}= std_velocity_2{kFiles,iEvents,ind}/sqrt(counter);
                std_acceleration_2{kFiles,iEvents,ind} =sqrt(acceleration_2/(counter-1));
                error_acceleration_2{kFiles,iEvents,ind}= std_acceleration_2{kFiles,iEvents,ind}/sqrt(counter);
                
               std_acceleration_pf2{kFiles,iEvents,ind} =sqrt(acceleration_pf2/(counter-1));
               error_acceleration_pf2{kFiles,iEvents,ind}= std_acceleration_pf2{kFiles,iEvents,ind}/sqrt(counter); 
                 
               data_velocity_pf1{kFiles,iEvents,ind}=velocity_pf1/counter;
               
            else
                std_velocity_2{kFiles,iEvents,ind} = zeros(1,3);
                error_velocity_2{kFiles,iEvents,ind}= zeros(1,3);
                std_acceleration_2{kFiles,iEvents,ind} = zeros(1,3); 
                error_acceleration_2{kFiles,iEvents,ind}= zeros(1,3);
                
                std_acceleration_pf2{kFiles,iEvents,ind} = zeros(1,3); 
                error_acceleration_pf2{kFiles,iEvents,ind}= zeros(1,3);
                
                data_velocity_pf1{kFiles,iEvents,ind}=zeros(1,3);
             end
             
             data_velocity_p{kFiles,iEvents,ind} = velocity_p;
            data_acceleration_p{kFiles,iEvents,ind} = acceleration_p;
            time_data{kFiles, iEvents,ind} = data.t(ind);
            
        end %ind
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%
           % basset error _part1
           
            dt = 1/160;
            tau = 6:-1:1;
            k_v=1E-6;

            density_l = 0.998 * 1000; % kg/m^3 ; % density of water ,temperature= 22
            D = 550E-6; %diameter of the solid particle [m]
            r = D/2;
            er_p=1E-8;
           

            for i = 1:length(data.t)
            %basset_force1{kFiles,iEvents,i} = zeros(1,3);    
            basset_force2{kFiles,iEvents,i} = zeros(1,3);
            basset_force3{kFiles,iEvents,i} = zeros(1,3);
            end
            
        
        for ind = 6:length(data.t)
            
            urel = [data_velocity_p_f1{kFiles,iEvents,ind-5:ind}]; % 6 x [U,V,W]
            
            urel_1 = [data_velocity_p_f2{kFiles,iEvents,ind-5:ind}]; % 6 x [U,V,W] 
            
            for k = 1:3
                % dt in seconds, urel in m/sec
                du = gradient(urel(k:3:end),dt);
                
                du1 = gradient(urel_1(k:3:end),dt);
                
                y1=du./sqrt(tau);
                y2 = du./ sqrt(tau * dt);
                y3=du1./sqrt(tau*dt);
                
                tmp1=trapz(y1);
                %tmp2 =  12* r*er_p*sqrt(pi*k_v)* density_l * dt * trapz(y2);
                tmp2= 6*r^2 * sqrt(pi*k_v)* density_l*trapz(y2);
                tmp3= 6*r^2 * sqrt(pi*k_v)* density_l*trapz(y3);
                
                %basset_force1{kFiles,iEvents,ind}(k) = tmp1;
                basset_force2{kFiles,iEvents,ind}(k) = tmp2;
                basset_force3{kFiles,iEvents,ind}(k)= tmp3;
                
            end % k
  
        end %basset/ind
end
end

    save measurments_uncertainty_80
    
%   %U_p  
% figure, hold on,
% for kFiles = 1:size(data_velocity_p,1)
%     for iEvents = 1:size(data_velocity_p,2)
%         if ~isempty([time_of_event{kFiles,iEvents}]);
%             tmp = cat(1,Up{kFiles,iEvents,:});
%             norm_Up = sum(tmp.^2,2).^(0.5);
%             plot([time_data{kFiles,iEvents,1:end}]-[time_data{kFiles,iEvents,1}],norm_Up);
%            
%         end
%     end
% end
% xlabel('time from the first event [frame]')
% ylabel('data U_p [m/s]');
% 
% %V_p
% figure, hold on,
% for kFiles = 1:size(data_velocity_p,1)
%     for iEvents = 1:size(data_velocity_p,2)
%         if ~isempty([time_of_event{kFiles,iEvents}]);
%             tmp = cat(1,Vp{kFiles,iEvents,:});
%             norm_Vp = sum(tmp.^2,2).^(0.5);
%             plot([time_data{kFiles,iEvents,1:end}]-[time_data{kFiles,iEvents,1}],norm_Vp);
%            
%         end
%     end
% end
% xlabel('time from the first event [frame]')
% ylabel('data V_p [m/s]');
% 
% %W_p
% figure, hold on,
% for kFiles = 1:size(data_velocity_p,1)
%     for iEvents = 1:size(data_velocity_p,2)
%         if ~isempty([time_of_event{kFiles,iEvents}]);
%             tmp = cat(1,Wp{kFiles,iEvents,:});
%             norm_Wp = sum(tmp.^2,2).^(0.5);
%             plot([time_data{kFiles,iEvents,1:end}]-[time_data{kFiles,iEvents,1}],norm_Wp);
%            
%         end
%     end
% end
% xlabel('time from the first event [frame]')
% ylabel('data W_p [m/s]');
% 
% 
% %velocity_p
% figure, hold on,
% for kFiles = 1:size(data_velocity_p,1)
%     for iEvents = 1:size(data_velocity_p,2)
%         if ~isempty([time_of_event{kFiles,iEvents}]);
%             tmp = cat(1,data_velocity_p{kFiles,iEvents,:});
%             norm_data_velocity_p= sum(tmp.^2,2).^(0.5);
%             plot([time_data{kFiles,iEvents,1:end}]-[time_data{kFiles,iEvents,1}],norm_data_velocity_p);
%            
%         end
%     end
% end
% xlabel('time from the first event [frame]')
% ylabel('data velocity_p [m/s]');
% 
% %axp  
% figure, hold on,
% for kFiles = 1:size(data_velocity_p,1)
%     for iEvents = 1:size(data_velocity_p,2)
%         if ~isempty([time_of_event{kFiles,iEvents}]);
%             tmp = cat(1,Axp{kFiles,iEvents,:});
%             norm_Axp = sum(tmp.^2,2).^(0.5);
%             plot([time_data{kFiles,iEvents,1:end}]-[time_data{kFiles,iEvents,1}],norm_Axp);
%            
%         end
%     end
% end
% xlabel('time from the first event [frame]')
% ylabel('data ax_p [m/s]');
% 
% %ayp
% figure, hold on,
% for kFiles = 1:size(data_velocity_p,1)
%     for iEvents = 1:size(data_velocity_p,2)
%         if ~isempty([time_of_event{kFiles,iEvents}]);
%             tmp = cat(1,Ayp{kFiles,iEvents,:});
%             norm_Ayp = sum(tmp.^2,2).^(0.5);
%             plot([time_data{kFiles,iEvents,1:end}]-[time_data{kFiles,iEvents,1}],norm_Ayp);
%            
%         end
%     end
% end
% xlabel('time from the first event [frame]')
% ylabel('data ay_p [m/s]');
% 
% %azp
% figure, hold on,
% for kFiles = 1:size(data_velocity_p,1)
%     for iEvents = 1:size(data_velocity_p,2)
%         if ~isempty([time_of_event{kFiles,iEvents}]);
%             tmp = cat(1,Azp{kFiles,iEvents,:});
%             norm_Azp= sum(tmp.^2,2).^(0.5);
%             plot([time_data{kFiles,iEvents,1:end}]-[time_data{kFiles,iEvents,1}],norm_Azp);
%            
%         end
%     end
% end
% xlabel('time from the first event [frame]')
% ylabel('data az_p [m/s]');
% 
% %acceleration_p
% figure, hold on,
% for kFiles = 1:size(data_acceleration_p,1)
%     for iEvents = 1:size(data_acceleration_p,2)
%         if ~isempty([time_of_event{kFiles,iEvents}]);
%             tmp = cat(1,data_acceleration_p{kFiles,iEvents,:});
%             norm_data_acceleration_p= sum(tmp.^2,2).^(0.5);
%             plot([time_data{kFiles,iEvents,1:end}]-[time_data{kFiles,iEvents,1}],norm_data_acceleration_p);
%            
%         end
%     end
% end
% xlabel('time from the first event [frame]')
% ylabel('data acceleration_p [m/s]');
% 
% 
% %velocity_f-average of velocity tracers
% figure, hold on,
% for kFiles = 1:size(data_velocity_f1,1)
%     for iEvents = 1:size(data_velocity_f1,2)
%         if ~isempty([time_of_event{kFiles,iEvents}]);
%             tmp = cat(1,data_velocity_f1{kFiles,iEvents,:});
%             norm_data_velocity_f1 = sum(tmp.^2,2).^(0.5);
%             plot([time_data{kFiles,iEvents,1:end}]-[time_data{kFiles,iEvents,1}],norm_data_velocity_f1);  
%         end
%     end
% end
% xlabel('time from the first event [frame]')
% ylabel('velocity_f [m/s]');
% 
% %error_velocity-f
% figure, hold on,
% for kFiles = 1:size( error_velocity_2,1)
%     for iEvents = 1:size( error_velocity_2,2)
%         if ~isempty([time_of_event{kFiles,iEvents}]);
%             tmp = cat(1, error_velocity_2{kFiles,iEvents,:});  
%             norm_error_velocity_2 = sum(tmp.^2,2).^(0.5);
%             plot([time_data{kFiles,iEvents,1:end}]-[time_data{kFiles,iEvents,1}],norm_error_velocity_2);  
%         end
%     end
% end
% xlabel('time from the first event [frame]')
% ylabel(' error velocity_f [m/s]');
% 
% 
% 
% %acceleration_f- average of acceleration tracers
% figure, hold on,
% for kFiles = 1:size(data_acceleration_f1,1)
%     for iEvents = 1:size(data_acceleration_f1,2)
%         if ~isempty([time_of_event{kFiles,iEvents}]);
%             tmp = cat(1,data_acceleration_f1{kFiles,iEvents,:});
%             norm_data_acceleration_f1 = sum(tmp.^2,2).^(0.5);
%             plot([time_data{kFiles,iEvents,1:end}]-[time_data{kFiles,iEvents,1}],norm_data_acceleration_f1);  
%         end
%     end
% end
% 
% xlabel('time from the first event [frame]')
% ylabel('acceleration_f [m/s^2]');
% 
% %error acceleration_f
% figure, hold on,
% for kFiles = 1:size(error_acceleration_2,1)
%     for iEvents = 1:size(error_acceleration_2,2)
%         if ~isempty([time_of_event{kFiles,iEvents}]);
%             tmp = cat(1,error_acceleration_2{kFiles,iEvents,:});
%             norm_error_acceleration_2 = sum(tmp.^2,2).^(0.5);
%             plot([time_data{kFiles,iEvents,1:end}]-[time_data{kFiles,iEvents,1}],norm_error_acceleration_2);  
%         end
%     end
% end
% xlabel('time from the first event [frame]')
% ylabel(' error acceleration_f [m/s^2]');
% 
% %acceleration_pf
% figure, hold on,
% for kFiles = 1:size(data_acceleration_p_f_1,1)
%     for iEvents = 1:size(data_acceleration_p_f_1,2)
%         if ~isempty([time_of_event{kFiles,iEvents}]);
%             tmp = cat(1,data_acceleration_p_f_1{kFiles,iEvents,:});
%             norm_data_acceleration_p_f_1 = sum(tmp.^2,2).^(0.5);
%             plot([time_data{kFiles,iEvents,1:end}]-[time_data{kFiles,iEvents,1}],norm_data_acceleration_p_f_1 );  
%         end
%     end
% end
% xlabel('time from the first event [frame]')
% ylabel('  acceleration_pf [m/s^2]');
% 
% %error acceleration_f
% figure, hold on,
% for kFiles = 1:size(error_acceleration_pf2,1)
%     for iEvents = 1:size(error_acceleration_pf2,2)
%         if ~isempty([time_of_event{kFiles,iEvents}]);
%             tmp = cat(1,error_acceleration_pf2{kFiles,iEvents,:});
%             norm_error_acceleration_pf2 = sum(tmp.^2,2).^(0.5);
%             plot([time_data{kFiles,iEvents,1:end}]-[time_data{kFiles,iEvents,1}],norm_error_acceleration_pf2);  
%         end
%     end
% end
% xlabel('time from the first event [frame]')
% ylabel(' error acceleration_pf [m/s^2]');

%basset_force_2
figure, hold on,
for kFiles = 1:size(basset_force2,1)
    for iEvents = 1:size(basset_force2,2)
        if ~isempty([time_of_event{kFiles,iEvents}]);
            tmp2 = cat(1,basset_force2{kFiles,iEvents,:});
            norm_basset_force2 = sum(tmp2.^2,2).^(0.5); 
            plot([time_data{kFiles,iEvents,1:end}]-[time_data{kFiles,iEvents,1}],norm_basset_force2);  
        end
    end
end
xlabel('time from the first event [frame]')
ylabel('basset force2 [N]');

%basset_force_3
figure, hold on,
for kFiles = 1:size(basset_force3,1)
    for iEvents = 1:size(basset_force3,2)
        if ~isempty([time_of_event{kFiles,iEvents}]);
            tmp2 = cat(1,basset_force3{kFiles,iEvents,:});
            norm_basset_force3 = sum(tmp2.^2,2).^(0.5); 
            plot([time_data{kFiles,iEvents,1:end}]-[time_data{kFiles,iEvents,1}],norm_basset_force3);  
        end
    end
end
xlabel('time from the first event [frame]')
ylabel('basset force3 [N]');



% %error basset
% figure, hold on,
% for kFiles = 1:size(error_basset1,1)
%     for iEvents = 1:size(error_basset1,2)
%         if ~isempty([time_of_event{kFiles,iEvents}]);
%             tmp = cat(1,error_basset1{kFiles,iEvents,:});
%             norm_error_basset1 = sum(tmp.^2,2).^(0.5);
%             plot([time_data{kFiles,iEvents,6:end}]-[time_data{kFiles,iEvents,6}],norm_error_basset1);  
%         end
%     end
% end
% xlabel('time from the first event [frame]')
% ylabel(' error_basset1 [N]')
% 











    
            