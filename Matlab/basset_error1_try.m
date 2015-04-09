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
     std_basset1,error_basset1, velocity_p_f2, data_velocity_p_f2] = deal({});



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
 
            axp = data.axf(ind)/1000; %[m/s^2]
            ayp = data.ayf(ind)/1000;
            azp = data.azf(ind)/1000;

            velocity_p = [up,vp,wp]; %[m/s]

            velocity_f1 = zeros(1,3);
            
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
                
                %velocity_p_f = [up-uf,vp-vf,wp-wf];
                %velocity_p_f1= velocity_p_f1 + velocity_p_f;
               

            end % smalltraj
          
            
            if counter > 0
                % Averaging in the sphere of radius R
                data_velocity_f1{kFiles,iEvents,ind} = velocity_f1/counter;
                %data_velocity_p_f1{kFiles,iEvents,ind} = velocity_p_f1/counter;

                %data_velocity_p_f1{kFiles,iEvents,ind} = velocity_p_f1/counter;
                     
            else
                data_velocity_f1{kFiles,iEvents,ind} = zeros(1,3);
               %data_velocity_p_f1{kFiles,iEvents,ind} = zeros(1,3);

               % data_velocity_p_f1{kFiles,iEvents,ind} = zeros(1,3);
            end
            
          
       
               
                velocity_p_f= (velocity_p - data_velocity_f1{kFiles,iEvents,ind});
               
        end % ind
      
         
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
            basset_force2{kFiles,iEvents,i} = zeros(1,3);
            basset_force3{kFiles,iEvents,i} = zeros(1,3);
            end
            
           for ind = 6:length(data.t)
            
%             xp = data.xf(ind);
%             yp = data.yf(ind);
%             zp = data.zf(ind);
%             
%             up = data.uf(ind)/1000; %[m/s]
%             vp = data.vf(ind)/1000;
%             wp = data.wf(ind)/1000;
%             
%             velocity_p = [up,vp,wp]; %[m/s]

%             velocity_f1 = zeros(1,3);
%             velocity_p_f1 = zeros(1,3);
%            
            counter = 0;
            
            for j = 1:numSmallTraj % for all small ones
                
                distance = dist([small(j).xf,small(j).yf,small(j).zf]',[xp,yp,zp])';
                
                relevant = find(small(j).t == data.t(ind) &  distance <= R);
                
                if isempty(relevant), continue, end
                
                tmp = small(j);
                
                if tmp.yf(relevant) <= -21.5 && abs(tmp.yf(relevant)-yp )>= 10  ,continue, end
                counter = counter + 1;
%                 
%                 uf = tmp.uf(relevant)/1000; %[m/s]
%                 vf = tmp.vf(relevant)/1000;
%                 wf = tmp.wf(relevant)/1000;
%                 
% 
%                 velocity_f = [uf,vf,wf];
% 
%                 velocity_p_f = [up-uf,vp-vf,wp-wf];
                %data_velocity_p_f2{kFiles,iEvents,ind-5:ind,j}= velocity_p_f ;
    
            
            urel = [data_velocity_p_f{kFiles,iEvents,ind-5:ind}]; % 6 x [U,V,W]
            %data_velocity_p_f{kFiles,iEvents,ind,j} = velocity_p_f1;
            %urel_1 = [data_velocity_p_f2{kFiles,iEvents,ind-5:ind,j}];
            
            for k = 1:3
                % dt in seconds, urel in m/sec
                du = gradient(urel(k:3:end),dt);
                %du1 = gradient(urel_1(k:3:end),dt);
                
               % y = du./ sqrt(tau * dt);
                %y1=du./sqrt(tau);
                %y2=du1./sqrt(tau);
                
              %  tmp =   trapz(y1);
                tmp1= trapz(y2);
                
                %basset_force2{kFiles,iEvents,ind}(k) = tmp;
                basset_force3{kFiles,iEvents,ind}(k)= tmp1;
                
            end % k
  

        
            end % smalltraj 
           end  % ind
    end % iEevent
end %Kfiles

    save measurments_uncertainty_80
    
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

% %basset_force_3
% figure, hold on,
% for kFiles = 1:size(basset_force3,1)
%     for iEvents = 1:size(basset_force3,2)
%         if ~isempty([time_of_event{kFiles,iEvents}]);
%             tmp2 = cat(1,basset_force3{kFiles,iEvents,:});
%             norm_basset_force3 = sum(tmp2.^2,2).^(0.5); 
%             plot([time_data{kFiles,iEvents,1:end}]-[time_data{kFiles,iEvents,1}],norm_basset_force3);  
%         end
%     end
% end
% xlabel('time from the first event [frame]')
% ylabel('basset force3 [N]');
% 
