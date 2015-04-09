clc,clear all

dist = @(x,y) sqrt((x(1,:)-y(1)).^2 + (x(2,:)-y(2)).^2 + (x(3,:)-y(3)).^2);

if ismac
    matdirectory = '/Users/alex/Dropbox/resuspension/2011/trajectories(186-194)';
else
    matdirectory ='C:\Users\hadar\Dropbox\resuspension\Matlab\trajectories_alex_2012(186-194)\events_90';
end

large_ones = dir(fullfile(matdirectory,'large*'));

data_large =[];
data_small = [];

D = 550 * 10^(-6); %diameter of the solid particle [m]
r = D/2;
A= pi * r^2 ;% the projected grain area perpendicular to flow direction

density_l = 0.998 *1000; % density of water ,temperature= 22 
density_s = 1.065 *1000; % density of solid patricles 

CD = 4; % drag coefficient

%ws = (density_s - density_l)* 1000 * pi  * 9.81 *4 / 3 * r^3 ;%submerged particle weight 
%fv = (1+0.5*(density_l/(density_s-density_l)));  
% f_drag_critical =  fv * ws  ; 

V = 4* pi* r^3 / 3; 
mass_particle = density_s * V; % [Kg]
g = 9.81;

R = 25;

counter1= 0;
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
        
        time_of_event = events.events(iEvents,2);
        trajind = trajIds == events.events(iEvents,1);
        data = large(trajind);
        ind = find(data.t == time_of_event);
            
            xp = data.xf(ind);
            yp = data.yf(ind);
            zp = data.zf(ind);
            
            axp = data.axf(ind)/1000^2;
            ayp = data.ayf(ind)/1000^2;
            azp = data.azf(ind)/1000^2;
            
            up = data.uf(ind)/1000;
            vp = data.vf(ind)/1000;
            wp = data.wf(ind)/1000;
            
            counter1 = counter1 + 1;
            
            %force_y = mass_particle * ayp;
            force_x_w = mass_particle * sqrt(axp.^2+azp.^2);
            data_drag_force (counter1) = force_x_w;
            
            counter = 0;
            drag_force = 0;
            
            for j = 1:numSmallTraj % for all small ones
                
                distance = dist([small(j).xf,small(j).yf,small(j).zf]',[xp,yp,zp])';
                
                relevant = find(small(j).t == data.t(ind) &  distance <= R);
                
                if ~isempty(relevant)
                    
                    tmp = small(j);
                    
                    tmp.uf = tmp.uf(relevant)/1000;
                    tmp.vf = tmp.vf(relevant)/1000;
                    tmp.wf = tmp.wf(relevant)/1000;
                    
                    tmp.yf= tmp.yf (relevant);
                    
                    if tmp.yf <= -21 ,continue, end
                    
            %velocity_p_f = [tmp.uf-up,tmp.vf-vp,tmp.wf-wp];
            velocity_p_f = [up-tmp.uf,0,wp-tmp.wf];
            
           % drag_force = drag_force + 0.5* density_l*1000* CD * A * norm(velocity_p_f) * velocity_p_f; % sum of all the relevant small particle       
            drag_force = drag_force + 0.5* density_l* A * norm(velocity_p_f) * velocity_p_f; % sum of all the relevant small particles      
                  
            counter = counter + 1;
                   
                  
                end % relevant
            end % smalltraj
            
            data_drag_liquid(counter1) = norm (drag_force) / counter; % average the small particles kinetic energy
         
    end
end


plot(data_drag_liquid,data_drag_force,'r.','MarkerSize',10);
