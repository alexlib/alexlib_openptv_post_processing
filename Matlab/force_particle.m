function [hf1] = force_particle(filename)
% Copyright (c) 2011, Turbulence Structure Laboratory

if ismac
   matdirectory = '/Users/alex/Dropbox/resuspension/2011/trajectories(186-194)';
else
    matdirectory ='C:\Users\hadar\Dropbox\resuspension\2011\trajectories(186-194)';
end

[path,filename,ext] = fileparts(filename);

% Load the largeXXXX.mat
large = load(fullfile(matdirectory,filename));
% and rename it to 'large'
large = large.(filename);

% Load the eventsXXXXX.mat
events = load(fullfile(matdirectory,strrep(filename,'large','events')));

numEvents = size(events.events,1);
%numTraj = length(small);

trajIds = unique(cat(1,large.trajid));

for i = 1:numEvents
    % Select the trajectory and the time of event
    time_of_event = events.events(i,2);
    trajind = trajIds == events.events(i,1);
    
    data = large(trajind);
    ind = find(data.t == time_of_event);
    
     axp = data.axf;
     ayp = data.ayf;
     azp = data.azf; 
     
    
     D = 550 * 10^(-6); %diameter of the solid particle [m]
     r = D/2; %[m]
     density_s = 1.065 ; %density of solid patricles [gr/cm^3]
     volume_s = 4* pi* (r*100)^3 / 3; % [cm^3]
     mass_particle = density_s * volume_s/1000; % [Kg]

    %force_y = mass_particle * ayp;
    force_x_w = mass_particle * sqrt(axp.^2+azp.^2+ayp.^2);

    
    
     %Plot the large particles
  
   % plot((data.t(1:end-2)-data.t(1)(/160,force_y(1:end-2),'LineWidth',1,'DisplayName',sprintf('%d',data.trajid(1)));
   % plot((data.t(ind)-data.t(1))/160,force_y(ind),'ro','markersize',10);%marks the time of event.
   

    plot((data.t(3:end-2)-data.t(1))/160,force_x_w(3:end-2),'LineWidth',1,'DisplayName',sprintf('%d',data.trajid(1)));
    plot((data.t(ind)-data.t(1))/160,force_x_w(ind),'ro','markersize',10);%marks the time of event.
end
 


