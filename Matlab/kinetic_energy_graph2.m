clear all

dist = @(x,y) sqrt((x(1,:)-y(1)).^2 + (x(2,:)-y(2)).^2 + (x(3,:)-y(3)).^2);

if ismac
    matdirectory = '/Users/alex/Dropbox/resuspension/2011/trajectories(186-194)';
else
    matdirectory ='C:\Users\hadar\Dropbox\resuspension\2011\trajectories(186-194)\events_90';
end

large_ones = dir(fullfile(matdirectory,'large*'));

data_large =[];
data_small = [];


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
            up = data.uf(ind);
            vp = data.vf(ind);
            wp = data.wf(ind);
            kinetic_energy_p = 1/2* (up^2 + vp^2 + wp^2 );
            counter1 = counter1 + 1;
            
            % storage
            data_large(counter1) = kinetic_energy_p;
            
            counter = 0;
            tmp_kinetic_energy = 0;
            
            for j = 1:numSmallTraj % for all small ones
                
                distance = dist([small(j).xf,small(j).yf,small(j).zf]',[xp,yp,zp])';
                
                relevant = find(small(j).t == data.t(ind) &  distance <= R);
                
                if ~isempty(relevant)
                    
                    tmp = small(j);
                    tmp.uf = tmp.uf(relevant);
                    tmp.vf = tmp.vf(relevant);
                    tmp.wf = tmp.wf(relevant);
                    
                    tmp_kinetic_energy = tmp_kinetic_energy +1/2 *(tmp.uf^2+ tmp.vf^2 + tmp.wf^2); % sum of all the relevant small particles
                    %tmp_kinetic_energy = tmp_kinetic_energy +1/2 *(tmp.uf^2+tmp.wf^2);
                    %tmp_kinetic_energy = tmp_kinetic_energy +1/2 *( tmp.vf^2 );
                    
                    counter = counter + 1;
                    %correlation = tmp_kinetic_energy * kinetic_energy_p
                    % if verbose,plot_long_trajectories(tmp(relevant),1,hf);   end
                    
                end % relevant
            end % smalltraj
            
            data_small(counter1) =  tmp_kinetic_energy / counter; % average the small particles kinetic energy
          %correlation = data_small * kinetic_energy_p;
       
    end
end


%figure, 
%scatter(data_large,data_small,'b');
%xlabel(sprintf(' kinetic energy (particle)'));
%ylabel(sprintf('kinetic enery (flow)'));
%set(gca,'xtick',0); set(gca,'ytick',0);
%box on; grid on; 

plot(data_large,data_small,'r.','MarkerSize',10);
