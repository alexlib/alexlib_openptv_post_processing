clc,clear

if ismac
    matdirectory = '/Users/alex/Dropbox/resuspension/2011/trajectories(186-194)';
else
    matdirectory ='C:\Users\hadar\Dropbox\resuspension\2011/trajectories(186-194)/events_110';
end 

large_ones = dir(fullfile(matdirectory,'large*'));

    hf1 = figure; hold on
    ylabel(' Nominal  drag force [N]');
    xlabel('Time [s]')
    
   % hf2 = figure; hold on
   %ylabel('force_x_w');
   % xlabel('Time [frames]')
    
for i = 1:length(large_ones)
    force_particle(large_ones(i).name);
end
