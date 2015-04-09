% run_collect_statistics_type2_event_small

if ismac
    matdirectory = '/Users/alex/Dropbox/resuspension/2011/trajectories_2012(186-194)';
else
    matdirectory ='C:\Users\dell\Dropbox\resuspension\Matlab\trajectories_alex_2012(186-194)';
end
large_ones = dir(fullfile(matdirectory,'large*'));

 quantity1 = 'ax';
 quantity2 = 'ay';
 quantity3= 'az';
 
R = 25; %mm
%Ws = 12; % settlling velocity [mm/s]

data = [];

for i = 1:length(large_ones)
    disp(large_ones(i).name)
    collected_data = collect_statistics_event_small_nhist(large_ones(i).name,R); 
    data = cat(1,data,collected_data);
end

subplot(311),
nhist(data(:,1),21)
title(sprintf('small quantity: %s',quantity1))

hold on
subplot(312),
nhist(data(:,2),21)
title(sprintf('small quantity: %s',quantity2))

hold on
subplot(313),
nhist(data(:,3),21)
title(sprintf('small quantity: %s',quantity3))
