% run_collect_statistics_type2_event_larg_and_small

if ismac
    matdirectory = '/Users/alex/Dropbox/resuspension/2011/trajectories(186-194)';
else
    matdirectory ='C:\Users\hadar\Dropbox\resuspension\2011\trajectories(186-194)';
end
large_ones = dir(fullfile(matdirectory,'large*'));

quantity1 = 'vf';
quantity2 = 'vf';
R = 40; %mm

data = [];

for i = 1:length(large_ones)
    disp(large_ones(i).name)
    % collected_data = collect_statistics_type2_event_small(large_ones(i).name,quantity1,quantity2,R); 
    collected_data = collect_statistics_type2_event_large_and_small(large_ones(i).name,quantity1,quantity2,R); 
    data = cat(1,data,collected_data);
end

figure, 
scatter(data(:,1),data(:,2),'r.');
xlabel(sprintf('large quantity %s',quantity1));
ylabel(sprintf('small quantity %s',quantity2));
set(gca,'xtick',0); set(gca,'ytick',0);
box on; grid on; axis equal

% Instead of scatter, we can make histogram
 hf1 = figure;
    subplot(211), hold on
    hist(data(:,1),10)
    title(sprintf('large quantity: %s',quantity1))
    subplot(212),
    hist(data(:,2),10)
title(sprintf('small quantity: %s',quantity2))
 
