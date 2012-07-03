% run_collect_statistics_type2_event_small
clc, clear

if ismac
    matdirectory = '/Users/alex/Dropbox/resuspension/2011/trajectories(186-194)';
else
    matdirectory ='C:\Users\hadar\Dropbox\resuspension\2011\trajectories(186-194)';
end
large_ones = dir(fullfile(matdirectory,'large*'));

quantity1 = 'uf';
quantity2 = 'vf';
quantity3= 'wf';
R = 25; %mm

data = [];

for i = 1:length(large_ones)
    disp(large_ones(i).name)
    collected_data = collect_statistics_type2_event_small(large_ones(i).name,quantity1,quantity2,quantity3,R); 
    data = cat(1,data,collected_data);
end

figure, 
scatter(data(:,1),data(:,2),'b.');
xlabel(sprintf('small quantity %s',quantity1));
ylabel(sprintf('small quantity %s',quantity2));
set(gca,'xtick',0); set(gca,'ytick',0);
box on; grid on; axis equal


figure, 
scatter(data(:,1),data(:,2),'b.');
xlabel(sprintf('small quantity %s',quantity3));
ylabel(sprintf('small quantity %s',quantity2));
set(gca,'xtick',0); set(gca,'ytick',0);
box on; grid on; axis equal


hf1=figure;
subplot(311), hold on
nhist(data(:,1),41),hold on
title(sprintf('small quantity: %s',quantity1))

hold on
subplot(312),
nhist(data(:,2),41)
title(sprintf('small quantity: %s',quantity2))

hold on
subplot(313),
nhist(data(:,3),41)
title(sprintf('small quantity: %s',quantity3))
