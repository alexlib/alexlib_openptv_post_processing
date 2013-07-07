% run_velocity_correlation_large_small for all events
clc,clear
if ismac
    matdirectory = '/Users/alex/Dropbox/resuspension/2011/trajectories_2012(186-194)';
else
    matdirectory ='C:\Users\hadar\Dropbox\resuspension\2011\trajectories_2012(186-194)';
end

large_ones = dir(fullfile(matdirectory,'large*'));

data_correlation = [];

for i = 1:length(large_ones)
    disp(large_ones(i).name)
    collected_data = velocity_correlation_large_small(large_ones(i).name); 
    data_correlation = cat(1,data_correlation,collected_data);
end

temp  = sortrows(data_correlation , 1); % organizes tha data from small to large 
 figure
plot (temp(:,1), temp(:,2),'.')
figure
plot (temp(:,1), temp(:,3),'.')
 

figure

hold on
bins = {[0:2:70],[0:5:70],[5:2:70],[5:5:70],[8:3:70],[10:1:70]};
for i = 1:length(bins)
  [bin_centers,average_in_bins,n] = plot_average_in_bins(temp,bins{i});
 plot(bin_centers,average_in_bins(1:end-1)/max(average_in_bins(1:end-1)),'r.')
end



figure

hold on
bins = {[0:2:70],[0:5:70],[5:2:70],[5:5:70],[5:1:70]};
for i = 1:length(bins)
  [bin_centers,average_in_bins,n] = plot_average_in_bins(temp,bins{i});
 plot(bin_centers,average_in_bins(1:end-1),'b.')
end 



figure

bar(bin_centers,n(1:end-1))

figure
plot(bin_centers,cumsum(n(1:end-1))



