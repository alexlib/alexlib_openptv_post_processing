% run_velocity_correlation_large_small for all events

if ismac
    matdirectory = '/Users/alex/Dropbox/resuspension/2011/trajectories(186-194)';
else
    matdirectory ='C:\Users\hadar\Dropbox\resuspension\2011\trajectories(186-194)';
end

large_ones = dir(fullfile(matdirectory,'large*'));

data_correlation = [];

for i = 1:length(large_ones)
    disp(large_ones(i).name)
    collected_data = velocity_correlation_large_small(large_ones(i).name); 
    data_correlation = cat(1,data_correlation,collected_data);
end

temp  = sortrows(data_correlation , 1); % organizes tha data from small to large 
 %plot (temp(:,1), temp(:,2)./temp(1,2),'o')
 
  [bin_centers,average_in_bins,n] = plot_average_in_bins(temp,0:5:90);
  
  
  figure
  plot(bin_centers,average_in_bins(1:end-1)/max(average_in_bins(1:end-1)),'o')
   
  
  figure
  bar(bin_centers,n(1:end-1))