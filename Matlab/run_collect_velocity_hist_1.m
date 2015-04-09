
if ismac
    matdirectory = '/Users/alex/Dropbox/resuspension/2011/trajectories(186-194)';
else
    matdirectory ='C:\Users\hadar\Desktop\trajectories_alex_2012(186-194)';
end
large_ones = dir(fullfile(matdirectory,'large*'));

data = [];

for i = 1:length(large_ones)
    disp(large_ones(i).name)
    collected_data = velocity_hist_1(large_ones(i).name); 
    data = cat(1,data,collected_data);
end


hf1=figure;
subplot(311), hold on
nhist(data(:,1),21),hold on

hold on
subplot(312),
nhist(data(:,2),21)

hold on
subplot(313),
nhist(data(:,3),21)

% hold on
% subplot(414),
% nhist(data(:,3),21)
