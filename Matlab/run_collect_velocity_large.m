% run_.... is the outer loop for any function ...
% that performs the same operation on all the data
% collected from 3D-PTV in the resuspension experiment
% 


if ismac
    matdirectory = '/Users/alex/Dropbox/resuspension/2011/trajectories(186-194)';
else
    matdirectory ='C:\Users\hadar\Dropbox\resuspension\2011\trajectories(186-194)';
end
large_ones = dir(fullfile(matdirectory,'large*'));


styles = {'ro','bs','cd','g*','k+','r<','b>'};

figure;
hold on;
counter = 0;
for k = [-8:2:2]
    data = [];
    for i = 1:length(large_ones)
        disp(large_ones(i).name)
        collected_data = collect_velocity_large(fullfile(matdirectory,large_ones(i).name),k); 
        data = cat(1,data,collected_data);
    end
    counter = counter + 1;
   
    scatter(data(:,1),data(:,2),styles{counter},'DisplayName',int2str(k));
%     scatter(data(:,3),data(:,2),styles{counter},'DisplayName',int2str(k));

end





