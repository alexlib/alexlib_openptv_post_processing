clc, clear

if ismac
    matdirectory = '/Users/alex/Dropbox/resuspension/2011/trajectories(186-194)';
else
    matdirectory ='C:\Users\hadar\Dropbox\resuspension\Matlab\trajectories_2012(186-194)\events_80';
end
large_ones = dir(fullfile(matdirectory,'large*'));

quantity1 = 'xf';
quantity2 = 'zf';
data = [];

for i = 1:length(large_ones)
    disp(large_ones(i).name)
    collected_data = type2_location_events(large_ones(i).name,quantity1,quantity2); 
    data = cat(1,data,collected_data);
end

%figure, 
%scatter(data(:,1),data(:,2),'b.');
%xlabel(sprintf('x'));
%ylabel(sprintf('z'));
%set(gca,'xtick',0); set(gca,'ytick',0);
%box on; grid on; 

figure
plot(data(:,1)./295,data(:,2)./295,'r.')
