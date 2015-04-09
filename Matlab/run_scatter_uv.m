clc, clear

if ismac
    matdirectory = '/Users/alex/Dropbox/resuspension/2011/trajectories_2012(186-194)';
else
    matdirectory ='C:\Users\hadar\Dropbox\resuspension\Matlab\trajectories_alex_2012(186-194)\events_80';
end
large_ones = dir(fullfile(matdirectory,'large*'));

 quantity1 = 'vf';
 quantity2 = 'uf';
 quantity3= 'wf';
 
R = 25; %mm
Ws = 12; % settlling velocity [mm/s]

data = [];

for i = 1:length(large_ones)
    disp(large_ones(i).name)
    collected_data = scatter_uv(large_ones(i).name,R); 
    data = cat(1,data,collected_data);
end

figure, 
scatter(data(:,1)./(Ws)^2,data(:,2)./(Ws)^2,'b.');
xlabel(sprintf('u_fw_f/W_s^2'));
ylabel(sprintf('u_fv_f/W_s^2'));
set(gca,'xtick',0); set(gca,'ytick',0);
box on; grid on; 

figure, 
scatter(data(:,1)./(Ws)^2,data(:,3)./(Ws)^2,'b.');
xlabel(sprintf('u_fw_f/W_s^2'));
ylabel(sprintf('v_fw_f/W_s^2'));
set(gca,'xtick',0); set(gca,'ytick',0);
box on; grid on; 

figure, 
scatter(data(:,2)./(Ws)^2,data(:,3)./(Ws)^2,'b.');
xlabel(sprintf('u_fv_f/W_s^2'));
ylabel(sprintf('w_fv_f/W_s^2'));
set(gca,'xtick',0); set(gca,'ytick',0);
box on; grid on; 



%figure;
%hist(data(:,2)./Ws,41)
% title(sprintf('small quantity: %s',quantity1))
% 
% hold on
% subplot(312),
% nhist(data(:,2),21)
% title(sprintf('small quantity: %s',quantity2))
% 
% hold on
% subplot(313),
% nhist(data(:,3),21)
% title(sprintf('small quantity: %s',quantity3))
