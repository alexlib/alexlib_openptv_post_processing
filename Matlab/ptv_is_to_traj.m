function [varargout] = ptv_is_to_traj(directory,minlength,dt,first,last)
% directory = 'C:\Documents and Settings\user\My Documents\People\ETH\Beat\colloidbreakage\Colloid\Colloid\ptv_is';
% directory = 'c:\PTV\Experiments\res_128_Re1435ptvis\'
% directory = 'c:\PTV\Experiments\res_scene126_Re_890\'
% directory = 'c:\PTV\Experiments\res_scene130_Re1920\';
% directory = '/Users/alex/Desktop/res25'; % resuspension project, Aug 2, 2010

disp('Reading ...')
if nargin == 1 % only directory name
    minlength = 5; 
    dt = 1;
end

if nargin < 5
    data = read_ptv_is_files(directory);
else
    data = read_ptv_is_files(directory,first,last);
end
disp('Done .')
save tmp data
disp('Building trajectories')
newdata = building_trajectories(data);
save tmp newdata
disp('Done')
disp('PTV to Traj')
[traj,trajLen] = ptv2traj(newdata,minlength,dt);
save tmp traj trajLen
disp('Done')
traj = ensure_order_traj(traj);
save tmp traj

% If manual cleaning is needed
%traj = graphical_cleaning_traj(traj,'xy');

varargout{1} = traj;



% plot_long_trajectories(newtraj,5)
% 
% [traj,removelist] = haitao_linking_criteria(newtraj,10,1.5);
% 
% figure, plot_long_trajectories(traj,5)




% save(fullfile(directory,'traj.mat'),'traj')
% load(fullfile(directory,'traj.mat'),'traj')
% 
% s = 100;
% for i = 1:length(traj)
%     newtraj(i) = link_trajectories_smoothn(traj(i),s);
% end
% % 
% for i = 1:length(traj)
%     newtraj(i) = link_trajectories_rbf(traj(i),0.1);
% end
% [traj,removelist] = haitao_linking_criteria(newtraj,10,.1);
% save colloid_traj traj


% plot_long_trajectories(traj,5)
% hold on
% plot_long_trajectories(newtraj(removelist),1,1)