function [traj] = ptv_is_to_traj(directory,first,last,minlength,dt)
% converts a directory full of ptv_is.* files into a structure 
% 'traj' that contains Lagrangian trajectories


disp('Reading ...')
if nargin == 1
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

    
[traj,~] = ptv2traj_v2(newdata,minlength,dt);
% save tmp traj
% disp('Done')
traj = ensure_order_traj(traj);
% save tmp traj

% If manual cleaning is needed
%traj = graphical_cleaning_traj(traj,'xy');

% varargout{1} = traj;



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