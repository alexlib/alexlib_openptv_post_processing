function [varargout] = ptv_is_to_traj(directory,first,last,minlength,dt)
% TRAJ = PTV_IS_TO_TRAJ(DIRECTORY,FIRST,LAST,MINLENGTH,DT)
% reads the ptv_is.* files from the directory and 
% converts it to trajectories into a structure TRAJ
% Optional arguments:
% FIRST,LAST - numbers of the frames from which to which to read
% MINLENGTH - minimum length of trajectory
% DT - 1/frame rate of the recording, time interval between two frames
%
% Example:
%{
    directory = '../test_data'
    data = ptv_is_to_traj(directory,10000,10004,3,1/100)
    plot_long_trajectories(data,5)
%}


if nargin < 5
    dt = 1; % default dt is not known
end
if nargin < 4
    minlength = 5; % shortest trajectory is 5 frames
end

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

    
[traj,trajLen] = ptv2traj_v2(newdata,minlength,dt);
save tmp traj
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