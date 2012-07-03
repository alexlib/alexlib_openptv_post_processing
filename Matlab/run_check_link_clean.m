% run_check_link_clean is the script that operates the various functions
% that use the xuap.* files and create the trajectory structure with the
% clean and long trajectories
%

% Author: Alex Liberzon
% Copyright (c) 2009 OpenPIV http://www.openpiv.net


% directory = 'c:\PTV\Experiments\Mark3DPTVSingleParticle\res';
% directory = '.';


% initial values
smoothness = 0.1; % for the Radial Basis Functions Splines
jump_thresh = 0.1; % empirical value
distThresh = 0.001; % 3 mm
maxDt = 5;
xmax = 0.05; % 42 x 2 = 84 mm is the size of our field of view.
dt = 1/125; % 50 fps

directory = '../exp_031110/res/';
% dt =  0.005; % sec
% xmax = Inf;

%{
if ~exist('xuap','var')
    if exist(fullfile(directory,'xuap.mat'),'file')
        load(fullfile(directory,'xuap.mat'));
    else
        xuap = readxuapfiles(directory);
    end
end

xuaptraj = building_trajectories(xuap(1:end));
xuaptraj = xuap_cleaning_from_zeros(xuaptraj,xmax);
traj = xuap2traj(xuaptraj);
%}

% Load data
% if ~exist('traj','var')
%     load(fullfile(directory,'traj.mat'));
% end

% clean the "jumps", compare the actual velocity and the one from the
% post_process.exe - something is fishy there.
tmp = zeros(1,length(traj));
for i = 1:length(traj)
    tmp(i) = norm([gradient(traj(i).xf,traj(i).t*dt) - traj(i).uf,...
        gradient(traj(i).yf,traj(i).t*dt) - traj(i).vf,...
        gradient(traj(i).zf,traj(i).t*dt) - traj(i).wf]);
    if any(diff(traj(i).t)>50), tmp(i) = 2*jump_thresh; end
end
%
figure, plot(1:length(traj),tmp)
%
% ind = (tmp > jump_thresh);
% tmptraj = traj(~ind); % all those above the threshold are "jumps"


% check what we can link



traj = haitao_linking_criteria(traj,maxDt,distThresh);

% k = 0;

newtraj = traj; % backup
trajid_toremove = [];
links = find([linkid{:,3}] < distThresh);
trajid_to_interpolate = zeros(length(newtraj),1);
for i = 1:length(links)
    [k,m,dist] = linkid{links(i),:};
    newtraj(m) = combine_trajectories(newtraj(k),newtraj(m));
    trajid_toremove = cat(1,trajid_toremove,k);
    trajid_to_interpolate(m) = 1;
    trajid_to_interpolate(k) = 0;
end

trajid_to_interpolate = find(trajid_to_interpolate);


for i = 1:length(trajid_to_interpolate)
    tmp = link_trajectories_cubicspline(newtraj(trajid_to_interpolate(i)),dt);
    fields = fieldnames(tmp);
    for j = 1:length(fields)
        newtraj(trajid_to_interpolate(i)).(fields{j}) = tmp.(fields{j});
    end
end
newtraj(trajid_toremove) = [];
newtraj = ensure_order_traj(newtraj);


for i = 1:length(newtraj)
    tmp  = link_trajectories_rbf(newtraj(i),.5,dt);
       fields = fieldnames(tmp);
    for j = 1:length(fields)
        newtraj(i).(fields{j}) = tmp.(fields{j});
    end

end
newtraj = ensure_order_traj(newtraj);


%% iterations
%{
for kIter = 1:5
    linkid = haitao_linking_criteria(newtraj,maxDt);
    trajid_toremove = [];
    links = find([linkid{:,3}] < distThresh*kIter);
    trajid_to_interpolate = zeros(length(newtraj),1);
    for i = 1:length(links)
        [k,m,dist] = linkid{links(i),:};
        newtraj(m) = combine_trajectories(newtraj(k),newtraj(m));
        trajid_toremove = cat(1,trajid_toremove,k);
        trajid_to_interpolate(m) = 1;
        trajid_to_interpolate(k) = 0;
    end
    trajid_to_interpolate = find(trajid_to_interpolate)
    for i = 1:length(trajid_to_interpolate);
        newtraj(trajid_to_interpolate(i)) = link_trajectories_rbf(newtraj(trajid_to_interpolate(i)),smoothness,dt);
    end
    newtraj(trajid_toremove) = [];
    newtraj = ensure_order_traj(newtraj);
    plot_trajectory_properties_vs_t(newtraj(1:10),0)
end
%}

%% Manual 
%{
newtraj = visual_manual_linking(newtraj,24,27,smoothness,dt)
%}


