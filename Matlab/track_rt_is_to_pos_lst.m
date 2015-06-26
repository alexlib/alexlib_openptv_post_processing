% function traj = rt_is_to_pos_lst(directory,basename)
% DATA = RT_IS_TO_POS_LST(DIRECTORY,FILES)
% reads rt_is.* files in the directory and converts into a single
% ASCII file that is pos_lst for the 2D Matlab based tracking
% http://www.physics.georgetown.edu/matlab/tutorial.html
% BASENAME could be 'rt_is.35*', default = 'rt_is.*'
%
%

% See also: read_rt_is_files

%if ~nargin,
% directory = '/Users/alex/Documents/PTV/ptv_alex/res';


% directory = '/Volumes/MY PASSPORT/PTV/19/res_19_try';
% basename = 'rt_is.19*';

% directory = '/Volumes/MY PASSPORT/PTV/res_35_try';
% basename = 'rt_is.35*';



%elseif nargin == 1
%    basename = 'rt_is.35*';
% end
filelist = dir(fullfile(directory,basename));


rtis = read_rt_is_files(directory,filelist);

for i = 1:length(rtis)
    rtis(i).t = repmat(rtis(i).t,length(rtis(i).xr),1);
end
pos_lst = cat(2,cat(1,rtis.xr),cat(1,rtis.yr),cat(1,rtis.zr),cat(1,rtis.t));


% If storage is required
%
% fid = fopen('pos_list.txt','w');
%
% fprintf(fid,'%6.4f %6.4f %6.4f %d\n',  pos_lst');
%
% fclose(fid);


% Track using Matlab tracking
%

param.mem = 2;
param.dim = 3;
param.good = 5;
param.quiet = 0;


result = track(pos_lst, 5, param);
nTraj = result(end,5);


traj = repmat(struct('xf',[],'yf',[],'zf',[],'uf',[],'vf',[],'wf',[],...
    'axf',[],'ayf',[],'azf',[],'t',[],'trajid',[]),nTraj,1);

minLength = param.good;

i = 0;
for k = 1:nTraj
    id = find(result(:,5) == k);
    if length(id) > minLength
        i = i + 1;
        trajLen(i) = length(id);
        traj(i).xf = result(id,1);
        traj(i).yf = result(id,2);
        traj(i).zf = result(id,3);
        traj(i).t  = result(id,4);
        traj(i).trajid = result(id,5);
    end
end
plot_long_trajectories(traj,param.good);



k = 0;
for i = 1:length(traj)
    if length(traj(i).xf > 1)
        k = k + 1;
        newtraj(k) = link_trajectories_smoothn(traj(i));
    end
end

traj = newtraj;
k = 0;
for i = 1:length(traj)
    if length(traj(i).xf > 1)
        k = k + 1;
        newtraj(k) = link_trajectories_smoothn(traj(i));
    end
end
plot_long_trajectories(newtraj,param.good);

plot_long_trajectories(newtraj,param.good);

% manually identified pieces of a single trajectory:
% for the scene 38
listoftraj = [13,16,46];
% tmptraj = manual_linking(traj,listoftraj);
% plot_long_trajectories(tmptraj)
% plot_trajectory_properties_vs_t(tmptraj)





