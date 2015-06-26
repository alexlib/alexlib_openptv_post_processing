% Due to large gaps:
% two stage tracking:
% 1. use clustering to get geometrically closed points without time issues
% 2. store in trajectories, then relink with gaps.
%

% data is in the form of the particle clouds
% X Y Z T
% Option 1: read from the file
% data = dlmread('tmp.txt');

% Option 2: use the xuap or traj structure:
% data = cat(2,cat(1,data.xf),cat(1,data.yf),cat(1,data.zf),cat(1,data.t));

% Option 3: use ptv_is files:
% for example:
% directory = 'c:\PTV\Experiments\res_128_Re1435ptvis\';
% ptvis = read_ptv_is_files(directory);
% for i = 1:length(ptvis)
%     ptvis(i).t = repmat(ptvis(i).t,length(ptvis(i).xr),1);
% end
% data = cat(2,cat(1,ptvis.xr),cat(1,ptvis.yr),cat(1,ptvis.zr),cat(1,ptvis.t));

% Option 4: use rt_is files:
directory = '/Volumes/MY PASSPORT/PTV/19/res_19_try'; % '/Users/alex/Documents/PTV/ptv_alex/res';
filelist = dir(fullfile(directory,'rt_is.19*'));
rtis = read_rt_is_files(directory,filelist);

for i = 1:length(rtis)
    rtis(i).t = repmat(rtis(i).t,length(rtis(i).xr),1);
end
data = cat(2,cat(1,rtis.xr),cat(1,rtis.yr),cat(1,rtis.zr),cat(1,rtis.t));

% chunk_size = 100;
% maxNumParticles = 100;

nump = zeros(length(rtis),1);
for i = 1:length(rtis)
    nump(i) = length(rtis(i).xr);
end
maxNumParticles = 2*max(nump);
chunk_size = length(data);




k = 0;
for i = 1:floor(length(data)/chunk_size)
    X = data(chunk_size*(i-1)+1:chunk_size*i,1:3); % only geometry
    Y = data(chunk_size*(i-1)+1:chunk_size*i,4);
    T = clusterdata(X,'maxclust',maxNumParticles,'linkage','single'); %,'single');
    % scatter3(X(:,1),X(:,3),X(:,2),20,T,'filled'); keyboard;
    for j = 1:max(T)
        k = k + 1;
        ind = find(T == j);
        traj(k).xf = X(ind,1);
        traj(k).yf = X(ind,2);
        traj(k).zf = X(ind,3);
        traj(k).t =  Y(ind,1);
        traj(k).trajid = repmat(k,length(ind),1);
    end
end
figure, plot_long_trajectories(traj);
pause

% Smooth and get the velocities
% initial values
smoothness = 0.1; % for the Radial Basis Functions Splines
distThresh = 0.01; % 3 mm
maxDt = 3;
dt = 1/100; % 100 fps


for i = 1:length(traj)
    newtraj(i) = link_trajectories_rbf(traj(i),smoothness);
end

% Stage 2

newtraj = haitao_linking_criteria(newtraj,maxDt,distThresh);

plot_long_trajectories(newtraj);

save traj_scene32 newtraj



