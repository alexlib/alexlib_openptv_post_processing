function [traj,numTraj,lenTraj] = building_trajectories_ptc(directory,varargin)
% BUILDING_TRAJECTORIES_PTC "tracks" the particles using the PTCXXX.dat
% files. Since every particle has its unique ID the tracking is simple
% book keeping. 
% [TRAJ, NUMTRAJ, LENTRAJ] = BUILDING_TRAJECTORIES(DIRECTORY,MINLENGTH)
% creates Matlab structure data set of trajectories:
% traj.x, traj.y, traj.z of the 1:numTraj trajectories, each of the
% lenTraj(id) length. 
% 
% Example:
% >> [traj,numTraj,lenTraj]= building_trajectories_ptc('./case352/ptc',1);
% >> figure, [maxLen,i] = max(lenTraj); scatter3(traj(i).x,traj(i).y,traj(i).z)
% >> figure, hist(lenTraj);
%
% See also: HELP BUILDING_TRAJECTORIES

d = dir(fullfile(directory,'*.dat')); 
d(end) = [];

tmp = load(fullfile(directory,d(end).name));
numTraj = max(tmp(:,1)); % last file, largest id


traj = struct('x',[],'y',[],'z',[]);
traj = repmat(traj,numTraj,1);

for k = 1:length(d)
    tmp = load(fullfile(directory,d(k).name));
    for j = 1:length(tmp(:,1))
        traj(tmp(j,1)).x = cat(1,traj(tmp(j,1)).x,tmp(j,2));
        traj(tmp(j,1)).y = cat(1,traj(tmp(j,1)).y,tmp(j,3));
        traj(tmp(j,1)).z = cat(1,traj(tmp(j,1)).z,tmp(j,4));
    end
end
% clean empty ones
numTraj = length(traj);

if nargin > 1
    minLength = varargin{1};
else
    minLength = 5;
end

emptyId = [];
lenTraj(numTraj) = 0;
for i = 1:numTraj
    len = length(traj(i).x);
    if len < minLength
        emptyId = cat(1,emptyId,i);
    else
        lenTraj(i) = len;
    end
end

traj(emptyId) = [];
numTraj = length(traj);
lenTraj = lenTraj(lenTraj > 0);



    

            