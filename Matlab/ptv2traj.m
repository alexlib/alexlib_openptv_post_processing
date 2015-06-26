function [traj,trajLen] = ptv2traj(ptv,minLength)
% ptv2TRAJ converts structure from ptv where every element is the particles
% in a specific time (i.e. per frame, with trajectory id's)
% to the structure TRAJ in which every element of the structure is a
% trajectory in time and space.

% Author: Alex Liberzon, alex.dot.liberzon.at.gmail.dot.com
% Last modified: 25-Feb-2010

if nargin < 2
    minLength = 2;
end

trajid = cat(2,ptv.trajid);
trajLen = trajid*0;

xf = cat(1,ptv.xr);
yf = cat(1,ptv.yr);
zf = cat(1,ptv.zr);

for i = 1:length(ptv)
    if length(ptv(i).xr) > 1
        ptv(i).t = repmat(ptv(i).t,[length(ptv(i).xr),1]);
    end
end

t = cat(1,ptv.t);
dt = min(diff(t));


unique_trajid = unique(trajid);
unique_trajid(unique_trajid < 0)= [];

traj = repmat(struct('xf',[],'yf',[],'zf',[],'uf',[],'vf',[],'wf',[],...
    'axf',[],'ayf',[],'azf',[],'t',[],'trajid',[]),length(unique_trajid),1);

i = 0;
for k = 1:length(unique_trajid)
    id = find(trajid == unique_trajid(k));
    if length(id) > minLength
        i = i + 1;
        trajLen(i) = length(id);
        
        [b,xi] = spaps(t(id),xf(id),0);
        [b,yi] = spaps(t(id),yf(id),0);
        [b,zi] = spaps(t(id),zf(id),0);
        
        traj(i).xf = xi;
        traj(i).yf = yi;
        traj(i).zf = zi;
        
        
        traj(i).uf = gradient5(traj(i).xf,t(id)*dt);
        traj(i).vf = gradient5(traj(i).yf,t(id)*dt);
        traj(i).wf = gradient5(traj(i).zf,t(id)*dt);
        
        traj(i).axf = gradient5(traj(i).uf,t(id)*dt);
        traj(i).ayf = gradient5(traj(i).vf,t(id)*dt);
        traj(i).azf = gradient5(traj(i).wf,t(id)*dt);
        
        
        traj(i).t = t(id);
        traj(i).trajid = trajid(id);
        
        xf(id) = [];
        yf(id) = [];
        zf(id) = [];
        trajid(id) = [];
    end
end
traj(i+1:end) = [];
trajLen(i+1:end) = [];



