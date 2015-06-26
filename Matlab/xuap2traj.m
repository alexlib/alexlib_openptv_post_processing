function [traj,trajLen] = xuap2traj(xuap,minLength)
% XUAP2TRAJ converts structure from XUAP where every element is the particles
% in a specific time (i.e. per frame, with trajectory id's)
% to the structure TRAJ in which every element of the structure is a
% trajectory in time and space. 

% Author: Alex Liberzon, alex.dot.liberzon.at.gmail.dot.com
% Last modified: 25-Feb-2010

if nargin < 2
    minLength = 2;
end

trajid = cat(2,xuap.trajid);
trajLen = trajid*0;

xf = cat(1,xuap.xf);
yf = cat(1,xuap.yf);
zf = cat(1,xuap.zf);

% xf = cat(1,xuap.xr);
% yf = cat(1,xuap.yr);
% zf = cat(1,xuap.zr);

uf = cat(1,xuap.uf);
vf = cat(1,xuap.vf);
wf = cat(1,xuap.wf);
axf = cat(1,xuap.axf);
ayf = cat(1,xuap.ayf);
azf = cat(1,xuap.azf);

for i = 1:length(xuap)
    if length(xuap(i).xf) > 1
        xuap(i).t = repmat(xuap(i).t,[length(xuap(i).xf),1]);
    end
end

t = cat(1,xuap.t);

unique_trajid = unique(trajid);
unique_trajid(unique_trajid < 1)= [];

traj = repmat(struct('xf',[],'yf',[],'zf',[],'uf',[],'vf',[],'wf',[],...
    'axf',[],'ayf',[],'azf',[],'t',[],'trajid',[]),length(unique_trajid),1);

i = 0;
for k = 1:length(unique_trajid)
    id = find(trajid == unique_trajid(k));
    if length(id) > minLength
        i = i + 1;
        trajLen(i) = length(id);
        traj(i).xf = xf(id);
        traj(i).yf = yf(id);
        traj(i).zf = zf(id);
        traj(i).uf = uf(id);
        traj(i).vf = vf(id);
        traj(i).wf = wf(id);
        traj(i).axf = axf(id);
        traj(i).ayf = ayf(id);
        traj(i).azf = azf(id);
        traj(i).t = t(id);
        traj(i).trajid = trajid(id);
    end
end
traj(i+1:end) = [];
trajLen(i+1:end) = [];



