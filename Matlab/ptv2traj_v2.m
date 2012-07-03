function [traj,trajLen] = ptv2traj_v2(ptv,minLength,dt)
% ptv2TRAJ_v2 converts structure from ptv where every element is the particles
% in a specific time (i.e. per frame, with trajectory id's)
% to the structure TRAJ in which every element of the structure is a
% trajectory in time and space.

% Author: Alex Liberzon, alex.dot.liberzon.at.gmail.dot.com
% Created: 2-Apr-2011
% Last modified 2-Apr-2011
% - should work using sortrows and diff of the last column approach
%

if nargin < 2
    minLength = 3;
    dt = NaN;
elseif nargin < 3
    dt = NaN;
end



for i = 1:length(ptv)
    if length(ptv(i).xr) > 1
        ptv(i).t = repmat(ptv(i).t,[length(ptv(i).xr),1]);
    end
end

if isnan(dt)
    dt = min(diff(t));
end
try
    tmp = cat(2,cat(1,ptv.xr),cat(1,ptv.yr),cat(1,ptv.zr),cat(2,ptv.trajid).',cat(1,ptv.t));
catch
        tmp = cat(2,cat(1,ptv.xr),cat(1,ptv.yr),cat(1,ptv.zr),cat(1,ptv.trajid),cat(1,ptv.t));
end
tmp = sortrows(tmp,4);
tmp = tmp(find(tmp(:,4)==-999,1,'last')+1:end,:); % cut out -999
[len,val,first,last] = mttrlencode(tmp(:,4));
idl = find(len > minLength); % only long ones

traj = repmat(struct('xf',[],'yf',[],'zf',[],'uf',[],'vf',[],'wf',[],...
    'axf',[],'ayf',[],'azf',[],'t',[],'trajid',[]),length(idl),1);

%h = waitbar(0,'Please wait ...');
for k = 1:length(idl)
    
    i = idl(k);
    %waitbar(k/length(idl),h);
    
    traj(k).xf = tmp(first(i):last(i),1);
    traj(k).yf = tmp(first(i):last(i),2);
    traj(k).zf = tmp(first(i):last(i),3);
    traj(k).trajid  = tmp(first(i):last(i),4);
    traj(k).t  = tmp(first(i):last(i),5);
    
    [b,traj(k).xf] = spaps(traj(k).t,traj(k).xf,0);
    [b,traj(k).yf] = spaps(traj(k).t,traj(k).yf,0);
    [b,traj(k).zf] = spaps(traj(k).t,traj(k).zf,0);
    
    traj(k).uf = gradient5(traj(k).xf,traj(k).t*dt);
    traj(k).vf = gradient5(traj(k).yf,traj(k).t*dt);
    traj(k).wf = gradient5(traj(k).zf,traj(k).t*dt);
    
    traj(k).axf = gradient5(traj(k).uf,traj(k).t*dt);
    traj(k).ayf = gradient5(traj(k).vf,traj(k).t*dt);
    traj(k).azf = gradient5(traj(k).wf,traj(k).t*dt);
    
    trajLen(k) = len(i);
end



