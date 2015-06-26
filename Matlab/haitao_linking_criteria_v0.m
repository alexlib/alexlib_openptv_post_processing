function [linkid] = haitao_linking_criteria(traj,xdist)
% HAITAO_LINKNG_CRITERIA builds containers of trajectories that can be
% linked together. Simple criteria of distance between the end of the
% previous trajectory and the beginning of the next trajectory is used.
% no check on velocity, acceleration, etc.

% dist = inline('(x-y)^2','x','y');

% first container is the first trajectory
clear linkid
dt = 1/50;
k = 1;
linkid{k} = 1;

%threshold of the distance
if nargin < 2
    xdist = 3e-3; %default
end
figure, 
hold on

%%
for i = 2:length(traj)
    dt = traj(i).t(1) - traj(i-1).t(end);
    dt = dt/50;
    predicted_x = traj(i-1).xf(end) + dt*traj(i-1).uf(end);
    predicted_y = traj(i-1).yf(end) + dt*traj(i-1).vf(end);
    predicted_z = traj(i-1).zf(end) + dt*traj(i-1).wf(end);
%     predicted_x = traj(i-1).xf(end);
%     predicted_y = traj(i-1).yf(end);
%     predicted_z = traj(i-1).zf(end);

    disp(i) 
    disp([predicted_x,predicted_y,predicted_z]);
    disp([traj(i).xf(1),traj(i).yf(1),traj(i).zf(1)]);
    
    plot(i,norm([predicted_x,predicted_y,predicted_z]-[traj(i).xf(1),traj(i).yf(1),traj(i).zf(1)]),'.')
    
    if norm([predicted_x,predicted_y,predicted_z]- [traj(i).xf(1),traj(i).yf(1),traj(i).zf(1)])  < xdist
        linkid{k} = cat(1,linkid{k},i); % continue filling this container
    else
        k = k + 1; % next container
        linkid{k} = i;
    end
end
        
        