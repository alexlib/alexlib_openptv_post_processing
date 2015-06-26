function varargout = show_length_of_trajectories(data)

trajID = cat(2,data.traj);
trajLen = trajID*0;

for i = 1:max(trajID)
    % traj = find(trajID == i);
    % trajLen(i) = length(traj);
    trajLen(i) = sum(trajID == i);
 
end
varargout{1} = trajLen;