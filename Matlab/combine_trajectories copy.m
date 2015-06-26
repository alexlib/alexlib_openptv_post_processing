function traj1 = combine_trajectories(traj1,traj2)
% TRAJ12 = COMBINE_TRAJECTORIES(TRAJ1,TRAJ2) combines the trajectories
fields = fieldnames(traj1);
    for j = 1:length(fields)
        traj1.(fields{j}) = cat(1,traj1.(fields{j})(:),traj2.(fields{j})(:));
    end
end