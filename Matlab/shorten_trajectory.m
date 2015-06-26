function traj1 = trim_trajectory(traj1,breakpoint)
fields = fieldnames(traj1);
    for j = 1:length(fields)
        traj1.(fields{j}) = traj1.(fields{j})(1:breakpoint);
    end
end