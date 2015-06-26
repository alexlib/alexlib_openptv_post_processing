function traj2 = split_trajectory(traj1,trajNum,breakpoint)
fields = fieldnames(traj1);
traj2(1:trajNum-1) = traj1(1:trajNum-1);
for j = 1:length(fields)
    traj2(trajNum).(fields{j}) = traj1(trajNum).(fields{j})(1:breakpoint);
    traj2(trajNum+1).(fields{j}) = traj1(trajNum).(fields{j})(breakpoint+1:end);
end
traj2(trajNum+2:length(traj1)+1) = traj1(trajNum+1:end);