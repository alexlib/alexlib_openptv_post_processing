function traj2 = cutout_trajectory(traj1,trajNum,cutout)
fields = fieldnames(traj1);
traj2(1:trajNum-1) = traj1(1:trajNum-1);
for j = 1:length(fields)
    traj2(trajNum).(fields{j}) = traj1(trajNum).(fields{j})(1:cutout(1)-1);
    traj2(trajNum+1).(fields{j}) = traj1(trajNum).(fields{j})(cutout(end)+1:end);
end
traj2(trajNum+2:length(traj1)+1) = traj1(trajNum+1:end);