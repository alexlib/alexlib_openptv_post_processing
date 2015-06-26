function traj1 = trim_trajectory(traj1,fromFirstEnd,fromLastEnd)
% trim_trajectory(traj1,fromFirstEnd,fromLastEnd)
fields = fieldnames(traj1);
for k = 1:length(traj1)
    len = length(traj1(k).(fields{1}));
    for j = 1:length(fields)
        traj1(k).(fields{j}) = traj1(k).(fields{j})(fromFirstEnd:len-fromLastEnd+1);
    end
end
