function out = traj_length(traj)
out = zeros(length(traj),1);
fields = fieldnames(traj);
for i = 1:length(traj)
    out(i) = length(traj(i).(fields{1}));
end