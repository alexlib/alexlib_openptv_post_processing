function traj = ensure_order_traj(traj)
fields = fieldnames(traj);

for i = 1:length(traj)
    for j = 1:length(fields)
        traj(i).(fields{j}) = traj(i).(fields{j})(:);
    end
end

t = zeros(length(traj),1);
for i = 1:length(traj)
    t(i) = traj(i).t(1);
end
[t,ind] = sort(t);
traj = traj(ind);