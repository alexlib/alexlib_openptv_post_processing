function traj = clean_traj(traj)
% remove zeros

flds = fieldnames(traj);
% flds(strcmp(flds,'t')) = [];




id = [];
frameid = [];
for i = 1:length(traj)
    for j = 1:length(traj(i).xf)
        if (traj(i).xf(j) == 0 || traj(i).yf(j) == 0 || traj(i).zf(j) == 0 )
            id = cat(1,id,j);
        end
    end
    if ~isempty(id)
        for k = 1:length(flds)
%             if length(traj(i).(flds{k})) > 1
                traj(i).(flds{k})(id) = [];
%             end
        end
        id = [];
    end
    if isempty(traj(i).xf)
        frameid = cat(1,frameid,i);
    end
end
traj(frameid) = [];
