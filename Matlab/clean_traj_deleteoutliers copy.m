flds = fieldnames(traj);
frameid = zeros(length(traj),1);
% flds(strcmp(flds,'t')) = []; % t is not touched
for i = 1:length(traj)
    [b,idx,outliers] = deleteoutliers(traj(i).xf,.6);
     [b,idy,outliers] = deleteoutliers(traj(i).yf,.6);
     [b,idz,outliers] = deleteoutliers(traj(i).zf,0.66);

    id = union(union(idx,idy),idz);

    if ~isempty(id)
        for k = 1:length(flds)
            %              if length(xuap(i).(flds{k})) > 1
            traj(i).(flds{k})(id) = [];
            %              end
        end
    end
    if isempty(traj(i).xf)
        frameid(i) = 1;
    end
end
traj(find(frameid)) = [];