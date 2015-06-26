function xuap = clean_xuap(xuap)
% added time to xuap



% remove zeros

% tmp = xuap;
flds = fieldnames(xuap);
flds(strcmp(flds,'t')) = [];




id = [];
frameid = [];
for i = 1:length(xuap)
    for j = 1:length(xuap(i).xf)
        if (abs(xuap(i).xf(j)) < eps || abs(xuap(i).yf(j)) < eps || abs(xuap(i).zf(j)) < eps )
            id = cat(1,id,j);
        end
    end
    if ~isempty(id)
        for k = 1:length(flds)
%             if length(xuap(i).(flds{k})) > 1
                xuap(i).(flds{k})(id) = [];
%             end
        end
        id = [];
    end
    if isempty(xuap(i).xf)
        frameid = cat(1,frameid,i);
    end
end
xuap(frameid) = [];
