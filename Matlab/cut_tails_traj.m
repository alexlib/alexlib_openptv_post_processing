function traj = cut_tails_traj(traj,tail)
% TRAJ = CUT_TAILS_TRAJ(TRAJ,LENGTH_OF_TAIL) cuts all objects in the 
% structure TRAJ. Default length is 2 (if not provided). Objects which are
% shorter than the 2 x LENGTH_OF_TAIL then the object is removed. 

% Author: Alex Liberzon
% Copyright (c) 2010, TAU

if ~nargin
    help cut_tails_traj
    return
end
if nargin < 2
    tail = 2;
end
fields = fieldnames(traj);
remove = [];
for i = 1:length(traj)
    for j = 1:length(fields)
        traj(i).(fields{j}) = traj(i).(fields{j})(3:end-2);
    end
    if isempty(traj(i).(fields{1}))
        remove = cat(1,remove,i);
    end
end
traj(remove) = [];
