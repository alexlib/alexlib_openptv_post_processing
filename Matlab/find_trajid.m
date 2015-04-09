function id = find_trajid(traj,list_of_id)
% FIND_TRAJID(TRAJECTORY_DATA,LIST_OF_ID_TO_FIND) finds the trajectories by their ID
% 
% 
% Author: Alex Liberzon
% Created:
% 14-12-2012
%
% Last modified: 
%14-12-2012

all_id_ind = zeros(length(traj),1);

for i = 1:length(traj)
    all_id_ind(i) = traj(i).trajid(1);
end



id = find(ismember(all_id_ind,list_of_id));

