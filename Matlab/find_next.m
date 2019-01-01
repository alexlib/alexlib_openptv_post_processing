function [traj, num_in_traj] = find_next(traj, ptv, num_in_traj,columns,first,last)

% global traj;
% global ptv;

i=traj(num_in_traj,1);
j=traj(num_in_traj,2);
if i<last-first+1
    if ptv(i,j,2)>-1 & ptv(i+1,ptv(i,j,2)+1,columns+1)==0
        num_in_traj=num_in_traj+1;
        traj(num_in_traj,1)=i+1;
        traj(num_in_traj,2)=ptv(i,j,2)+1;%% C -> Matlab
    end
end