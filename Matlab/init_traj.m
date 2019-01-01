function [traj,num_in_traj] = init_traj(traj, i,j)

% global traj;

traj=[0 0];
traj(1,1)=i;
traj(1,2)=j;
num_in_traj=1;