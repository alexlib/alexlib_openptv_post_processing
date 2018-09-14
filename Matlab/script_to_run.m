%script_to_run_
directory = '../myexp_alex_run/res';
first = 220;
last = 240;
minlength = 1;
dt = 1/100;

[traj] = ptv_is_to_traj(directory,first,last,minlength,dt);

plot_long_trajectories(traj,5);


% Eulerian approach
% read it the processed data in traj (otherwise we do not have velocity)
% and store in a frame by frame type.
% ! Warning - the present function is not optimal as it grows memory
% dynamically
% we could possibly improve it by providing a maximum likely number of
% particles and initializing a bit wasted memory rather than growing it
% during the run. #TODO

[ scatter_t ] = traj2scatter_t( traj );

% gridding from min to max, copied from scatter_to_grid
grid = scatter_to_grid (scatter_t,5);

figure
quiver3(grid.x,grid.y,grid.z,mean(grid.u,4),mean(grid.v,4),mean(grid.w,4),5)
