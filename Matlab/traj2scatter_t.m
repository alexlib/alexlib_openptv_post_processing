function [ scatter_t ] = traj2scatter_t( traj ) 
%Temporal transformation 
%This function is used to create a scatter data of the 3D position and velocities 
%of the particles from the trajectory data obtained by ptv is to traj .m 
%Its required inputs are the traj data . 
%Input:[ traj] 
%Output:[scatter data respect to time] 
%Jin Tae Kim 2014 University of Illinois at Urbana Champaign 
% Modified 2018 Alex Liberzon, Tel Aviv University

t_min = traj(1).t(1) ; 
t_max = traj(end).t(end) ; 
idl = t_max - t_min+1;

scatter_t = repmat(struct( 'x' ,[] , 'y' ,[] , 'z' ,[] , 'u' ,[] , 'v' ,[] , 'w' ,[] ,'ax',[],'ay',[],'az',[],'trajid',[]),idl,1); 

for i=1:length(traj) 
    for j=1:length(traj(i).xf)
        k = traj(i).t(j) - t_min + 1; % this line was changed
        scatter_t(k).x = [scatter_t(k).x; traj(i).xf(j)]; % dynamically growing?
        scatter_t(k).y = [scatter_t(k).y; traj(i).yf(j)]; 
        scatter_t(k).z  = [scatter_t(k).z; traj(i).zf(j)]; 
        scatter_t(k).u  = [scatter_t(k).u; traj(i).uf(j)]; 
        scatter_t(k).v  = [scatter_t(k).v; traj(i).vf(j)]; 
        scatter_t(k).w  = [scatter_t(k).w; traj(i).wf(j)]; 
        scatter_t(k).ax = [scatter_t(k).ax; traj(i).axf(j)]; 
        scatter_t(k).ay = [scatter_t(k).ay; traj(i).ayf(j)]; 
        scatter_t(k).az = [scatter_t(k).az; traj(i).azf(j)]; 
        scatter_t(k).trajid = [scatter_t(k).trajid; traj(i).trajid(j)]; 
    end 
end
