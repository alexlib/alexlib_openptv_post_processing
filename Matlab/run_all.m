% the procedure to get the traj data
xuap = readXUAPFiles('../exp_031110/res_scene93/');
xuap = graphical_cleaning_xuap(xuap,'xy');
xuap = graphical_cleaning_xuap(xuap,'yz');
xuap = graphical_cleaning_xuap(xuap,'xz');
traj = xuap2traj(building_trajectories(xuap));
% save traj traj
trajLen = plot_long_trajectories(traj,3);
plot_xf_vs_t(traj)

return
%% 
id = [];
colors = {'r','g','b','k','m','c','y','r','y','r','g','b'};
figure, hold on
for i = 1:length(traj)
    if length(traj(i).xf)>3
        id = cat(1,id,i);
        for j = 1:length(traj(i).xf) 
            % plot(traj(i).t,traj(i).xf(j),'.','Color',colors{j});
            plot3(traj(i).xf(j),traj(i).yf(j),traj(i).zf(j),'o','Color',colors{j});
        end
    end
end
