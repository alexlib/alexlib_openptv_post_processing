function plot_xf_vs_t(traj)
% plot vs time shows that there are many zeros and non-single particle
% frames

colors = {'r','y','m','k','c','g','b'};
figure, hold on, 
for i = 1:length(traj)
    for j = 1:length(traj(i).xf), 
        plot(traj(i).t/50,traj(i).xf(j),'.','Color',colors{j}); 
    end
end
hold off
xlabel('t [sec]');
ylabel('x [mm]');
