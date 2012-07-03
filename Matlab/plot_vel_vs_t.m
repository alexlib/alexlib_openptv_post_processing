function plot_vel_vs_t(traj)
% plot velocity vs time 

% colors = {'r','y','m','k','c','g','b'};
figure, hold on,
for i = 1:length(traj)
    %     for j = 1:length(traj(i).xf),
    %         plot(traj(i).t(j)/50,traj(i).xf(j),'r.');
    %         plot(traj(i).t(j)/50,traj(i).yf(j),'g.');
    %         plot(traj(i).t(j)/50,traj(i).zf(j),'b.');
    %     end
    plot(traj(i).t(1)/50,traj(i).uf(1),'ro');
    plot(traj(i).t(1)/50,traj(i).vf(1),'go');
    plot(traj(i).t(1)/50,traj(i).wf(1),'bo');
    
    l1=plot(traj(i).t/50,traj(i).uf,'r-','DisplayName',num2str(i));
    l2=plot(traj(i).t/50,traj(i).vf,'g--','DisplayName',num2str(i));
    l3=plot(traj(i).t/50,traj(i).wf,'b-.','DisplayName',num2str(i));
    
    
    plot(traj(i).t(end)/50,traj(i).uf(end),'rs');
    plot(traj(i).t(end)/50,traj(i).vf(end),'gs');
    plot(traj(i).t(end)/50,traj(i).wf(end),'bs');
    

end
hold off
xlabel('$t$ [sec]');
ylabel('$u$ [m/sec]');
legend([l1,l2,l3],{'$u$','$v$','$w$'})

