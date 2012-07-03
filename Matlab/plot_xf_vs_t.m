function plot_xf_vs_t(traj)
% plot vs time shows that there are many zeros and non-single particle
% frames

% colors = {'r','y','m','k','c','g','b'};
figure, hold on,
for i = 1:length(traj)
    %     for j = 1:length(traj(i).xf),
    %         plot(traj(i).t(j)/50,traj(i).xf(j),'r.');
    %         plot(traj(i).t(j)/50,traj(i).yf(j),'g.');
    %         plot(traj(i).t(j)/50,traj(i).zf(j),'b.');
    %     end
    plot(traj(i).t(1)/50,traj(i).xf(1),'r.');
    plot(traj(i).t(1)/50,traj(i).yf(1),'g.');
    plot(traj(i).t(1)/50,traj(i).zf(1),'b.');
    
    l1=plot(traj(i).t/50,traj(i).xf,'r-','LineWidth',2,'DisplayName',num2str(i));
    l2=plot(traj(i).t/50,traj(i).yf,'g--','LineWidth',2,'DisplayName',num2str(i));
    l3=plot(traj(i).t/50,traj(i).zf,'b-.','LineWidth',2,'DisplayName',num2str(i));
    
    
    plot(traj(i).t(end)/50,traj(i).xf(end),'rx');
    plot(traj(i).t(end)/50,traj(i).yf(end),'gx');
    plot(traj(i).t(end)/50,traj(i).zf(end),'bx');
    

end
hold off
xlabel('$t$ [sec]','Interpreter','Latex');
ylabel('$x$ [mm]','Interpreter','Latex');
legend([l1,l2,l3],{'$x$','$y$','$z$'})
