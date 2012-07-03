function plot_acc_vs_t(traj)
% plot acceleration vs time 

% colors = {'r','y','m','k','c','g','b'};
figure, hold on,
for i = 1:length(traj)
    %     for j = 1:length(traj(i).xf),
    %         plot(traj(i).t(j)/50,traj(i).xf(j),'r.');
    %         plot(traj(i).t(j)/50,traj(i).yf(j),'g.');
    %         plot(traj(i).t(j)/50,traj(i).zf(j),'b.');
    %     end
    plot(traj(i).t(3)/50,traj(i).axf(3),'ro');
    plot(traj(i).t(3)/50,traj(i).ayf(3),'go');
    plot(traj(i).t(3)/50,traj(i).azf(3),'bo');
    
    l1=plot(traj(i).t(3:end-2)/50,traj(i).axf(3:end-2),'r-','DisplayName',num2str(i));
    l2=plot(traj(i).t(3:end-2)/50,traj(i).ayf(3:end-2),'g--','DisplayName',num2str(i));
    l3=plot(traj(i).t(3:end-2)/50,traj(i).azf(3:end-2),'b-.','DisplayName',num2str(i));
    
    
    plot(traj(i).t(end-2)/50,traj(i).axf(end-2),'rs');
    plot(traj(i).t(end-2)/50,traj(i).ayf(end-2),'gs');
    plot(traj(i).t(end-2)/50,traj(i).azf(end-2),'bs');
    

end
hold off
xlabel('$t$ [sec]','Interpreter','latex');
ylabel('$a$ [m/s$^2$]','Interpreter','latex');
legend([l1,l2,l3],{'$a_x$','$a_y$','$a_z$'})
