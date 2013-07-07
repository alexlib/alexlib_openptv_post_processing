function plot_trajectory_properties_vs_t(traj,smooth)
if nargin < 2
    smooth = 0;
end
% plot acceleration vs time
%%
if smooth == 1
figure,
subplot(311)
box on, grid on, hold on,
for i = 1:length(traj)
    if length(traj(i).xf) > 5
        traj1(i).t = traj(i).t - traj(i).t(3);
        plot(traj1(i).t(3)/160,traj(i).xf(3),'ro');
        plot(traj1(i).t(3)/160,traj(i).yf(3),'go');
        plot(traj1(i).t(3)/160,traj(i).zf(3),'bo');

        l1=plot(traj1(i).t(3:end-2)/160,traj(i).xf(3:end-2),'r-','DisplayName',num2str(i));
        l2=plot(traj1(i).t(3:end-2)/160,traj(i).yf(3:end-2),'g--','DisplayName',num2str(i));
        l3=plot(traj1(i).t(3:end-2)/160,traj(i).zf(3:end-2),'b-.','DisplayName',num2str(i));

        plot(traj1(i).t(end-2)/160,traj(i).xf(end-2),'rs');
        plot(traj1(i).t(end-2)/160,traj(i).yf(end-2),'gs');
        plot(traj1(i).t(end-2)/160,traj(i).zf(end-2),'bs');
    end
end
hold off
xlabel('t [sec]','Interpreter','latex');
ylabel('x [mm]','Interpreter','latex');
legend([l1,l2,l3],{'x','y','z'})
%
subplot(312)
box on, grid on, hold on,
for i = 1:length(traj)
    if length(traj(i).xf) > 5
        traj1(i).t(i) = traj(i).t - traj(i).t(3);
        plot(traj1(i).t(3)/160,traj(i).uf(3),'ro');
        plot(traj1(i).t(3)/160,traj(i).vf(3),'go');
        plot(traj1(i).t(3)/160,traj(i).wf(3),'bo');

        l1=plot(traj1(i).t(3:end-2)/160,traj(i).uf(3:end-2),'r-','DisplayName',num2str(i));
        l2=plot(traj1(i).t(3:end-2)/160,traj(i).vf(3:end-2),'g--','DisplayName',num2str(i));
        l3=plot(traj1(i).t(3:end-2)/160,traj(i).wf(3:end-2),'b-.','DisplayName',num2str(i));

        plot(traj1(i).t(end-2)/160,traj(i).uf(end-2),'rs');
        plot(traj1(i).t(end-2)/160,traj(i).vf(end-2),'gs');
        plot(traj1(i).t(end-2)/160,traj(i).wf(end-2),'bs');
    end
end
hold off
xlabel('t [sec]','Interpreter','latex');
ylabel('u [mm/s]','Interpreter','latex');
legend([l1,l2,l3],{'u','v','w'})
%
subplot(313)
box on, grid on, hold on,
for i = 1:length(traj)
    if length(traj(i).xf) > 5
        traj1(i).t = traj(i).t - traj(i).t(3);
        plot(traj1(i).t(3)/160,traj(i).axf(3),'ro');
        plot(traj1(i).t(3)/160,traj(i).ayf(3),'go');
        plot(traj1(i).t(3)/160,traj(i).azf(3),'bo');

        l1=plot(traj1(i).t(3:end-2)/160,traj(i).axf(3:end-2),'r-','DisplayName',num2str(i));
        l2=plot(traj1(i).t(3:end-2)/160,traj(i).ayf(3:end-2),'g--','DisplayName',num2str(i));
        l3=plot(traj1(i).t(3:end-2)/160,traj(i).azf(3:end-2),'b-.','DisplayName',num2str(i));

        plot(traj1(i).t(end-2)/160,traj(i).axf(end-2),'rs');
        plot(traj1(i).t(end-2)/160,traj(i).ayf(end-2),'gs');
        plot(traj1(i).t(end-2)/160,traj(i).azf(end-2),'bs');
    end
end
hold off
xlabel('t [sec]','Interpreter','latex');
ylabel('a [mm/s$^2$]','Interpreter','latex');
legend([l1,l2,l3],{'a_x','a_y','a_z'})

else
%%
figure,
subplot(311)
box on, grid on, hold on,
for i = 1:length(traj)
        traj1(i).t = traj(i).t - traj(i).t(1);
        plot(traj1(i).t(1)/160,traj(i).xf(1),'r.','MarkerSize',4);
        plot(traj1(i).t(1)/160,traj(i).yf(1),'g.','MarkerSize',4);
        plot(traj1(i).t(1)/160,traj(i).zf(1),'b.','MarkerSize',4);

        l1=plot(traj1(i).t/160,traj(i).xf,'r-','DisplayName',num2str(i));
        l2=plot(traj1(i).t/160,traj(i).yf,'g--','DisplayName',num2str(i));
        l3=plot(traj1(i).t/160,traj(i).zf,'b-.','DisplayName',num2str(i));

        plot(traj1(i).t(end)/160,traj(i).xf(end),'r.','MarkerSize',4);
        plot(traj1(i).t(end)/160,traj(i).yf(end),'g.','MarkerSize',4);
        plot(traj1(i).t(end)/160,traj(i).zf(end),'b.','MarkerSize',4);
end
hold off
xlabel('t [sec]','Interpreter','latex');
ylabel('x [mm]','Interpreter','latex');
legend([l1,l2,l3],{'x','y','z'})
%
subplot(312)
box on, grid on, hold on,
for i = 1:length(traj)
    if length(traj(i).xf) > 5
        traj1(i).t = traj(i).t - traj(i).t(1);
        plot(traj1(i).t(1)/160,traj(i).uf(1),'r.','MarkerSize',4);
        plot(traj1(i).t(1)/160,traj(i).vf(1),'g.','MarkerSize',4);
        plot(traj1(i).t(1)/160,traj(i).wf(1),'b.','MarkerSize',4);

        l1=plot(traj1(i).t/160,traj(i).uf,'r-','DisplayName',num2str(i));
        l2=plot(traj1(i).t/160,traj(i).vf,'g--','DisplayName',num2str(i));
        l3=plot(traj1(i).t/160,traj(i).wf,'b-.','DisplayName',num2str(i));

        plot(traj1(i).t(end)/160,traj(i).uf(end),'r.','MarkerSize',4);
        plot(traj1(i).t(end)/160,traj(i).vf(end),'g.','MarkerSize',4);
        plot(traj1(i).t(end)/160,traj(i).wf(end),'b.','MarkerSize',4);
    end
end
hold off
xlabel('t [sec]','Interpreter','latex');
ylabel('u [mm/s]','Interpreter','latex');
legend([l1,l2,l3],{'u','v','w'})
%
subplot(313)
box on, grid on, hold on,
for i = 1:length(traj)
    if length(traj(i).xf) > 5
        traj1(i).t = traj(i).t - traj(i).t(1);
        plot(traj1(i).t(1)/160,traj(i).axf(1),'r.','MarkerSize',4);
        plot(traj1(i).t(1)/160,traj(i).ayf(1),'g.','MarkerSize',4);
        plot(traj1(i).t(1)/160,traj(i).azf(1),'b.','MarkerSize',4);

        l1=plot(traj1(i).t/160,traj(i).axf,'r-','DisplayName',num2str(i));
        l2=plot(traj1(i).t/160,traj(i).ayf,'g--','DisplayName',num2str(i));
        l3=plot(traj1(i).t/160,traj(i).azf,'b-.','DisplayName',num2str(i));

        plot(traj1(i).t(end)/160,traj(i).axf(end),'r.','MarkerSize',4);
        plot(traj1(i).t(end)/160,traj(i).ayf(end),'g.','MarkerSize',4);
        plot(traj1(i).t(end)/160,traj(i).azf(end),'b.','MarkerSize',4);
    end
end
hold off
xlabel('t [sec]','Interpreter','latex');
ylabel('a [mm/s$^2$]','Interpreter','latex');
legend([l1,l2,l3],{'a_x','a_y','a_z'})
end