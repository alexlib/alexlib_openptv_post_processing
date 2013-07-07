function plot_trajectory_properties_vs_t_normalized(traj,dt,scene_num)
% if nargin < 2
%     smooth = 0;
% end
% plot acceleration vs time
%%
% if smooth == 1
figure,
subplot(211)
box on, grid on, hold on,

xf=(-1)*cat(1,traj.xf)/10;

for i = 1:length(traj)
    if length(traj(i).xf) > 5
        x = traj(i).xf(1)/10;
        y = traj(i).yf(1)/10;
        %r = sqrt(traj(i).yf.^2 + traj(i).zf.^2/10)/10;
        %r0 = r(1);
        
        plot(traj(i).t(2)/dt,(-1)*(traj(i).xf(2)/10-x)/abs(mean(xf)),'r.')
        plot(traj(i).t(2)/dt,(-10)*(traj(i).yf(2)/10-y)/abs(mean(xf)),'g.');
     

       l1=plot(traj(i).t(3:end-1)/dt,(-1)*(traj(i).xf(3:end-1)/10-x)/abs(mean(xf)),'r--','DisplayName',num2str(i));
       l2=plot(traj(i).t(3:end-1)/dt,(-10)*(traj(i).yf(3:end-1)/10-y)/abs(mean(xf)),'g--','DisplayName',num2str(i));
        
        %l2 = plot(traj(i).t(3:end-1)/dt,r(3:end-1)-r0,'g--','DisplayName',num2str(i));
         %l2=plot(traj(i).t(3:end-1)/dt,(-1)*(traj(i).yf(3:end-1)/10-y),'g--','DisplayName',num2str(i));
        %l3=plot(traj(i).t(3:end-2)/dt,traj(i).zf(3:end-2)/10,'b-.','DisplayName',num2str(i));

        plot(traj(i).t(end)/dt,(-1)*(traj(i).xf(end)/10-x)/abs(mean(xf)),'r*');
        plot(traj(i).t(end)/dt,(-10)*(traj(i).yf(end)/10-y)/abs(mean(xf)),'g*','MarkerSize',3);
      
        

        
        %plot(traj(i).t(end)/dt,(-1)*(traj(i).yf(end)/10-y),'g*','MarkerSize',3);
      %  plot(traj(i).t(end-2)/dt,traj(i).zf(end-2)/10,'b.');
    end
end
hold off
xlabel(' [sec]');
legend([l1,l2],{'x/x_mean','10*y/x_mean'});
title(sprintf('scene %g',scene_num));

%
subplot(212)
box on, grid on, hold on,
ux=(-1)*cat(1,traj.xf)*dt/10;
for i = 1:length(traj)
    if length(traj(i).xf) > 5
        plot(traj(i).t(3)/dt,(-1)*(traj(i).uf(3)*dt/10)/abs(mean(ux)),'r.');
        plot(traj(i).t(3)/dt,(-10)*(traj(i).vf(3)*dt/10)/abs(mean(ux)),'g.');
       
      %  plot(traj(i).t(3)/dt,traj(i).wf(3)*dt/10,'b.');

        l1=plot(traj(i).t(3:end-2)/dt,(-1)*(traj(i).uf(3:end-2)*dt/10)/abs(mean(ux)),'r-','DisplayName',num2str(i));
        l2=plot(traj(i).t(3:end-2)/dt,(-10)*(traj(i).vf(3:end-2)*dt/10)/abs(mean(ux)),'g--','DisplayName',num2str(i));
      %  l3=plot(traj(i).t(3:end-2)/dt,traj(i).wf(3:end-2)*dt/10,'b-.','DisplayName',num2str(i));

        plot(traj(i).t(end-2)/dt,(-1)*(traj(i).uf(end-2)*dt/10)/abs(mean(ux)),'r*');
        plot(traj(i).t(end-2)/dt,(-10)*(traj(i).vf(end-2)*dt/10)/abs(mean(ux)),'g*');
       % plot(traj(i).t(end-2)/dt,traj(i).wf(end-2)*dt/10,'b.');
    end
end
hold off
xlabel('[sec]');
legend([l1,l2],{'u_x/u_mean','10*u_y/u_mean'});
% 
% %
% subplot(313)
% box on, grid on, hold on,
% 
% ax=(-1)*cat(1,traj.xf)*(dt^2)/10;
% for i = 1:length(traj)
%     if length(traj(i).xf) > 5
%         plot(traj(i).t(3)/dt,(-1)*(traj(i).axf(3)*dt^2/10)/mean(ax),'r.');
%         plot(traj(i).t(3)/dt,(-1)*(traj(i).ayf(3)*dt^2/10)/mean(ax),'g.');
% % %         plot(traj(i).t(3)/dt,traj(i).azf(3)/10,'b.')*dt^2;
% % % 
%          l1=plot(traj(i).t(3:end-2)/dt,(-1)*(traj(i).axf(3:end-2)*dt^2/10)/mean(ax),'r-','DisplayName',num2str(i));
%          l2=plot(traj(i).t(3:end-2)/dt,(-1)*(traj(i).ayf(3:end-2)*dt^2/10)/mean(ax),'g--','DisplayName',num2str(i));
% % %         l3=plot(traj(i).t(3:end-2)/dt,traj(i).azf(3:end-2)*dt^2/10,'b-.','DisplayName',num2str(i));
% % % 
%           plot(traj(i).t(end-2)/dt,(-1)*(traj(i).axf(end-2)*dt^2/10)/mean(ax),'r*');
%           plot(traj(i).t(end-2)/dt,(-1)*(traj(i).ayf(end-2)*dt^2/10)/mean(ax),'g*');
% % %         plot(traj(i).t(end-2)/dt,traj(i).azf(end-2),'b.')*dt^2/10;
%      end
%  end
%  hold off
%  xlabel('[sec]');
% legend([l1,l2],{'ax/ax_max','ay/ax_max'});
% 
% % else
% % %%
% % figure,
% % subplot(311)
% % box on, grid on, hold on,
% % for i = 1:length(traj)
% %         plot(traj(i).t(1)/dt,traj(i).xf(1),'r.','MarkerSize',4);
% %         plot(traj(i).t(1)/dt,traj(i).yf(1),'g.','MarkerSize',4);
% %         plot(traj(i).t(1)/dt,traj(i).zf(1),'b.','MarkerSize',4);
% % 
% %         l1=plot(traj(i).t/dt,traj(i).xf,'r-','DisplayName',num2str(i));
% %         l2=plot(traj(i).t/dt,traj(i).yf,'g--','DisplayName',num2str(i));
% %         l3=plot(traj(i).t/dt,traj(i).zf,'b-.','DisplayName',num2str(i));
% % 
% %         plot(traj(i).t(end)/dt,traj(i).xf(end),'r.','MarkerSize',4);
% %         plot(traj(i).t(end)/dt,traj(i).yf(end),'g.','MarkerSize',4);
% %         plot(traj(i).t(end)/dt,traj(i).zf(end),'b.','MarkerSize',4);
% % end
% % hold off
% % xlabel('$t$ [sec]','Interpreter','latex');
% % ylabel('$x$ [m]','Interpreter','latex');
% % legend([l1,l2,l3],{'$x$','$y$','$z$'})
% % %
% % subplot(312)
% % box on, grid on, hold on,
% % for i = 1:length(traj)
% %     if length(traj(i).xf) > 5
% %         plot(traj(i).t(1)/dt,traj(i).uf(1)*dt,'r.','MarkerSize',4);
% %         plot(traj(i).t(1)/dt,traj(i).vf(1)*dt,'g.','MarkerSize',4);
% %         plot(traj(i).t(1)/dt,traj(i).wf(1)*dt,'b.','MarkerSize',4);
% % 
% %         l1=plot(traj(i).t/dt,traj(i).uf*dt,'r-','DisplayName',num2str(i));
% %         l2=plot(traj(i).t/dt,traj(i).vf*dt,'g--','DisplayName',num2str(i));
% %         l3=plot(traj(i).t/dt,traj(i).wf*dt,'b-.','DisplayName',num2str(i));
% % 
% %         plot(traj(i).t(end)/dt,traj(i).uf(end)*dt,'r.','MarkerSize',4);
% %         plot(traj(i).t(end)/dt,traj(i).vf(end)*dt,'g.','MarkerSize',4);
% %         plot(traj(i).t(end)/dt,traj(i).wf(end)*dt,'b.','MarkerSize',4);
% %     end
% % end
% % hold off
% % xlabel('$t$ [sec]','Interpreter','latex');
% % ylabel('$u$ [m/s]','Interpreter','latex');
% % legend([l1,l2,l3],{'$u$','$v$','$w$'})
% % %
% % subplot(313)
% % box on, grid on, hold on,
% % for i = 1:length(traj)
% %     if length(traj(i).xf) > 5
% %         plot(traj(i).t(1)/dt,traj(i).axf(1)*dt^2,'r.','MarkerSize',4);
% %         plot(traj(i).t(1)/dt,traj(i).ayf(1)*dt^2,'g.','MarkerSize',4);
% %         plot(traj(i).t(1)/dt,traj(i).azf(1)*dt^2,'b.','MarkerSize',4);
% % 
% %         l1=plot(traj(i).t/dt,traj(i).axf*dt^2,'r-','DisplayName',num2str(i));
% %         l2=plot(traj(i).t/dt,traj(i).ayf*dt^2,'g--','DisplayName',num2str(i));
% %         l3=plot(traj(i).t/dt,traj(i).azf*dt^2,'b-.','DisplayName',num2str(i));
% % 
% %         plot(traj(i).t(end)/dt,traj(i).axf(end)*dt^2,'r.','MarkerSize',4);
% %         plot(traj(i).t(end)/dt,traj(i).ayf(end)*dt^2,'g.','MarkerSize',4);
% %         plot(traj(i).t(end)/dt,traj(i).azf(end)*dt^2,'b.','MarkerSize',4);
% %     end
% % end
% % hold off
% % xlabel('$t$ [sec]','Interpreter','latex');
% % ylabel('$a$ [m/s$^2$]','Interpreter','latex');
% % legend([l1,l2,l3],{'$a_x$','$a_y$','$a_z$'})
% % end