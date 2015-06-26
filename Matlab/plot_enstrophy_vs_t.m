function hf = plot_enstrophy_vs_t(traj,fps)
% plot enstrophy and separate vorticity components
% per trajectory in time
% Usage:
%       hf = plot_enstrophy_vs_t(traj,fps)
%
% default 'fps' = 100
% the first 2 points are cropped
if nargin == 1
    fps = 100; % default frames-per-second
end


%%
hf = figure;
subplot(211)
box on, grid on, hold on,
for i = 1:length(traj)
    if length(traj(i).xf) > 5
        traj(i).t = traj(i).t - traj(i).t(3);
%         plot(traj(i).t(3)/fps,traj(i).omegax(3),'ro');
%         plot(traj(i).t(3)/fps,traj(i).omegay(3),'go');
%         plot(traj(i).t(3)/fps,traj(i).omegaz(3),'bo');
        
        l1=plot(traj(i).t(3:end-2)/fps,traj(i).omegax(3:end-2),'r-','DisplayName',num2str(i));
        l2=plot(traj(i).t(3:end-2)/fps,traj(i).omegay(3:end-2),'g--','DisplayName',num2str(i));
        l3=plot(traj(i).t(3:end-2)/fps,traj(i).omegaz(3:end-2),'b-.','DisplayName',num2str(i));
        
        plot(traj(i).t(end-2)/fps,traj(i).omegax(end-2),'rs');
        plot(traj(i).t(end-2)/fps,traj(i).omegay(end-2),'gs');
        plot(traj(i).t(end-2)/fps,traj(i).omegaz(end-2),'bs');
    end
end
hold off
xlabel('t [sec]');
ylabel('x [mm]');
legend([l1,l2,l3],{'\omega_x','\omega_y','\omega_z'})

subplot(212)
box on, grid on, hold on,
for i = 1:length(traj)
    if length(traj(i).omegax) > 5
        c = rand(3,1);
        ens =  traj(i).omegax.^2 + traj(i).omegay.^2 + traj(i).omegaz.^2 ;
        traj(i).t = traj(i).t - traj(i).t(3);
        l1=plot(traj(i).t(3:end-2)/fps,ens(3:end-2),'LineStyle','-','Color',c,'DisplayName',num2str(i));
        plot(traj(i).t(end-2)/fps,ens(end-2));
    end
end
hold off
xlabel('t [sec]');
ylabel('u [mm/s]');
legend(l1,'\omega^2')
%

