function plot_trajectory_properties_vs_t_2(traj,dt)
col=['r';'b';'k';'g';'m';'c'];
    
% if nargin < 2
%     smooth = 0;
% end
% plot acceleration vs time
%%
% if smooth == 1
figure,
subplot(221)
box on, grid on, hold on,
for i = 1:length(traj)
     randcolor = ceil(length(col)*rand(1));
    if length(traj(i).xf) > 5
        plot(traj(i).t(3)/dt,(-1)*traj(i).xf(3)/10,'Color',col(randcolor));
        plot(traj(i).t(end-2)/dt,(-1)*traj(i).xf(end-2)/10,'Color',col(randcolor));
      
        
        l1=plot(traj(i).t(3:end-2)/dt,traj(i).xf(3:end-2)/10,'Color',col(randcolor),'DisplayName',num2str(i));
    end
end
hold off
xlabel('$t$ [sec]','Interpreter','latex');
ylabel('$x$ [cm]','Interpreter','latex');

%
subplot(222)

box on, grid on, hold on,
for i = 1:length(traj)
    if length(traj(i).xf) > 5
      randcolor = ceil(length(col)*rand(1));
        plot(traj(i).t(3)/dt,(-1)*traj(i).yf(3)/10,'Color',col(randcolor));
        plot(traj(i).t(end-2)/dt,(-1)*traj(i).yf(end-2)/10,'Color',col(randcolor));
       
          l2=plot(traj(i).t(3:end-2)/dt,traj(i).yf(3:end-2)/10,'g--','DisplayName',num2str(i));
    end
end
hold off
xlabel('$t$ [sec]','Interpreter','latex');
ylabel('$y$ [cm]','Interpreter','latex');


subplot(223)
box on, grid on, hold on,
for i = 1:length(traj)
    if length(traj(i).xf) > 5
        
         randcolor = ceil(length(col)*rand(1));
        
        plot(traj(i).t(3)/dt,(-1)*traj(i).uf(3)*dt/10,'Color',col(randcolor));
        plot(traj(i).t(end-2)/dt,(-1)*traj(i).uf(end-2)*dt/10,'Color',col(randcolor));
       
           l1=plot(traj(i).t(3:end-2)/dt,traj(i).uf(3:end-2)*dt/10,'r-','DisplayName',num2str(i));
    end
end
hold off
xlabel('$t$ [sec]','Interpreter','latex');
ylabel('$ux$ [cm/s]','Interpreter','latex');


subplot(224)
box on, grid on, hold on,
for i = 1:length(traj)
    if length(traj(i).xf) > 5
            randcolor = ceil(length(col)*rand(1));
            
        plot(traj(i).t(3)/dt,traj(i).vf(3)*dt/10,'g.');
      
        plot(traj(i).t(end-2)/dt,traj(i).vf(end-2)*dt/10,'g.');
       
         l2=plot(traj(i).t(3:end-2)/dt,traj(i).vf(3:end-2)*dt/10,'g--','DisplayName',num2str(i));
    end
end
hold off
xlabel('$t$ [sec]','Interpreter','latex');
ylabel('$ur$ [cm/s]','Interpreter','latex');


