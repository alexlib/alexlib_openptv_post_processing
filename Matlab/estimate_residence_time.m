cd('/Users/alex/Dropbox/PTV (1)/trajectories/')
lst = dir('*.mat')
dt = 1/100;
colors = {'r','g','b','k','m'};
markers = {'o','s','^','*'};

dist = @(x,y,z) sqrt((x - x(1)).^2 + (y - y(1)).^2 + (z - z(1)).^2);



for f = 1:length(lst)
    disp(lst(f).name)
    clear('t','res_time','mean_res_time','traj')
    tmp = load(lst(f).name);
    fld = fieldnames(tmp);
    traj = tmp.(fld{1});
    
    figure(f), hold on
    
    k = 0;
    
    
    for r0 = 10e-3:50e-3:250e-3  % millimeter, i.e. 200 micron
        
        
        res_time = zeros(length(traj),1);
        t = zeros(length(traj),1);
        
        for i = 1:length(traj)
            d = dist(traj(i).xf,traj(i).yf,traj(i).zf);
            rt = traj(i).t(find(d > r0,1,'first'));
            t(i) = traj(i).t(1);
            if ~isempty(rt)
                res_time(i) = rt - traj(i).t(1);
            else
                res_time(i) = NaN;
            end
        end
        %
        % t = cat(1,traj.t);
        % t = t - t(1);
        % dt = 1/100 % sec
        t = t*dt;
        
        k = k + 1;
        plot(t,res_time*dt,colors{k});
        mean_res_time(k) = nanmean(res_time);
    end
    xlabel('Time [sec]')
    ylabel('Residence time [sec]')
    legend(num2str(mean_res_time'))
    
    figure(length(lst)+1), hold on
    plot(10e-3:50e-3:250e-3,mean_res_time,markers{f})
    
end

figure(length(lst)+1); 
xlabel('r_0 [mm]')
ylabel('Mean residence time [sec]')
tmp = {};
[tmp{1:length(lst),1}] = deal(lst.name);
legend(tmp);




% maxdisp = zeros(length(traj),1);
%
% for i = 1:length(traj)
%     mxf = max((traj(i).xf - traj(i).xf(1)).^2);
%     myf = max((traj(i).yf - traj(i).yf(1)).^2);
%     maxdisp(i) = sqrt(mxf + myf);
% end
