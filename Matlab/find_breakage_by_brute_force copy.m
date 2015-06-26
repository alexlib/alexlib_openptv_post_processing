function [minDist] = find_breakage_by_brute_force(traj,maxDt,distThresh)
% HAITAO_LINKNG_CRITERIA builds containers of trajectories that can be
% linked together. Simple criteria of distance between the end of the
% previous trajectory and the beginning of the next trajectory is used.
% no check on velocity, acceleration, etc.


% first container is the first trajectory
k = 0;

% tmp = [];

trajNum = length(traj);

% for the test
% trajNum = 300;

% linkid = cell(trajNum*100,3);

minDist = Inf;


%%
for i = 1:trajNum - 1
    %     if traj(i).t(end) > 21040 && length(traj(i).xf) >= 5
    %     for j = i+numTrajFrame:min(trajNum,i+numTrajFrame*2) % <--- hard coded number, need something more
    k = 0;
    candidate = [];
    candidateDist = [];
    for j = i+1:trajNum
        dt = traj(j).t(1) - traj(i).t(end);
%         if dt < 0
%             tmp = traj(i);
%             traj(i) = traj(j);
%             traj(j) = tmp;
%             dt = traj(j).t(1) - traj(i).t(end);
%         end
        if dt > 0 && dt <= maxDt && length(traj(i).xf) > 5 && length(traj(j).xf) > 5
            predicted_x = interp1(traj(i).t,traj(i).xf,traj(i).t(end)+dt,'pchip','extrap');
            predicted_y = interp1(traj(i).t,traj(i).yf,traj(i).t(end)+dt,'pchip','extrap');
            predicted_z = interp1(traj(i).t,traj(i).zf,traj(i).t(end)+dt,'pchip','extrap');
            predicted_x1 = interp1(traj(j).t,traj(j).xf,traj(j).t(1)-dt,'pchip','extrap');
            predicted_y1 = interp1(traj(j).t,traj(j).yf,traj(j).t(1)-dt,'pchip','extrap');
            predicted_z1 = interp1(traj(j).t,traj(j).zf,traj(j).t(1)-dt,'pchip','extrap');
            dist1 = norm([predicted_x,predicted_y,predicted_z] - [traj(j).xf(1),traj(j).yf(1),traj(j).zf(1)]);
            dist2 = norm([predicted_x1,predicted_y1,predicted_z1]- [traj(i).xf(end),traj(i).yf(end),traj(i).zf(end)]);
            if dist1 < minDist
                minDist = dist1;
            end
            %             tmp = cat(1,tmp,dist);
            if dist1 < distThresh && dist2 < distThresh
                k = k + 1;
                candidate(k) = j;
                candidateDist(k) = max([dist1,dist2]);
            end
        end
    end

    if k > 1 % at least 2 candidates
        % check if they come from the same place:
        for m = 1:k-1
            x1 = interp1(traj(candidate(m)).t,traj(candidate(m)).xf,traj(i).t(end),'pchip','extrap');
            y1 = interp1(traj(candidate(m)).t,traj(candidate(m)).yf,traj(i).t(end),'pchip','extrap');
            z1 = interp1(traj(candidate(m)).t,traj(candidate(m)).zf,traj(i).t(end),'pchip','extrap');
            for n = 2:k
                x2 = interp1(traj(candidate(n)).t,traj(candidate(n)).xf,traj(i).t(end),'pchip','extrap');
                y2 = interp1(traj(candidate(n)).t,traj(candidate(n)).yf,traj(i).t(end),'pchip','extrap');
                z2 = interp1(traj(candidate(n)).t,traj(candidate(n)).zf,traj(i).t(end),'pchip','extrap');
                if norm([x1,y1,z1]-[x2,y2,z2]) < distThresh
                    plot_long_trajectories(traj([i,candidate(m),candidate(n)]),5);
                    hold on
                    plot3(x1,z1,y1,'rs');
                    plot3(x2,z2,y2,'ms');                  
                    title(sprintf('%5d %5d %5d',i,candidate(m),candidate(n)));
                end
            end
        end

    end
end
