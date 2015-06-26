function [minDist] = find_breakage_by_brute_force_spline(traj,maxDt,distThresh)


trajNum = length(traj);

minDist = Inf;

%% 
% i = 72
% f = csaps(traj(i).t,[traj(i).xf,traj(i).yf,traj(i).zf]')
% g = fnxtr(f)
% fnplt(g,[traj(i).t(1)-3,traj(i).t(end)+3],'r')
% fnval(g,traj(i).t(end)+dt)

%%
for i = 1:trajNum - 1
    %     if traj(i).t(end) > 21040 && length(traj(i).xf) >= 5
    %     for j = i+numTrajFrame:min(trajNum,i+numTrajFrame*2) % <--- hard coded number, need something more
    k = 0;
    candidate = [];
    candidateDist = [];
    for j = i+1:trajNum
        dt = traj(j).t(1) - traj(i).t(end);
        if dt > 0 && dt <= maxDt && length(traj(i).xf) > 5 && length(traj(j).xf) > 5
            predicted_x = extrapolate(traj(i),traj(i).t(end)+dt);
            predicted_x1 = extrapolate(traj(j),traj(j).t(1)-dt);
            
            dist1 = norm(predicted_x' - [traj(j).xf(1),traj(j).yf(1),traj(j).zf(1)]);
            dist2 = norm(predicted_x1'- [traj(i).xf(end),traj(i).yf(end),traj(i).zf(end)]);
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
            x1 = extrapolate(traj(candidate(m)),traj(i).t(end));
            for n = 2:k
                x2 = extrapolate(traj(candidate(n)),traj(i).t(end));
                if norm([x1]-[x2]) < distThresh
                    plot_long_trajectories(traj([i,candidate(m),candidate(n)]),5);
                    hold on
                    plot3(x1(1),x1(2),x1(3),'rs');
                    plot3(x2(1),x2(2),x2(3),'ms');                  
                    title(sprintf('%5d %5d %5d',i,candidate(m),candidate(n)));
                end
            end
        end

    end
end

function x = extrapolate(traj,t)
x = fnval(fnxtr(csaps(traj.t,[traj.xf,traj.yf,traj.zf]')),t);

