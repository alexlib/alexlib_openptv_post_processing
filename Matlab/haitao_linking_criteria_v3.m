function [traj,removelist] = haitao_linking_criteria_v3(traj,maxDt,distThresh)
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

removelist = [];


%%
for i = 1:trajNum - 1
    %     for j = i+numTrajFrame:min(trajNum,i+numTrajFrame*2) % <--- hard coded number, need something more
    minDist = distThresh;
    for j = i+1:trajNum
        dt = traj(j).t(1) - traj(i).t(end);
        if dt > 0 && dt <= maxDt
            if length(traj(i).xf) < 3
                predicted_x = traj(i).xf(end) + dt*traj(i).uf(end);
                predicted_y = traj(i).yf(end) + dt*traj(i).vf(end);
                predicted_z = traj(i).zf(end) + dt*traj(i).wf(end);

            else
                predicted_x = interp1(traj(i).t,traj(i).xf,traj(i).t(end)+dt,'cubic','extrap');
                predicted_y = interp1(traj(i).t,traj(i).yf,traj(i).t(end)+dt,'cubic','extrap');
                predicted_z = interp1(traj(i).t,traj(i).zf,traj(i).t(end)+dt,'cubic','extrap');
            end

            if length(traj(j).xf) < 3
                predicted_x1 = traj(j).xf(1) - dt*traj(j).uf(1);
                predicted_y1 = traj(j).xf(1) - dt*traj(j).vf(1);
                predicted_z1 = traj(j).xf(1) - dt*traj(j).wf(1);
            else
                predicted_x1 = interp1(traj(j).t,traj(j).xf,traj(j).t(1)-dt,'cubic','extrap');
                predicted_y1 = interp1(traj(j).t,traj(j).yf,traj(j).t(1)-dt,'cubic','extrap');
                predicted_z1 = interp1(traj(j).t,traj(j).zf,traj(j).t(1)-dt,'cubic','extrap');

            end
            dist = norm([predicted_x,predicted_y,predicted_z]- [traj(j).xf(1),traj(j).yf(1),traj(j).zf(1)]) + ...
                norm([predicted_x1,predicted_y1,predicted_z1]- [traj(i).xf(end),traj(i).yf(end),traj(i).zf(end)]);
%             tmp = cat(1,tmp,dist);
            if dist < minDist
                minDist = dist;
                linkid = j;
            end
        end
    end

    
    
    if minDist < distThresh
%          hold on;
%          plot(traj(i).t,traj(i).xf,'r-o',traj(i).t(end)+dt,predicted_x,'rs',traj(linkid).t(1)-dt,predicted_x1,'b^',traj(linkid).t,traj(linkid).xf,'b--^')
%          title(norm([predicted_x,predicted_y,predicted_z]- [traj(linkid).xf(1),traj(linkid).yf(1),traj(linkid).zf(1)]))

         traj(linkid) = link_trajectories_rbf(combine_trajectories(traj(i),traj(linkid)));
        
%         
%          plot(traj(linkid).t,traj(linkid).xf,':');
% % 
%          drawnow
%          pause
%          hold off

        removelist = cat(1,removelist,i);
    end

    if ~mod(i,100), fprintf(1,'.'); end; % fprintf(1,'\n'); end

end
traj(removelist) = [];
