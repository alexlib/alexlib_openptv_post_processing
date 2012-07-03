function [linkid] = haitao_linking_criteria(traj,xdist)
% HAITAO_LINKNG_CRITERIA builds containers of trajectories that can be
% linked together. Simple criteria of distance between the end of the
% previous trajectory and the beginning of the next trajectory is used.
% no check on velocity, acceleration, etc.

% dist = inline('(x-y)^2','x','y');

% first container is the first trajectory
k = 1;
% linkid{k} = 1;

%threshold of the distance
if nargin < 2
    xdist = 3e-3; %default
end
% figure,
% hold on

%%
for i = 1:length(traj)-1
        linkid{k} = i;
        len = length(traj(i).xf);
        
    for j = i+1:min(length(traj),i+50)
        dt = traj(j).t(1) - traj(i).t(end);
        len1 = length(traj(j).xf);
        % dt = dt/50;
        %% replace linear prediction with cubic spline prediction using all
        %% or last 3-5 points. maybe not the last points which are usually
        %% wrong.
%         predicted_x = traj(i).xf(end) + dt*traj(i).uf(end);
%         predicted_y = traj(i).yf(end) + dt*traj(i).vf(end);
%         predicted_z = traj(i).zf(end) + dt*traj(i).wf(end);
%         predicted_x = traj(i).xf(end) + dt*(traj(i).xf(end)-traj(i).xf(end-1));
%         predicted_y = traj(i).yf(end) + dt*(traj(i).yf(end)-traj(i).yf(end-1));
%         predicted_z = traj(i).zf(end) + dt*(traj(i).zf(end)-traj(i).zf(end-1));
        
predicted_x = interp1(traj(i).t,traj(i).xf,traj(i).t(end)+dt,'cubic','extrap');
predicted_y = interp1(traj(i).t,traj(i).yf,traj(i).t(end)+dt,'cubic','extrap');
predicted_z = interp1(traj(i).t,traj(i).zf,traj(i).t(end)+dt,'cubic','extrap');

predicted_x1 = interp1(traj(j).t,traj(j).xf,traj(j).t(1)-dt,'cubic','extrap');
predicted_y1 = interp1(traj(j).t,traj(j).yf,traj(j).t(1)-dt,'cubic','extrap');
predicted_z1 = interp1(traj(j).t,traj(j).zf,traj(j).t(1)-dt,'cubic','extrap');



%         plot([1:length(traj(i).xf)+1,length(traj(i).xf)+1],[traj(i).xf;predicted_x;traj(j).xf(1)],'o')
%         plot([1:length(traj(i).xf)+1,length(traj(i).xf)+1],[traj(i).yf;predicted_y;traj(j).yf(1)],'o')
%         plot([1:length(traj(i).xf)+1,length(traj(i).xf)+1],[traj(i).zf;predicted_z;traj(j).zf(1)],'o')
%       
        
        %     predicted_x = traj(i-1).xf(end);
        %     predicted_y = traj(i-1).yf(end);
        %     predicted_z = traj(i-1).zf(end);

%         disp(i)
%          disp([predicted_x,predicted_y,predicted_z]);
%          disp([traj(j).xf(1),traj(j).yf(1),traj(j).zf(1)]);

%           plot(k,norm([predicted_x,predicted_y,predicted_z]-[traj(j).xf(1),traj(j).yf(1),traj(j).zf(1)]),'.')

        if (norm([predicted_x,predicted_y,predicted_z]- [traj(j).xf(1),traj(j).yf(1),traj(j).zf(1)]) + ...
            norm([predicted_x1,predicted_y1,predicted_z1]- [traj(i).xf(end),traj(i).yf(end),traj(i).zf(end)])) < xdist
             linkid{k} = cat(1,linkid{k},j); % continue filling this container
%              plot(traj(i).t,traj(i).xf,'-o',traj(i).t(end)+dt,predicted_x,'s',traj(j).t(1)-dt,predicted_x1,'^',traj(j).t,traj(j).xf,'--^')
%              title(norm([predicted_x,predicted_y,predicted_z]- [traj(j).xf(1),traj(j).yf(1),traj(j).zf(1)]))
%              figure, subplot(211), plot(traj(i).t,traj(i).xf,'-o'), hold on, plot(traj(j).t,traj(j).xf,'-s'),
%              plot(traj(i).t+dt,predicted_x,'^')
%              subplot(212), plot(traj(i).t,traj(i).uf,'-o'), hold on, plot(traj(j).t,traj(j).uf,'-s')
%              plot(traj(i).t, gradient5(traj(i).xf,traj(i).t),'--x');
%          else
%              k = k + 1; % next container
%              linkid{k} = j;
        end
         
    end
    k = k + 1;
end

for i = 1:length(linkid)-1
    for j = i+1:length(linkid)
        if any(ismember(linkid{i},linkid{j}))
            linkid{i} = union(linkid{i},linkid{j});
            linkid{j} = [];
        end
    end
end

