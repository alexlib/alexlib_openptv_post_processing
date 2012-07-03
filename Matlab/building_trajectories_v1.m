function xuap = building_trajectories_v1(xuap)

% if isfield(xuap,'xf')
%
%     % xuap = struct('prev','next','x','y','z'....
%
%     id = 0;
%     xuap(1).trajid = -999; % only initialization, overwritten immediately
%
%     for j = 1:length(xuap(1).prev) % all points in the FIRST xuap file
%         %    if xuap(1).prev(j) == 0 % check
%         if xuap(1).next(j) ~= -1 && xuap(1).xf(j) ~= 0
%             id = id + 1; % add new trajectory ID
%             xuap(1).trajid(j) = id; % assign to the line
%         else
%             xuap(1).trajid(j) = -999;
%         end
%     end
%
%
%     for i = 2:length(xuap) % all the following XUAP files
%         if xuap(i).t - xuap(i-1).t ~= 1 || id == 0
%             for j = 1:length(xuap(i).prev) % all points in the FIRST xuap file
%                 if xuap(i).next(j) ~= -1 && xuap(i).xf(j) ~= 0
%                     id = id + 1; % add new trajectory ID
%                     xuap(i).trajid(j) = id; % assign to the line
%                 else
%                     xuap(i).trajid(j) = -999;
%                 end
%             end
%         else
%             for j = 1:length(xuap(i).prev) % all points in the xuap files
%                 if xuap(i).prev(j) == 0 || (length(xuap(i-1).trajid) < xuap(i).prev(j)) % open new trajectory
%                     id = id + 1; % add new trajectory ID
%                     xuap(i).trajid(j) = id; % assign to the line
%                 else %if xuap(i).prev(j) <= length(xuap(i-1).trajid) % && j == xuap(i-1).next(xuap(i).prev(j)) % double check
%                     xuap(i).trajid(j) = xuap(i-1).trajid(xuap(i).prev(j));
%                 end
%             end
%         end
%     end
% else % ADDED FILES


% first file

trajid = 1; % last trajid

ind = xuap(1).next > -2; % find all relevant
xuap(1).trajid(1:length(xuap(1).next),1) = NaN;
xuap(1).trajid(ind) = [trajid:trajid+sum(ind)-1];

trajid = trajid+sum(ind);


for k = 2:length(xuap)
    
    xuap(k).trajid(1:length(xuap(k).next),1) = NaN;
    
    old = xuap(k).prev > -1;
    
    xuap(k).trajid(old) = xuap(k-1).trajid(xuap(k).prev(old)+1);
    
    new = xuap(k).prev < 0 & xuap(k).next > -2;
    
    xuap(k).trajid(new) = [trajid:trajid+sum(new)-1];
    
    trajid = trajid + sum(new);
end

rmfield(xuap,'prev');
rmfield(xuap,'next');

for i = 1:length(xuap)
    ind = isnan(xuap(i).trajid);
    xuap(i).xr(ind) = [];
    xuap(i).yr(ind) = [];
    xuap(i).zr(ind) = [];
    xuap(i).trajid(ind) = [];
    xuap(i).t = repmat(xuap(i).t,length(xuap(i).xr),1);
end



%
%
%
%
%
%
% %     traj = struct('x',[],'y',[],'z',[],'t',[],'trajid',[]);
% %     traj = repmat(traj,length(xuap)*sum(ind)*2,1); %rough estimate of the number of trajectories
%
%
%
%
% %     ind = find(ind);
%
% %     for i = 1:length(ind)
% %         traj(i).x = xuap(1).xr(ind(i));
% %         traj(i).y = xuap(1).yr(ind(i));
% %         traj(i).z = xuap(1).zr(ind(i));
% %         traj(i).t = xuap(1).t;
% %         traj(i).trajid = i;
% %     end
%
%     trajid = i; % last trajectory
%
%     % loop over files
%
%     tmp = xuap(2);
%
%     % old trajectories
%     old = find(tmp.prev > -1);
%
%     for i = 1:length(old)
%
%         traj(i).x = xuap(1).xr(ind(i));
%         traj(i).y = xuap(1).yr(ind(i));
%         traj(i).z = xuap(1).zr(ind(i));
%         traj(i).t = xuap(1).t;
%         traj(i).trajid = i;
%
%
%
%     for j = 1:length(xuap(1).prev) % all points in the FIRST xuap file
%         %    if xuap(1).prev(j) == 0 % check
%         if xuap(1).next(j) ~= -2 && xuap(1).xr(j) ~= 0
%             id = id + 1; % add new trajectory ID
%             xuap(1).trajid(j) = id; % assign to the line
%         else
%             xuap(1).trajid(j) = -999;
%         end
%     end
%
%
%     for i = 2:length(xuap) % all the following XUAP files
%         if xuap(i).t - xuap(i-1).t ~= 1 || id == 0  % IF THERE WAS A JUMP IN TIME
%             for j = 1:length(xuap(i).prev) % all points in the FIRST xuap file
%                 if xuap(i).next(j) ~= -2 && xuap(i).xr(j) ~= 0
%                     id = id + 1; % add new trajectory ID
%                     xuap(i).trajid(j) = id; % assign to the line
%                 else
%                     xuap(i).trajid(j) = -999;
%                 end
%             end
%         else
%             for j = 1:length(xuap(i).prev) % all points in the xuap files
%                 if xuap(i).prev(j) == -1 || (length(xuap(i-1).trajid) < xuap(i).prev(j)) % open new trajectory
%                     id = id + 1; % add new trajectory ID
%                     xuap(i).trajid(j) = id; % assign to the line
%                 else %if xuap(i).prev(j) <= length(xuap(i-1).trajid) % && j == xuap(i-1).next(xuap(i).prev(j)) % double check
%                     xuap(i).trajid(j) = xuap(i-1).trajid(xuap(i).prev(j)+1); % +1 is because added is numbered from 0 to N-1, not from 1 to N
%                 end
%             end
%         end
%     end
% end
