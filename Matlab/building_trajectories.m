function xuap = building_trajectories(xuap)

if isfield(xuap,'xf')

    % xuap = struct('prev','next','x','y','z'....

    id = 0;
    xuap(1).trajid = -999; % only initialization, overwritten immediately

    for j = 1:length(xuap(1).prev) % all points in the FIRST xuap file
        %    if xuap(1).prev(j) == 0 % check
        if xuap(1).next(j) ~= -1 && xuap(1).xf(j) ~= 0
            id = id + 1; % add new trajectory ID
            xuap(1).trajid(j) = id; % assign to the line
        else
            xuap(1).trajid(j) = -999;
        end
    end


    for i = 2:length(xuap) % all the following XUAP files
        if xuap(i).t - xuap(i-1).t ~= 1 || id == 0
            for j = 1:length(xuap(i).prev) % all points in the FIRST xuap file
                if xuap(i).next(j) ~= -1 && xuap(i).xf(j) ~= 0
                    id = id + 1; % add new trajectory ID
                    xuap(i).trajid(j) = id; % assign to the line
                else
                    xuap(i).trajid(j) = -999;
                end
            end
        else
            for j = 1:length(xuap(i).prev) % all points in the xuap files
                if xuap(i).prev(j) == 0 || (length(xuap(i-1).trajid) < xuap(i).prev(j)) % open new trajectory
                    id = id + 1; % add new trajectory ID
                    xuap(i).trajid(j) = id; % assign to the line
                else %if xuap(i).prev(j) <= length(xuap(i-1).trajid) % && j == xuap(i-1).next(xuap(i).prev(j)) % double check
                    xuap(i).trajid(j) = xuap(i-1).trajid(xuap(i).prev(j));
                end
            end
        end
    end
else % ADDED FILES 

    id = 0;
    xuap(1).trajid = -999;

    for j = 1:length(xuap(1).prev) % all points in the FIRST xuap file
        %    if xuap(1).prev(j) == 0 % check
        if xuap(1).next(j) ~= -2 && xuap(1).xr(j) ~= 0
            id = id + 1; % add new trajectory ID
            xuap(1).trajid(j) = id; % assign to the line
        else
            xuap(1).trajid(j) = -999;
        end
    end


    for i = 2:length(xuap) % all the following XUAP files
        if xuap(i).t - xuap(i-1).t ~= 1 || id == 0  % IF THERE WAS A JUMP IN TIME
            for j = 1:length(xuap(i).prev) % all points in the FIRST xuap file
                if xuap(i).next(j) ~= -2 && xuap(i).xr(j) ~= 0
                    id = id + 1; % add new trajectory ID
                    xuap(i).trajid(j) = id; % assign to the line
                else
                    xuap(i).trajid(j) = -999;
                end
            end
        else
            for j = 1:length(xuap(i).prev) % all points in the xuap files
                if xuap(i).prev(j) == -1 || (length(xuap(i-1).trajid) < xuap(i).prev(j)) % open new trajectory
                    id = id + 1; % add new trajectory ID
                    xuap(i).trajid(j) = id; % assign to the line
                else %if xuap(i).prev(j) <= length(xuap(i-1).trajid) % && j == xuap(i-1).next(xuap(i).prev(j)) % double check
                    xuap(i).trajid(j) = xuap(i-1).trajid(xuap(i).prev(j)+1); % +1 is because added is numbered from 0 to N-1, not from 1 to N
                end
            end
        end
    end
end
