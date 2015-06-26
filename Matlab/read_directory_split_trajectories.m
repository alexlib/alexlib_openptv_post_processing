function traj = read_directory_split_trajectories(filepath,minlength)

if ~nargin
    filepath = '/Volumes/MY_PASSPORT/Folding/datasets/april26_trajAcc_Beat/';
    minlength = 10;
end

dr = dir(fullfile(filepath,'trajAcc.*'));
id = 1;
N = 0;

traj = struct('x',[],'y',[],'z',[],'u',[],'v',[],'w',[],'t',[],'id',[]);

for i = 1:length(dr)
    if ~mod(i,200), disp(i); end
    fid = fopen(fullfile(filepath,dr(i).name),'rb');
    data = textscan(fid,'%f');
    fclose(fid);
    try 
        d = reshape(data{1},34,[]).';
    catch
        disp(dr(i).name);
        save traj traj
        continue
    end
    
    x = d(:,1);
    y = d(:,2);
    z = d(:,3);
    u = d(:,4);
    v = d(:,5);
    w = d(:,6);
    age = d(:,34);
    
    [s,e,l] = splitTraj(age,minlength);
    
    N = N + length(s);
    
    for j = 1:length(s)
        id = id + 1;
        traj(id).x = x(s(j):e(j));
        traj(id).y = y(s(j):e(j));
        traj(id).z = z(s(j):e(j));
        traj(id).u = u(s(j):e(j));
        traj(id).v = v(s(j):e(j));
        traj(id).w = w(s(j):e(j));
        traj(id).t = age(s(j):e(j))+i;
        traj(id).id = id;
    end
end






end

% ---------------------------------------------------------------------------------------------------------
function [startTraj, endTraj, lenTraj] = splitTraj(age,minlength)
% Split the file into single trajectories, according to the
% AGE vector, column 29 of the trajPoint.* files
startTraj = find(age == 0);      % start points of the trajectories
endTraj = startTraj(2:end)-1;    % end points of the trajectories
endTraj = [endTraj;length(age)];
lenTraj = endTraj - startTraj + 1;

if nargin < 2, minlength = 10; end

longTraj = find(lenTraj >= minlength);
if isempty(longTraj), return, end; % if there are no long trajectories, move to the next file


%    numLongTraj = length(longTraj);     % number of long trajectories
startTraj = startTraj(longTraj);    % starting
endTraj = endTraj(longTraj);        % ending points
lenTraj = lenTraj(longTraj);

end % of function