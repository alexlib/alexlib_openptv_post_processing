function traj = read_new_traj_files(directory,first,last)
% READ NEW_TRAJ files
if ~nargin,
    directory = 'c:\ptv\experiments\re1000'; % '.';
end

% wd = cd;
% cd(directory);
d = dir(fullfile(directory,'new_traj.*'));
if nargin < 2
    first = 1;
end
if nargin < 3
    last = length(d);
else
    last = min(last,length(d));
end

m = 0;

for k = first:last
    tmp = textread(fullfile(directory,d(k).name));

    if ~isempty(tmp)
        age = tmp(:,1);
        starts = find(age == 0);
        ends = [starts(2:end)-1;length(age)];
        for i = 1:length(starts) % number of traj
            m  = m + 1;
            traj(m).xf = tmp(starts(i):ends(i),4);
            traj(m).yf = tmp(starts(i):ends(i),5);
            traj(m).zf = tmp(starts(i):ends(i),6);
            traj(m).t =  tmp(starts(i):ends(i),2);
            traj(m).trajid = m(:,ones(length(traj(m).xf),1));
        end
    end
end
