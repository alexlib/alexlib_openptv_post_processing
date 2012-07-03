function data = read_ptv_is_files(directory,first,last)
% to be added - READ_PTV_IS_FILES(DIRECTORY,FIRST,LAST)

if ~nargin,
    directory = '.';
end

% wd = cd;
% cd(directory);
d = dir(fullfile(directory,'ptv_is.*'));
if nargin < 2
    first = 1;
end
if nargin < 3
    last = length(d);
else
    last = min(last,first + length(d));
end


tmp = textread(fullfile(directory,d(1).name));
numRows = tmp(1);
tmp = tmp(2:end,1:5);

% ind = (tmp(:,1)~=-1 & tmp(:,2)~=-2);
% tmp = tmp(ind,:);

% if ~isempty(tmp)
data(1).prev = tmp(:,1);
data(1).next = tmp(:,2);

data(1).xr = tmp(:,3);
data(1).yr = tmp(:,4);
data(1).zr = tmp(:,5);

[junk,junk,ext,junk] = fileparts(fullfile(directory,d(1).name));
data(1).t = str2double(ext(2:end));
% end

last-first+1
data(last-first+1) = data(1);

for k = first+1:last
    % for i = 1:40000
    i = k - first + 1;
    tmp = textread(fullfile(directory,d(i).name));
    numRows = tmp(1);
    tmp = tmp(2:end,1:5);

    %     ind = (tmp(:,1)~=0 & tmp(:,2)~=-1);
    %     tmp = tmp(ind,:);

    
    data(i).prev = tmp(:,1);
    data(i).next = tmp(:,2);


    data(i).xr = tmp(:,3);
    data(i).yr = tmp(:,4);
    data(i).zr = tmp(:,5);

    data(i).t = data(i-1).t + 1;
    clear tmp
    %
end