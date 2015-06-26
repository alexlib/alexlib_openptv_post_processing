function data = xuap2pos_lst(directory,first,last)
% reads xuap files in the directory and converts into a single
% ASCII file that is pos_lst for the 2D Matlab based tracking
% http://www.physics.georgetown.edu/matlab/tutorial.html
%
%

% See also: readxuapfiles

if ~nargin,
    directory = '.';
end

% wd = cd;
% cd(directory);
d = dir(fullfile(directory,'xuap.*'));
if nargin < 2
    first = 1;
end
if nargin < 3
    last = length(d);
else
    last = min(last,length(d));
end

fid = fopen('pos_list.txt','w');

%% first file
tmp = textread(fullfile(directory,d(first).name));
ind = (tmp(:,1)~=0 & tmp(:,2)~=-1);
tmp = tmp(ind,:);

[junk,junk,ext,junk] = fileparts(fullfile(directory,d(first).name));
t = str2double(ext(2:end));
fprintf(fid,'%6.4f %6.4f %6.4f %d\n',  [tmp(:,3), tmp(:,4), tmp(:,5), t*ones(size(tmp(:,1)))]');

%% second to last
for k = first+1:last

    tmp = textread(fullfile(directory,d(k).name));
    ind = (tmp(:,1)~=0 & tmp(:,2)~=-1);
tmp = tmp(ind,:);

    t = t + 1;
    fprintf(fid,'%6.4f %6.4f %6.4f %d\n',  [tmp(:,3), tmp(:,4), tmp(:,5), t*ones(size(tmp(:,1)))]');
end
fclose(fid);
