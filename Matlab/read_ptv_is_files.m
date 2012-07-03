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


[data(1).prev, data(1).next, data(1).xr, data(1).yr, data(1).zr] = ... 
	textread(fullfile(directory,d(1).name), '%d %d %f %f %f', 'headerlines', 1);

[junk,junk,ext,junk] = fileparts(fullfile(directory,d(1).name));
data(1).t = str2double(ext(2:end));

for k = first+1:last
    i = k - first + 1;
	[data(i).prev, data(i).next, data(i).xr, data(i).yr, data(i).zr] = ... 
		textread(fullfile(directory,d(i).name), '%d %d %f %f %f', 'headerlines', 1);

    data(i).t = data(i-1).t + 1;
end