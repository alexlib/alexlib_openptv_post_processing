function data = read_rt_is_files(directory,filelist)
% to be added - READ_RT_IS_FILES(DIRECTORY,FILELIST)
% Example:
% directory = '/Users/alex/Documents/PTV/ptv_alex/res';
% rt_is_list = dir('/Users/alex/Documents/PTV/ptv_alex/res/rt_is.36*')
% data = read_rt_is_files(directory,rt_is_list);

% assuming filelist is the structure obtained by dir function:
tmp = textread(fullfile(directory,filelist(1).name));
numRows = tmp(1);

data(1).xr = tmp(2:end,2);
data(1).yr = tmp(2:end,3);
data(1).zr = tmp(2:end,4);

extpos = find(filelist(1).name == '.',1,'last');
extension = filelist(1).name(extpos+1:end);
data(1).t = str2double(extension);
% end


data(length(filelist)) = data(1);

for k = 2:length(filelist)
    % for i = 1:40000
    tmp = textread(fullfile(directory,filelist(k).name));
    numRows = tmp(1);
    
    data(k).xr = tmp(2:end,2);
    data(k).yr = tmp(2:end,3);
    data(k).zr = tmp(2:end,4);
    
    data(k).t = data(k-1).t + 1;
    clear tmp
    %
end