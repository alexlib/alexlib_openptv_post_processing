function pos_lst = rt_is_to_pos_lst(directory,basename)
% DATA = RT_IS_TO_POS_LST(DIRECTORY,FILES)
% reads rt_is.* files in the directory and converts into a single
% ASCII file that is pos_lst for the 2D Matlab based tracking
% http://www.physics.georgetown.edu/matlab/tutorial.html
% BASENAME could be 'rt_is.35*', default = 'rt_is.*'
%
%

% See also: read_rt_is_files

if ~nargin,
    directory = '/Users/alex/Documents/PTV/ptv_alex/res';
    basename = 'rt_is.35*';
elseif nargin == 1
    basename = 'rt_is.35*';
end
        filelist = dir(fullfile(directory,basename));


rtis = read_rt_is_files(directory,filelist);

for i = 1:length(rtis)
    rtis(i).t = repmat(rtis(i).t,length(rtis(i).xr),1);
end
pos_lst = cat(2,cat(1,rtis.xr),cat(1,rtis.yr),cat(1,rtis.zr),cat(1,rtis.t));


fid = fopen('pos_list.txt','w');

fprintf(fid,'%6.4f %6.4f %6.4f %d\n',  pos_lst');

fclose(fid);
