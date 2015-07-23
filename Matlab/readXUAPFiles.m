function data = readXUAPFiles(directory,first,last)
% DATA = READXUAPFILES(DIRECTORY) reads all the data from the xuap.* files
% in the directory (DIRECTORY) into a structure DATA.
% DATA has the following fields:
% DATA.prev, DATA.next, DATA.xr, DATA.yr, DATA.zr,
% DATA.xf,yf,zf,uf,vf,wf,axf,ayf,azf
%{
The output is in two kinds of files. For each ptv_is.* there will be a xuap.*
containing column by column the following (r=raw, f=filtered):
link_past, link_future, x_r,y_r,z_r,x_f,y_f,z_f,u_f,v_f,w_f,ax_f,ay_f,az_f,sucessfull
%}

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


tmp = textread(fullfile(directory,d(first).name));

% ind = (tmp(:,1)>0 & tmp(:,2)>-1);
ind = (tmp(:,6)~=0);


tmp = tmp(ind,:);


data(1).prev = tmp(:,1);
data(1).next = tmp(:,2);

data(1).xr = tmp(:,3);
data(1).yr = tmp(:,4);
data(1).zr = tmp(:,5);
%
data(1).xf = tmp(:,6);
data(1).yf = tmp(:,7);
data(1).zf = tmp(:,8);
%
data(1).uf = tmp(:,9);
data(1).vf = tmp(:,10);
data(1).wf = tmp(:,11);
%
data(1).axf = tmp(:,12);
data(1).ayf = tmp(:,13);
data(1).azf = tmp(:,14);
[junk,junk,ext,junk] = fileparts(fullfile(directory,d(first).name));
data(1).t = str2double(ext(2:end));


data(last-first+1) = data(1);

for k = first+1:last
    % for i = 1:40000
    tmp = textread(fullfile(directory,d(k).name));
%     ind = (tmp(:,1)>0 & tmp(:,2)>-1);
    ind = (tmp(:,6)~=0);

    if ~isempty(tmp)

        tmp = tmp(ind,:);



        i = k - first + 1;
        data(i).prev = tmp(:,1);
        data(i).next = tmp(:,2);


        data(i).xr = tmp(:,3);
        data(i).yr = tmp(:,4);
        data(i).zr = tmp(:,5);
        %
        data(i).xf = tmp(:,6);
        data(i).yf = tmp(:,7);
        data(i).zf = tmp(:,8);
        %
        data(i).uf = tmp(:,9);
        data(i).vf = tmp(:,10);
        data(i).wf = tmp(:,11);
        %
        data(i).axf = tmp(:,12);
        data(i).ayf = tmp(:,13);
        data(i).azf = tmp(:,14);
        data(i).t = data(i-1).t + 1;
    end
    %
end
return