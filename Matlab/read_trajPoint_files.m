function trajPoint_to_traj(directory,first,last)
% TRAJPOINT_TO_TRAJ is reading all the trajPoint.xxxx files
% from the FIRST to LAST or all the files in the directory
% and stores it in a TRAJ structure
%
% Example:
% 
%   traj = trajPoint_to_traj('/Volumes/ALEX/SAFL/sep12/res',1,10);
%


% Initialization
traj = repmat(struct('xf',[],'yf',[],'zf',[],...
    'uf',[],'vf',[],'wf',[],...
    'axf',[],'ayf',[],'azf',[],...
    't',[],'trajid',[]),length(unique_trajid),1);


xr,yr,zr,...
    u,v,w,...
    ax,ay,az,...
    alx,aly,alz,...
    w1,w2,w3,...
    s11,s12,s13,s22,s23,s33,...
    axx,axy,axz,ayx,ayy,ayz,azx,azy,azz,...
    reldiv, age


i = 0;
for k = 1:length(unique_trajid)
    id = find(trajid == unique_trajid(k));
    if length(id) > minLength
        i = i + 1;
        trajLen(i) = length(id);
        traj(i).xf = xf(id);
        traj(i).yf = yf(id);
        traj(i).zf = zf(id);
        traj(i).uf = uf(id);
        traj(i).vf = vf(id);
        traj(i).wf = wf(id);
        traj(i).axf = axf(id);
        traj(i).ayf = ayf(id);
        traj(i).azf = azf(id);
        traj(i).t = t(id);
        traj(i).trajid = trajid(id);
    end
end
traj(i+1:end) = [];
trajLen(i+1:end) = [];


if ~nargin, 
    directory = '/Volumes/ALEX/SAFL/sep12/res/'; % default for the test
end
fname = fullfile(directory,'trajPoint.');

if nargin < 2 % we need to check what is first, what is last
    d = dir([fname,'*']);
    first = d(1).name;
    last = d(end).name;
    [~,first] = strtok(first,'.');
    [~,last] = strtok(last,'.');
    first = str2double(first(2:end));
    last = str2double(last(2:end));
end

count = 1;
data  = zeros(2e6,32);
for i = first:last
    fprintf(1,'%d\n',i);
    tmp = load(sprintf('%s%d',fname,i));
    s = size(tmp,1);
    if s > 0
        [xr,yr,zr,...
    u,v,w,...
    ax,ay,az,...
    alx,aly,alz,...
    w1,w2,w3,...
    s11,s12,s13,s22,s23,s33,...
    axx,axy,axz,ayx,ayy,ayz,azx,azy,azz,...
    reldiv, age] = parse_trajPoint(tmp);

        data(count:count+s-1,:) = tmp;
        count = count + s;
    end
end


data = data(1:count-1,:);

[xr,yr,zr,...
    u,v,w,...
    ax,ay,az,...
    alx,aly,alz,...
    w1,w2,w3,...
    s11,s12,s13,s22,s23,s33,...
    axx,axy,axz,ayx,ayy,ayz,azx,azy,azz,...
    reldiv, age] = parse_f(data);


[startTraj,endTraj,lenTraj] = split_traj(age);



traj = data_to_traj(data);







