% function readPTV2HDF(filename)
% reads in the HDF5 file, formatted as
% it's explained on a ptvwiki.netcipia.net
% frame[i].x[particleJ]
% each Particles has x,y,z,u,v,w,ID,TRAJID,prev,next
%
% August 12, 2007
%
% Programmer: Alex Liberzon
%
%


% Get info:
hinfo = hdf5info('tmp.h5')
disp(hinfo.GroupHierarchy)
disp(hinfo.GroupHierarchy.Datasets)
disp(hinfo.GroupHierarchy.Datasets(1).Datatype.Elements)

% Prepare the structure
numDataSets = length(hinfo.GroupHierarchy.Datasets);
numElements = length(hinfo.GroupHierarchy.Datasets(1).Datatype.Elements);
tmp = hinfo.GroupHierarchy.Datasets(1).Datatype.Elements;

fieldNames = {};
data = struct();
for i = 1:numElements
    fieldNames{i} = tmp{i,1};
    data.(fieldNames{i}) = [];
end


% read the data:
tmp = hdf5read(hinfo.GroupHierarchy.Datasets(1));

% two hard-coded values (later will appear in the attributes, etc.)
numFrames = 100;
numParticles = 100;

data = repmat(data,numFrames,numParticles);

for i = 1:numFrames
    for k = 1:numElements
        for j = 1:numParticles        
        data(i,j).(fieldNames{k}) = tmp(i).Data{k}.Data(j);
        end
    end
end

% end of a reading part

% example of showing trajectories
trajs = [data(:).trajid];
x = [data(:).x];
y = [data(:).y];
z = [data(:).z];
tr = unique(trajs(trajs > 0 & trajs ~= 999));
minLength = 5;
figure, hold on, view(3)
for i = 1:length(tr)
    ind = find(trajs == tr(i));
    if length(ind) > minLength
    plot3(x(ind),y(ind),z(ind),'DisplayName',int2str(tr(i)));
   drawnow
    end
end


