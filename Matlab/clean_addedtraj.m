function added = clean_addedtraj(added,xmax)
% CLEANING_added cleans the numerious zeros values in xr,yr,zr
%

if nargin < 2
    % default is no check on the limits
    xmax = Inf;
end

flds = fieldnames(added);
flds(strcmp(flds,'t')) = []; % t is not touched
% flds(strcmp(flds,'trajid')) = []; % t is not touched


frameid = zeros(length(added),1);
for i = 1:length(added)
    j = 1:length(added(i).xr);
    id = (added(i).xr(j) == 0 | added(i).yr(j) == 0 | added(i).zr(j) == 0 ...
                | abs(added(i).xr(j)) > xmax | abs(added(i).yr(j)) > xmax ...
                | abs(added(i).zr(j)) > xmax);
    if any(id)
        for k = 1:length(flds)
%              if length(added(i).(flds{k})) > 1
                added(i).(flds{k})(id) = [];
%              end
        end
%         id = [];
    end
    if isempty(added(i).xr)
        frameid(i) = 1;
    end
end
added(find(frameid)) = [];


