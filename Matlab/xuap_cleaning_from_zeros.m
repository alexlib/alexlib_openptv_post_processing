function xuap = xuap_cleaning_from_zeros(xuap,xmax)
% CLEANING_XUAP cleans the numerious zeros values in xf,yf,zf
%

if nargin < 2
    % default is no check on the limits
    xmax = Inf;
end

flds = fieldnames(xuap);
flds(strcmp(flds,'t')) = []; % t is not touched
% flds(strcmp(flds,'trajid')) = []; % t is not touched


frameid = zeros(length(xuap),1);
for i = 1:length(xuap)
    j = 1:length(xuap(i).xf);
    id = (xuap(i).xf(j) == 0 | xuap(i).yf(j) == 0 | xuap(i).zf(j) == 0 ...
                | abs(xuap(i).xf(j)) > xmax | abs(xuap(i).yf(j)) > xmax ...
                | abs(xuap(i).zf(j)) > xmax);
            
    if any(id)
        for k = 1:length(flds)
%              if length(xuap(i).(flds{k})) > 1
                xuap(i).(flds{k})(id) = [];
%              end
        end
%         id = [];
    end
    if isempty(xuap(i).xf)
        frameid(i) = 1;
    end
end
xuap(find(frameid)) = [];


