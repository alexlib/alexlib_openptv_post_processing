function xuap = graphical_cleaning_xuap(xuap,plane)
% GRAPHICAL_CLEANING_XUAP allows for the graphical input of a polygon.
% the points found in the polygon will be removed from the trajectory
% database. Works in 2D only, provide 'xy' or 'xz' or 'yz' input. Default
% is 'xy';




if nargin < 2
    error('wrong number of inputs');
end

save tmp xuap


flds = fieldnames(xuap);
flds(strcmp(flds,'t')) = [];

switch plane
    case{'xy','yx'}
        x = 'xf';
        y = 'yf';
    case{'xz','zx'}
        x = 'xf';
        y = 'zf';
    case{'yz','zy'}
        x = 'yf';
        y = 'zf';
    otherwise
        error('wrong second input: xy,yz,xz')
end


figure, hold on
plot(cat(1,xuap.(x)),cat(1,xuap.(y)),'.');
[xv,yv] = ginput;

id = [];
frameid = [];
for i = 1:length(xuap)
    for j = 1:length(xuap(i).xf)
        if inpolygon(xuap(i).(x)(j),xuap(i).(y)(j),xv,yv)
            id = cat(1,id,j);
            plot(xuap(i).(x)(j),xuap(i).(y)(j),'ro');
        end
    end
    if ~isempty(id)
        for k = 1:length(flds)
%             if length(xuap(i).(flds{k})) > 1
                xuap(i).(flds{k})(id) = [];
%             end
        end
        id = [];
    end
    if isempty(xuap(i).xf)
        frameid = cat(1,frameid,i);
    end
end
xuap(frameid) = [];


