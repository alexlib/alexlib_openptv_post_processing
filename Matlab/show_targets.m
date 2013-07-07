function show_targets(fname,varargin)
% SHOW_TARGETS(IMAGE_NAME) shows the image and overlays the
% identified centroids of the particles in the associated image_name_targets file
%
% SHOW_TARGETS(IMAGE_NAME_TARGETS) shows only the
% identified centroids of the particles in the associated
% image_name_targets file

%
% Example:
% imname = 'E:\resuspension\img\cam1.320269';
% find_balls_normxcorr(imname);
% show_targets(imname);

%
% ver 1.0
% Date: July 11, 2010
% Copyright (c) 2010, Alex dot Liberzon at gmail dot com

%
if isempty(findstr(fname,'_targets'))
    data = textread([fname,'_targets']);
    
    if nargin > 1
        ax = varargin{1};
    else
        ax = gca;
    end
    
    imshow(fname,'Parent',ax); hold on
    if length(data) > 1
        scatter(data(2:end,2),data(2:end,3),'r+');
    end
else
    data = textread(fname);
    if nargin > 1
        ax = varargin{1};
    else
        ax = gca;
    end
    
    if length(data) > 1
        scatter(data(2:end,2),data(2:end,3),'r+');
    end
end
title(fname)
drawnow




