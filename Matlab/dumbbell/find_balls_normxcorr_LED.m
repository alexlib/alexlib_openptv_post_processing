function [stats] = find_balls_normxcorr_LED(imname, dumbbels,threshold)
%FIND_BALLS_NORMXCORR_LED (IMAGE_NAME, BALLS, THRESHOLD)
%
% Example:
% directory = 'E:\resuspension\img';
% first = 320209;
% last =  320217;
% for cameraNum = 1:4
%   ball_detection_v1(directory,cameraNum,first,last)
% end
%
%
% Example:
% imname = 'E:\resuspension\img\cam1.320269';
% find_balls_normxcorr(imname);
% show_targets(imname);


% ver 1.0
% Date: July 11, 2010
% Copyright (c) 2010, Alex dot Liberzon at gmail dot com

% if nargin  < 2 % load the saved ball.mat
%     % load('ball.mat');
%     load led1
%     load led2
%     balls{1} = led1;
%     balls{2} = led2;
%     threshold = 0.4;
%     append = true; % for LEDs where we do two loops
% end

if nargin ~= 3
    error('wrong number of inputs');
end


verbose = 1; % 0 or 1


% load the image:
im = imread(imname);

if length(dumbbels) > 1
    append = true;
else
    append = false;
end


for nLed = 1:length(dumbbels)

    ball = dumbbels{nLed};
    
    % image registration using normalized cross-correlation
    % C = conv2(single(db2_1a),single(ball1(end:-1:1,end:-1:1)),'same');
    c = normxcorr2(ball,im);
    
    shiftx = ceil(size(ball,1)/2);
    shifty = ceil(size(ball,2)/2);
    c = c(shiftx:end-shiftx,shifty:end-shifty);
    
    
    
    if append % multiple LEDs, leave only one
        d = zeros(size(c));
        [xc,yc] = find(c == max(c(:)));
        d(xc-shiftx:xc+shiftx,yc-shifty:yc+shifty) = ...
            c(xc-shiftx:xc+shiftx,yc-shifty:yc+shifty);
        c = d;
    end
    
    % threshold
    c( c < threshold) = 0;
    
    
    % binarization
    bw  = im2bw(c,threshold);
    % bw = logical(c);
    if verbose
        figure, imagesc(c); colorbar;
        figure, imshow(bw);
    end
    
    
    
    % labeling:
    L = logical(bw);
    
    % feature properties:
    stats = regionprops(L,{'Centroid','Area','MajorAxisLength','MinorAxisLength','PixelIdxList'});
    for i = 1:length(stats)
        stats(i).sumg = sum(im(stats(i).PixelIdxList));
    end
    
    % fname = [directory,'\firstview\db',int2str(n_img),'.',int2str(filenum),'_targets'];
    % fname = tempname('.');
    fname = [imname,'_targets'];
    
    
    write_targets(fname,stats,append);
    
    
    
    if verbose
        figure,imshow(im); hold on;
        for i = 1:length(stats)
            scatter(stats(i).Centroid(1),stats(i).Centroid(2),'r+');
            % scatter(stats(i).Centroid(2),50,'go','fill');
        end
        hold off
    end
    
end

% post-analysis show:
% data = textread('E:\resuspension\Matlab\tpd4e99c19_d90b_49c8_ac48_54fd3e65006c_targets')
% figure, imshow(im)
% scatter(data(2:end,2),data(2:end,3),'ro')
