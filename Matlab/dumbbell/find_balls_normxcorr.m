function find_balls_normxcorr(imname, ball,threshold)
%FIND_BALLS_NORMXCORR (IMAGE_NAME, BALL, THRESHOLD)
%
% Example: 
% imname = 'E:\resuspension\img\cam1.320269';
% find_balls_normxcorr(imname);
% show_targets(imname);


% ver 1.0
% Date: July 11, 2010
% Copyright (c) 2010, Alex dot Liberzon at gmail dot com



verbose = 0; % 0 or 1


% load the image:
im = imread(imname);
% if ~isempty(findstr(imname,'cam2')), 
im = medfilt2(im,[1 10]); % remove white strips in camera 2
% end
% remove background
% im = imsubtract(im,medfilt2(im,[31 31]));


if nargin  < 2
    threshold = 0.1;
    disp('Select one of the particles and double click the selection ...');
    figure, imshow(im,[]); 
    ball = imcrop;
    [fpath,fname,fext] = fileparts(imname);
    save(fullfile(fpath,fname),'ball');
end

    

% image registration using normalized cross-correlation
% C = conv2(single(db2_1a),single(ball1(end:-1:1,end:-1:1)),'same');
c = normxcorr2(ball,im);

shiftx = ceil(size(ball,1)/2);
shifty = ceil(size(ball,2)/2);
c = c(shiftx:end-shiftx,shifty:end-shifty);

% threshold
c( c < threshold) = 0;


% binarization
bw  = im2bw(c);
% bw = logical(c);
if verbose
    figure, imagesc(c); colorbar;
    figure, imshow(bw);
end



% labeling:
L = bwlabel(bw);

% feature properties:
stats = regionprops(L,{'Centroid','Area','MajorAxisLength','MinorAxisLength','PixelIdxList'});
for i = 1:length(stats)
    stats(i).sumg = sum(im(stats(i).PixelIdxList));
end

% fname = [directory,'\firstview\db',int2str(n_img),'.',int2str(filenum),'_targets'];
% fname = tempname('.');
fname = [imname,'_targets'];

fprintf(1,'Writing %s \n',fname); 
write_targets(fname,stats);

if verbose
    figure,imshow(im); hold on;
    for i = 1:length(stats)
        scatter(stats(i).Centroid(1),stats(i).Centroid(2),'r+');
        % scatter(stats(i).Centroid(2),50,'go','fill');
    end
    hold off
end

% post-analysis show:
% data = textread('E:\resuspension\Matlab\tpd4e99c19_d90b_49c8_ac48_54fd3e65006c_targets')
% figure, imshow(im)
% scatter(data(2:end,2),data(2:end,3),'ro')
 