function background_motion_estimate(directory)
if nargin < 1
    directory = '/Users/alex/Documents/PTV/ptv_alex/img_35_try';
    % directory = '/Users/alex/Documents/PTV/ptv_alex/img34_10';
end
for ncam = 1:4
    sub_background_motion_estimate(directory,ncam,20)
end

function sub_background_motion_estimate(varargin)
% BACKGROUND_MOTION_ESTIMATION(DIRECTORY,N_CAMERA,N_IMAGES_TO_AVERAGE)
%
% Default:
%   d = '/Users/alex/Documents/PTV/ptv_alex/img_35_try';
%   ncam = 3;
%	nimg = 25;
%

d = varargin{1};
ncam = varargin{2};
nimg = varargin{3};


% change to 1 if you want to see the running image
verbose = 1;

filenames = getfiles(d,ncam);

% let's take nimg images average and std

tmp = repmat(readim(filenames{1}),[1 1 nimg]);
for i = 2:nimg
    tmp(:,:,i) = readim(filenames{i});
end
imavg = uint8(mean(single(tmp),3));
% imstd = uint8(std(single(tmp),1,3));


% imnew = repmat(imavg,[1,1,nimg]);

if verbose, figure, end

for i = 1:length(filenames)
    im = readim(filenames{i});
    tmp(:,:,1:end-1) = tmp(:,:,2:end);
    tmp(:,:,end) = im;
    imavg = uint8(mean(single(tmp),3));
    % imavg = imlincomb((nimg-1)/nimg, imavg, 1/nimg, im);
%     imstd = imlincomb((nimg-1)/nimg, imstd, 1/nimg, im.^2);
%     if verbose,
%         imshow(im);
%         hold on
%     end
    % im(imabsdiff(imavg,im) < 1*sqrt(single(imstd))) = 0;
    % bw = im2bw(im);
    % bw = bwmorph(im2bw(imfilter(imadjust(imlincomb(1,imabsdiff(im,imavg),-.5,imstd)),fspecial('gaussian'))),'clean');
    im1 = imopen(imsubtract(im,imavg),ones(2,2)); 
    bw = im2bw(im1,.1);
    if verbose,
        imshow(im,[]);
        hold on
    end
    % im1 = adapthisteq(imsubtract(im,imavg));

    % bw = im2bw(im1,.2);
    % imshow(bw);
    stats  = regionprops(logical(bw), {'Centroid','Area','MajorAxisLength','MinorAxisLength','PixelIdxList'});
%     area = cat(1,stats.Area);
%     valid = area > 10 & area < 300;
%     stats = stats(valid);    
    
    for j = 1:length(stats)
        stats(j).sumg = sum(im(stats(j).PixelIdxList));
    end
    
    
    fprintf(1,'Writing %s \n',[filenames{i},'_targets']);
    write_targets([filenames{i},'_targets'],stats);
    
    
    if verbose
        centroids = cat(1, stats.Centroid);
        plot(centroids(:,1), centroids(:,2), 'gs');
        title(sprintf('%d',i));
        drawnow
        hold off
    end
end
return

function filenames = getfiles(d,ncam)
% GETFILES constructs a cell array of file names

% read the files
imlist = dir(sprintf('%s/cam%d.*',d,ncam));
filenames = cell(1,length(imlist));
k = 1;
for i = 1:length(imlist)
    tmp = fullfile(d,imlist(i).name);
    if isempty(findstr(tmp,'_targets')) && isempty(findstr(tmp,'.mat'))
        filenames{k} = tmp;
        k = k + 1;
    end
end
filenames = filenames(1:k-1);
return


function im = readim(filename)
% reading and filtering the image
% if ~isempty(findstr(filename,'cam3.')) || ~isempty(findstr(filename,'cam2.'))
%     im = imerode(imread(filename),strel('disk',3));
% else
    im = imread(filename);
% end
% im = uint8(imfilter(im,fspecial('gaussian',[3 3],0.5)));
im = imclose(imopen(im,ones(3,3)),ones(3,3));

return

