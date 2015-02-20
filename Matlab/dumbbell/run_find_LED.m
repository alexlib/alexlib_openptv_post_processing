% run find_LED

directory = '../Dumbell';
first = 10000;
last =  10005;

% load the preloaded images or create those first:
load dumbbels

% or
% prepare dumbbels first if these do not exist:
%{
for nCam = 1:4
    % we must prepare the LED images for each camera first
    imname = fullfile(directory,sprintf('cam%d.%05d',nCam,first));
    im = imread(imname);
    for nLED = 1:2
        dumbbels{nCam,nLED} = imcrop(im);
    end
end
save dumbbels
%}

% try one image
verbose = true;
nCam = 4;
imname = fullfile(directory,sprintf('cam%d.%05d',nCam,first));
find_LED(imname,dumbbels(nCam,:),verbose)




% now run the batch:
load dumbbels
for nCam = 1:4
    for nFile = first:last
        imname = fullfile(directory,sprintf('cam%d.%05d',nCam,nFile));
        find_LED(imname,dumbbels(nCam,:),false);
    end
end