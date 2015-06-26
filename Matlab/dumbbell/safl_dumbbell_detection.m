function safl_dumbbell_detection(directory,n_img,first,last,withbg,verbose)
% SAFL_dumbbell_detection(DIRECTORY,N_IMG,FIRST,LAST)
%
% The user selects manually the two balls of the dumbell and then uses
% those as the template matching for the rest of the files in the directory
%
% Example:
%{
directory = 'C:\PTV\Experiments\Sep13\db2400'
first = 1;
last = 399;
withbg = 0; % 1 if you have a background image to subtract, named cam[N]_bg.tif, N = 1..4 
verbose = 1; % 1 if you want to see the images, slow 
% Initial step - create "ball" images by cropping the first image
  for cameraNum = 1:4
    safl_dumbbell_detection(directory,cameraNum,first,first,withbg,verbose);
  end
 
verbose = 0; % too slow for large sets
for cameraNum = 1:4
  safl_dumbbell_detection(directory,cameraNum,first,last,withbg,verbose)
end
%}

% History:
% - Sep 13, 2012 - based on the tau_dumbbell_detection_db_v3.m, updated for
% SAFL



% Authors: Alex Liberzon and modified by Beat L?thi
% Copyright (c) 2010, Alex Liberzon
% Copyright (c) 2010, Beat L?thi
% All rights reserved.
%

% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
%
% * Redistributions of source code must retain the above copyright
% notice, this list of conditions and the following disclaimer.
% * Redistributions in binary form must reproduce the above copyright
% notice, this list of conditions and the following disclaimer in
% the documentation and/or other materials provided with the distribution
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.

if nargin < 6
    verbose = 1; % or change to 0 if you want silent work
end

if nargin < 5, withbg = false; end;

if withbg
% load the background image
    bg = imread(sprintf('%s%scam%d_bg.tif',directory,filesep,n_img));
end

% load the 'balls' images or prepare them manually
try
    load(sprintf('%s%scam%d.mat',directory,filesep,n_img));
catch

    im = imread([directory,filesep,'cam',int2str(n_img),'.',sprintf('1%06d',first)]);
    if withbg
        im1 = (imlincomb(1,im,-.8,bg));
    else
        im1 = im;  % SAFL camera is dark
    end
    
    figure, imshow(imadjust(im1),[]);
    [~,rect1]  = imcrop;
    ball1 = imcrop(im1,rect1);
    imshow(imadjust(im1),[]);
    [~,rect2]  = imcrop;
    ball2 = imcrop(im1,rect2);

    save(sprintf('%s%scam%d.mat',directory,filesep,n_img),'ball1','ball2');
end



for filenum = first:last

    %     try % error-proofing - writes empty files

    % read the image
    % im = imread([directory,filesep,'db2_',int2str(n_img),'.',sprintf('%02d',filenum)]);
    imname = [directory,filesep,'cam',int2str(n_img),'.',sprintf('1%06d',filenum)];
    im = imread(imname);
    if withbg
        im1 = (imlincomb(1,im,-1,bg));
    else
        im1 = im; 
    end
    
%     if verbose
%         figure, imshow([im,bg,im1]);
%     end

    [i1,j1,i2,j2] = find_targets_normxcorr(im, ball1, ball2);

    stats(1).Centroid(1) = j1;
    stats(1).Centroid(2) = i1;

    stats(1).MajorAxisLength = size(ball1,1);
    stats(1).MinorAxisLength = size(ball1,2);
    stats(1).sumg = sum(ball1(:));
    stats(1).Area = numel(ball1);


    stats(2).Centroid(1) = j2;
    stats(2).Centroid(2) = i2;
    stats(2).MajorAxisLength = size(ball2,1);
    stats(2).MinorAxisLength = size(ball2,2);
    stats(2).sumg = sum(ball2(:));
    stats(2).Area = numel(ball2);
    fname = [imname,'_targets'];
    write_targets(fname,stats);
    if verbose
        show_targets(fname);
    end

end
