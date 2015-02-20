function tau_dumbbell_detection_db_y1(directory, filename_pattern, n_cam, ...
    first, last, withbg, verbose)
% Based on tau_dumbbell_detection_db_v3b. Modified by Yosef Meller.
% This isn't under VCS, so let's try for some sane metadata in the comments.
% Changelog: 
%    21/6/2012. More flexible file selection, improved docs.
%
% tau_dumbbell_detection_db_y1(filename_pattern, n_img, first, last, withbg, verbose)
%
% this routine is written to binary filter round objects from an image
% the objects should have a radius that is always slightly larger than
% what is defined in 'min_pixel_radius'
% In TAU case we modified the original detection_proc_matlab_db.m to
% tune it for the 'poor' images taken with our dumbbell.
% in principle we do not use imextendedmax, but rather the simple threshold
% and the size filter.
%
% Arguments:
% directory - path to the dumbell files.
% filename_pattern - A sprintf pattern of image files. Should contain two 
%     %d placeholders, for the camera number and image number, in this
%     order. Assumes one and only one dot in the name, because MATLAB can't
%     process strings like a grown-up.
% n_cam - camera numbers to substitute into filename_pattern, a vector!
% first, last - numbers of respective frames in the image directory.
% withbg - subtract a background image from each image. The background name
%     is cam%d_bg.tif. The %d is replaced using n_cam.
%     Optional. Not well-tested!
% verbose - print information about what was found. Optional, default yes.
%
% Example:
%     directory = 'C:\PTV\PyPTV\test\img_16_db\'
%     filename_pattern = 'cam%d.%d';
%     first = 160145;
%     last = 160260;
%     withbg = 0; 
%     verbose = 1; %
%     tau_dumbbell_detection_db_y1(directory, filename_pattern, 1:4, first,
%        last, withbg, verbose);
%
% History before the y versions:
% - v1.0 is updated in TAU on 24.06.10 to use our dumbbell, for Tracey
% - v2.0 is updated for the data in Scene 45, 07.07.10
% - v3.0 is updated for the data in Scene 45 and 23, using normxcorr, 07.07.10
% - v3.0 is updated for the data of Dikla, 07.11.10

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

if nargin < 7, verbose = 0; end
if nargin < 6, withbg = false; end;

if withbg
% load the background image
    for cam = n_cam
        bg(cam) = imread(sprintf('%s%scam%d_bg.tif', directory, filesep, n_cam));
    end
end

for cam = n_cam
    % Prepares a cropped image around the balls at their starting position. If
    % a file with the results exists, loads it instead.
    cache_file = sprintf('%s%scam%d.mat', directory, filesep, cam);
    try
        load(cache_file);
    catch
        im = imread(sprintf(['%s%s', filename_pattern], directory, ...
            filesep, cam, first));
        
        if withbg
            im1 = (imlincomb(1, im, -.8, bg(cam)));
        else
            im1 = im; 
        end
        
        [ball1,rect1]  = imcrop(im1);
        [ball2,rect2]  = imcrop(im1);

        save(cache_file, 'ball1', 'ball2');
    end
    
    % Now iterate over frames, and find the targets using the cropped part.
    for filenum = first:last
        % read the image
        imname = sprintf(['%s%s', filename_pattern], directory, filesep,...
            cam, filenum);
        im = imread(imname);
        
        % Possibly subtract background:
        if withbg
            im1 = (imlincomb(1, im, -1, bg(cam)));
        else
            im1 = im; 
        end

        % Data extraction part starts here:
        [i1,j1,i2,j2] = find_targets_normxcorr(im, ball1, ball2);

        stats(1).Centroid(1) = j1;
        stats(1).Centroid(2) = i1;

        stats(1).MajorAxisLength = size(ball1,1);
        stats(1).MinorAxisLength = size(ball1,2);
        stats(1).sumg = sum(ball1(:));
        stats(1).Area = prod(size(ball1));

        stats(2).Centroid(1) = j2;
        stats(2).Centroid(2) = i2;
        stats(2).MajorAxisLength = size(ball2,1);
        stats(2).MinorAxisLength = size(ball2,2);
        stats(2).sumg = sum(ball2(:));
        stats(2).Area = prod(size(ball2));

        % Save results:
        fname = strrep(imname, '.', '_targets.');
        write_targets(fname,stats);
        if verbose
            show_targets(fname);
        end
    end
end
