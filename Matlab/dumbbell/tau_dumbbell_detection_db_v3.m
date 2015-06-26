function tau_dumbbell_detection_db_v3(directory,n_img,first,last,withbg,verbose)
% tau_dumbbell_detection_db_v2(DIRECTORY,N_IMG,FIRST,LAST)
% this routine is written to binary filter round objects from an image
% the objects should have a radius that is always slightly larger than
% what is defined in 'min_pixel_radius'
% In TAU case we modified the original detection_proc_matlab_db.m to
% tune it for the 'poor' images taken with our dumbbell.
% in principle we do not use imextendedmax, but rather the simple threshold
% and the size filter.
%
% Example:
% directory = 'e:\Resuspension\Calibration Files\';
% first = 10674;
% last = 10976;
% withbg = 0; % 1 if you have a background image to subtract, named cam[N]_bg.tif, N = 1..4 
% verbose = 0; % 1 if you want to see the image after marking the
% dumbbells, slow 
% % Initial step - create "ball" images by cropping the first image
%   for cameraNum = 1:4
%     tau_dumbbell_detection_db_v3(directory,cameraNum,first,first,withbg,verbose);
%   end
%  
% for cameraNum = 1:4
%   tau_dumbbell_detection_db_v3(directory,cameraNum,first,last,withbg,verbose)
% end

% History:
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

    im = imread([directory,filesep,'cam',int2str(n_img),'.',sprintf('%05d',first)]);
    if withbg
        im1 = (imlincomb(1,im,-.8,bg));
    else
        im1 = im; 
    end
        
    [ball1,rect1]  = imcrop(im1);
    [ball2,rect2]  = imcrop(im1);

    save(sprintf('%s%scam%d.mat',directory,filesep,n_img),'ball1','ball2');
end



for filenum = first:last

    %     try % error-proofing - writes empty files

    % read the image
    % im = imread([directory,filesep,'db2_',int2str(n_img),'.',sprintf('%02d',filenum)]);
    imname = [directory,filesep,'cam',int2str(n_img),'.',sprintf('%05d',filenum)];
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
    stats(1).Area = prod(size(ball1));


    stats(2).Centroid(1) = j2;
    stats(2).Centroid(2) = i2;
    stats(2).MajorAxisLength = size(ball2,1);
    stats(2).MinorAxisLength = size(ball2,2);
    stats(2).sumg = sum(ball2(:));
    stats(2).Area = prod(size(ball2));
    %
    %     catch
    %         disp('error in file')
    %         disp([n_img filenum])
    %
    %         stats(1).Centroid(1) = 0;
    %         stats(1).Centroid(2) = 0;
    %
    %         stats(1).MajorAxisLength = 999;
    %         stats(1).MinorAxisLength = 999;
    %         stats(1).sumg = 999;
    %         stats(1).Area = 999;
    %
    %         stats(2).Centroid(1) = 0;
    %         stats(2).Centroid(2) = 0;
    %         stats(2).MajorAxisLength = 999;
    %         stats(2).MinorAxisLength = 999;
    %         stats(2).sumg = 999;
    %         stats(2).Area = 999;
    %
    %
    %     end
    fname = [imname,'_targets'];
    write_targets(fname,stats);
    if verbose
        show_targets(fname);
    end

end
