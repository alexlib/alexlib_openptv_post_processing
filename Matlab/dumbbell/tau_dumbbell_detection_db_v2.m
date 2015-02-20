function tau_dumbbell_detection_db_v2(directory,n_img,first,last,min_pixel_radius)
% tau_dumbbell_detection_db_v2(DIRECTORY,N_IMG,FIRST,LAST,MIN_PIXEL_RADIUS)
% this routine is written to binary filter round objects from an image
% the objects should have a radius that is always slightly larger than
% what is defined in 'min_pixel_radius'
% In TAU case we modified the original detection_proc_matlab_db.m to
% tune it for the 'poor' images taken with our dumbbell.
% in principle we do not use imextendedmax, but rather the simple threshold
% and the size filter.
%
% Example:
% directory = 'H:\PTV';
% first = 338;
% last = 452;
% min_pixel_radius = 35;
% for cameraNum = 1:4
% tau_dumbbell_detection_db(directory,cameraNum,first,last,min_pixel_radius)
% end
% History:
% - v1.0 is updated in TAU on 24.06.10 to use our dumbbell, for Tracey
% - v2.0 is updated for the data in Scene 45, 07.07.10

% Authors: Alex Liberzon and modified by Beat Lüthi
% Copyright (c) 2010, Alex Liberzon
% Copyright (c) 2010, Beat Lüthi
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
% min_pixel_radius=40; % moved to inputs, 24.06.10

se = strel('disk',min_pixel_radius);
for img = 1:n_img

for filenum = first:last
     disp([n_img filenum])

    % read the file
     try % error-proofing - writes empty files
        f1 = imread([directory,'\firstview\Movie1_Scene19_',int2str(n_img),'-DVR Express CLFC_0',int2str(filenum),'.TIF']);
        
        
        % found iterative way of getting rid of the background:
        % imshow([db2_2,imsubtract(db2_2,imadjust(imsubtract(db2_1,imabsdif
        % f(db2_1,db2_2))))])
        
        % f1 = imadjust(medfilt2(f1,[15 15]));
        % background = imopen(f1,strel('disk',2*min_pixel_radius));
        
        % f1 = db2_1;
        % f1a = imsubtract(f1,imadjust(imsubtract(f1,imabsdiff(f1,db2_2))));
        f1 = imsubtract(f1,db2_2);

        % f1 = f1 - background;
        f1(1:min_pixel_radius,:) = 0;
        f1(:,1:min_pixel_radius) = 0;
        f1(end-min_pixel_radius:end,:) = 0;
        f1(:,end-min_pixel_radius:end) = 0;
        
        f1 = imadjust(f1);

        f2 = imopen(f1,se);

        %
        stats = [];

        factor = 2;
        n_iter = 0;

        while ((length(stats) < 1 || length(stats) > 2 ) && n_iter < 10)
            % 24.06.10 - change from previous version
            f3 = im2bw(f2,factor*graythresh(f2));
            
            imshow([f3]); drawnow

            % f3 = imextendedmax(f2,10);
            % imshow(f3);
            stats = regionprops(bwlabel(f3),{'Centroid','Area','MajorAxisLength','MinorAxisLength','PixelIdxList'});
            if length(stats) > 2
                factor = factor *1.1;
            else
                factor = factor *0.9;
            end
            n_iter = n_iter + 1;
        end
        
        for i = 1:length(stats)
            stats(i).sumg = sum(f1(stats(i).PixelIdxList));
        end
%     catch
%         disp('error in file')
%         disp([n_img filenum])

%         stats = [];

%     end
%     fname = [directory,'\firstview\db',int2str(n_img),'.',int2str(filenum),'_targets'];
%     write_targets(fname,stats);
% 
% end
% 
% end