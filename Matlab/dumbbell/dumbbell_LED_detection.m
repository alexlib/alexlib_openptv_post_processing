function dumbbell_LED_detection(directory,n_img,first,last,dumbbels)
% dumbbell_LED_detection is based on
% BALL_DETECTION_V1(DIRECTORY,N_IMG,FIRST,LAST)
% similar to TAU_DUMBBELL_DETECTION, using normxcorr2
% but it is for TWO DIFFERENT LED images
%
% Example:
%{
% directory = 'C:\PTV\PyPTV\Mark_Ron_Alex_18_02_15_spherical_grid\Dumbell';
directory = '../Dumbell';
first = 10000;
last =  10005;
load led1
load led2
dumbbels{1} = led1; dumbbels{2} = led2;
% try one image:
last = first;

dumbbell_LED_detection(directory,1,first,first,dumbbels)

% or a set
for cameraNum = 1:4
  dumbbell_LED_detection(directory,cameraNum,first,last,dumbbels,0)
end
%}

% % If first/last are not given, all the files in the directory for the same cam%d number will be processed.
% Example:
%

% History:
% - v1.0 is written in TAU on 11.07.10 to deal with the relfective balls
% which are not circular objects.

% Authors: Alex Liberzon
% Copyright (c) 2010-2014, Alex Liberzon
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

% Prepare led1.mat and led2.mat manually by imcrop(image)
% store it in the same directory

% if nargin == 2
%     d = dir([directory,filesep,sprintf('cam%d.*',n_img)]);
%     for i = 1:length(d)
%         imname = fullfile(directory,d(i).name);
%         find_balls_normxcorr_LED(imname,[],[],[],dumbbels);
%     end
% elseif nargin == 1
% for n_img = 1:4
%     d = dir([directory,filesep,sprintf('cam%d.*',n_img)]);
%     for i = 1:length(d)
%         imname = fullfile(directory,d(i).name);
%         find_balls_normxcorr_LED(imname);
%     end
% end
% elseif nargin == 4
for n_img 
    for filenum = first:last
        % read the image
        imname = ([directory,filesep,'cam',int2str(n_img),'.',sprintf('%05d',filenum)]);
        find_balls_normxcorr_LED(imname,dumbbels,0.4);
    end
end
% end
