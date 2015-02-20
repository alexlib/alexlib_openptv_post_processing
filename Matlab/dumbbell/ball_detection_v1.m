function ball_detection_v1(directory,n_cam,first,last,varargin)
% BALL_DETECTION_V1(DIRECTORY,N_IMG,FIRST,LAST,VARARGIN)
% similar to TAU_DUMBBELL_DETECTION, using normxcorr2
%

% Usage:

%
% Example:
% % all the parameters given manually:
% % otherwise, skipping either of those will take a default
% % first = 1; last = last image
% % ball = load('ball.mat')
% % threshold = 0.4;
%
% directory = 'E:\resuspension\img';
% first = 320209;
% last =  320217;
% ball = load('small_ball.mat');
% thresh = 0.1;
%
% for cameraNum = 1:4
%   ball_detection_v1(directory,cameraNum,first,last,ball,thresh)
% end
%
%
% % If first/last are not given, all the files in the directory for the same cam%d number will be processed.
% Example:
%

% History:
% - v1.0 is written in TAU on 11.07.10 to deal with the relfective balls
% which are not circular objects.

% Authors: Alex Liberzon
% Copyright (c) 2010, Alex Liberzon
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


if nargin < 1
    error('Usage: ball_detection_v1(directory)')
end

d = dir_no_targets(directory,1);

if nargin < 4 || isempty(first)
    [pathstr,name,ext] = fileparts(fullfile(directory,d(1).name));
    first = str2num(ext(2:end));
    [pathstr,name,ext] = fileparts(fullfile(directory,d(end).name));
    last = str2num(ext(2:end));
end

if nargin < 2
    n_cam = 1:4;
end


for n_img = n_cam
    
    d = dir_no_targets(directory,n_img);
    
    for i = first:last
        imname = ([directory,filesep,'cam',int2str(n_img),'.',sprintf('%d',i)]);
        find_balls_normxcorr(imname,varargin{:});
    end
end

function d = dir_no_targets(directory,n_img)

d = dir([directory,filesep,sprintf('cam%d.*',n_img)]);

k = 0;
ind = [];

for i = 1:length(d)
    if ~isempty(findstr(d(i).name,'targets')) || ~isempty(findstr(d(i).name,'mat'))
        k = k + 1;
        ind(k) = i;
    end
end
if k > 0
    d(ind) = [];
end
return

