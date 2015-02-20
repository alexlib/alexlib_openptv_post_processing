


% provide the path to the images
directory = '\PTV\origo\working_folder_Dumbbell_b_10_08\img_35_try\';

% we assume there are images in the format of cam1.xxxx , cam2.xxxx
% first create the "ball" - crop manually one image to get the picture of
% the transparent ball illuminated from a side


%% this part we do once for each image directory if the balls look different in this set
%{
for n_cam = 1:4
   ball_detection_v1(directory,n_cam);
end
close all
%}

% if the mat files are in the "directory", we can use it
% Test on the first files from all cameras, show results
% test on all the 4 cameras, using only first file:

d = dir(fullfile(directory,'cam1.*'));
first = str2num(d(1).name(findstr(d(1).name,'.')+1:end));

%% this part we do once for each image directory if the balls look different in this set
% to run it using the saved 'ball.mat' on all the first files and show it
  for n_cam = 3
      ball_detection_v1(directory,n_cam,first,first);
      show_targets(fullfile(directory,sprintf('cam%d.%d',n_cam,first)));
  end

% to run it using the saved 'ball.mat' on all the first files and show it
threshold = 0.1;
  figure, hold on
  for n_cam = 3
      load(fullfile(directory,sprintf('cam%d.mat',n_cam)));
      ball_detection_v1(directory,n_cam,first,first,ball,threshold);
      ax = subplot(2,2,n_cam);
     show_targets(fullfile(directory,sprintf('cam%d.%d',n_cam,first)),ax);
  end
  
  %% Run all the files of the same camera in the directory, manually chosen
  % 
  threshold = 0.1;
  figure, hold on
  for n_cam = 3
      for img = 351160:351180
          load(fullfile(directory,sprintf('cam%d.mat',n_cam)));
          ball_detection_v1(directory,n_cam,img,img,ball,threshold);
          % ax = subplot(2,2,n_cam);
          ax = axes;
          show_targets(fullfile(directory,sprintf('cam%d.%d',n_cam,img)),ax);
          drawnow
      end
  end

%% Run all the files in the directory
% to run on all the files, no graphics
for n_cam = 1:4
   load(fullfile(directory,sprintf('cam%d.mat',n_cam)));
   ball_detection_v1(directory,n_cam,[],[],ball,threshold);
    % ax = subplot(2,2,n_cam);
    % show_targets(fullfile(directory,sprintf('cam%d.%d',n_cam,first)),ax);
end