close all
clear all
clc
figure;




% first=19729;
% last=19729;

% first=21000;%------One breaks into two
% last=21250;

% first=19600;
% last=19750;

% first=19200;
% last=19450;
% first=19700;
% last=19750;


% first=8400;
% last=8450;

% first=19735;
% last=19870;

% 
% first=19720;%----One breaks into 4
% last=19750;
% 

% first=5902;
% last=5902;


% first=1270;
% last=1350;

first=2649;
last=2649;

dummy=0;

%winsize(1:2) = [0 0];%%-----------------------------------
%colordef white

for kamera=1:1
    for i=first:last
        dummy=dummy+1;
        %name=['E:\PTV\Working_folder\WF_1\img_for_full_processing\cam',num2str(kamera),'.',num2Str(i)];
        name=['E:\PTV\Working_folder\Exp25b_09112010\image\filt_Cam',num2str(kamera),'.',num2Str(i)];
        f2=imread(name);
        %-----addition-----
        %f2=imcrop(f2,[427 8 350 900]);
        %-------------------%
        imshow(f2);
        tit=[ 'Camera:',num2str(kamera) ' ' 'Frame number: ',num2Str(i)];
        title(tit);
        A(dummy)=getframe;
    end
end
%close all
%movie(A,1,50)
%movie(A,-1,5)% what about other input i.e., window size ???
%save movie.mat A %--------This is computationaly expensive, i.e., .mat
%                          file. Its a Matlab movie that runs only in Matlab platform
%mpgwrite(A,gray,'movie.mpg');

% we can also specify the color map as jet / gray in may case (after A and without any string)
%mpeg_play movie.mpg



%[A,rect]=imcrop(f2);----This is the way to extract the croped area










