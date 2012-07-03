cam{1} = 'E:\resuspension\photo_video\Movie1_Scene22_1-DVR Express CLFC.avi';
cam{2} = 'E:\resuspension\photo_video\Movie1_Scene22_2-DVR Express CLFC.avi';
cam{3} = 'E:\resuspension\photo_video\Movie1_Scene22_3-DVR Express CLFC.avi';
cam{4} = 'E:\resuspension\photo_video\Movie1_Scene22_4-DVR Express CLFC.avi';

finfo = aviinfo(cam{1});

hf = figure;
% set(hf,'NextPlot','replacechildren');
% for i = 1:4
%     ha{i} = subplot(2,2,i);
% end
ha = tight_subplot(2,2);
% for i = 1:4
%     axes(ha(i));
%     xlabel(sprintf('Cam %d',i));
% end


for n = 1:finfo.NumFrames
    for i = 1:4
        tmp = aviread(cam{i},n);
        subplot(ha(i));
        subimage(imadjust(tmp.cdata,[]));
        axis off
    end
    suptitle(sprintf('%05d%',n));
    drawnow
end
        