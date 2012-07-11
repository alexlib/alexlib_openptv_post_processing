function tif_rename(first,width,scene)
% source=['Y:/img/Movie1_Scene',num2str(scene),'_'];
% source=['E:/TRANSPORT/Movie1_Scene',num2str(scene),'_'];
source=['E:\beat708/img/Movie1_Scene',num2str(scene),'_'];
target=['D:/Krug/PTV/working_folder/img_',num2str(scene),'/cam'];
% target=['D:/Krug/PTV/working_folder/img_',num2str(scene),'/cam'];
mkdir(['D:/Krug/PTV/working_folder/img_',num2str(scene)]);

for i=first:first+width
    i
    for cam=1:4
        if i<10
            movefile([source,num2str(cam),'-DVR Express CLFC_00',num2str(i),'.TIF'],[target,num2str(cam),'.',num2str(i)])
        elseif i<100
            movefile([source,num2str(cam),'-DVR Express CLFC_0',num2str(i),'.TIF'],[target,num2str(cam),'.',num2str(i)])
        elseif i<1000
            movefile([source,num2str(cam),'-DVR Express CLFC_',num2str(i),'.TIF'],[target,num2str(cam),'.',num2str(i)])
        end
        
    end
end


end
