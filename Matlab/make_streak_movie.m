first=31;
last=100;
n=30;

count=0;
for i=first-n:first-1
    count=count+1;
    if i<10
        name_load=['\\Ifu-gwh-disk\Hydromechanik3\Beat_Colloid\Mean_flow_streak_visualization/tb_mean00000',num2Str(i),'.tif'];
    elseif i<100
        name_load=['\\Ifu-gwh-disk\Hydromechanik3\Beat_Colloid\Mean_flow_streak_visualization/tb_mean0000',num2Str(i),'.tif'];
    elseif i<1000
        name_load=['\\Ifu-gwh-disk\Hydromechanik3\Beat_Colloid\Mean_flow_streak_visualization/tb_mean000',num2Str(i),'.tif'];
    else
        name_load=['\\Ifu-gwh-disk\Hydromechanik3\Beat_Colloid\Mean_flow_streak_visualization/tb_mean00',num2Str(i),'.tif'];
    end
    A(count,1:1024,1:1024)=imread(name_load);
end

for i=first:last
    if mod(i,10)==0
        i
    end
    A(1,:,:)=[];
    if i<10
        name_load=['\\Ifu-gwh-disk\Hydromechanik3\Beat_Colloid\Mean_flow_streak_visualization/tb_mean00000',num2Str(i),'.tif'];
    elseif i<100
        name_load=['\\Ifu-gwh-disk\Hydromechanik3\Beat_Colloid\Mean_flow_streak_visualization/tb_mean0000',num2Str(i),'.tif'];
    elseif i<1000
        name_load=['\\Ifu-gwh-disk\Hydromechanik3\Beat_Colloid\Mean_flow_streak_visualization/tb_mean000',num2Str(i),'.tif'];
    else
        name_load=['\\Ifu-gwh-disk\Hydromechanik3\Beat_Colloid\Mean_flow_streak_visualization/tb_mean00',num2Str(i),'.tif'];
    end
    A(n,1:1024,1:1024)=imread(name_load);
    B=squeeze(max(A,[],1));
    name_save=['\\Ifu-gwh-disk\Hydromechanik3\Beat_Colloid\Mean_flow_streak_visualization/streak_',num2Str(100000+i),'.tif'];
    
    imwrite(B,name_save,'tif');
end

