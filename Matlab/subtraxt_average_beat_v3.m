

for m=1:1
    
    direc=['E:\Exp_12112010\Exp',num2str(m),'_12112010'];
    d = [direc,'/Cam1*'];
    
    
    dd=dir(d)
    filename1=dd(1).name
    %filename2=dd(end).name
    
    
    
    first = eval(filename1(6:end));
    if first==1
        last = size(dd,1);
    else
        filename2=dd(end).name
        last=eval(filename2(6:end))
    end
    
    
    
    
    %first =100001;
    %last =101799;
    
    %av=zeros(509-118+1,314-1+1);
    av=zeros(512,512);
    
    %make average
    for k=1:2
        for i = first:last
            i=i;
            %tmp=imread(['E:\PTV\Working_folder\WF_3\image/Cam',num2str(k),'.',num2Str(i)]);
            tmp=imread([direc,'\Cam',num2str(k),'.',num2Str(i)]);
            %A=double(tmp(118:509,1:314,1));
            A=double(tmp(1:512,1:512,1));
            av=av+A;
            
            averaging=[m k i last];
            display(averaging)
        end
        
        av=av/(last-first+1);
        %     figure;
        %imshow(uint8(av))
        %     title('average')
        
        for i = first:1:last
            %i=i;
            %tmp=imread(['E:\PTV\Working_folder\WF_3\image/Cam',num2str(k),'.',num2Str(i)]);
            tmp=imread([direc,'\Cam',num2str(k),'.',num2Str(i)]);
            %A=double(tmp(118:509,1:314,1));
            A=double(tmp(1:512,1:512,1));
            B=A-av;%
            ind=find(B<0);%anything below avg px value will be eliminated
            B(ind)=0;
            %figure;
            %imshow(uint8(B))
            %     h=surf(log(B+1));
            %     set(h,'EdgeColor','none');
            %tit=['frame: ',num2Str(i)];
            %title(tit);
            %imwrite(uint8(B),['E:\PTV\Working_folder\WF_3\image/filt_Cam',num2str(k),'.',num2Str(i),'.tif'],'tiff','compression','none');
            imwrite(uint8(B),[direc,'\filt_Cam',num2str(k),'.',num2Str(i)],'tiff','compression','none');
            writing=[m k i last];
            display(writing)
        end
    end
    
    
    
    recycle off;
    for k=1:4
        for i=first:last
            delete([direc,'\Cam',num2str(k),'.',num2Str(i)]);
            deleting=[m k i last]
        end
    end
end

% figure;
% i=700;
% k=1;
% tmpp=imread(['E:\PTV\Working_folder\Exp23multib_15112010\Copy_of_Exp23multib_15112010/filt_Cam',num2str(k),'.',num2Str(i)]);
%
% imshow(tmpp)
