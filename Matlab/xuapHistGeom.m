%------xuap histogram------%
%n=19233;
first=1;
last=21405;

xmin=10/1000;
xmax=15/1000;
ymin=20/1000;

histN_1_LS=zeros(1,10000);
histN_2_LS=zeros(1,10000);
histN_3_LS=zeros(1,10000);
histN_4_LS=zeros(1,10000);

histN_1_L=zeros(1,10000);
histN_2_L=zeros(1,10000);
histN_3_L=zeros(1,10000);
histN_4_L=zeros(1,10000);

for n=first:last;

    if mod(n,2)==0
        n
    end
    f_name_xuap=['E:\PTV\Working_folder\WF_1\res\xuap.',num2Str(n)];


    xuap=load(f_name_xuap);
    sz=size(xuap);

    if sz(1,1)>0
        past=xuap(:,1);
        future=xuap(:,2);
        x=xuap(:,3);
        y=xuap(:,4);
        spline=xuap(:,15);

        link=find (future>0 & x> xmin & x< xmax & y>ymin);
        link_splined=find(future>0 & x> xmin & x< xmax & y>ymin & spline==1);

    end

    f1=['E:\PTV\Working_folder\WF_1\f_files\f_',num2Str(n),'.mat'];
    clear f
    f=load(f1);

    
    npx_1=f.a1_a;
    npx_2=f.a2_a;
    npx_3=f.a3_a;
    npx_4=f.a4_a;




   %------------ Camera 1------------------------------
    
        % for only linked for camera 1
    
        px_L_1=npx_1(link);
        px_L_1(end+1)=0;% sets up an artificial zero in front
        [a1 b1]=hist(px_L_1,max(px_L_1));
        histN_1_L(1:max(px_L_1))=histN_1_L(1:max(px_L_1))+a1;



        % for linked and splined for camera 1

        px_LS_1=npx_1(link_splined);
        
        px_LS_1(end+1)=0;% sets up an artificial zero in front
        [a1 b1]=hist(px_LS_1,max(px_LS_1));
        histN_1_LS(1:max(px_LS_1))=histN_1_LS(1:max(px_LS_1))+a1;  

        
        
        
        %-------------- Camera 2---------------------------
        
        
        % for only linked for camera 2
    
        px_L_2=npx_2(link);
        px_L_2(end+1)=0;% sets up an artificial zero in front
        [a2 b2]=hist(px_L_2,max(px_L_2));
        histN_2_L(1:max(px_L_2))=histN_2_L(1:max(px_L_2))+a2;



        % for linked and splined for camera 2

        px_LS_2=npx_2(link_splined);
        
        px_LS_2(end+1)=0;% sets up an artificial zero in front
        [a2 b2]=hist(px_LS_2,max(px_LS_2));
        histN_2_LS(1:max(px_LS_2))=histN_2_LS(1:max(px_LS_2))+a2;  

        
        %-------Camera 3----------
        
        % for only linked for camera 3
    
        px_L_3=npx_3(link);
        px_L_3(end+1)=0;% sets up an artificial zero in front
        [a3 b3]=hist(px_L_3,max(px_L_3));
        histN_3_L(1:max(px_L_3))=histN_3_L(1:max(px_L_3))+a3;



        % for linked and splined for camera 3

        px_LS_3=npx_3(link_splined);
        
        px_LS_3(end+1)=0;% sets up an artificial zero in front
        [a3 b3]=hist(px_LS_3,max(px_LS_3));
        histN_3_LS(1:max(px_LS_3))=histN_3_LS(1:max(px_LS_3))+a3;  

        %----------------Camera 4---------
        
        
        % for only linked for camera 4
    
        px_L_4=npx_4(link);
        px_L_4(end+1)=0;% sets up an artificial zero in front
        [a4 b4]=hist(px_L_4,max(px_L_4));
        histN_4_L(1:max(px_L_4))=histN_4_L(1:max(px_L_4))+a4;



        % for linked and splined for camera 4

        px_LS_4=npx_4(link_splined);
        
        px_LS_4(end+1)=0;% sets up an artificial zero in front
        [a4 b4]=hist(px_LS_4,max(px_LS_4));
        histN_4_LS(1:max(px_LS_4))=histN_4_LS(1:max(px_LS_4))+a4;  

              
end

% This figure shows bin width == 1 px
% figure;
% plot(histN_1_LS)
% figure;
% plot(histN_1_lns)


%----------- figure for camera 1----------------

bin_width=100;
bin_range=sort(linspace(length(histN_1_LS),bin_width));
%bin_range=bin_width:bin_width:length(histN1);

ss=size(bin_range,2);
sm_LS1(1)=sum(histN_1_LS(1:bin_range(2)));
sm_L1(1)=sum(histN_1_L(1:bin_range(2)));

for jj=2:ss
    sm_LS1(jj)=sum(histN_1_LS(bin_range(jj-1):bin_range(jj)));
    sm_L1(jj)=sum(histN_1_L(bin_range(jj-1):bin_range(jj)));
end
%figure,plot(histN_1_l)
subplot(2,2,1)
loglog(bin_range,sm_L1,'r')
%title('Linked')
%text(200,3000,'\leftarrow Linked (red line)')

hold on
%figure,plot(histN_1_LS)
loglog(bin_range,sm_LS1,'b')
%title('Linked and splined')
xlabel('bin range (px)');
ylabel('frequency')
%text(150,500,'\Rightarrow Linked and Splined(blue line)')

hold off

title('Camera 1')



%------figure for camera 2----

bin_width=100;
bin_range=sort(linspace(length(histN_2_LS),bin_width));
%bin_range=bin_width:bin_width:length(histN2);

ss=size(bin_range,2);
sm_LS2(1)=sum(histN_2_LS(1:bin_range(2)));
sm_L2(1)=sum(histN_2_L(1:bin_range(2)));

for jj=2:ss
    sm_LS2(jj)=sum(histN_2_LS(bin_range(jj-1):bin_range(jj)));
    sm_L2(jj)=sum(histN_2_L(bin_range(jj-1):bin_range(jj)));
end
%figure,plot(histN_2_l)
subplot(2,2,2)
loglog(bin_range,sm_L2,'r')
%title('Linked')
%text(200,3000,'\leftarrow Linked (red line)')

hold on
%figure,plot(histN_2_LS)
loglog(bin_range,sm_LS2,'b')
%title('Linked and splined')
xlabel('bin range (px)');
ylabel('frequency')
%text(150,500,'\Rightarrow Linked and Splined(blue line)')

hold off

title('Camera 2')

%------figure for camera 3------


bin_width=100;
bin_range=sort(linspace(length(histN_3_LS),bin_width));
%bin_range=bin_width:bin_width:length(histN3);

ss=size(bin_range,2);
sm_LS3(1)=sum(histN_3_LS(1:bin_range(2)));
sm_L3(1)=sum(histN_3_L(1:bin_range(2)));

for jj=2:ss
    sm_LS3(jj)=sum(histN_3_LS(bin_range(jj-1):bin_range(jj)));
    sm_L3(jj)=sum(histN_3_L(bin_range(jj-1):bin_range(jj)));
end
%figure,plot(histN_3_l)
subplot(2,2,4)
loglog(bin_range,sm_L3,'r')
%title('Linked')
%text(200,3000,'\leftarrow Linked (red line)')

hold on
%figure,plot(histN_3_LS)
loglog(bin_range,sm_LS3,'b')
%title('Linked and splined')
xlabel('bin range (px)');
ylabel('frequency')
%text(150,500,'\Rightarrow Linked and Splined(blue line)')

hold off

title('Camera 3')

%----figure for camera 4----

bin_width=100;
bin_range=sort(linspace(length(histN_4_LS),bin_width));
%bin_range=bin_width:bin_width:length(histN4);

ss=size(bin_range,2);
sm_LS4(1)=sum(histN_4_LS(1:bin_range(2)));
sm_L4(1)=sum(histN_4_L(1:bin_range(2)));

for jj=2:ss
    sm_LS4(jj)=sum(histN_4_LS(bin_range(jj-1):bin_range(jj)));
    sm_L4(jj)=sum(histN_4_L(bin_range(jj-1):bin_range(jj)));
end
%figure,plot(histN_4_l)
subplot(2,2,3)
loglog(bin_range,sm_L4,'r')
%title('Linked')
%text(200,3000,'\leftarrow Linked (red line)')

hold on
%figure,plot(histN_4_LS)
loglog(bin_range,sm_LS4,'b')
%title('Linked and splined')
xlabel('bin range (px)');
ylabel('frequency')
%text(150,500,'\Rightarrow Linked and Splined(blue line)')

hold off

title('Camera 4')





