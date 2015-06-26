%------xuap histogram------%
%n=19233;
first=7550;
Last=8000;

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

for n=first:Last;

    if mod(n,2)==0
        n
    end
    f_name_xuap=['E:\PTV\Working_foLder\WF_1\res\xuap.',num2Str(n)];
    xuap=load(f_name_xuap);
    sz=size(xuap);

    if sz(1,1)>0
        past=xuap(:,1);
        future=xuap(:,2);
        x=xuap(:,3);
        y=xuap(:,4);
        spLine=xuap(:,15);

        Link=find (future>0 & x> xmin & x< xmax & y>ymin);
        Link_spLined=find(future>0 & x> xmin & x< xmax & y>ymin & spLine==1);

    end

    f1=['E:\PTV\Working_foLder\WF_1\f_fiLes\f_',num2Str(n),'.mat'];
    clear f
    f=load(f1);

    for k=1:4;
        switch k
            case 1
                npx_1=f.a1_a;

                % for onLy Linked for camera 1
                px_L_1=npx_1(Link);
                px_L_1(end+1)=0;% sets up an artificiaL zero in front
                [a1 b1]=hist(px_L_1,max(px_L_1));
                histN_1_L(1:max(px_L_1))=histN_1_L(1:max(px_L_1))+a1;
                % for Linked and spLined for camera 1
                px_LS_1=npx_1(Link_spLined);
                px_LS_1(end+1)=0;% sets up an artificiaL zero in front
                [a1 b1]=hist(px_LS_1,max(px_LS_1));
                histN_1_LS(1:max(px_LS_1))=histN_1_LS(1:max(px_LS_1))+a1;

            case 2
                npx_2=f.a2_a;

                % for onLy Linked for camera 1
                px_L_2=npx_2(Link);
                px_L_2(end+1)=0;% sets up an artificiaL zero in front
                [a2 b2]=hist(px_L_2,max(px_L_2));
                histN_2_L(1:max(px_L_2))=histN_2_L(1:max(px_L_2))+a2;
                % for Linked and spLined for camera 1
                px_LS_2=npx_2(Link_spLined);
                px_LS_2(end+1)=0;% sets up an artificiaL zero in front
                [a2 b2]=hist(px_LS_2,max(px_LS_2));
                histN_2_LS(1:max(px_LS_2))=histN_2_LS(1:max(px_LS_2))+a2;



            case 3
                npx_3=f.a3_a;

                % for onLy Linked for camera 1
                px_L_3=npx_3(Link);
                px_L_3(end+1)=0;% sets up an artificiaL zero in front
                [a3 b3]=hist(px_L_3,max(px_L_3));
                histN_3_L(1:max(px_L_3))=histN_3_L(1:max(px_L_3))+a3;
                % for Linked and spLined for camera 1
                px_LS_3=npx_3(Link_spLined);
                px_LS_3(end+1)=0;% sets up an artificiaL zero in front
                [a3 b3]=hist(px_LS_3,max(px_LS_3));
                histN_3_LS(1:max(px_LS_3))=histN_3_LS(1:max(px_LS_3))+a3;

            case 4
                npx_4=f.a4_a;

                % for onLy Linked for camera 1
                px_L_4=npx_4(Link);
                px_L_4(end+1)=0;% sets up an artificiaL zero in front
                [a4 b4]=hist(px_L_4,max(px_L_4));
                histN_4_L(1:max(px_L_4))=histN_4_L(1:max(px_L_4))+a4;
                % for Linked and spLined for camera 1
                px_LS_4=npx_4(Link_spLined);
                px_LS_4(end+1)=0;% sets up an artificiaL zero in front
                [a4 b4]=hist(px_LS_4,max(px_LS_4));
                histN_4_LS(1:max(px_LS_4))=histN_4_LS(1:max(px_LS_4))+a4;

        end

    end
end



%%%-------------figures


% This figure shows bin width == 1 px
% figure;
% loglog(histN_1_LS)
% figure;
% loglog(histN_1_Lns)


bin_width=100;
bin_range=sort(linspace(length(histN_1_LS),bin_width));
%bin_range=bin_width:bin_width:Length(histN1);

ss=size(bin_range,2);
sm_LS(1)=sum(histN_1_LS(1:bin_range(2)));
sm_L(1)=sum(histN_1_L(1:bin_range(2)));

for jj=2:ss
    sm_LS(jj)=sum(histN_1_LS(bin_range(jj-1):bin_range(jj)));
    sm_L(jj)=sum(histN_1_L(bin_range(jj-1):bin_range(jj)));
end
%figure,loglog(histN_1_L)
subplot(2,2,1)

loglog(bin_range,sm_L,'r',bin_range,sm_LS,'b')
%titLe('Linked')
%text(200,3000,'\Leftarrow Linked (red Line)')

%figure,loglog(histN_1_LS)

%titLe('Linked and spLined')
xlabel('bin range (px)');
ylabel('frequency')
%text(150,500,'\Rightarrow Linked and SpLined(bLue Line)')
title('Camera 1')




%--------------


bin_width=100;
bin_range=sort(linspace(length(histN_2_LS),bin_width));
%bin_range=bin_width:bin_width:Length(histN1);

ss=size(bin_range,2);
sm_LS(1)=sum(histN_2_LS(1:bin_range(2)));
sm_L(1)=sum(histN_2_L(1:bin_range(2)));

for jj=2:ss
    sm_LS(jj)=sum(histN_2_LS(bin_range(jj-1):bin_range(jj)));
    sm_L(jj)=sum(histN_2_L(bin_range(jj-1):bin_range(jj)));
end
%figure,loglog(histN_1_L)
subplot(2,2,2)

loglog(bin_range,sm_L,'r',bin_range,sm_LS,'b')

%titLe('Linked')
%text(200,3000,'\Leftarrow Linked (red Line)')

%figure,loglog(histN_1_LS)

%titLe('Linked and spLined')
xlabel('bin range (px)');
ylabel('frequency')
%text(150,500,'\Rightarrow Linked and SpLined(bLue Line)')

title('Camera 2')




%------------


bin_width=100;
bin_range=sort(linspace(length(histN_3_LS),bin_width));
%bin_range=bin_width:bin_width:Length(histN1);

ss=size(bin_range,2);
sm_LS(1)=sum(histN_3_LS(1:bin_range(2)));
sm_L(1)=sum(histN_3_L(1:bin_range(2)));

for jj=2:ss
    sm_LS(jj)=sum(histN_3_LS(bin_range(jj-1):bin_range(jj)));
    sm_L(jj)=sum(histN_3_L(bin_range(jj-1):bin_range(jj)));
end
%figure,loglog(histN_1_L)
subplot(2,2,4)

loglog(bin_range,sm_L,'r',bin_range,sm_LS,'b')

%titLe('Linked')
%text(200,3000,'\Leftarrow Linked (red Line)')

xlabel('bin range (px)');
ylabel('frequency')
%figure,loglog(histN_1_LS)

%titLe('Linked and spLined')

%text(150,500,'\Rightarrow Linked and SpLined(bLue Line)')

title('Camera 3')
%---------


bin_width=100;
bin_range=sort(linspace(length(histN_4_LS),bin_width));
%bin_range=bin_width:bin_width:Length(histN1);

ss=size(bin_range,2);
sm_LS(1)=sum(histN_4_LS(1:bin_range(2)));
sm_L(1)=sum(histN_4_L(1:bin_range(2)));

for jj=2:ss
    sm_LS(jj)=sum(histN_4_LS(bin_range(jj-1):bin_range(jj)));
    sm_L(jj)=sum(histN_4_L(bin_range(jj-1):bin_range(jj)));
end
%figure,loglog(histN_1_L)
subplot(2,2,3)

loglog(bin_range,sm_L,'r',bin_range,sm_LS,'b')

%titLe('Linked')
%text(200,3000,'\Leftarrow Linked (red Line)')

xlabel('bin range (px)');
ylabel('frequency')
%figure,loglog(histN_1_LS)

%titLe('Linked and spLined')

%text(150,500,'\Rightarrow Linked and SpLined(bLue Line)')

title('Camera 4')
