%------xuap histogram------%
%n=19233;
first=7500;
last=8000;

histN_1_ls=zeros(1,10000);
histN_2_ls=zeros(1,10000);
histN_3_ls=zeros(1,10000);
histN_4_ls=zeros(1,10000);

histN_1_l=zeros(1,10000);
histN_2_l=zeros(1,10000);
histN_3_l=zeros(1,10000);
histN_4_l=zeros(1,10000);

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
        spline=xuap(:,15);

        link=find(future>0);
        link_splined=find(future>0 & spline==1);

    end

    f1=['E:\PTV\Working_folder\WF_1\f_files\f_',num2Str(n),'.mat'];
    clear f
    f=load(f1);

    % for linked and splined for camera 1
    npx_1=f.a1_a;
    px_ls_1=npx_1(link_splined);
    px_ls_1(end+1)=0;% sets up an artificial zero in front
    [a b]=hist(px_ls_1,max(px_ls_1));
    histN_1_ls(1:max(px_ls_1))=histN_1_ls(1:max(px_ls_1))+a;

    % for linked but not splined for camera 1

    npx_1=f.a1_a;
    px_l_1=npx_1(link);
    px_l_1(end+1)=0;% sets up an artificial zero in front
    [a b]=hist(px_l_1,max(px_l_1));
    histN_1_l(1:max(px_l_1))=histN_1_l(1:max(px_l_1))+a;


end

% This figure shows bin width == 1 px 
% figure;
% plot(histN_1_ls)
% figure;
% plot(histN_1_lns)

bin_width=100;
bin_range=sort(linspace(length(histN_1_ls),bin_width));
%bin_range=bin_width:bin_width:length(histN1);

ss=size(bin_range,2);
sm_ls(1)=sum(histN_1_ls(1:bin_range(2)));
sm_l(1)=sum(histN_1_l(1:bin_range(2)));

for jj=2:ss
    sm_ls(jj)=sum(histN_1_ls(bin_range(jj-1):bin_range(jj)));
    sm_l(jj)=sum(histN_1_l(bin_range(jj-1):bin_range(jj)));
end

figure,plot(histN_1_ls)
figure,plot(bin_range,sm_ls)
title('Linked and splined')

figure,plot(histN_1_l)
figure,plot(bin_range,sm_l)
title('Linked')

% npx_2=f.a2_a;
% 
% npx_3=f.a3_a;
% 
% npx_4=f.a4_a;












