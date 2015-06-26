
clear all
clear all
clear all
%close all
clc
clc
count=0;
eps=9e-5;
r_0=1.0e-3;
tau_b=((r_0)^(2/3))/(eps^(1/3))*250;
pp=0;
%%% WF 70 and WF 114 are invalid and WF 105,11
WF=[1 2 6 10 12 13 15 16 18 21 22 23 25 26 27 28 29 30 31 32 33 34 35 36 37 38 41 42 45 46 50 51 52 53 55 56 58 59 61 62 63 64 65 66 67 68 69  71 72 73 74 75 76 77 79 80 81 82 83 84 86 87 88 89 90 92 93 95 97 98  99 100 101 102 103 104 106 107 108 109 110 111 112 113 115 116 117 118 119 120 121 122 123];
fold='D:\Turbulent_data_backup\flow_induced_aggregates_breakage_at_100rpm\WF';
prefactor=[0.01:0.05:2];
for pre=prefactor;
    pp=pp+1;
int_window=round(prefactor*tau_b);



axis=zeros(length(WF),2);
for j=1:length(WF)
   [pre j]
    par_var=[fold,num2str(WF(j)),'/parent_variables'];
    child_var=[fold,num2str(WF(j)),'/children_variables'];
    break_var=[fold,num2str(WF(j)),'/break_frame'];
    
    parent=load(par_var);
    children=load(child_var);
    break_frame=load(break_var);
    
    break_fr_ptv=break_frame.break_frame_ptv;
    
    
    Aij=parent.parent_Aij;
    sij=parent.parent_Sij;
    
    v1=parent.parent_v1;
    v2=parent.parent_v2;
    v3=parent.parent_v3;
    e1=(parent.parent_e1);
    e2=(parent.parent_e2);
    e3=(parent.parent_e3);
    
    time=parent.T;
    x_=parent.X;
    y_=parent.Y;
    z_=parent.Z;
    
    
    
    index=find(break_fr_ptv(1)==time);
    
    minor_axis=nanmean(parent.parent_max_shrt_ax_mm(1:index));
    major_axis=nanmean(parent.parent_max_lng_ax_mm(1:index));
    
    axis(j,1)=minor_axis;
    axis(j,2)=major_axis;
    
    c_time=1;
    dt=1/250;
    B=[];
    lam_1=[];
    lam_2=[];
    lam_3=[];
    w1=[];
    w2=[];
    w3=[];
    
    
    B(c_time,1:3,1:3)=[1 0 0;0 1 0;0 0 1];
    B(c_time,1:3,1:3)=B(c_time,1:3,1:3);%*1e-3;
    
    
        int_t0_tb=find(break_fr_ptv(1)==time);
        if isempty(int_t0_tb)
            continue;
        else
            frame_break=int_t0_tb;
        end
        count=count+1;
        
        if int_window>int_t0_tb
            int_window=int_t0_tb;
            frame_0=1;
        end
        frame_0=int_t0_tb-int_window+1;
        
        
        %%%% Integration over time
        
        
        %for k=1:length(bin)
            for frame=frame_0+5:frame_break
                h=squeeze(Aij(frame,:,:));
                c_time=c_time+1;
                dbdt=h*squeeze(B(c_time-1,1:3,1:3));
                B(c_time,1:3,1:3)=squeeze(B(c_time-1,1:3,1:3))+dt*dbdt;
                lam_1(c_time,1:3)=v1(frame,:)';
                lam_2(c_time,1:3)=v2(frame,:)';
                lam_3(c_time,1:3)=v3(frame,:)';
            end
        
        % Eigen analysis of Cauchy Green
        for p=1:size(B,1)
            W=squeeze(B(p,1:3,1:3))*squeeze(B(p,1:3,1:3))';
            [V,D] = eig(W);
            [D,k] = sort(diag(D)); 	% ascending order
            D = diag(D(end:-1:1));  % descending order
            V = V(:,k(end:-1:1)); 	% the same order for eigenvectors
            w1(p)=D(1,1)^0.5;       % eigen value
            w2(p)=D(2,2)^0.5;
            w3(p)=D(3,3)^0.5;
            [w1 w2 w3];
            W1(p,1:3)=V(:,1);
            W2(p,1:3)=V(:,2);
            W3(p,1:3)=V(:,3);
        end
        %end
        
        cX1=children.cx{1}(1); %%% Children coordinate
        cY1=children.cy{1}(1);
        cZ1=children.cz{1}(1);
        
        cX2=parent.X(frame_break);
        cY2=parent.Y(frame_break);
        cZ2=parent.Z(frame_break);
        
        r=[cX1-cX2;cY1-cY2;cZ1-cZ2];
        r=r/sum(r.^2)^0.5; %%% Separation @ breakage between parent and fragment
        
        a=W1(end,:);
        cos_W1_r(pp,j)=dot(a',r);%%% Angle between separation vetor and W1
        cos_lambda1_r(j)=dot(lam_1(end,:)',r);%%% Angle between separation vetor and lambda_1
        cos_W1_l1(j)=dot(a',lam_1(end,:)');%%% Angle between W1 and lambda_1
        
        b=W2(end,:);
        cos_W2_r(pp,j)=dot(b',r);
        cos_lambda2_r(j)=dot(lam_2(end,:)',r);
        cos_W2_l2(j)=dot(b',lam_2(end,:)');
        
        c=W3(end,:);
        cos_W3_r(pp,j)=dot(c',r);
        cos_lambda3_r(j)=dot(lam_3(end,:)',r);
        cos_W3_l3(j)=dot(c',lam_3(end,:)');
        
        w1_stat(j)=w1(end);
        w2_stat(j)=w2(end);
        w3_stat(j)=w3(end);
    %end
end
end

figure('Name','pdf of abs(cos(W1..3,r))','NumberTitle','on');
hold on;
[nout1,xout1]=nhist(abs(cos_W1_r),10);
m1=sum(nout1.*xout1)/sum(nout1);
plot(xout1,nout1,'r')

[nout2,xout2]=nhist(abs(cos_W2_r),10);
plot(xout2,nout2,'g')

[nout3,xout3]=nhist(abs(cos_W3_r),10);
m3=sum(nout3.*xout3)/sum(nout3);
plot(xout3,nout3,'b')
ylim([0 3])
tit=['pdf of abs(cos(W1..3,r)) integrated over ',num2str(prefactor),' \tau_b', ' and r_0=',num2str(r_0)];
title(tit)



for jj=1:20
    var1=abs(cos_W1_r(jj,:));
    var1=var1(var1>0);
    [nout11,xout11]=nhist(var1);
m1(jj)=sum(nout11.*xout11)/sum(nout11);
end
    
for jj=1:20
    var3=abs(cos_W3_r(jj,:));
    var3=var3(var3>0);
    [nout33,xout33]=nhist(var3);
m3(jj)=sum(nout33.*xout33)/sum(nout33);
end

figure;
plot([0.1:0.1:2],m1,'r')
hold on
plot([0.1:0.1:2],m3,'b')



% figure('Name','pdf of abs(cos(lambda..3,r))','NumberTitle','on');
% hold on;
% [nout,xout]=nhist(abs(cos_lambda1_r),10);
% plot(xout,nout,'r')
% [nout,xout]=nhist(abs(cos_lambda2_r),10);
% plot(xout,nout,'g')
% [nout,xout]=nhist(abs(cos_lambda3_r),10);
% plot(xout,nout,'b')
% ylim([0 3])
% tit=['pdf of abs(cos(lambda..3,r)) integrated over ',num2str(prefactor),' \tau_b'];
% title(tit)


% figure('Name','pdf of log(w2)/log(w1) & log(w1)/log(w3)@breakage','NumberTitle','on');
% [n1,x1]=nhist(log(w2_stat)./log(w1_stat));
% plot(x1,n1,'r')
% hold on;
% [n1,x1]=nhist(log(w1_stat)./log(w3_stat));
% plot(x1,n1,'g')
% title('pdf of log(w2)/log(w1) & log(w1)/log(w3)@breakage')



% figure('Name','pdf of w1..3^{0.5}@breakage','NumberTitle','on');
% hold on;
% [n1,x1]=nhist(w1_stat.^0.5);
% plot(x1,n1,'r')
% 
% % figure('Name','pdf of w2^{0.5}@breakage','NumberTitle','on');
% [n2,x2]=nhist(w2_stat.^0.5);
% plot(x2,n2,'g')
% 
% % figure('Name','pdf of w3^{0.5}@breakage','NumberTitle','on');
% [n3,x3]=nhist(w3_stat.^0.5);
% plot(x3,n3,'b')
% title('pdf of w1..3^{0.5}@breakage')
% 
% 
% figure;
% hist(abs(cos_W1_l1));
% title('cos(W_1,\lambda_1)')
% 
% figure('Name','Eigen values of Cauchy Green ','NumberTitle','on');
% plot(log(w1),'r')
% hold on
% plot(log(w2),'g')
% plot(log(w3),'b')
% xlabel('Frame')
% ylabel('\bf ln (l_{i})')
% legend('w1','w2','w3')
% title('Eigen values of Cauchy Green')
% 
% figure('Name','log(w1_stat),log(w2_stat) ','NumberTitle','on');
% scatter(log(w1_stat),log(w2_stat))
% title('log(w1-stat),log(w2-stat)')



















