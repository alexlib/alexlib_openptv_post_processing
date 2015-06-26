

clear all
clear all
clear all
%close all
clc
WF=[1 2 6 10  12 13 15 16 18 21 22 23 25 26 27 28 29 30 31 32 33 34 35 36 37 38 41 42 45 46 50 51 52 53 55 56 58 59 61 62 63 64 65 66 67 68 69  71 72 73 74 75 76 77 79 80 81 82 83 84 86 87 88 89  92 93 95 97 98  99 100 101 102  104  106 107 108 109 110 111 112 113];

fold='D:\Turbulent_data_backup\flow_induced_aggregates_breakage_at_100rpm\WF';
dt=1/250;
eps=9e-5;
tau_b=((1.0e-3)^(2/3))/(eps^(1/3))*250;
nu=1e-6;



    count=0;
    
    int_window=round(1*tau_b);
     
    for j=1:length(WF)
  
        clear B  w1 w2 w3 W1 W2 W3 time t0_tb frame_break...
          frame_0 x y z x_ y_ z_ Aij Sij rij v1 v2 v3 e1 e2 e3...
          cX1 cY1 cZ1 cX2 cY2 cZ2  
        
        
        
        par=load([fold,num2str(WF(j)),'/parent_variables_v3']);
        time=par.pTIME;
        
        break_frame=load([fold,num2str(WF(j)),'/break_frame']);
        break_frame_manual=load([fold,num2str(WF(j)),'/breakage_manual']);
        
        break_fr_ptv=break_frame.break_frame_ptv;
        break_fr_manual=break_frame_manual.break_fr_manual;
        
  
        t0_tb=find(break_fr_ptv==time);%%% Breakage detected by PTV
        
        if isempty(t0_tb)==1
            min_search_ptv=abs(time-break_fr_ptv);
            index_ptv=find(min_search_ptv==min(min_search_ptv));
            t0_tb=index_ptv;
        end
        
        
        %     t0_tb=find(break_fr_manual==time);%%% Breakage detected manually
        %     if isempty(t0_tb)==1
        %         min_search_manual=abs(time-break_fr_manual);
        %         index_manual=find(min_search_manual==min(min_search_manual));
        %         t0_tb=index_manual;
        %     end
        
        if  int_window<=t0_tb
            count=count+1;
            %WF_(pf,count)=WF(j);%%% Surviving WF for a certain int_window
            
        children=load([fold,num2str(WF(j)),'/children_variables_v3']);
        
        minor_axis(count)=nanmean(par.parent_max_shrt_ax_mm(1:t0_tb));
        major_axis(count)=nanmean(par.parent_max_lng_ax_mm(1:t0_tb));
            
       
        
        
        frame_break=t0_tb;
        frame_0=t0_tb-int_window+1;
        
        x=par.X;
        y=par.Y;
        z=par.Z;
        x_=x(frame_0:frame_break);
        y_=y(frame_0:frame_break);
        z_=z(frame_0:frame_break);
        
        Aij=par.parent_Aij;
        sij=par.parent_Sij;
        rij=par.parent_rij;
        
        v1=par.parent_v1;
        v2=par.parent_v2;
        v3=par.parent_v3;
        e1=par.parent_e1;
        e2=par.parent_e2;
        e3=par.parent_e3;
        
        
        c_time=1;
        
        B(c_time,1:3,1:3)=[1 0 0;0 1 0;0 0 1];
        B(c_time,1:3,1:3)=B(c_time,1:3,1:3);%*1e-3;
        
   
            
            % Integration over time
            for frame=frame_0:frame_break
                h=squeeze(Aij(frame,:,:));
                c_time=c_time+1;
                dbdt=h*squeeze(B(c_time-1,1:3,1:3));
                B(c_time,1:3,1:3)=squeeze(B(c_time-1,1:3,1:3))+dt*dbdt;
                Sij_v1(count,c_time-1,1:3)=v1(frame,:);
                Sij_v2(count,c_time-1,1:3)=v2(frame,:);
                Sij_v3(count,c_time-1,1:3)=v3(frame,:);
                Sij_e1(count,c_time-1)=e1(frame);
                Sij_e2(count,c_time-1)=e2(frame);
                Sij_e3(count,c_time-1)=e3(frame);
                epsilon(count,c_time-1)=par.parent_strain(frame).*2*nu;
                
                Sij_int_win(count,c_time-1,1:3,1:3)=sij(frame,1:3,1:3);
                rij_int_win(count,c_time-1,1:3,1:3)=rij(frame,1:3,1:3);
                Aij_int_win(count,c_time-1,1:3,1:3)=Aij(frame,1:3,1:3);
                
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
                
                W1(p,1:3)=V(:,1);
                W2(p,1:3)=V(:,2);
                W3(p,1:3)=V(:,3);
            end
            
            w1_all(count,1:size(B,1))=w1; %% along the int window
            w2_all(count,1:size(B,1))=w2;
            w3_all(count,1:size(B,1))=w3;
            
            
            
            
            
            cX1=children.c_X(1); %%% Children coordinate
            cY1=children.c_Y(1);
            cZ1=children.c_Z(1);
            
            cX2=par.X(frame_break);
            cY2=par.Y(frame_break);
            cZ2=par.Z(frame_break);
            
            r1(count,1:3)=[cX1-cX2;cY1-cY2;cZ1-cZ2];
            r(count,1:3)=r1(count,1:3)./sum(r1(count,1:3).^2)^0.5; %%% Separation @ breakage between parent and fragment
            
            a=W1(end,:);
            cos_W1_r(count)=dot(a',r(count,1:3));%%% Angle between separation vetor and W1
            cos_lambda1_r(count)=dot(squeeze(Sij_v1(count,end,:))',r(count,1:3));%%% Angle between separation vetor and lambda_1
            cos_W1_l1(count)=dot(a',squeeze(Sij_v1(count,end,:))');%%% Angle between W1 and lambda_1
            
            b=W2(end,:);
            cos_W2_r(count)=dot(b',r(count,1:3));
            cos_lambda2_r(count)=dot(squeeze(Sij_v2(count,end,:))',r(count,1:3));
            cos_W2_l2(count)=dot(b',squeeze(Sij_v2(count,end,:))');
            
            c=W3(end,:);
            cos_W3_r(count)=dot(c',r(count,1:3));
            cos_lambda3_r(count)=dot(squeeze(Sij_v3(count,end,:))',r(count,1:3));
            cos_W3_l3(count)=dot(c',squeeze(Sij_v3(count,end,:))');
            
            w1_at_break(count)=w1(end);
            w2_at_break(count)=w2(end);
            w3_at_break(count)=w3(end);
            
            Sij_e1_at_break(count)=Sij_e1(end);
            Sij_e2_at_break(count)=Sij_e2(end);
            Sij_e3_at_break(count)=Sij_e3(end);
            
            
            Sij_at_break(count)=Sij_e1(end)^2+Sij_e2(end)^2+Sij_e3(end)^2;
            %Sij_at_break(count)=parent.parent_strain(frame_break);
            int_eps(count)=sum(epsilon(count,:));
            
            
        end
        
    end


%%%-----projection of strain(full tensor) along separation vector

for exp=1:count % no of exp
    
    r2(exp,1:3) = cross(rand(1,3),r(exp,:));
    r2(exp,1:3) = r2(exp,1:3)./(sum(r2(exp,1:3).^2).^0.5);
    r3(exp,1:3) = cross(r(exp,1:3),r2(exp,1:3));
    sep_cord(exp,1:3,1:3)=[r(exp,1:3); r2(exp,1:3); r3(exp,1:3)];
    
    %line([-sep_cord(1,1) sep_cord(1,1)], [-sep_cord(1,2) sep_cord(1,2)], [-sep_cord(1,3) sep_cord(1,3)],'LineWidth',2,'color','r','LineStyle','-')
    %line([-sep_cord(2,1) sep_cord(2,1)], [-sep_cord(2,2) sep_cord(2,2)], [-sep_cord(2,3) sep_cord(2,3)],'LineWidth',2,'color','g','LineStyle','-')
    %line([-sep_cord(3,1) sep_cord(3,1)], [-sep_cord(3,2) sep_cord(3,2)], [-sep_cord(3,3) sep_cord(3,3)],'LineWidth',2,'color','b','LineStyle','-')
    
   %c=[squeeze(Sij_v1(exp,end,:))';squeeze(Sij_v2(exp,end,:))';squeeze(Sij_v3(exp,end,:))'];
   %d=c*squeeze(Sij_int_win(exp,end,:,:))*c';
    
   for t=1:int_window
        rot_Sij(exp,t,1:3,1:3)=squeeze(sep_cord(exp,1:3,1:3))*squeeze(Sij_int_win(exp,t,:,:))*squeeze(sep_cord(exp,1:3,1:3))';
  
    end
end

for exp=1:count
    %hold on
    figure;
    %plot(squeeze(Sij_int_win(exp,:,1,1)))
    plot(squeeze(rot_Sij(exp,:,1,1)))
    title(['WF=', num2str(WF(exp))])
end

for exp=1:count
    s11_break(exp)=squeeze(rot_Sij(exp,end,1,1));
end
figure;
scatter(major_axis,s11_break,'.')
xlabel('major axis')
ylabel('S11 at break')


for eps=1:counter
    int_eps(eps)=sum(epsilon(eps,:))*(1/250);
end
figure;
scatter(major_axis,int_eps)

% figure(1);hold on;
% for jj=1:length(prefactor)
%     var1=abs(cos_W1_r(jj,1:counter(1,jj)));
%     [nout1,xout1]=nhist(var1);
%     m1(jj)=sum(nout1.*xout1)/sum(nout1);
%     plot(xout1,nout1,'r')
% end
% title(' pdf of cos(W_{i=1},r)')


figure;
[nout1,xout1]=nhist(abs(cos_W1_r),10);
plot(xout1,nout1,'r')
title(' pdf of cos(W_{i=1},r)')




% figure(2);hold on;
% for jj=1:length(prefactor)
%     var3=abs(cos_W3_r(jj,1:counter(1,jj)));
%     [nout3,xout3]=nhist(var3);
%     m3(jj)=sum(nout3.*xout3)/sum(nout3);
%     plot(xout3,nout3,'b')
% end
% title(' pdf of cos(W_{i=3},r)')


figure
[nout3,xout3]=nhist(abs(cos_W3_r));
plot(xout3,nout3,'b')
title(' pdf of cos(W_{i=3},r)')



% figure(3);
% plot(prefactor,m1,'r')
% hold on
% plot(prefactor,m3,'b')
% plot(prefactor,m1-m3,'k')
% legend('\langlecos(W_{i=1},r)\rangle','\langlecos(W_{i=3},r)\rangle','W1-W3')
% xlabel('\tau_b')
% ylabel('\langlecos(W_{i=1,3},r)\rangle')
% 


% figure(4);hold on;
% for jj=1:length(prefactor)
%     var5=abs(cos_lambda1_r(jj,1:counter(1,jj)));
%     [nout5,xout5]=nhist(var5);
%     m5(jj)=sum(nout5.*xout5)/sum(nout5);
%     plot(xout5,nout5,'r')
% end
% title(' pdf of cos(\lambda_{i=1},r)')



figure
[nout5,xout5]=nhist(abs(cos_lambda1_r),10);
plot(xout5,nout5,'r')
title(' pdf of cos(\lambda_{i=1},r)')



% figure(5);hold on;
% for jj=1:length(prefactor)
%     var7=abs(cos_lambda3_r(jj,1:counter(1,jj)));
%     [nout7,xout7]=nhist(var7);
%     m7(jj)=sum(nout7.*xout7)/sum(nout7);
%     plot(xout7,nout7,'b')
% end
% title(' pdf of cos(\lambda_{i=3},r)')


figure
[nout7,xout7]=nhist(abs(cos_lambda3_r),10);
plot(xout7,nout7,'b')
title(' pdf of cos(\lambda_{i=3},r)')



% figure(6);
% plot(prefactor,m5,'r.')
% hold on
% plot(prefactor,m7,'b.')
% plot(prefactor,m5-m7,'k')
% legend('\lambda_1','\lambda_3','\lambda_1-\lambda_3')
% 
% xlabel('\tau_b')
% ylabel('\langlecos(\lambda_{i=1,3},r)\rangle')



% figure(7);
% plot(prefactor,m1,'r')
% hold on
% plot(prefactor,m3,'b')
% %plot(prefactor,m1-m3,'k')
% 
% plot(prefactor,m5,'r.')
% plot(prefactor,m7,'b.')
% %plot(prefactor,m5-m7,'k')

% legend('W1','W3','\lambda_1','\lambda_3')
% 
% xlabel('\tau_b')
% ylabel('\langlecos(W_{i=1,3},r)\rangle,\langlecos(\lambda_{i=1,3},r)\rangle')



% figure(8);
% scatter(Sij_e1_at_break(4,1:counter(1,4)),Sij_e2_at_break(4,1:counter(1,4)),100,major_axis(4,1:counter(1,4)),'filled')
% xlabel('\Lambda_1')
% ylabel('\Lambda_2')
% caxis([min(major_axis(4,1:counter(1,4))), max(major_axis(4,1:counter(1,4)))])
% colorbar 


figure;
scatter(Sij_e1_at_break,Sij_e2_at_break,100,major_axis,'filled')
xlabel('\Lambda_1')
ylabel('\Lambda_2')
caxis([min(major_axis), max(major_axis)])
colorbar 



% figure(9);hold on;
% for jj=1:length(prefactor)
%     var8=abs(w1_at_break(jj,1:counter(1,jj)));
%     [nout8,xout8]=nhist(var8);
%     m8(jj)=sum(nout8.*xout8)/sum(nout8);
%     plot(xout8,nout8,'r')
% end

figure
[nout8,xout8]=nhist(abs(w1_at_break));
plot(xout8,nout8,'r')




% figure(10);hold on;
% for jj=1:length(prefactor)
%     var9=abs(w2_at_break(jj,1:counter(1,jj)));
%     [nout9,xout9]=nhist(var9);
%     m9(jj)=sum(nout9.*xout9)/sum(nout9);
%     plot(xout9,nout9,'g')
% end

figure
[nout9,xout9]=nhist(abs(w2_at_break));
plot(xout9,nout9,'g')



% figure(11);
% scatter(log(w1_at_break(4,1:counter(1,4))),log(w2_at_break(4,1:counter(1,4))),100,major_axis(4,1:counter(1,4)),'filled')
% xlabel('log(w1)')
% ylabel('log(w2)')
% title(['\tau_b= ',num2str(prefactor)])
% caxis([min(major_axis(4,1:counter(1,4))), max(major_axis(4,1:counter(1,4)))])
% colorbar


figure;
scatter(log(w1_at_break),log(w2_at_break),100,major_axis,'filled')
xlabel('log(w1)')
ylabel('log(w2)')
title(['\tau_b= ',num2str(prefactor)])
caxis([min(major_axis), max(major_axis)])
colorbar




% figure(12);
% [nout10,xout10]=nhist(Sij_at_break(4,1:counter(1,4)));
% plot(xout10,nout10)

figure;
[nout10,xout10]=nhist(Sij_at_break);
plot(xout10,nout10)


% figure(13);
% scatter(major_axis(4,1:counter(1,4)),Sij_at_break(4,1:counter(1,4)))

figure;
scatter(major_axis,Sij_at_break)



% figure(14);
% [nout11,xout11]=nhist(log(w2_at_break(4,1:counter(1,4)))./log(w1_at_break(4,1:counter(1,4))));
% plot(xout11,nout11,'r')
% hold on;
% [nout12,xout12]=nhist(log(w1_at_break(4,1:counter(1,4)))./log(w3_at_break(4,1:counter(1,4))));
% plot(xout12,nout12,'g')
% title('pdf of log(w2)/log(w1) & log(w1)/log(w3)@breakage')
% legend('w2/w1','w1/w3')


figure;
[nout11,xout11]=nhist(log(w2_at_break)./log(w1_at_break));
plot(xout11,nout11,'r')
hold on;
[nout12,xout12]=nhist(log(w1_at_break)./log(w3_at_break));
plot(xout12,nout12,'g')
title('pdf of log(w2)/log(w1) & log(w1)/log(w3)@breakage')
legend('w2/w1','w1/w3')




% figure('Name','pdf of w2^{0.5}@breakage','NumberTitle','on');
[nout16,xout16]=nhist(w2_at_break(4,1:counter(1,4)));
plot(xout16,nout16,'g')

% figure('Name','pdf of w3^{0.5}@breakage','NumberTitle','on');
[nout17,xout17]=nhist(w3_at_break(4,1:counter(1,4)));
plot(xout17,nout17,'b')
title('pdf of w_{1,2,3}@breakage')


figure;
hold on;
[nout15,xout15]=nhist(w1_at_break);
plot(xout15,nout15,'r')

% figure('Name','pdf of w2^{0.5}@breakage','NumberTitle','on');
[nout16,xout16]=nhist(w2_at_break);
plot(xout16,nout16,'g')

% figure('Name','pdf of w3^{0.5}@breakage','NumberTitle','on');
[nout17,xout17]=nhist(w3_at_break);
plot(xout17,nout17,'b')
title('pdf of w_{1,2,3}@breakage')






% figure(16);
% [nout18,xout18]=nhist(abs(cos_W1_l1(4,1:counter(1,4))));
% plot(xout18,nout17,'r')
% title('cos(W_1,\lambda_1)')


figure;
[nout18,xout18]=nhist(abs(cos_W1_l1));
plot(xout18,nout17,'r')
title('cos(W_1,\lambda_1)')



% w1_alng_traj=[];
% w2_alng_traj=[];
% 
% for i=4%:length(counter)
%     j=0;
%     while j<counter(i)
%         j=j+1;
% w1_All=squeeze(w1_all(i,j,:)); %% (pf,number of wf survived that pf, all w1's along the int_window)
% w1_All=w1_All(w1_All>0);
% w1_alng_traj=[w1_alng_traj;w1_All];
%   
% w2_All=squeeze(w2_all(i,j,:));
% w2_All=w2_All(w2_All>0);
% w2_alng_traj=[w2_alng_traj;w1_All];
% 
%     end
% end
% figure(17);
% scatter(log(w1_alng_traj),log(w2_alng_traj),'.')
% xlabel('log(w1 alng traj )')
% ylabel('log(w2 alng traj)')



w1_alng_traj=[];
w2_alng_traj=[];


    j=0;
    while j<counter
        j=j+1;
w1_All=w1_all(j,:)'; %% (pf,number of wf survived that pf, all w1's along the int_window)
w1_All=w1_All(w1_All>0);
w1_alng_traj=[w1_alng_traj;w1_All];
  
w2_All=squeeze(w2_all(j,:))';
w2_All=w2_All(w2_All>0);
w2_alng_traj=[w2_alng_traj;w1_All];

    end

figure;
scatter(log(w1_alng_traj),log(w2_alng_traj),'.')
xlabel('log(w1 alng traj )')
ylabel('log(w2 alng traj)')







%%% rij----------------

% Intial alignment with the strain eigen frame


% l_lambda(1,1:3,1:3)=[v1(1,1) v1(1,2) v1(1,3); v2(1,1) v2(1,2) v2(1,3); v3(1,1) v3(1,2) v3(1,3)];
% counter_l_lambda=1;
% for frame=frame_break-int_window+1:frame_break%size(time)%175
%     h=squeeze(rij(frame,:,:));
%     counter_l_lambda=counter_l_lambda+1;
%     dbdt_l_lambda=h*squeeze(l_lambda(counter_l_lambda-1,1:3,1:3));
%     l_lambda(counter_l_lambda,1:3,1:3)=squeeze(l_lambda(counter_l_lambda-1,1:3,1:3))+dt*dbdt_l_lambda;
% end




% Intial alignment with the strain eigen frame
l_lambda(int_window,1:3,1:3)=[Sij_v1(end,1) Sij_v1(end,2) Sij_v1(end,3); Sij_v2(end,1) Sij_v2(end,2) Sij_v2(end,3); Sij_v3(end,1) Sij_v3(end,2) Sij_v3(end,3)];
counter_l_lambda=0;
for frame=frame_break:-1:frame_break-int_window+2 %size(time)%175
    r_ij=squeeze(rij(frame,:,:));
    counter_l_lambda=counter_l_lambda+1;
    dbdt_l_lambda=r_ij*squeeze(l_lambda(int_window-counter_l_lambda+1,1:3,1:3));
    l_lambda(int_window-counter_l_lambda,1:3,1:3)=squeeze(l_lambda(int_window-counter_l_lambda+1,1:3,1:3))-dt*dbdt_l_lambda;
end


% Initial alignment with separation vector; the other two directions are
% arbitrary (but normal to each other of course)
r2 = cross(rand(1,3),r);
r2 = r2./(sum(r2.^2).^0.5);
r3 = cross(r,r2);
l_r(int_window,1:3,1:3)=[r(1) r(2) r(3); r2; r3];
counter_l_r=0;

for frame=frame_break:-1:frame_break-int_window+2 %int_window+1:-1:2%size(time)
    counter_l_r=counter_l_r+1;
    
    r_ij=squeeze(rij(frame,:,:));
    dbdt_l_r=r_ij*squeeze(l_r(int_window-counter_l_r+1,1:3,1:3));
    l_r(int_window-counter_l_r,1:3,1:3)=squeeze(l_r(int_window-counter_l_r+1,1:3,1:3))-dt*dbdt_l_r;
end



% Initial alignment with eigenvectors of Cauchy green;
counter_l_cg=0;
l_cg(int_window,1:3,1:3)=[W1(end,:); W2(end,:); W3(end,:)];
for frame=frame_break:-1:frame_break-int_window+2%size(time)%175
    counter_l_cg=counter_l_cg+1;
    r_ij=squeeze(rij(frame,:,:));
    dbdt_l_cg=r_ij*squeeze(l_cg(int_window-counter_l_cg+1,1:3,1:3));
    l_cg(int_window-counter_l_cg,1:3,1:3)=squeeze(l_cg(int_window-counter_l_cg+1,1:3,1:3))-dt*dbdt_l_cg;
end



% figure;
% hold on
% for j=round(linspace(1,int_window,10))%1:10:int_window
%     line_cg=squeeze(l_cg(j,1:3,1:3));
%
%     line([-line_cg(1,1) line_cg(1,1)]+x_(int_window-j+1)*1000, [-line_cg(1,2) line_cg(1,2)]+y_(int_window-j+1)*1000, [-line_cg(1,3) line_cg(1,3)]+z_(int_window-j+1)*1000,'LineWidth',2,'color','r','LineStyle','-')
%     line([-line_cg(2,1) line_cg(2,1)]+x_(int_window-j+1)*1000, [-line_cg(2,2) line_cg(2,2)]+y_(int_window-j+1)*1000, [-line_cg(2,3) line_cg(2,3)]+z_(int_window-j+1)*1000,'LineWidth',2,'color','g','LineStyle','-')
%     line([-line_cg(3,1) line_cg(3,1)]+x_(int_window-j+1)*1000, [-line_cg(3,2) line_cg(3,2)]+y_(int_window-j+1)*1000, [-line_cg(3,3) line_cg(3,3)]+z_(int_window-j+1)*1000,'LineWidth',2,'color','b','LineStyle','-')
% end
% axis equal
% title('Initial alignment with the CG')
%
% figure;
% hold on
% for j=round(linspace(1,int_window,10))%1:10:int_window
%     line_r=squeeze(l_r(j,1:3,1:3));
%
%     line([-line_r(1,1) line_r(1,1)]+x_(int_window-j+1)*1000, [-line_r(1,2) line_r(1,2)]+y_(int_window-j+1)*1000, [-line_r(1,3) line_r(1,3)]+z_(int_window-j+1)*1000,'LineWidth',2,'color','r','LineStyle','--')
%     line([-line_r(2,1) line_r(2,1)]+x_(int_window-j+1)*1000, [-line_r(2,2) line_r(2,2)]+y_(int_window-j+1)*1000, [-line_r(2,3) line_r(2,3)]+z_(int_window-j+1)*1000,'LineWidth',2,'color','g','LineStyle','--')
%     line([-line_r(3,1) line_r(3,1)]+x_(int_window-j+1)*1000, [-line_r(3,2) line_r(3,2)]+y_(int_window-j+1)*1000, [-line_r(3,3) line_r(3,3)]+z_(int_window-j+1)*1000,'LineWidth',2,'color','b','LineStyle','--')
% end
% title('Initial alignment with sep vector')
% axis equal
%
%
% figure;
% hold on
% for j=round(linspace(1,int_window,10))%1:10:int_window
%     line_lambda=squeeze(l_lambda(j,1:3,1:3));
%
%     line([-line_lambda(1,1) line_lambda(1,1)]+x_(int_window-j+1)*1000, [-line_lambda(1,2) line_lambda(1,2)]+y_(int_window-j+1)*1000, [-line_lambda(1,3) line_lambda(1,3)]+z_(int_window-j+1)*1000,'color','r','LineStyle','--')
%     line([-line_lambda(2,1) line_lambda(2,1)]+x_(int_window-j+1)*1000, [-line_lambda(2,2) line_lambda(2,2)]+y_(int_window-j+1)*1000, [-line_lambda(2,3) line_lambda(2,3)]+z_(int_window-j+1)*1000,'color','g','LineStyle','--')
%     line([-line_lambda(3,1) line_lambda(3,1)]+x_(int_window-j+1)*1000, [-line_lambda(3,2) line_lambda(3,2)]+y_(int_window-j+1)*1000, [-line_lambda(3,3) line_lambda(3,3)]+z_(int_window-j+1)*1000,'color','b','LineStyle','--')
% end
% title('Initial alignment with strain eigen frame')
% axis equal
%
%
% figure;
% hold on
% for j=round(linspace(1,int_window,10))%1:10:int_window
%     line_cg=squeeze(l_cg(j,1:3,1:3)); % initial alignment with cg
%
%     line([-line_cg(1,1) line_cg(1,1)]+x_(int_window-j+1)*1000, [-line_cg(1,2) line_cg(1,2)]+y_(int_window-j+1)*1000, [-line_cg(1,3) line_cg(1,3)]+z_(int_window-j+1)*1000,'LineWidth',2,'color','r','LineStyle','-')
%     line([-line_cg(2,1) line_cg(2,1)]+x_(int_window-j+1)*1000, [-line_cg(2,2) line_cg(2,2)]+y_(int_window-j+1)*1000, [-line_cg(2,3) line_cg(2,3)]+z_(int_window-j+1)*1000,'LineWidth',2,'color','g','LineStyle','-')
%     line([-line_cg(3,1) line_cg(3,1)]+x_(int_window-j+1)*1000, [-line_cg(3,2) line_cg(3,2)]+y_(int_window-j+1)*1000, [-line_cg(3,3) line_cg(3,3)]+z_(int_window-j+1)*1000,'LineWidth',2,'color','b','LineStyle','-')
%
%     line_r=squeeze(l_r(j,1:3,1:3)); % initial alignment with separation vector
%
%     line([-line_r(1,1) line_r(1,1)]+x_(int_window-j+1)*1000, [-line_r(1,2) line_r(1,2)]+y_(int_window-j+1)*1000, [-line_r(1,3) line_r(1,3)]+z_(int_window-j+1)*1000,'LineWidth',2,'color','r','LineStyle','--')
%     line([-line_r(2,1) line_r(2,1)]+x_(int_window-j+1)*1000, [-line_r(2,2) line_r(2,2)]+y_(int_window-j+1)*1000, [-line_r(2,3) line_r(2,3)]+z_(int_window-j+1)*1000,'LineWidth',2,'color','g','LineStyle','--')
%     line([-line_r(3,1) line_r(3,1)]+x_(int_window-j+1)*1000, [-line_r(3,2) line_r(3,2)]+y_(int_window-j+1)*1000, [-line_r(3,3) line_r(3,3)]+z_(int_window-j+1)*1000,'LineWidth',2,'color','b','LineStyle','--')
% end
% title('particle coordinates, one is initialy alinged with cauchy green and another with sep vector')
% axis equal
%
%
%
%
% figure(10);
% hold on
% for j=round(linspace(1,int_window,10))%1:10:int_window
%     line_lambda=squeeze(l_lambda(j,1:3,1:3));
%
%     line([-line_lambda(1,1) line_lambda(1,1)]+x_(int_window-j+1)*1000, [-line_lambda(1,2) line_lambda(1,2)]+y_(int_window-j+1)*1000, [-line_lambda(1,3) line_lambda(1,3)]+z_(int_window-j+1)*1000,'LineWidth',1,'color','r','LineStyle','--')
%     line([-line_lambda(2,1) line_lambda(2,1)]+x_(int_window-j+1)*1000, [-line_lambda(2,2) line_lambda(2,2)]+y_(int_window-j+1)*1000, [-line_lambda(2,3) line_lambda(2,3)]+z_(int_window-j+1)*1000,'LineWidth',1,'color','g','LineStyle','--')
%     line([-line_lambda(3,1) line_lambda(3,1)]+x_(int_window-j+1)*1000, [-line_lambda(3,2) line_lambda(3,2)]+y_(int_window-j+1)*1000, [-line_lambda(3,3) line_lambda(3,3)]+z_(int_window-j+1)*1000,'LineWidth',1,'color','b','LineStyle','--')
%
%     line([-Sij_v1(j,1)*Sij_e1(j)/2 Sij_v1(j,1)*Sij_e1(j)/2]+x_(int_window-j+1)*1000,[-Sij_v1(j,2)*Sij_e1(j)/2 Sij_v1(j,2)*Sij_e1(j)/2]+y_(int_window-j+1)*1000,[-Sij_v1(j,3)*Sij_e1(j)/2 Sij_v1(j,3)*Sij_e1(j)/2]+z_(int_window-j+1)*1000,'LineWidth',4,'Color','r')
%     line([-Sij_v2(j,1)*Sij_e2(j)/2 Sij_v2(j,1)*Sij_e2(j)/2]+x_(int_window-j+1)*1000,[-Sij_v2(j,2)*Sij_e2(j)/2 Sij_v2(j,2)*Sij_e2(j)/2]+y_(int_window-j+1)*1000,[-Sij_v2(j,3)*Sij_e2(j)/2 Sij_v2(j,3)*Sij_e2(j)/2]+z_(int_window-j+1)*1000,'LineWidth',4,'Color','g')
%     line([-Sij_v3(j,1)*Sij_e3(j)/2 Sij_v3(j,1)*Sij_e3(j)/2]+x_(int_window-j+1)*1000,[-Sij_v3(j,2)*Sij_e3(j)/2 Sij_v3(j,2)*Sij_e3(j)/2]+y_(int_window-j+1)*1000,[-Sij_v3(j,3)*Sij_e3(j)/2 Sij_v3(j,3)*Sij_e3(j)/2]+z_(int_window-j+1)*1000,'LineWidth',4,'Color','b')
% %
%     %line(([-W1(j,1)*w1(j)/2 W1(j,1)*w1(j)/2]+x_(j))*1000,([-W1(j,2)*w1(j)/2 W1(j,2)*w1(j)/2]+y_(j))*1000,([-W1(j,3)*w1(j)/2 W1(j,3)*w1(j)/2]+z_(j))*1000,'LineWidth',4,'Color','r')
%     %line(([-W2(j,1)*w2(j)/2 W2(j,1)*w2(j)/2]+x_(j))*1000,([-W2(j,2)*w2(j)/2 W2(j,2)*w2(j)/2]+y_(j))*1000,([-W2(j,3)*w2(j)/2 W2(j,3)*w2(j)/2]+z_(j))*1000,'LineWidth',4,'Color','g')
%     %line(([-W3(j,1)*w3(j)/2 W3(j,1)*w3(j)/2]+x_(j))*1000,([-W3(j,2)*w3(j)/2 W3(j,2)*w3(j)/2]+y_(j))*1000,([-W3(j,3)*w3(j)/2 W3(j,3)*w3(j)/2]+z_(j))*1000,'LineWidth',4,'Color','b')
% end
% title('particle coordinate initialy aligned with strain eigen frame vs strain eigen vector')
% axis equal
%
%
for n=1:int_window
    cos_theta_1(n)=abs(dot(Sij_v1(n,:),squeeze(l_r(n,1,:))'));
    cos_theta_2(n)=abs(dot(Sij_v2(n,:),squeeze(l_r(n,2,:))'));
    cos_theta_3(n)=abs(dot(Sij_v3(n,:),squeeze(l_r(n,3,:))'));
end
figure;
plot(cos_theta_1,'r')
hold on
plot(cos_theta_2,'g')
plot(cos_theta_3,'b')

legend('cos(v1,l_{r}_1)','cos(v2,l_{r}_2)','cos(v3,l_{r}_3)')



for n=1:int_window
    cos_gamma_1(n)=abs(dot(W1(int_window-n+2,:),squeeze(l_cg(n,1,:))'));
    cos_gamma_2(n)=abs(dot(W2(int_window-n+2,:),squeeze(l_cg(n,2,:))'));
    cos_gamma_3(n)=abs(dot(W3(int_window-n+2,:),squeeze(l_cg(n,3,:))'));
end
figure;
plot(cos_gamma_1,'r')
hold on
plot(cos_gamma_2,'g')
plot(cos_gamma_3,'b')

legend('cos(W1,l_{cg}_1)','cos(W2,l_{cg}_2)','cos(W3,l_{cg}_3)')



for n=1:int_window
    cos_alpha_1(n)=abs(dot(squeeze(l_cg(n,1,:)),squeeze(l_r(n,1,:))));
    cos_alpha_2(n)=abs(dot(squeeze(l_cg(n,1,:)),squeeze(l_r(n,2,:))));
    cos_alpha_3(n)=abs(dot(squeeze(l_cg(n,1,:)),squeeze(l_r(n,3,:))));
end
figure;plot(cos_alpha_1,'r')
hold on
plot(cos_alpha_2,'g')
plot(cos_alpha_3,'b')

legend('cos(l_{cg}_1,l_{r}_1)','cos(2_{cg}_2,l_{r}_2)','cos(l_{cg}_3,l_{r}_3)')

% figure;hold on;
% plot(w1,'r');
% plot(w2,'g');
% plot(w3,'b');
%
%
% figure;hold on;
% plot(Sij_e1,'r');
% plot(Sij_e2,'g');
% plot(Sij_e3,'b');
% legend('\Lambda_1','\Lambda_2','\Lambda_3')









% END OF rij PART
% Ellipsoid visualization part




for kk=round(linspace(2,size(W1,1),5))%1:10:size(W1,1)%round(linspace(2,size(W1,1),10))
    ii=kk-1;
    R11_cg=W1(kk,1);%l1x;
    R12_cg=W1(kk,2);%l1y;
    R13_cg=W1(kk,3);%l1z;
    
    R21_cg=W2(kk,1);%l2x;
    R22_cg=W2(kk,2);%l2y;
    R23_cg=W2(kk,3);%l2z;
    
    R31_cg=W3(kk,1);%l3x;
    R32_cg=W3(kk,2);%l3y;
    R33_cg=W3(kk,3);%l3z;
    
    % rij Matrix
    R_cg=[R11_cg R12_cg R13_cg; R21_cg R22_cg R23_cg;R31_cg R32_cg R33_cg]';
    % mesh gridpoints of elipsoid
    [x_cg,y_cg,z_cg]=ellipsoid(0,0,0,w1(kk)/2,w2(kk)/2,w3(kk)/2);
    
    % rotate coordinates
    for i=1:size(x_cg,1)
        for j=1:size(x_cg,1)
            dummy_cg = R_cg*[x_cg(i,j) y_cg(i,j) z_cg(i,j)]';
            x_rcg(i,j) = dummy_cg(1)+x_(ii);
            y_rcg(i,j) = dummy_cg(2)+y_(ii);
            z_rcg(i,j) = dummy_cg(3)+z_(ii);
        end
    end
    
    
    
    
    R11_sij=Sij_v1(ii,1);%l1x;
    R12_sij=Sij_v1(ii,2);%l1y;
    R13_sij=Sij_v1(ii,3);%l1z;
    
    R21_sij=Sij_v2(ii,1);%l2x;
    R22_sij=Sij_v2(ii,2);%l2y;
    R23_sij=Sij_v2(ii,3);%l2z;
    
    R31_sij=Sij_v3(ii,1);%l3x;
    R32_sij=Sij_v3(ii,2);%l3y;
    R33_sij=Sij_v3(ii,3);%l3z;
    
    % rij Matrix
    R_sij=[R11_sij R12_sij R13_sij; R21_sij R22_sij R23_sij;R31_sij R32_sij R33_sij]';
    % mesh gridpoints of elipsoid
    [x_s,y_s,z_s]=ellipsoid(0,0,0,Sij_e1(ii)/2*1E-3,Sij_e2(ii)/2*1E-3,Sij_e3(ii)/2*1E-3);  % HOW NORMALIZE???
    
    % rotate coordinates
    for i=1:size(x_s,1)
        for j=1:size(x_s,1)
            dummy_sij = R_sij*[x_s(i,j) y_s(i,j) z_s(i,j)]';
            x_rs(i,j) = dummy_sij(1)+x_(ii);
            y_rs(i,j) = dummy_sij(2)+y_(ii);
            z_rs(i,j) = dummy_sij(3)+z_(ii);
        end
    end
    
    
    R11_l_r=l_r(size(W1,1)-ii,1,1);%l1x;
    R12_l_r=l_r(size(W1,1)-ii,1,2);%l1y;
    R13_l_r=l_r(size(W1,1)-ii,1,3);%l1z;
    
    R21_l_r=l_r(size(W1,1)-ii,2,1);%l2x;
    R22_l_r=l_r(size(W1,1)-ii,2,2);%l2y;
    R23_l_r=l_r(size(W1,1)-ii,2,3);%l2z;
    
    R31_l_r=l_r(size(W1,1)-ii,3,1);%l3x;
    R32_l_r=l_r(size(W1,1)-ii,3,2);%l3y;
    R33_l_r=l_r(size(W1,1)-ii,3,3);%l3z;
    
    % rij Matrix
    R_l_r=[R11_l_r R12_l_r R13_l_r; R21_l_r R22_l_r R23_l_r;R31_l_r R32_l_r R33_l_r]';
    % mesh gridpoints of elipsoid
    [x_l_r,y_l_r,z_l_r]=ellipsoid(0,0,0,0.001,0.0002,0.0002);
    
    
    % rotate coordinates
    for i=1:size(x_l_r,1)
        for j=1:size(x_l_r,1)
            dummy_l_r = R_l_r*[x_l_r(i,j) y_l_r(i,j) z_l_r(i,j)]';
            x_l_r(i,j) = dummy_l_r(1)+x_(ii);
            y_l_r(i,j) = dummy_l_r(2)+y_(ii);
            z_l_r(i,j) = dummy_l_r(3)+z_(ii);
        end
    end
    
    
    ll=squeeze(l_cg(size(W1,1)-ii,1:3,1:3));
    
    figure(100);
    
    hold on
    box on
    line([-ll(1,1) ll(1,1)]/1000+x_(ii), [-ll(1,2) ll(1,2)]/1000+y_(ii), [-ll(1,3) ll(1,3)]/1000+z_(ii),'color','r','LineStyle','--')
    line([-ll(2,1) ll(2,1)]/1000+x_(ii), [-ll(2,2) ll(2,2)]/1000+y_(ii), [-ll(2,3) ll(2,3)]/1000+z_(ii),'color','g','LineStyle','--')
    line([-ll(3,1) ll(3,1)]/1000+x_(ii), [-ll(3,2) ll(3,2)]/1000+y_(ii), [-ll(3,3) ll(3,3)]/1000+z_(ii),'color','b','LineStyle','--')
    surf(x_l_r,y_l_r,z_l_r);
    surf(x_rcg,y_rcg,z_rcg);
    
    %      line([-ll(1,1) ll(1,1)]/1000+x_(ii), [-ll(1,2) ll(1,2)]/1000+y_(ii)+0.005, [-ll(1,3) ll(1,3)]/1000+z_(ii),'color','r','LineStyle','--')
    %      line([-ll(2,1) ll(2,1)]/1000+x_(ii), [-ll(2,2) ll(2,2)]/1000+y_(ii)+0.005, [-ll(2,3) ll(2,3)]/1000+z_(ii),'color','g','LineStyle','--')
    %      line([-ll(3,1) ll(3,1)]/1000+x_(ii), [-ll(3,2) ll(3,2)]/1000+y_(ii)+0.005, [-ll(3,3) ll(3,3)]/1000+z_(ii),'color','b','LineStyle','--')
    %      %surf(x_l_r,y_l_r+0.005,z_l_r);
    %      surf(x_rs,y_rs+0.005,z_rs);
    %
    
    
    %     scatter3(x_(kk),y_(kk),z_(kk),100,'filled')
    %     plot3(x_(kk),y_(kk),z_(kk),'r')
    shading flat
    %  alpha(0.2)
    %     tit=['frame=',num2str(kk+frame_0-1)];
    %     title(tit)
    xlabel('\bf x')
    ylabel('\bf y')
    zlabel('\bf z')
    %hold on;
    %quiver3([0 0 0]',[0 0 0]',[0 0 0]',V(:,1)*w1(kk)/2,V(:,2)*w2(kk)/2,V(:,3)*w3(kk)/2,0)
    %         line([-V(1,1)*w1(kk)/2 V(1,1)*w1(kk)/2],[-V(2,1)*w1(kk)/2 V(2,1)*w1(kk)/2],[-V(3,1)*w1(kk)/2 V(3,1)*w1(kk)/2],'Color','r')
    %         line([-V(1,2)*w2(kk)/2 V(1,2)*w2(kk)/2],[-V(2,2)*w2(kk)/2 V(2,2)*w2(kk)/2],[-V(3,2)*w2(kk)/2 V(3,2)*w2(kk)/2],'Color','g')
    %         line([-V(1,3)*w3(kk)/2 V(1,3)*w3(kk)/2],[-V(2,3)*w3(kk)/2 V(2,3)*w3(kk)/2],[-V(3,3)*w3(kk)/2 V(3,3)*w3(kk)/2],'Color','b')
    %quiver3([0 0 0]',[0 0 0]',[0 0 0]',V(:,1),V(:,2),V(:,3),0)
    axis equal
    %[V D.^0.5]
    
    
    
end




















