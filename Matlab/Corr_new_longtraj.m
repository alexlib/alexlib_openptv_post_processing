% Creating correlation of the xuap data of longest trajectories
%  from reading the data from directory with XUAP files ,creating the
%  trajectories measuring the length after that chosing the trajectories to
%  calculate the correlation and the calculation itself.

trajNum = length(newtraj);
trajectories = zeros(trajNum,1);

for i = 1:trajNum
    trajectories(i)  = length(newtraj(i).uf);
end

%creating arrays for the data needed to process
% xf = cat(1,data.xf);
% yf = cat(1,data.yf);
% zf = cat(1,data.zf);
% uf = cat(1,data.uf);
% vf = cat(1,data.vf);
% wf = cat(1,data.wf);

%choosing  the trajectories to correlate
% trajID = cat(2,data.traj);

%choosing trajectory to work on
% example trajID2Corr=177723;
trajID2Corr = find(trajectories>=5);



% Lagrangian  time correlation 
Traj_t=zeros(max(trajectories),6); % Creating maximum length array
% Traj_t_f=zeros(max(trajectories),6); % correlation containing array
% creating Traj array to process the data into more proper form
for l = 1:length(trajID2Corr)
    traj = trajID2Corr(l);
%     Traj=zeros(length(traj)-2,6);
%     for i=1:(size(Traj,1))
%         Traj(i,1)=uf(traj(i+1));Traj(i,2)=vf(traj(i+1));Traj(i,3)=wf(traj(i+1));
%     end 
    %full processing
    for i=1:length(newtraj(traj).uf)-1 % index of jump
        for j=1:length(newtraj(traj).uf)-i % index that run through the Traj array
            ua = [newtraj(traj).uf(j),newtraj(traj).vf(j),newtraj(traj).wf(j)];
%             ua = [newtraj(traj).uf(1),newtraj(traj).vf(1),newtraj(traj).wf(1)];
            ub = [newtraj(traj).uf(i+j-1),newtraj(traj).vf(i+j-1),newtraj(traj).wf(i+j-1)];
            % inserting the data into the array
            Traj_t(i,1)=norm(ua)*norm(ub)+Traj_t(i,1);  %  ua_norm*ub_norm
            Traj_t(i,2)=ua*ub'+Traj_t(i,2);  % ua_full*ub_full
            Traj_t(i,3)=ua(1)*ub(1)+Traj_t(i,3);  % ua_x*ub_x
            Traj_t(i,4)=ua(2)*ub(2)+Traj_t(i,4);  % ua_y*ub_y
            Traj_t(i,5)=ua(3)*ub(3)+Traj_t(i,5);  % ua_z*ua_z
            Traj_t(i,6)=1+Traj_t(i,6);  % % counting array
        end
    end
end


% averaging by counting array
% for i=1:length(Traj_t)
    Traj_t(:,1)=Traj_t(:,1)./Traj_t(:,6);
    Traj_t(:,2)=Traj_t(:,2)./Traj_t(:,6);
    Traj_t(:,3)=Traj_t(:,3)./Traj_t(:,6);
    Traj_t(:,4)=Traj_t(:,4)./Traj_t(:,6);
    Traj_t(:,5)=Traj_t(:,5)./Traj_t(:,6);
% end
% Plotting the Correlation results of time correlation
% normalization of the Traj_t vector
cor_norm=[Traj_t(:,1)];
cor_full=[Traj_t(:,2)];
cor_x=[ Traj_t(:,3)];
cor_y=[ Traj_t(:,4)];
cor_z=[ Traj_t(:,5)];

del_t=1/150; % here is the frame rate 
x=del_t:del_t:del_t*(length(cor_norm));
plot(x,cor_norm/cor_norm(1),'r','MarkerSize',3,'DisplayName','norm')
xlabel('time [s]')
ylabel('normal correlation')
hold on
plot(x,cor_full/cor_full(1),'g','MarkerSize',3,'DisplayName','full')
plot(x,cor_x/cor_x(1),'b','MarkerSize',3,'DisplayName','x')
plot(x,cor_y/cor_y(1),'y','MarkerSize',3,'DisplayName','y')
plot(x,cor_z/cor_z(1),'p','MarkerSize',3,'DisplayName','z')
hold off

%{
% pdf velocity increments for tau=1,2,4,20,50

hef_y1=[];hef_x1=[];hef_full1=[];
hef_y2=[];hef_x2=[];hef_full2=[];
hef_y4=[];hef_x4=[];hef_full4=[];
hef_y20=[];hef_x20=[];hef_full20=[];
hef_y50=[];hef_x50=[];hef_full50=[];

for l=1:length(trajID2Corr)
    traj = find(trajID == trajID2Corr(l));
    Traj=zeros(length(traj)-2,6);
    for i=1:(size(Traj,1))
        Traj(i,1)=xf(traj(i+1));Traj(i,2)=yf(traj(i+1));Traj(i,3)=zf(traj(i+1));
        Traj(i,4)=uf(traj(i+1));Traj(i,5)=vf(traj(i+1));Traj(i,6)=wf(traj(i+1));
    end
    for i=1:(size(Traj,1)-1)
        hef_full1=[hef_full1;norm([Traj(i+1,4),Traj(i+1,5),Traj(i+1,6)])-norm([Traj(i,4),Traj(i,5),Traj(i,6)])];
        hef_x1=[hef_x1;abs(Traj(i+1,4))-abs(Traj(i,4))];
        hef_y1=[hef_y1;abs(Traj(i+1,5))-abs(Traj(i,5))];
    end
    for i=1:(size(Traj,1)-2)
        hef_full2=[hef_full2;norm([Traj(i+2,4),Traj(i+2,5),Traj(i+2,6)])-norm([Traj(i,4),Traj(i,5),Traj(i,6)])];
        hef_y2=[hef_y2;abs(Traj(i+2,5))-abs(Traj(i,5))];
        hef_x2=[hef_x2;abs(Traj(i+2,4))-abs(Traj(i,4))];
    end
    if size(Traj,1)>4
        for i=1:(size(Traj,1)-4)
            hef_full4=[hef_full4;norm([Traj(i+4,4),Traj(i+4,5),Traj(i+4,6)])-norm([Traj(i,4),Traj(i,5),Traj(i,6)])];
            hef_y4=[hef_y4;abs(Traj(i+4,5))-abs(Traj(i,5))];
            hef_x4=[hef_x4;abs(Traj(i+4,4))-abs(Traj(i,4))];
        end
    end
    if size(Traj,1)>20
        for i=1:(size(Traj,1)-20)
            hef_full20=[hef_full20;norm([Traj(i+20,4),Traj(i+20,5),Traj(i+20,6)])-norm([Traj(i,4),Traj(i,5),Traj(i,6)])];
            hef_x20=[hef_x20;abs(Traj(i+20,4))-abs(Traj(i,4))];
            hef_y20=[hef_y20;abs(Traj(i+20,5))-abs(Traj(i,5))];
        end
    end
    if size(Traj,1)>50
        for i=1:(size(Traj,1)-50)
            hef_full50=[hef_full50;norm([Traj(i+50,4),Traj(i+50,5),Traj(i+50,6)])-norm([Traj(i,4),Traj(i,5),Traj(i,6)])];
            hef_x50=[hef_x50;abs(Traj(i+50,4))-abs(Traj(i,4))];
            hef_y50=[hef_y50;abs(Traj(i+50,5))-abs(Traj(i,5))];
        end
    end
end
% normalizing by rms values
hef_y1=hef_y1./rms2(hef_y1);hef_x1=hef_x1./rms2(hef_x1);hef_full1=hef_full1./rms2(hef_full1);
hef_y2=hef_y2./rms2(hef_y2);hef_x2=hef_x2./rms2(hef_x2);hef_full2=hef_full2./rms2(hef_full2);
hef_y4=hef_y4./rms2(hef_y4);hef_x4=hef_x4./rms2(hef_x4);hef_full4=hef_full4./rms2(hef_full4);
hef_y20=hef_y20./rms2(hef_y20);hef_x20=hef_x20./rms2(hef_x20);hef_full20=hef_full20./rms2(hef_full20);
hef_y50=hef_y50./rms2(hef_y50);hef_x50=hef_x50./rms2(hef_x50);hef_full50=hef_full50./rms2(hef_full50);

[f_f_1,x_f_1] = ksdensity(hef_full1);[f_x_1,x_x_1] = ksdensity(hef_x1);[f_y_1,x_y_1] = ksdensity(hef_y1);
[f_f_2,x_f_2] = ksdensity(hef_full2);[f_x_2,x_x_2] = ksdensity(hef_x2);[f_y_2,x_y_2] = ksdensity(hef_y2);
[f_f_4,x_f_4] = ksdensity(hef_full4);[f_x_4,x_x_4] = ksdensity(hef_x4);[f_y_4,x_y_4] = ksdensity(hef_y4);
[f_f_20,x_f_20] = ksdensity(hef_full20);[f_x_20,x_x_20] = ksdensity(hef_x20);[f_y_20,x_y_20] = ksdensity(hef_y20);
[f_f_50,x_f_50] = ksdensity(hef_full50);[f_x_50,x_x_50] = ksdensity(hef_x50);[f_y_50,x_y_50] = ksdensity(hef_y50);

figure;
plot(x_y_1,f_y_1,'r')
hold on
plot(x_y_2,f_y_2,'b')
plot(x_y_4,f_y_4,'g')
plot(x_y_20,f_y_20,'y')
plot(x_y_50,f_y_50,'p')
hold off


%calculation of structure function
SF2=zeros(max(trajectories),3);
SF4=zeros(max(trajectories),3);
SF6=zeros(max(trajectories),3);
SF_c=zeros(max(trajectories),1); % counting vector
% correlation containing array
% creating Traj array to process the data into more proper form
for l=1:length(trajID2Corr)
    traj = find(trajID == trajID2Corr(l));
    Traj=zeros(length(traj)-2,6);
    for i=1:(size(Traj,1))
        Traj(i,1)=xf(traj(i+1));Traj(i,2)=yf(traj(i+1));Traj(i,3)=zf(traj(i+1));
        Traj(i,4)=uf(traj(i+1));Traj(i,5)=vf(traj(i+1));Traj(i,6)=wf(traj(i+1));
    end
    % processing
    Traj_t2=zeros(max(trajectories),3);
    Traj_t4=zeros(max(trajectories),3);
    Traj_t6=zeros(max(trajectories),3);
    Traj_t_c=zeros(max(trajectories),1);
    for i=1:size(Traj,1) % index of jump
        for j=1:(size(Traj,1)-i) % index that run through the Traj array
            ua=[Traj(j,4),Traj(j,5),Traj(j,6)];% velocity of point a
            ub=[Traj(i+j,4),Traj(i+j,5),Traj(i+j,6)]; % velocity of point b
            % inserting the data into the array
            Traj_t2(i,1)=(norm(ua)-norm(ub))^2+Traj_t2(i,1);  %  ua_norm*ub_norm^2
            Traj_t2(i,2)=(ua(1)-ub(1))^2+Traj_t2(i,2);  % ua_x*ub_x^2
            Traj_t2(i,3)=(ua(2)-ub(2))^2+Traj_t2(i,3);  % ua_y*ub_y^2
            Traj_t4(i,1)=(norm(ua)-norm(ub))^4+Traj_t4(i,1);  %  ua_norm*ub_norm^2
            Traj_t4(i,2)=(ua(1)-ub(1))^4+Traj_t4(i,2);  % ua_x*ub_x^2
            Traj_t4(i,3)=(ua(2)-ub(2))^4+Traj_t4(i,3);  % ua_y*ub_y^2
            Traj_t6(i,1)=(norm(ua)-norm(ub))^6+Traj_t6(i,1);  %  ua_norm*ub_norm^2
            Traj_t6(i,2)=(ua(1)-ub(1))^6+Traj_t6(i,2);  % ua_x*ub_x^2
            Traj_t6(i,3)=(ua(2)-ub(2))^6+Traj_t6(i,3);  % ua_y*ub_y^2
            Traj_t_c(i)=1+Traj_t_c(i);  % % counting array
        end
    end
    % normalization by counting
    for i=1:(size(Traj,1)-1)
        Traj_t2(i,1)= Traj_t2(i,1)./Traj_t_c(i);
        Traj_t2(i,2)=Traj_t2(i,2)./Traj_t_c(i);
        Traj_t2(i,3)=Traj_t2(i,3)./Traj_t_c(i);
        Traj_t4(i,1)= Traj_t4(i,1)./Traj_t_c(i);
        Traj_t4(i,2)=Traj_t4(i,2)./Traj_t_c(i);
        Traj_t4(i,3)=Traj_t4(i,3)./Traj_t_c(i);
        Traj_t6(i,1)= Traj_t6(i,1)./Traj_t_c(i);
        Traj_t6(i,2)=Traj_t6(i,2)./Traj_t_c(i);
        Traj_t6(i,3)=Traj_t6(i,3)./Traj_t_c(i);
    end
    for i=1:(size(Traj,1)-1)
        SF2(i,1)=Traj_t2(i,1)+SF2(i,1);
        SF2(i,2)=Traj_t2(i,2)+SF2(i,2);
        SF2(i,3)=Traj_t2(i,3)+SF2(i,3);
        SF4(i,1)=Traj_t4(i,1)+SF4(i,1);
        SF4(i,2)=Traj_t4(i,2)+SF4(i,2);
        SF4(i,3)=Traj_t4(i,3)+SF4(i,3);
        SF6(i,1)=Traj_t6(i,1)+SF6(i,1);
        SF6(i,2)=Traj_t6(i,2)+SF6(i,2);
        SF6(i,3)=Traj_t6(i,3)+SF6(i,3);
        SF_c(i)=1+SF_c(i);
    end
end

SF_p2_x=SF2(:,2)./SF_c(:); SF_p2_y=SF2(:,3)./SF_c(:); SF_p2_f=SF2(:,1)./SF_c(:);
SF_p4_x=SF4(:,2)./SF_c(:); SF_p4_y=SF4(:,3)./SF_c(:); SF_p4_f=SF4(:,1)./SF_c(:);
SF_p6_x=SF6(:,2)./SF_c(:); SF_p6_y=SF6(:,3)./SF_c(:); SF_p6_f=SF6(:,1)./SF_c(:);

%plotting the results
figure;
plot(SF_p2_f,'r')
hold on
plot(SF_p4_x,'y')
plot(SF_p6_y,'g')
hold off

%eps
eps_4_f=gradient(log(SF_p4_f),abs(log(SF_p2_f)));
eps_4_x=gradient(log(SF_p4_x),log(SF_p2_x));
eps_4_y=gradient(log(SF_p4_y),log(SF_p2_y));
eps_6_f=gradient(log(SF_p6_f),log(SF_p2_f));
eps_6_x=gradient(log(SF_p6_x),log(SF_p2_x));
eps_6_y=gradient(log(SF_p6_y),log(SF_p2_y));

figure;
plot(eps_4_f,'r')
hold on
plot(eps_4_x,'y')
plot(eps_4_y,'g')
hold off

figure;
plot(eps_6_f,'r')
hold on
plot(eps_6_x,'y')
plot(eps_6_y,'g')
hold off

%}
