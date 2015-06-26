clear all
clear all
%close all
clc
tol=0.800e-3;
radius_sphere=0.006;
offset=12;
long_traj_count=0;
max_length=10;
tau_eta=0.14;
frame_rate=250;
factor=tau_eta*frame_rate;
load map_turbulent_box

already_there=1;
%for WF=1:102;
WF=6;
    fold='D:\Turbulent_data_backup\flow_induced_aggregates_breakage_at_100rpm\WF';

    ff=[];
    d = [fold,num2str(WF),'\res_agg','/rt_is_v2*'];
    dd=dir(d);
    for k=1:size(dd,1)
        f=str2double(dd(k).name(10:end));
        ff=[ff; f];
    end
    first=min(ff);
    last=max(ff);
    clear ff

      %first=1;
      %last=169;
    
    fold2=[fold,num2str(WF),'\'];
    

   
    agg=repmat(zeros,[last-first+1,200,19]);

    name_root=[fold,num2str(WF),'\res_agg'];
    
    
    if already_there==1
        load([fold,num2str(WF),'/variables'])
    else
        for i=first:last
            ['act1: Reading agg info from WF = ',num2str(WF),' from ',num2str(i),' of ',num2str(last)]
            name=[name_root,'\rt_is_v2.',num2Str(i)];
            fid = fopen(name, 'r');
            num_points = fscanf(fid, '%i', [1 1]);  % It has two rows now
            tmp = fscanf(fid, '%i %f %f %f %i %i %i %i %f %f %f %f ', [12 num_points]);
            A=tmp';
            A(:,2:4)=A(:,2:4)*0.001;
            agg(i-first+1,1:num_points,1:7)=A(:,[2 3 4 9 10 11 12]);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            agg(i-first+1,1:num_points,8)=i;
            [mean_shortaxis max_shortaxis mean_longaxis max_longaxis orientation]=ellipse_measure(i,fold2);
            agg(i-first+1,1:num_points,10)=mean_shortaxis;
            agg(i-first+1,1:num_points,11)=max_shortaxis;
            agg(i-first+1,1:num_points,12)=mean_longaxis;
            agg(i-first+1,1:num_points,13)=max_longaxis;
            agg(i-first+1,1:num_points,14)=orientation;
            agg(i-first+1,1:num_points,15)=interp3(m_x,m_y,m_z,m_px,agg(i-first+1,1:num_points,1)*1000,agg(i-first+1,1:num_points,2)*1000,agg(i-first+1,1:num_points,3)*1000);
            agg(i-first+1,1:num_points,16)=agg(i-first+1,1:num_points,10)./agg(i-first+1,1:num_points,15);%mean_shortaxis in mm
            agg(i-first+1,1:num_points,17)=agg(i-first+1,1:num_points,11)./agg(i-first+1,1:num_points,15);%max_shortaxis in mm
            agg(i-first+1,1:num_points,18)=agg(i-first+1,1:num_points,12)./agg(i-first+1,1:num_points,15);%mean_longaxis in mm
            agg(i-first+1,1:num_points,19)=agg(i-first+1,1:num_points,13)./agg(i-first+1,1:num_points,15);%max_longaxis in mm
            
            
            fclose(fid);
            
        end
        
        wd= 'E:\OLD-Data-DISK\Matlab for PTV';
        cd([fold,num2str(WF)]);
        save  variables agg
        cd(wd) 

    end
    %clear agg first last d num_points tmp A
%end
  
    xa=[];
    ya=[];
    za=[];
    ta=[];
    for i=1:last-first+1
        
        xa=[xa;agg(i,:,1)'];
        ya=[ya;agg(i,:,2)'];
        za=[za;agg(i,:,3)'];
        t=repmat(i,length(agg(i,:,1)),1);
        ta=[ta;t];
        
    end
    
    %track aggregates using nearest neighbour, 'tol' tolerance
    tid=0;
    for i=1:last-first
       ['act2: linking agg in time =', num2str(i),' of ', num2str(last-first+1)]
        ind=find(abs(agg(i,:,1))>0 & abs(agg(i,:,2))>0);
        %ind=find(agg(i,:,4)>0);
        
        ind2=find(agg(i,ind,9)==0);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ind=ind(ind2);
        num_points=length(ind);
        for j=1:num_points
            tid=tid+1;%trajectory id
            pid=1;
            traj(tid,pid,1)=i;%frame id
            traj(tid,pid,2)=ind(j);%pos in frame
            agg(traj(tid,pid,1),traj(tid,pid,2),9)=1;%occupied as of now
            %find next ?
            end_of_traj=0;
            while end_of_traj==0
                pos_agg=[agg(traj(tid,pid,1),traj(tid,pid,2),1) agg(traj(tid,pid,1),traj(tid,pid,2),2) agg(traj(tid,pid,1),traj(tid,pid,2),3)];
                ind2=find(abs(agg(traj(tid,pid,1)+1,:,1))>0 & abs(agg(traj(tid,pid,1)+1,:,2))>0 & agg(traj(tid,pid,1)+1,:,9)<1);
                %ind2=find(agg(traj(tid,pid,1)+1,:,4)>0  & agg(traj(tid,pid,1)+1,:,9)<1);
                pos=repmat(pos_agg,length(ind2),1);
                if length(ind2>0)%%%%%
                    if length(ind2)>1
                        delta=squeeze(agg(traj(tid,pid,1)+1,ind2,1:3))-pos;
                        dist=sum(delta.^2,2).^0.5;
                    else
                        delta=squeeze(agg(traj(tid,pid,1)+1,ind2,1:3))-pos';
                        dist=sum(delta.^2,1).^0.5;
                    end
                    ind3=find(dist<tol);
                    if length(ind3)>0
                        ind2=ind2(ind3);
                    else
                        ind2=[];
                    end
                else
                    ind2=[];
                end
                
                if length(ind2)>0 & traj(tid,pid,1)<last-first 
                    pid=pid+1;
                    traj(tid,pid,1)=traj(tid,pid-1,1)+1;
                    traj(tid,pid,2)=ind2(1);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    agg(traj(tid,pid,1),traj(tid,pid,2),9)=1;%occupied as of now;
                else
                    end_of_traj=1;
                    
                    if pid>10%20 %% making this small can also underestimate the 
                        
                        long_traj_count=long_traj_count+1;
                        
                        traj_len(tid)=pid;
                        
                        if pid>max_length 
                            max_length=pid;
                        end
                        
                        for k=1:pid
                            tx(tid,k)=agg(traj(tid,k,1),traj(tid,k,2),1);
                            ty(tid,k)=agg(traj(tid,k,1),traj(tid,k,2),2);
                            tz(tid,k)=agg(traj(tid,k,1),traj(tid,k,2),3);
                            totalpix(tid,k)=agg(traj(tid,k,1),traj(tid,k,2),4);
                            xpix(tid,k)=agg(traj(tid,k,1),traj(tid,k,2),5);
                            ypix(tid,k)=agg(traj(tid,k,1),traj(tid,k,2),6);
                            grv(tid,k)=agg(traj(tid,k,1),traj(tid,k,2),7);
                            time(tid,k)=agg(traj(tid,k,1),traj(tid,k,2),8);
                            
                            mn_shrt_ax(tid,k)=agg(traj(tid,k,1),traj(tid,k,2),10);
                            max_shrt_ax(tid,k)=agg(traj(tid,k,1),traj(tid,k,2),11);
                            mn_lng_ax(tid,k)=agg(traj(tid,k,1),traj(tid,k,2),12);
                            max_lng_ax(tid,k)=agg(traj(tid,k,1),traj(tid,k,2),13);
                            mn_shrt_ax_mm(tid,k)=agg(traj(tid,k,1),traj(tid,k,2),16);
                            max_shrt_ax_mm(tid,k)=agg(traj(tid,k,1),traj(tid,k,2),17);
                            mn_lng_ax_mm(tid,k)=agg(traj(tid,k,1),traj(tid,k,2),18);
                            max_lng_ax_mm(tid,k)=agg(traj(tid,k,1),traj(tid,k,2),19);
                            
                           
            
                        end
                        
                    else
                        
                        tid=tid-1;
                    end
                end
            end
        end
    end
    long_traj_count=long_traj_count




name_root=[fold,num2str(WF),'/res/'];

for i=first:last
    if mod(i-first,100)==0
        ['act3: reading XUAP',' = ', num2str(i), ' of ', num2str(last)]
    end
    name=[name_root,'xuap.',num2Str(i)];
    clear A;
    A=load(name);
    si=size(A);
    num_points=si(1,1);
    xu(i-first+1,1:num_points,1:6)=A(:,6:11);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% [x y z vx vy vz]
    xu(i-first+1,1:num_points,7)=A(:,15);%%%%%%%%%%%%%%%%%%%%%% [whether cubic spline works or not]
end

x=zeros(tid,max_length,100);
y=zeros(tid,max_length,100);
z=zeros(tid,max_length,100);
Aij=zeros(tid,max_length,3,3);
co=zeros(tid,max_length,1);
A=zeros(tid,max_length,3,3);
%get Aij along all aggregate trajectories
for i=1:tid
    ['act4: Calculating grad around traj = ',num2str(i),' of ',num2str(tid)]
    for j=1:traj_len(i)
        %find swarm
        ind=find(xu(traj(i,j,1),:,7)==1);
        pos_agg=[agg(traj(i,j,1),traj(i,j,2),1) agg(traj(i,j,1),traj(i,j,2),2) agg(traj(i,j,1),traj(i,j,2),3)];
        pos=repmat(pos_agg,length(ind),1);
        delta=squeeze(xu(traj(i,j,1),ind,1:3))-pos;
        dist=sum(delta.^2,2).^0.5;
        ind2=find(dist<radius_sphere);
        ind=ind(ind2);
        co(i,j)=length(ind); % this tells me that how many points there are around each point on a trajectory

        if length(ind)>=4
            %get Aij around this point
            Gij=get_G_ij_linear(xu(traj(i,j,1),ind,1),xu(traj(i,j,1),ind,2),xu(traj(i,j,1),ind,3),xu(traj(i,j,1),ind,4),xu(traj(i,j,1),ind,5),xu(traj(i,j,1),ind,6));
            rel_div_1(i,j)=Gij(4,1);
            
            Aij(i,j,1:3,1:3)=Gij(1:3,1:3);
        else
            Aij(i,j,1:3,1:3)=zeros(3,3);
        end
    end
    %compute omega, sij and eigen-stuff
    for j=1:traj_len(i)
        w11(i,j)=0;
        w12(i,j)=0.5*(Aij(i,j,1,2)-Aij(i,j,2,1));
        w13(i,j)=0.5*(Aij(i,j,1,3)-Aij(i,j,3,1));%G(1,3)-G(3,1);%
        
        w21(i,j)=0.5*(Aij(i,j,2,1)-Aij(i,j,1,2));%G(2,1)-G(1,2);%
        w22(i,j)=0;
        w23(i,j)=0.5*(Aij(i,j,2,3)-Aij(i,j,3,2));
       
        w31(i,j)=0.5*(Aij(i,j,3,1)-Aij(i,j,1,3));
        w32(i,j)=0.5*(Aij(i,j,3,2)-Aij(i,j,2,3));%G(3,2)-G(2,3);%
        w33(i,j)=0;
        
        rot_M(i,j,1:3,1:3)=[w11(i,j)  w12(i,j) w13(i,j); w21(i,j) w22(i,j) w23(i,j); w31(i,j)  w32(i,j) w33(i,j)];
        
        w1(i,j)=2*w23(i,j);
        w2(i,j)=2*w31(i,j);
        w3(i,j)=2*w12(i,j);
        
        enstro(i,j)=w1(i,j)^2+w2(i,j)^2+w3(i,j)^2;
        
        % DEVIATORIC STRAIN COMPONENT  AS Aij IS TRACE FREE
        s11(i,j)=0.5*(Aij(i,j,1,1)+Aij(i,j,1,1));%(G(1,1)+G(1,1));%
        s12(i,j)=0.5*(Aij(i,j,1,2)+Aij(i,j,2,1));%(G(1,2)+G(2,1));%
        s13(i,j)=0.5*(Aij(i,j,1,3)+Aij(i,j,3,1));%(G(1,3)+G(3,1));%
        s22(i,j)=0.5*(Aij(i,j,2,2)+Aij(i,j,2,2));%(G(2,2)+G(2,2));%
        s23(i,j)=0.5*(Aij(i,j,2,3)+Aij(i,j,3,2));%(G(2,3)+G(3,2));%
        s33(i,j)=0.5*(Aij(i,j,3,3)+Aij(i,j,3,3));%(G(3,3)+G(3,3));%

        %A=[s11(i,j) s12(i,j) s13(i,j);s12(i,j) s22(i,j) s23(i,j);s13(i,j) s23(i,j) s33(i,j)];
        %[V,D] = eig(A);
        A(i,j,1:3,1:3)=[s11(i,j) s12(i,j) s13(i,j);s12(i,j) s22(i,j) s23(i,j);s13(i,j) s23(i,j) s33(i,j)];
        [V,D] = eig(squeeze(A(i,j,:,:)));
        [D,k] = sort(diag(D)); 	% ascending order
        D = diag(D(end:-1:1)); 			% descending order
        V = V(:,k(end:-1:1)); 	% the same order for eigenvectors
        v1(i,j,1:3)=V(:,1)';%eigen direction
        v2(i,j,1:3)=V(:,2)';
        v3(i,j,1:3)=V(:,3)';
        e1(i,j)=D(1,1);%eigen value
        e2(i,j)=D(2,2);
        e3(i,j)=D(3,3);
        strain(i,j)=e1(i,j)^2+e2(i,j)^2+e3(i,j)^2;
   end
    %remove crap
%     w1(i,:)=remove_outlayer(w1(i,:),traj_len(i));
%     w2(i,:)=remove_outlayer(w2(i,:),traj_len(i));
%     w3(i,:)=remove_outlayer(w3(i,:),traj_len(i));
%     
%     w11(i,:)=remove_outlayer(w11(i,:),traj_len(i));
%     w12(i,:)=remove_outlayer(w12(i,:),traj_len(i));
%     w13(i,:)=remove_outlayer(w13(i,:),traj_len(i));
%     w21(i,:)=remove_outlayer(w21(i,:),traj_len(i));
%     w22(i,:)=remove_outlayer(w22(i,:),traj_len(i));
%     w23(i,:)=remove_outlayer(w23(i,:),traj_len(i));
%     w31(i,:)=remove_outlayer(w31(i,:),traj_len(i));
%     w32(i,:)=remove_outlayer(w32(i,:),traj_len(i));
%     w33(i,:)=remove_outlayer(w33(i,:),traj_len(i));
%     
%     s11(i,:)=remove_outlayer(s11(i,:),traj_len(i));
%     s12(i,:)=remove_outlayer(s12(i,:),traj_len(i));
%     s13(i,:)=remove_outlayer(s13(i,:),traj_len(i));
%     s22(i,:)=remove_outlayer(s22(i,:),traj_len(i));
%     s23(i,:)=remove_outlayer(s23(i,:),traj_len(i));
%     s33(i,:)=remove_outlayer(s33(i,:),traj_len(i));
    
    %low pass filter
    width=41; % Beat =41
    order=2;
    w1(i,:) = smooth(squeeze(w1(i,:)),width,'sgolay',order);
    w2(i,:) = smooth(squeeze(w2(i,:)),width,'sgolay',order);
    w3(i,:) = smooth(squeeze(w3(i,:)),width,'sgolay',order);
    
    w11(i,:)=smooth(squeeze(w11(i,:)),width,'sgolay',order);
    w12(i,:)=smooth(squeeze(w12(i,:)),width,'sgolay',order);
    w13(i,:)=smooth(squeeze(w13(i,:)),width,'sgolay',order);
    w21(i,:)=smooth(squeeze(w21(i,:)),width,'sgolay',order);
    w22(i,:)=smooth(squeeze(w22(i,:)),width,'sgolay',order);
    w23(i,:)=smooth(squeeze(w23(i,:)),width,'sgolay',order);
    w31(i,:)=smooth(squeeze(w31(i,:)),width,'sgolay',order);
    w32(i,:)=smooth(squeeze(w32(i,:)),width,'sgolay',order);
    w33(i,:)=smooth(squeeze(w33(i,:)),width,'sgolay',order);
    
    s11(i,:) = smooth(squeeze(s11(i,:)),width,'sgolay',order);
    s12(i,:) = smooth(squeeze(s12(i,:)),width,'sgolay',order);
    s13(i,:) = smooth(squeeze(s13(i,:)),width,'sgolay',order);
    s22(i,:) = smooth(squeeze(s22(i,:)),width,'sgolay',order);
    s23(i,:) = smooth(squeeze(s23(i,:)),width,'sgolay',order);
    s33(i,:) = smooth(squeeze(s33(i,:)),width,'sgolay',order);
    
    

    for j=1:traj_len(i)
        
        rot_M(i,j,1:3,1:3)=[w11(i,j)  w12(i,j) w13(i,j); w21(i,j) w22(i,j) w23(i,j); w31(i,j)  w32(i,j) w33(i,j)];
        enstro(i,j)=w1(i,j)^2+w2(i,j)^2+w3(i,j)^2;
        %A=[s11(i,j) s12(i,j) s13(i,j);s12(i,j) s22(i,j) s23(i,j);s13(i,j) s23(i,j) s33(i,j)];
        %[V,D] = eig(A);
        
        A(i,j,1:3,1:3)=[s11(i,j) s12(i,j) s13(i,j);s12(i,j) s22(i,j) s23(i,j);s13(i,j) s23(i,j) s33(i,j)];
        [V,D] = eig(squeeze(A(i,j,:,:)));
        
        [D,k] = sort(diag(D)); 	% ascending order
        D = diag(D(end:-1:1));  % descending order
        V = V(:,k(end:-1:1)); 	% the same order for eigenvectors
        v1(i,j,1:3)=V(:,1)';
        v2(i,j,1:3)=V(:,2)';
        v3(i,j,1:3)=V(:,3)';
        e1(i,j)=D(1,1);
        e2(i,j)=D(2,2);
        e3(i,j)=D(3,3);
        strain(i,j)=e1(i,j)^2+e2(i,j)^2+e3(i,j)^2;
        % THIS IS DIVERGENCE FREE
        rel_div(i,j)=abs(s11(i,j)+s22(i,j)+s33(i,j))/(abs(s11(i,j))+abs(s22(i,j))+abs(s33(i,j)));
    end
end


wd= 'E:\OLD-Data-DISK\Matlab for PTV';
cd([fold,num2str(WF)]);
save  variables WF wd fold fold2 xa ya za ta tx ty tz...
                totalpix xpix ypix grv time tol tid traj traj_len co...
                Aij A rot_M w1 w2 w3 v1 v2 v3 e1 e2 e3 strain enstro...
                mn_shrt_ax    max_shrt_ax   mn_lng_ax   max_lng_ax...
                mn_shrt_ax_mm   max_shrt_ax_mm  mn_lng_ax_mm  max_lng_ax_mm...
                rel_div  radius_sphere  -append 

            %cd(wd)

% agg traj linked in time


figure;
for traj_id=1:tid
%for traj_id=[parent child]
    t=time(traj_id,:);
    t=t(t>0);
    x=tx(traj_id,:);
    x=x(1:numel(t));%%% This is a bug remover !!
    y=ty(traj_id,:);
    y=y(1:numel(t));
    z=tz(traj_id,:);
    z=z(1:numel(t));
    
    hold on
    scatter3(x,y,z,'k.')
    plot3(x,y,z,'k')
    scatter3(x(1),y(1),z(1),100,'g','filled')
    text(x(1),y(1),z(1),[' T (',num2str(traj_id),')', '   ', 'f',num2str(t(1))]);%num2str(traj(traj_id,1))
    scatter3(x(end),y(end),z(end),100,'r','filled')
    text(x(end),y(end),z(end),[' f',num2str(t(end))]);
    tit=['tol=',num2str(tol*1000),'mm','    ', 'Num\_traj=',num2str(tid)];
    title(tit)
    hold off
end

parent=[1  ];%parent traj ID
child=[3]; % children traj id

save variables parent child -append


X=[];
Y=[];
Z=[];
T=[];
parent_totalpix=[];
parent_xpix=[];
parent_ypix=[];
parent_grv=[];
parent_co=[];
parent_Aij=[];
parent_Sij=[];
parent_rij=[];
parent_w1=[];
parent_w2=[];
parent_w3=[];
parent_e1=[];
parent_e2=[];
parent_e3=[];
parent_v1=[];
parent_v2=[];
parent_v3=[];
parent_strain=[];
parent_enstro=[];
parent_mn_shrt_ax=[];
parent_mn_shrt_ax_mm=[];
parent_max_shrt_ax=[];
parent_max_shrt_ax_mm=[];
parent_mn_lng_ax=[];
parent_mn_lng_ax_mm=[];
parent_max_lng_ax=[];
parent_max_lng_ax_mm=[];
parent_rel_div=[];


for traj_id=1:length(parent)%%% If in case parent trajectory is fragmented


    t=time(parent(traj_id),:);
    index=find(t>0);
    T=[T;t(index)'];
    x=tx(parent(traj_id),:);
    X=[X;x(index)'];%%% This is a bug remover !!
    y=ty(parent(traj_id),:);
    Y=[Y;y(index)'];
    z=tz(parent(traj_id),:);
    Z=[Z;z(index)'];
    
    ptotalpix=totalpix(parent(traj_id),:);
    parent_totalpix=[parent_totalpix;ptotalpix(index)'];
    
    pxpix=xpix(parent(traj_id),:);
    parent_xpix=[parent_xpix;pxpix(index)'];
    
    pypix=ypix(parent(traj_id),:);
    parent_ypix=[parent_ypix;pypix(index)'];
    
    pgrv=grv(parent(traj_id),:);
    parent_grv=[parent_grv;pgrv(index)'];
    
    pco=co(parent(traj_id),:);
    parent_co=[parent_co;pco(index)'];
    
    pAij=squeeze(Aij(parent(traj_id),:,:,:));
    parent_Aij=[parent_Aij;pAij(index,:,:)];
    
    rotation=squeeze(rot_M(parent(traj_id),:,:,:));
    parent_rij=[parent_rij;rotation(index,:,:)];
    
    strain_tensor=squeeze(A(parent(traj_id),:,:,:));
    parent_Sij=[parent_Sij;strain_tensor(index,:,:)];
    
    pw1=w1(parent(traj_id),:);
    parent_w1=[parent_w1;pw1(index)'];
    
    pw2=w2(parent(traj_id),:);
    parent_w2=[parent_w2;pw2(index)'];
    
    pw3=w3(parent(traj_id),:);
    parent_w3=[parent_w3;pw3(index)'];
    
    pe1=e1(parent(traj_id),:);
    parent_e1=[parent_e1;pe1(index)'];
    
    pe2=e2(parent(traj_id),:);
    parent_e2=[parent_e2;pe2(index)'];
    
    pe3=e3(parent(traj_id),:);
    parent_e3=[parent_e3;pe3(index)'];
    
    pv1=squeeze(v1(parent(traj_id),:,:));
    parent_v1=[parent_v1;pv1(index,:)];
    
    pv2=squeeze(v2(parent(traj_id),:,:));
    parent_v2=[parent_v2;pv2(index,:)];
    
    pv3=squeeze(v3(parent(traj_id),:,:));
    parent_v3=[parent_v3;pv3(index,:)];
    
    pstrain=strain(parent(traj_id),:);
    parent_strain=[parent_strain;pstrain(:,index)'];
    
    penstro=enstro(parent(traj_id),:);
    parent_enstro=[parent_enstro;penstro(:,index)'];
    
    pmn_shrt_ax=mn_shrt_ax(parent(traj_id),:);
    parent_mn_shrt_ax=[parent_mn_shrt_ax;pmn_shrt_ax(:,index)'];
    
    pmax_shrt_ax=max_shrt_ax(parent(traj_id),:);
    parent_max_shrt_ax=[parent_max_shrt_ax;pmax_shrt_ax(:,index)'];
    
    pmn_shrt_ax_mm=mn_shrt_ax_mm(parent(traj_id),:);
    parent_mn_shrt_ax_mm=[parent_mn_shrt_ax_mm;pmn_shrt_ax_mm(:,index)'];
    
    pmax_shrt_ax_mm=max_shrt_ax_mm(parent(traj_id),:);
    parent_max_shrt_ax_mm=[parent_max_shrt_ax_mm;pmax_shrt_ax_mm(:,index)'];

    pmn_lng_ax=mn_lng_ax(parent(traj_id),:);
    parent_mn_lng_ax=[parent_mn_lng_ax;pmn_lng_ax(:,index)'];
    
    pmax_lng_ax=max_lng_ax(parent(traj_id),:);
    parent_max_lng_ax=[parent_max_lng_ax;pmax_lng_ax(:,index)'];
    
    pmn_lng_ax_mm=mn_lng_ax_mm(parent(traj_id),:);
    parent_mn_lng_ax_mm=[parent_mn_lng_ax_mm;pmn_lng_ax_mm(:,index)'];
    
    pmax_lng_ax_mm=max_lng_ax_mm(parent(traj_id),:);
    parent_max_lng_ax_mm=[parent_max_lng_ax_mm;pmax_lng_ax_mm(:,index)'];

    prel_div=rel_div(parent(traj_id),:);
    parent_rel_div=[parent_rel_div;prel_div(:,index)'];

end




% % % Breakage traj variables
% % % parent traj
% px=tx(par,:);py=ty(par,:);pz=tz(par,:);
% 
% ptotalpix=totalpix(par,:);pxpix=xpix(par,:);pypix=ypix(par,:);pgrv=grv(par,:);
% 
% ptime=time(par,:);ptraj=traj(par,:,:);
% 
% ptraj_len=traj_len(par);
% 
% 
% pco=co(par,:);
% pAij=Aij(par,:,:,:);
% pw1=w1(par,:);pw2=w2(par,:);pw3=w3(par,:);pe1=e1(par,:);pe2=e2(par,:);pe3=e3(par,:);pv1=v1(par,:,:);pv2=v2(par,:,:);pv3=v3(par,:,:);
% pstrain=strain(par,:); penstro=enstro(par,:);
% 
% pmn_shrt_ax=mn_shrt_ax(par,:);pmax_shrt_ax=max_shrt_ax(par,:);
% pmn_shrt_ax_mm=mn_shrt_ax_mm(par,:);pmax_shrt_ax_mm=max_shrt_ax_mm(par,:);
% 
% 
% pmn_lng_ax=mn_lng_ax(par,:);pmax_lng_ax=max_lng_ax(par,:);
% pmn_lng_ax_mm=mn_lng_ax_mm(par,:);pmax_lng_ax_mm=max_lng_ax_mm(par,:);
% 
% prel_div=rel_div(par,:); 



% In total 32 variables are stored

save parent_variables  X Y Z T...
     parent_totalpix  parent_xpix  parent_ypix  parent_grv  ...
     parent_co parent_Aij parent_Sij parent_rij parent_w1 parent_w2 parent_w3 parent_v1 parent_v2 parent_v3...
     parent_e1 parent_e2 parent_e3 parent_strain parent_enstro...
     parent_mn_shrt_ax    parent_max_shrt_ax   parent_mn_lng_ax   parent_max_lng_ax...
     parent_mn_shrt_ax_mm   parent_max_shrt_ax_mm  parent_mn_lng_ax_mm...
     parent_max_lng_ax_mm  parent_rel_div
 
%%% Children variables

for c=1:size(child,2)
    
    frag_time=time(child(c),:);
    indx=find(frag_time>0);
    frag_time=frag_time(:,indx);
    ctime{c}=frag_time;
    
    
    frag_x=tx(child(c),:);
    frag_x=frag_x(:,indx);
    cx{c}=frag_x;

    frag_y=ty(child(c),:);
    frag_y=frag_y(:,indx);
    cy{c}=frag_y;
    
    frag_z=tz(child(c),:);
    frag_z=frag_z(:,indx);
    cz{c}=frag_z;
    
    frag_totalpix=totalpix(child(c),:);
    frag_totalpix=frag_totalpix(:,indx);
    ctotalpix{c}=frag_totalpix;
    
    frag_xpix=xpix(child(c),:);
    frag_xpix=frag_xpix(:,indx);
    cxpix{c}=frag_xpix;

    frag_ypix=ypix(child(c),:);
    frag_ypix=frag_ypix(:,indx);
    cypix{c}=frag_ypix;
    
    frag_grv=grv(child(c),:);
    frag_grv=frag_grv(:,indx);
    cgrv{c}=frag_grv;
    
   
    
    frag_co=co(child(c),:);
    frag_co=frag_co(:,indx);
    cco{c}=frag_co;

    frag_Aij=squeeze(Aij(child(c),:,:,:));
    cAij{c}=frag_Aij(indx,:,:);
    
    frag_Sij=squeeze(A(child(c),:,:,:));
    cSij{c}=frag_Sij(indx,:,:);
    
    frag_rij=squeeze(rot_M(child(c),:,:,:));
    crij{c}=frag_rij(indx,:,:);
    
    frag_w1=w1(child(c),:);
    frag_w1=frag_w1(:,indx);
    cw1{c}=frag_w1;
    
    frag_w2=w2(child(c),:);
    frag_w2=frag_w2(:,indx);
    cw2{c}=frag_w2;
    
    frag_w3=w3(child(c),:);
    frag_w3=frag_w3(:,indx);
    cw3{c}=frag_w3;
    
    frag_e1=e1(child(c),:);
    frag_e1=frag_e1(:,indx);
    ce1{c}=frag_e1;
    
    frag_e2=e2(child(c),:);
    frag_e2=frag_e2(:,indx);
    ce2{c}=frag_e2;
    
    frag_e3=e3(child(c),:);
    frag_e3=frag_e3(:,indx);
    ce3{c}=frag_e3;
    
    frag_v1=v1(child(c),:);
    frag_v1=frag_v1(:,indx);
    cv1{c}=frag_v1;
    
    frag_v2=v2(child(c),:);
    frag_v2=frag_v2(:,indx);
    cv2{c}=frag_v2;
    
    frag_v3=v3(child(c),:);
    frag_v3=frag_v3(:,indx);
    cv3{c}=frag_v3;
    
    
    frag_strain=strain(child(c),:);
    frag_strain=frag_strain(:,indx);
    cstrain{c}=frag_strain;
    
    frag_enstro=enstro(child(c),:);
    frag_enstro=frag_enstro(:,indx);
    censtro{c}=frag_enstro;
    
    frag_mn_shrt_ax=mn_shrt_ax(child(c),:);
    frag_mn_shrt_ax=frag_mn_shrt_ax(:,indx);
    cmn_shrt_ax{c}=frag_mn_shrt_ax;
    
    frag_mn_shrt_ax_mm=mn_shrt_ax_mm(child(c),:);
    frag_mn_shrt_ax_mm=frag_mn_shrt_ax_mm(:,indx);
    cmn_shrt_ax_mm{c}=frag_mn_shrt_ax_mm;
    
    frag_max_shrt_ax=max_shrt_ax(child(c),:);
    frag_max_shrt_ax=frag_max_shrt_ax(:,indx);
    cmax_shrt_ax{c}=frag_max_shrt_ax;
    
    frag_max_shrt_ax_mm=max_shrt_ax_mm(child(c),:);
    frag_max_shrt_ax_mm=frag_max_shrt_ax_mm(:,indx);
    cmax_shrt_ax_mm{c}=frag_max_shrt_ax_mm;
    
    
    frag_mn_lng_ax=mn_lng_ax(child(c),:);
    frag_mn_lng_ax=frag_mn_lng_ax(:,indx);
    cmn_lng_ax{c}=frag_mn_lng_ax;
    
    frag_mn_lng_ax_mm=mn_lng_ax_mm(child(c),:);
    frag_mn_lng_ax_mm=frag_mn_lng_ax_mm(:,indx);
    cmn_lng_ax_mm{c}=frag_mn_lng_ax_mm;
    

    frag_max_lng_ax=max_lng_ax(child(c),:);
    frag_max_lng_ax=frag_max_lng_ax(:,indx);
    cmax_lng_ax{c}=frag_max_lng_ax;
    
    frag_max_lng_ax_mm=max_lng_ax_mm(child(c),:);
    frag_max_lng_ax_mm=frag_max_lng_ax_mm(:,indx);
    cmax_lng_ax_mm{c}=frag_max_lng_ax_mm;
    
    frag_rel_div=rel_div(child(c),:);
    frag_rel_div=frag_rel_div(:,indx);
    crel_div{c}=frag_rel_div;



clear indx frag_x frag_y frag_z frag_totalpix frag_xpix frag_ypix frag_grv...
      frag_traj_len   frag_time frag_co frag_Aij frag_Sij frag_rij frag_w1 frag_w2 frag_w3...
      frag_e1 frag_e2 frag_e3 frag_v1 frag_v2 frag_v3 frag_strain ...
      frag_enstro frag_mn_shrt_ax frag_mn_shrt_ax_mm frag_max_shrt_ax...
      frag_max_shrt_ax_mm frag_mn_lng_ax frag_mn_lng_ax_mm...
      frag_max_lng_ax frag_max_lng_ax_mm frag_rel_div;
end


save children_variables  cx cy cz...
    ctotalpix cxpix cypix cgrv ctime  child ...
    cco cAij cSij crij cw1 cw2 cw3 cv1 cv2 cv3 ce1 ce2 ce3 cstrain censtro...
    cmn_shrt_ax  cmax_shrt_ax   cmn_lng_ax   cmax_lng_ax...
    cmn_shrt_ax_mm   cmax_shrt_ax_mm  cmn_lng_ax_mm  cmax_lng_ax_mm  crel_div 
 




%%%%%%%%%%%%%%%%%%%%%%%%% Figure section %%%%%%%%%%%%%%%%%%%%%%% 

sep=zeros(size(child,2),500,7);
a=1;    
parent_time=T(T>0);
X=X(abs(X)>0);Y=Y(abs(Y)>0);Z=Z(abs(Z)>0);
for j=1:size(child,2) % j is/are the trajectory other than the parent one
    begi_fr=ctime{j}(1);
    if parent_time(end)>ctime{j}(end)
        end_fr=ctime{j}(end);
    else
        end_fr=parent_time(end);
    end
    frame_length=end_fr-begi_fr+1;
    cc=0;
    for gg=1:frame_length;
        %e=(begi_fr+gg-a)-parent_time(1)+1;% This was for ideal case when
                                           % parent traj was not fragmented anywhere
        e=find(ctime{j}(gg)==T);%This for when the parent line is fragmented either before or after breakage!
        if ~isempty(e)% This is if the parent line is fragmented after breakage frame
            cc=cc+1;
        sep(j,cc,1:7)=[X(e) cx{j}(cc) Y(e) cy{j}(cc) Z(e) cz{j}(cc)...
                       sqrt((X(e)-cx{j}(cc)).^2+(Y(e)-cy{j}(cc)).^2+(Z(e)-cz{j}(cc)).^2)];
        end
        %hold on
        %plot3([X(e) cx{j}(gg)], [Y(e) cy{j}(gg)] ,[Z(e) cz{j}(gg)],'b')
    end
end

%save sep



%figure;
%threshold=(mn_shrt_ax_mm(par,1)/2)*1e-3;
threshold=(parent_max_lng_ax_mm(1)/2)*1e-3;
for hh=1:size(child,2)

    sep_dist=sep(hh,:,7);
    sep_dist=sep_dist(sep_dist>0);
    if length(sep_dist)>35
        r=sep_dist(1:35);
        t=ctime{hh}(1:35);
    else
        r=sep_dist;
        t=ctime{hh}(1:length(r));
    end
    
    p=polyfit(t,r,1);
    break_frame_X1(hh)=round((threshold-p(2))/p(1));
   
    if  break_frame_X1(hh)>ctime{j}(1)  
    break_frame_X1(hh)=ctime{j}(1);
    end
    
%      hold on
%      plot(t,r)
%      clear sep_dist r t
end
% xlabel('Frame number')
% ylabel('r (m)')



%%% according to balastic regime.......
threshold=(parent_max_lng_ax_mm(1)/2)*1e-3;
epsi=9e-5;
K=11/3;
C2=2.13;
parent_time=T(T>0);
X=X(abs(X)>0);Y=Y(abs(Y)>0);Z=Z(abs(Z)>0);
for j=1:size(child,2) % j is/are the trajectory other than the parent one
    
    %begi_fr=ctime{j}(1);
    %e=begi_fr-parent_time(1)+1;
    e=find(ctime{j}(1)==T);
    if isempty(e)
        diference=abs(T-ctime{j}(1));
        a=find(diference==min(diference));
        e=a;
    end
        
    del_r(j)=sqrt((X(e)-cx{j}(1)).^2+(Y(e)-cy{j}(1)).^2+(Z(e)-cz{j}(1)).^2);
    tt(j)=ctime{j}(1)-round(sqrt((del_r(j)-threshold)^2/(K*C2*(epsi*threshold)^(2/3)))*250);
    if tt(j)>0
        index=find(T==tt(j));
        if~isempty(index)%% This is when the formula does not work !!!
                         %%% i.e., it goes back beyond the parent traj
                         %%% begining time(Example WF 10)
            break_frame_X2(j)=T(T==tt(j));
        
     else
            break_frame_X2(j)=ctime{j}(1);
        end
    end
    %hold on
    %plot3([X(e) cx{j}(1)], [Y(e) cy{j}(1)] ,[Z(e) cz{j}(1)],'r')
end



% or take as it is for breakage frame
for kk=1:size(child,2)
break_frame_ptv(kk)=time(child(kk),1);
end

save break_frame break_frame_ptv break_frame_X1 break_frame_X2



cd(wd)

for ee=1:size(child,2)
figure;
ptv_fr=break_frame_ptv(ee);
%ptv_fr=break_frame_X2(ee);

subplot(2,2,1)
imshow([fold2,'image_agg\','filt_Cam1.',num2str(ptv_fr)])
title(['breakage at frame ',num2str(ptv_fr)])
subplot(2,2,2)
imshow([fold2,'image_agg\','filt_Cam2.',num2str(ptv_fr)])
title(['breakage at frame ',num2str(ptv_fr)])
subplot(2,2,4)
imshow([fold2,'image_agg\','filt_Cam3.',num2str(ptv_fr)])
title(['breakage at frame ',num2str(ptv_fr)])
subplot(2,2,3)
imshow([fold2,'image_agg\','filt_Cam4.',num2str(ptv_fr)])
title(['breakage at frame ',num2str(ptv_fr)])
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


figure;% Time series of Strain
hold on;
break_ptv=find(T==break_frame_ptv);
break_manual=find(T==break_fr_manual);

time_axis=((1:length(T))-break_manual)/35;

%scatter(time_axis,parent_strain,50,parent_co,'filled')% parent_rel_div,
%parent_co
scatter(time_axis,parent_strain,'.k')% parent_rel_div

scatter(time_axis(break_ptv),parent_strain(break_ptv),75,'rs','filled')
scatter(time_axis(break_manual),parent_strain(break_manual),75,'gs','filled')

legend('strain','ptv','break\_fr\_manual' )
xlabel('t/\tau_\eta')
ylabel('S_{ij}^2')



figure;% Time series of enstro
hold on;
break_ptv=find(T==break_frame_ptv);
break_manual=find(T==break_fr_manual);

time_axis=((1:length(T))-break_manual)/35;

plot(time_axis,parent_enstro,'.k')

scatter(time_axis(break_ptv),parent_enstro(break_ptv),75,'rs','filled')
scatter(time_axis(break_manual),parent_enstro(break_manual),75,'gs','filled')

legend('enstrophy','ptv','break\_fr\_manual')
xlabel('t/\tau_\eta')
ylabel('\omega^2')












figure; % deformation 
c_time=1;
dt=1/250;
B(c_time)=[1 0 0;0 1 0;0 0 1];
for time=0:dt:1
    h=pAij;
    c_time=c_time+1;
    dbdt=h*B(c_time-1);
    B(c_time)=B(c_time-1)+dt*dbdt;
end

figure;% change in total pixel in time
for traj_id=1:length(tid)
    t=time(traj_id,:);
    t=t(t>0);
    total_pix=totalpix(traj_id,:);
    total_pix=total_pix(1:numel(t));%%% This is a bug remover !!
    hold on
    scatter(t,total_pix)
    hold off
end



figure;% how total pixel changes with depth
for traj_id=1:length(tid)
    t=time(traj_id,:);
    t=t(t>0);
    total_pix=totalpix(traj_id,:);
    total_pix=total_pix(1:numel(t));%%% This is a bug remover !!
   
    t_z=tz(traj_id,:);
    t_z=t_z(1:numel(t));%%% This is a bug remover !!
    hold on
    scatter(t_z,total_pix,'.')
    plot(t_z,total_pix,'r')
    hold off
end


figure;% How many particles are there for gradient
for traj_id=1:length(tid)
    t=time(traj_id,:);
    t=t(t>0);
    particle=co(traj_id,:);
    particle=particle(1:numel(t));%%% This is a bug remover !!
    hold on
    plot(t,particle)
    hold off
end

figure;% relative divergence
for traj_id=1:length(tid)
    t=time(traj_id,:);
    t=t(t>0);
    rldv=rel_div(traj_id,:);
    rldv=rldv(1:numel(t));%%% This is a bug remover !!
    hold on
    [nout,xout]=nhist(rldv);
    plot(xout,nout)
    hold off
    %length(rldv(rldv<0.2))/length(rldv)*100
end

































%%% This is Beat's figure part------------------
%----------------------------------------------------
a=146;b=-7; %looking down along omega 
%a=147;b=16; %visible separation

figure;hold on;
ind_a=find(xa<1e9);
ind_b=find(xa(ind_a)<0.0125);
ind_c=find(za(ind_a(ind_b))>0.002);
ind=ind_a(ind_b(ind_c));
scatter3(xa(ind),ya(ind),za(ind),30,'b','filled')
xlabel('x')
ylabel('y')
zlabel('z')
view(a,b)
axis equal;
box on;

for i=1:tid
    count=0;
    for j=1:traj_len(i)
        if 2*strain(i,j)<1000 & enstro(i,j)<1000 & strain(i,j)>0
            count=count+1;
            tx(count)=agg(traj(i,j,1),traj(i,j,2),1);
            ty(count)=agg(traj(i,j,1),traj(i,j,2),2);
            tz(count)=agg(traj(i,j,1),traj(i,j,2),3);
            ve1(count,1:3)=v1(i,j,1:3);
            ve2(count,1:3)=v2(i,j,1:3);
            ve3(count,1:3)=v3(i,j,1:3);
            ei1(count)=e1(i,j);
            ei2(count)=e2(i,j);
            ei3(count)=e3(i,j);
            wo1(count)=w1(i,j);
            wo2(count)=w2(i,j);
            wo3(count)=w3(i,j);
            st(count)=strain(i,j);
            en(count)=enstro(i,j);
        end
    end
    figure(10+offset);hold on;
    scatter3(tx(1:count),ty(1:count),st(1:count),15,'r','filled')
    xlabel('x')
    ylabel('y')
    zlabel('s^2')
    view(a,b)
    box on;
    figure(11+offset);hold on;
    scatter3(tx(1:count),ty(1:count),en(1:count),15,'r','filled')
    xlabel('x')
    ylabel('y')
    zlabel('\omega^2')
    view(a,b)
    box on;
    figure(12+offset);hold on;
    quiver3(tx(1:count),ty(1:count),tz(1:count),ve1(1:count,1),ve1(1:count,2),ve1(1:count,3),1,'r')
    quiver3(tx(1:count),ty(1:count),tz(1:count),-ve1(1:count,1),-ve1(1:count,2),-ve1(1:count,3),1,'r')
    quiver3(tx(1:count),ty(1:count),tz(1:count),ve2(1:count,1),ve2(1:count,2),ve2(1:count,3),1,'g')
    quiver3(tx(1:count),ty(1:count),tz(1:count),-ve2(1:count,1),-ve2(1:count,2),-ve2(1:count,3),1,'g')
    xlabel('x')
    ylabel('y')
    zlabel('z')
    title('\lambda_1 (red), \lambda_2 (green)')
    view(a,b)
    axis equal;
    box on;
    figure(13+offset);hold on;
    quiver3(tx(1:count),ty(1:count),tz(1:count),ve3(1:count,1),ve3(1:count,2),ve3(1:count,3),1,'b')
    quiver3(tx(1:count),ty(1:count),tz(1:count),-ve3(1:count,1),-ve3(1:count,2),-ve3(1:count,3),1,'b')
    xlabel('x')
    ylabel('y')
    zlabel('z')
    title('\lambda_3')
    view(a,b)
    axis equal;
    box on;
    figure(14+offset);hold on;
    quiver3(tx(1:count),ty(1:count),tz(1:count),wo1(1:count),wo2(1:count),wo3(1:count),2,'m')
    xlabel('x')
    ylabel('y')
    zlabel('z')
    title('\omega')
    view(a,b)
    axis equal;
    box on;
    figure(15+offset);hold on;
    scatter3(tx(1:count),ty(1:count),ei1(1:count),15,'k','filled')
    xlabel('x')
    ylabel('y')
    zlabel('z')
    title('\Lambda_1')
    view(a,b)
    box on;
end







