clear all
clear all
clear all
clear all

close all
clc
%tol=1.000e-3;
%radius_sphere=0.006;
offset=12;
long_traj_count=0;
max_length=10;
tau_eta=0.14;
frame_rate=250;
factor=tau_eta*frame_rate;
%load map_turbulent_box
order=3;width=51;
already_there=1;

wf=[1 2 6 10  12 13 15 16 18   22 23 25 26 27 28 29 30 31 32 33 34 35 36 37 38 41 42 45 46 50 51 52 53 55 56 58 59 61 62 63 64 65 66 67 68 69  71 72 73 74 75 76 77 79 80 81 82 83 84 86 87 88 89  92 93 95 97 98  99 100 101 102  104  106 107 108 109 110 111 112 113];
count=0;
for working_folder=wf
    count=count+1;
    WF=wf(count);
    true_WF=WF; 
    fold='D:\Turbulent_data_backup\flow_induced_aggregates_breakage_at_100rpm\WF';
    
    if already_there==1
        load([fold,num2str(WF),'/variables_v2'])
        WF=true_WF; % Reason: z.B., when loading variables_v2  
                    % it loads out WF=37 because this event was coming from
                    % WF37 and WF=37 was saved in variables_v2. Hence, in time of saving the variables,it
                    % over writes the variables in WF37 although processing
                    % WF104. That is why WF=37, as saved in variables_v2 is
                    % suppressed by current WF
    else
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

    fold2=[fold,num2str(WF),'\'];
    name_root=[fold,num2str(WF),'\res_agg'];
    
    agg=repmat(zeros,[last-first+1,200,19]);
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
        save  variables_v2 agg
        cd(wd) 

    end
    
  
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
    ['act4: Calculating grad around traj = ',num2str(i),' of ',num2str(tid),'in WF=',num2str(WF)]
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
        
        clear V D
    end
end 



for i=1:tid
    for j=1:traj_len(i)
    %remove crap
    w1_rc(i,:)=remove_outlayer(w1(i,:),traj_len(i));
    w2_rc(i,:)=remove_outlayer(w2(i,:),traj_len(i));
    w3_rc(i,:)=remove_outlayer(w3(i,:),traj_len(i));
    
    w11_rc(i,:)=remove_outlayer(w11(i,:),traj_len(i));
    w12_rc(i,:)=remove_outlayer(w12(i,:),traj_len(i));
    w13_rc(i,:)=remove_outlayer(w13(i,:),traj_len(i));
    w21_rc(i,:)=remove_outlayer(w21(i,:),traj_len(i));
    w22_rc(i,:)=remove_outlayer(w22(i,:),traj_len(i));
    w23_rc(i,:)=remove_outlayer(w23(i,:),traj_len(i));
    w31_rc(i,:)=remove_outlayer(w31(i,:),traj_len(i));
    w32_rc(i,:)=remove_outlayer(w32(i,:),traj_len(i));
    w33_rc(i,:)=remove_outlayer(w33(i,:),traj_len(i));
    
    s11_rc(i,:)=remove_outlayer(s11(i,:),traj_len(i));
    s12_rc(i,:)=remove_outlayer(s12(i,:),traj_len(i));
    s13_rc(i,:)=remove_outlayer(s13(i,:),traj_len(i));
    s22_rc(i,:)=remove_outlayer(s22(i,:),traj_len(i));
    s23_rc(i,:)=remove_outlayer(s23(i,:),traj_len(i));
    s33_rc(i,:)=remove_outlayer(s33(i,:),traj_len(i));
    end
end


    vari=load([fold,num2str(WF),'/parent_variables_v2']);
    parent=vari.parent;
    par_T=vari.T;
    par_difference=par_T(2:end)-par_T(1:end-1); % sometimes a traj begins immediately in the next frame, but still broken
    whatever=find(par_difference>1);            % in that case fill_vector is not needed and doing fill_vector actually showed error
if size(parent,2)>1 & ~isempty(whatever)        % that is why the contraint, i.e., the difference in time should be larger than 1 to use fill_vector
    
    X=[];
    Y=[];
    Z=[];
    T=[];
    
    parent_totalpix=[];
    parent_xpix=[];
    parent_ypix=[];
    parent_grv=[];
    parent_co=[];
    
    parent_mn_shrt_ax=[];
    parent_mn_shrt_ax_mm=[];
    parent_max_shrt_ax=[];
    parent_max_shrt_ax_mm=[];
    parent_mn_lng_ax=[];
    parent_mn_lng_ax_mm=[];
    parent_max_lng_ax=[];
    parent_max_lng_ax_mm=[];
    
    parent_rel_div=[];
    
    w1_p=[];
    w2_p=[];
    w3_p=[];
    w11_p=[];
    w12_p=[];
    w13_p=[];
    w21_p=[];
    w22_p=[];
    w23_p=[];
    w31_p=[];
    w32_p=[];
    w33_p=[];
    s11_p=[];
    s12_p=[];
    s13_p=[];
    s22_p=[];
    s23_p=[];
    s33_p=[];
    
    for traj_id=1:length(parent)
        
        t=time(parent(traj_id),:);
        index=find(t>0);
        
        T=[T;t(index)'];
        x=tx(parent(traj_id),:);
        X=[X;x(index)'];%%% This is a bug remover !!
        y=ty(parent(traj_id),:);
        Y=[Y;y(index)'];
        z=tz(parent(traj_id),:);
        Z=[Z;z(index)'];
        
        
        w1_p=[w1_p;w1_rc(parent(traj_id),index)'];
        w2_p=[w2_p;w2_rc(parent(traj_id),index)'];
        w3_p=[w3_p;w3_rc(parent(traj_id),index)'];
        
        w11_p=[w11_p;w11_rc(parent(traj_id),index)'];
        w12_p=[w12_p;w12_rc(parent(traj_id),index)'];
        w13_p=[w13_p;w13_rc(parent(traj_id),index)'];
        w21_p=[w21_p;w21_rc(parent(traj_id),index)'];
        w22_p=[w22_p;w22_rc(parent(traj_id),index)'];
        w23_p=[w23_p;w23_rc(parent(traj_id),index)'];
        w31_p=[w31_p;w31_rc(parent(traj_id),index)'];
        w32_p=[w32_p;w32_rc(parent(traj_id),index)'];
        w33_p=[w33_p;w33_rc(parent(traj_id),index)'];
        
        s11_p=[s11_p;s11_rc(parent(traj_id),index)'];
        s12_p=[s12_p;s12_rc(parent(traj_id),index)'];
        s13_p=[s13_p;s13_rc(parent(traj_id),index)'];
        s22_p=[s22_p;s22_rc(parent(traj_id),index)'];
        s23_p=[s23_p;s23_rc(parent(traj_id),index)'];
        s33_p=[s33_p;s33_rc(parent(traj_id),index)'];
        
        
        ptotalpix=totalpix(parent(traj_id),:);
        ptotalpix=ptotalpix(index)';
        ptotalpix=remove_outlayer(ptotalpix,length(ptotalpix));
        parent_totalpix=[parent_totalpix;ptotalpix];
        
        pxpix=xpix(parent(traj_id),:);
        pxpix=pxpix(index)';
        pxpix=remove_outlayer(pxpix,length(pxpix));
        parent_xpix=[parent_xpix;pxpix];
        
        pypix=ypix(parent(traj_id),:);
        pypix=pypix(index)';
        pypix=remove_outlayer(pypix,length(pypix));
        parent_ypix=[parent_ypix;pypix];
        
        pgrv=grv(parent(traj_id),:);
        pgrv=pgrv(index)';
        pgrv=remove_outlayer(pgrv,length(pgrv));
        parent_grv=[parent_grv;pgrv];
        
        
        pco=co(parent(traj_id),:);
        pco=pco(index)';
        pco=remove_outlayer(pco,length(pco));
        parent_co=[parent_co;pco];
        
        
        pmn_shrt_ax=mn_shrt_ax(parent(traj_id),:);
        pmn_shrt_ax=pmn_shrt_ax(index)';
        pmn_shrt_ax=remove_outlayer(pmn_shrt_ax, length(pmn_shrt_ax));
        parent_mn_shrt_ax=[parent_mn_shrt_ax;pmn_shrt_ax];
        
        pmax_shrt_ax=max_shrt_ax(parent(traj_id),:);
        pmax_shrt_ax=pmax_shrt_ax(index)';
        pmax_shrt_ax=remove_outlayer(pmax_shrt_ax,length(pmax_shrt_ax));
        parent_max_shrt_ax=[parent_max_shrt_ax;pmax_shrt_ax];
        
        pmn_shrt_ax_mm=mn_shrt_ax_mm(parent(traj_id),:);
        pmn_shrt_ax_mm=pmn_shrt_ax_mm(index)';
        pmn_shrt_ax_mm=remove_outlayer(pmn_shrt_ax_mm,length(pmn_shrt_ax_mm));
        parent_mn_shrt_ax_mm=[parent_mn_shrt_ax_mm;pmn_shrt_ax_mm];
        
        pmax_shrt_ax_mm=max_shrt_ax_mm(parent(traj_id),:);
        pmax_shrt_ax_mm=pmax_shrt_ax_mm(index)';
        pmax_shrt_ax_mm=remove_outlayer(pmax_shrt_ax_mm,length(pmax_shrt_ax_mm));
        parent_max_shrt_ax_mm=[parent_max_shrt_ax_mm;pmax_shrt_ax_mm];
        
        pmn_lng_ax=mn_lng_ax(parent(traj_id),:);
        pmn_lng_ax= pmn_lng_ax(index)';
        pmn_lng_ax=remove_outlayer(pmn_lng_ax,length(pmn_lng_ax));
        parent_mn_lng_ax=[parent_mn_lng_ax;pmn_lng_ax];
        
        pmax_lng_ax=max_lng_ax(parent(traj_id),:);
        pmax_lng_ax=pmax_lng_ax(index)';
        pmax_lng_ax=remove_outlayer(pmax_lng_ax,length(pmax_lng_ax));
        parent_max_lng_ax=[parent_max_lng_ax;pmax_lng_ax];
        
        pmn_lng_ax_mm=mn_lng_ax_mm(parent(traj_id),:);
        pmn_lng_ax_mm=pmn_lng_ax_mm(index)';
        pmn_lng_ax_mm=remove_outlayer(pmn_lng_ax_mm,length(pmn_lng_ax_mm));
        parent_mn_lng_ax_mm=[parent_mn_lng_ax_mm;pmn_lng_ax_mm];
        
        pmax_lng_ax_mm=max_lng_ax_mm(parent(traj_id),:);
        pmax_lng_ax_mm=pmax_lng_ax_mm(index)';
        pmax_lng_ax_mm=remove_outlayer(pmax_lng_ax_mm,length(pmax_lng_ax_mm));
        parent_max_lng_ax_mm=[parent_max_lng_ax_mm;pmax_lng_ax_mm];
        
        
    end
    
    %%%% Fill and smooth
    
    %width=51;%41;
    %order=3;%2;
    
    
    T_fill=fill_time(T);
    pTIME=T_fill;%smooth(T_fill,width,'sgolay',order);
    X=fill_vector(T,X);
    X=smooth(X,width,'sgolay',order);
    Y=fill_vector(T,Y);
    Y=smooth(Y,width,'sgolay',order);
    Z=fill_vector(T,Z);
    Z=smooth(Z,width,'sgolay',order);
    
    
    
    parent_totalpix=fill_vector(T,parent_totalpix);
    parent_totalpix=smooth(parent_totalpix,width,'sgolay',order);
    
    parent_xpix=fill_vector(T, parent_xpix);
    parent_xpix=smooth( parent_xpix,width,'sgolay',order);
    
    parent_ypix=fill_vector(T, parent_ypix);
    parent_ypix=smooth( parent_ypix,width,'sgolay',order);
    
    parent_grv=fill_vector(T, parent_grv);
    parent_grv=smooth( parent_grv,width,'sgolay',order);
    
    parent_co=fill_vector(T, parent_co);
    parent_co=smooth( parent_co,width,'sgolay',order);
    
    parent_mn_shrt_ax=fill_vector(T, parent_mn_shrt_ax);
    parent_mn_shrt_ax=smooth( parent_mn_shrt_ax,width,'sgolay',order);
    
    parent_max_shrt_ax=fill_vector(T, parent_max_shrt_ax);
    parent_max_shrt_ax=smooth( parent_max_shrt_ax,width,'sgolay',order);
    
    parent_mn_shrt_ax_mm=fill_vector(T, parent_mn_shrt_ax_mm);
    parent_mn_shrt_ax_mm=smooth( parent_mn_shrt_ax_mm,width,'sgolay',order);
    
    parent_max_shrt_ax_mm=fill_vector(T, parent_max_shrt_ax_mm);
    parent_max_shrt_ax_mm=smooth( parent_max_shrt_ax_mm,width,'sgolay',order);
    
    parent_mn_lng_ax=fill_vector(T, parent_mn_lng_ax);
    parent_mn_lng_ax=smooth( parent_mn_lng_ax,width,'sgolay',order);
    
    parent_max_lng_ax=fill_vector(T, parent_max_lng_ax);
    parent_max_lng_ax=smooth( parent_max_lng_ax,width,'sgolay',order);
    
    parent_mn_lng_ax_mm=fill_vector(T, parent_mn_lng_ax_mm);
    parent_mn_lng_ax_mm=smooth( parent_mn_lng_ax_mm,width,'sgolay',order);
    
    parent_max_lng_ax_mm=fill_vector(T, parent_max_lng_ax_mm);
    parent_max_lng_ax_mm=smooth( parent_max_lng_ax_mm,width,'sgolay',order);
    
    
    w1_p=fill_vector(T, w1_p);
    w2_p=fill_vector(T, w2_p);
    w3_p=fill_vector(T, w3_p);
    
    
    w11_p=fill_vector(T, w11_p);
    w12_p=fill_vector(T, w12_p);
    w13_p=fill_vector(T, w13_p);
    w21_p=fill_vector(T, w21_p);
    w22_p=fill_vector(T, w22_p);
    w23_p=fill_vector(T, w23_p);
    w31_p=fill_vector(T, w31_p);
    w32_p=fill_vector(T, w32_p);
    w33_p=fill_vector(T, w33_p);
    
    s11_p=fill_vector(T, s11_p);
    s12_p=fill_vector(T, s12_p);
    s13_p=fill_vector(T, s13_p);
    s22_p=fill_vector(T, s22_p);
    s23_p=fill_vector(T, s23_p);
    s33_p=fill_vector(T, s33_p);
    
    
    
    
    parent_w1 = smooth(w1_p,width,'sgolay',order);
    parent_w2 = smooth(w2_p,width,'sgolay',order);
    parent_w3 = smooth(w3_p,width,'sgolay',order);
    
    w11_sg=smooth(w11_p,width,'sgolay',order);
    w12_sg=smooth(w12_p,width,'sgolay',order);
    w13_sg=smooth(w13_p,width,'sgolay',order);
    w21_sg=smooth(w21_p,width,'sgolay',order);
    w22_sg=smooth(w22_p,width,'sgolay',order);
    w23_sg=smooth(w23_p,width,'sgolay',order);
    w31_sg=smooth(w31_p,width,'sgolay',order);
    w32_sg=smooth(w32_p,width,'sgolay',order);
    w33_sg=smooth(w33_p,width,'sgolay',order);
    
    s11_sg = smooth(s11_p,width,'sgolay',order);
    s12_sg = smooth(s12_p,width,'sgolay',order);
    s13_sg = smooth(s13_p,width,'sgolay',order);
    s22_sg = smooth(s22_p,width,'sgolay',order);
    s23_sg = smooth(s23_p,width,'sgolay',order);
    s33_sg = smooth(s33_p,width,'sgolay',order);
    
    
    
    for j=1:length(parent_w1)
        
        parent_rij(j,1:3,1:3)=[w11_sg(j)  w12_sg(j) w13_sg(j); w21_sg(j) w22_sg(j) w23_sg(j); w31_sg(j)  w32_sg(j) w33_sg(j)];
        parent_enstro(j)=parent_w1(j)^2+parent_w2(j)^2+parent_w3(j)^2;
        %A=[s11(i,j) s12(i,j) s13(i,j);s12(i,j) s22(i,j) s23(i,j);s13(i,j) s23(i,j) s33(i,j)];
        %[V,D] = eig(A);
        
        parent_Sij(j,1:3,1:3)=[s11_sg(j) s12_sg(j) s13_sg(j);s12_sg(j) s22_sg(j) s23_sg(j);s13_sg(j) s23_sg(j) s33_sg(j)];
        
        parent_Aij(j,1:3,1:3)=parent_Sij(j,1:3,1:3)+parent_rij(j,1:3,1:3);
        
        [V,D] = eig(squeeze(parent_Sij(j,1:3,1:3)));
        
        [D,k] = sort(diag(D)); 	% ascending order
        D = diag(D(end:-1:1));  % descending order
        V = V(:,k(end:-1:1)); 	% the same order for eigenvectors
        parent_v1(j,1:3)=V(:,1)';
        parent_v2(j,1:3)=V(:,2)';
        parent_v3(j,1:3)=V(:,3)';
        parent_e1(j)=D(1,1);
        parent_e2(j)=D(2,2);
        parent_e3(j)=D(3,3);
        parent_strain(j)=parent_e1(j)^2+parent_e2(j)^2+parent_e3(j)^2;
        % THIS IS DIVERGENCE FREE
        parent_rel_div (j)=abs(s11_sg(j)+s22_sg(j)+s33_sg(j))/(abs(s11_sg(j))+abs(s22_sg(j))+abs(s33_sg(j)));
   
        clear V D
    end
    
    
    wd= 'E:\OLD-Data-DISK\Matlab for PTV';
    cd([fold,num2str(WF)]);
    
    save parent_variables_v3 parent X Y Z pTIME...
        parent_totalpix  parent_xpix  parent_ypix  parent_grv   ...
        parent_co parent_Aij parent_Sij parent_rij parent_w1 parent_w2 parent_w3 parent_v1 parent_v2 parent_v3...
        parent_e1 parent_e2 parent_e3 parent_strain parent_enstro...
        parent_mn_shrt_ax    parent_max_shrt_ax   parent_mn_lng_ax   parent_max_lng_ax...
        parent_mn_shrt_ax_mm   parent_max_shrt_ax_mm  parent_mn_lng_ax_mm...
        parent_max_lng_ax_mm  parent_rel_div
    
    cd(wd);
    
    
    
else
    
    X=[];
    Y=[];
    Z=[];
    T=[];
    
    parent_totalpix=[];
    parent_xpix=[];
    parent_ypix=[];
    parent_grv=[];
    parent_co=[];
    
    parent_mn_shrt_ax=[];
    parent_mn_shrt_ax_mm=[];
    parent_max_shrt_ax=[];
    parent_max_shrt_ax_mm=[];
    parent_mn_lng_ax=[];
    parent_mn_lng_ax_mm=[];
    parent_max_lng_ax=[];
    parent_max_lng_ax_mm=[];
    
    parent_rel_div=[];
    
    w1_p=[];
    w2_p=[];
    w3_p=[];
    w11_p=[];
    w12_p=[];
    w13_p=[];
    w21_p=[];
    w22_p=[];
    w23_p=[];
    w31_p=[];
    w32_p=[];
    w33_p=[];
    s11_p=[];
    s12_p=[];
    s13_p=[];
    s22_p=[];
    s23_p=[];
    s33_p=[];
    
    for traj_id=1:length(parent)
        
        t=time(parent(traj_id),:);
        index=find(t>0);
        
        T=[T;t(index)'];
        pTIME=T;
        x=tx(parent(traj_id),:);
        X=[X;x(index)'];%%% This is a bug remover !!
        y=ty(parent(traj_id),:);
        Y=[Y;y(index)'];
        z=tz(parent(traj_id),:);
        Z=[Z;z(index)'];
        
        
        w1_p=[w1_p;w1_rc(parent(traj_id),index)'];
        w2_p=[w2_p;w2_rc(parent(traj_id),index)'];
        w3_p=[w3_p;w3_rc(parent(traj_id),index)'];
        
        w11_p=[w11_p;w11_rc(parent(traj_id),index)'];
        w12_p=[w12_p;w12_rc(parent(traj_id),index)'];
        w13_p=[w13_p;w13_rc(parent(traj_id),index)'];
        w21_p=[w21_p;w21_rc(parent(traj_id),index)'];
        w22_p=[w22_p;w22_rc(parent(traj_id),index)'];
        w23_p=[w23_p;w23_rc(parent(traj_id),index)'];
        w31_p=[w31_p;w31_rc(parent(traj_id),index)'];
        w32_p=[w32_p;w32_rc(parent(traj_id),index)'];
        w33_p=[w33_p;w33_rc(parent(traj_id),index)'];
        
        s11_p=[s11_p;s11_rc(parent(traj_id),index)'];
        s12_p=[s12_p;s12_rc(parent(traj_id),index)'];
        s13_p=[s13_p;s13_rc(parent(traj_id),index)'];
        s22_p=[s22_p;s22_rc(parent(traj_id),index)'];
        s23_p=[s23_p;s23_rc(parent(traj_id),index)'];
        s33_p=[s33_p;s33_rc(parent(traj_id),index)'];
        
        
        ptotalpix=totalpix(parent(traj_id),:);
        ptotalpix=ptotalpix(index)';
        ptotalpix=remove_outlayer(ptotalpix,length(ptotalpix));
        parent_totalpix=[parent_totalpix;ptotalpix];
        
        pxpix=xpix(parent(traj_id),:);
        pxpix=pxpix(index)';
        pxpix=remove_outlayer(pxpix,length(pxpix));
        parent_xpix=[parent_xpix;pxpix];
        
        pypix=ypix(parent(traj_id),:);
        pypix=pypix(index)';
        pypix=remove_outlayer(pypix,length(pypix));
        parent_ypix=[parent_ypix;pypix];
        
        pgrv=grv(parent(traj_id),:);
        pgrv=pgrv(index)';
        pgrv=remove_outlayer(pgrv,length(pgrv));
        parent_grv=[parent_grv;pgrv];
        
        
        pco=co(parent(traj_id),:);
        pco=pco(index)';
        pco=remove_outlayer(pco,length(pco));
        parent_co=[parent_co;pco];
        
        
        pmn_shrt_ax=mn_shrt_ax(parent(traj_id),:);
        pmn_shrt_ax=pmn_shrt_ax(index)';
        pmn_shrt_ax=remove_outlayer(pmn_shrt_ax, length(pmn_shrt_ax));
        parent_mn_shrt_ax=[parent_mn_shrt_ax;pmn_shrt_ax];
        
        pmax_shrt_ax=max_shrt_ax(parent(traj_id),:);
        pmax_shrt_ax=pmax_shrt_ax(index)';
        pmax_shrt_ax=remove_outlayer(pmax_shrt_ax,length(pmax_shrt_ax));
        parent_max_shrt_ax=[parent_max_shrt_ax;pmax_shrt_ax];
        
        pmn_shrt_ax_mm=mn_shrt_ax_mm(parent(traj_id),:);
        pmn_shrt_ax_mm=pmn_shrt_ax_mm(index)';
        pmn_shrt_ax_mm=remove_outlayer(pmn_shrt_ax_mm,length(pmn_shrt_ax_mm));
        parent_mn_shrt_ax_mm=[parent_mn_shrt_ax_mm;pmn_shrt_ax_mm];
        
        pmax_shrt_ax_mm=max_shrt_ax_mm(parent(traj_id),:);
        pmax_shrt_ax_mm=pmax_shrt_ax_mm(index)';
        pmax_shrt_ax_mm=remove_outlayer(pmax_shrt_ax_mm,length(pmax_shrt_ax_mm));
        parent_max_shrt_ax_mm=[parent_max_shrt_ax_mm;pmax_shrt_ax_mm];
        
        pmn_lng_ax=mn_lng_ax(parent(traj_id),:);
        pmn_lng_ax= pmn_lng_ax(index)';
        pmn_lng_ax=remove_outlayer(pmn_lng_ax,length(pmn_lng_ax));
        parent_mn_lng_ax=[parent_mn_lng_ax;pmn_lng_ax];
        
        pmax_lng_ax=max_lng_ax(parent(traj_id),:);
        pmax_lng_ax=pmax_lng_ax(index)';
        pmax_lng_ax=remove_outlayer(pmax_lng_ax,length(pmax_lng_ax));
        parent_max_lng_ax=[parent_max_lng_ax;pmax_lng_ax];
        
        pmn_lng_ax_mm=mn_lng_ax_mm(parent(traj_id),:);
        pmn_lng_ax_mm=pmn_lng_ax_mm(index)';
        pmn_lng_ax_mm=remove_outlayer(pmn_lng_ax_mm,length(pmn_lng_ax_mm));
        parent_mn_lng_ax_mm=[parent_mn_lng_ax_mm;pmn_lng_ax_mm];
        
        pmax_lng_ax_mm=max_lng_ax_mm(parent(traj_id),:);
        pmax_lng_ax_mm=pmax_lng_ax_mm(index)';
        pmax_lng_ax_mm=remove_outlayer(pmax_lng_ax_mm,length(pmax_lng_ax_mm));
        parent_max_lng_ax_mm=[parent_max_lng_ax_mm;pmax_lng_ax_mm];
        
        
    end
    
    %%%% Only smooth
    
    %width=51;%41;
    %order=3;%2;
    
    
    
    X=smooth(X,width,'sgolay',order);
    Y=smooth(Y,width,'sgolay',order);
    Z=smooth(Z,width,'sgolay',order);
    
    parent_totalpix=smooth(parent_totalpix,width,'sgolay',order);
    parent_xpix=smooth( parent_xpix,width,'sgolay',order);
    parent_ypix=smooth( parent_ypix,width,'sgolay',order);
    
    parent_grv=smooth( parent_grv,width,'sgolay',order);
    
    parent_co=smooth( parent_co,width,'sgolay',order);
    
    parent_mn_shrt_ax=smooth( parent_mn_shrt_ax,width,'sgolay',order);
    parent_max_shrt_ax=smooth( parent_max_shrt_ax,width,'sgolay',order);
    
    parent_mn_shrt_ax_mm=smooth( parent_mn_shrt_ax_mm,width,'sgolay',order);
    parent_max_shrt_ax_mm=smooth( parent_max_shrt_ax_mm,width,'sgolay',order);
    
    parent_mn_lng_ax=smooth( parent_mn_lng_ax,width,'sgolay',order);
    parent_max_lng_ax=smooth( parent_max_lng_ax,width,'sgolay',order);
    
    parent_mn_lng_ax_mm=smooth( parent_mn_lng_ax_mm,width,'sgolay',order);
    parent_max_lng_ax_mm=smooth( parent_max_lng_ax_mm,width,'sgolay',order);
    
    parent_w1 = smooth(w1_p,width,'sgolay',order);
    parent_w2 = smooth(w2_p,width,'sgolay',order);
    parent_w3 = smooth(w3_p,width,'sgolay',order);
    
    w11_sg=smooth(w11_p,width,'sgolay',order);
    w12_sg=smooth(w12_p,width,'sgolay',order);
    w13_sg=smooth(w13_p,width,'sgolay',order);
    w21_sg=smooth(w21_p,width,'sgolay',order);
    w22_sg=smooth(w22_p,width,'sgolay',order);
    w23_sg=smooth(w23_p,width,'sgolay',order);
    w31_sg=smooth(w31_p,width,'sgolay',order);
    w32_sg=smooth(w32_p,width,'sgolay',order);
    w33_sg=smooth(w33_p,width,'sgolay',order);
    
    s11_sg = smooth(s11_p,width,'sgolay',order);
    s12_sg = smooth(s12_p,width,'sgolay',order);
    s13_sg = smooth(s13_p,width,'sgolay',order);
    s22_sg = smooth(s22_p,width,'sgolay',order);
    s23_sg = smooth(s23_p,width,'sgolay',order);
    s33_sg = smooth(s33_p,width,'sgolay',order);
    
    
    
    for j=1:length(parent_w1)
        
        parent_rij(j,1:3,1:3)=[w11_sg(j)  w12_sg(j) w13_sg(j); w21_sg(j) w22_sg(j) w23_sg(j); w31_sg(j)  w32_sg(j) w33_sg(j)];
        parent_enstro(j)=parent_w1(j)^2+parent_w2(j)^2+parent_w3(j)^2;
        %A=[s11(i,j) s12(i,j) s13(i,j);s12(i,j) s22(i,j) s23(i,j);s13(i,j) s23(i,j) s33(i,j)];
        %[V,D] = eig(A);
        
        parent_Sij(j,1:3,1:3)=[s11_sg(j) s12_sg(j) s13_sg(j);s12_sg(j) s22_sg(j) s23_sg(j);s13_sg(j) s23_sg(j) s33_sg(j)];
        
        parent_Aij(j,1:3,1:3)=parent_Sij(j,1:3,1:3)+parent_rij(j,1:3,1:3);

        [V,D] = eig(squeeze(parent_Sij(j,1:3,1:3)));
        
        [D,k] = sort(diag(D)); 	% ascending order
        D = diag(D(end:-1:1));  % descending order
        V = V(:,k(end:-1:1)); 	% the same order for eigenvectors
        parent_v1(j,1:3)=V(:,1)';
        parent_v2(j,1:3)=V(:,2)';
        parent_v3(j,1:3)=V(:,3)';
        parent_e1(j)=D(1,1);
        parent_e2(j)=D(2,2);
        parent_e3(j)=D(3,3);
        parent_strain(j)=parent_e1(j)^2+parent_e2(j)^2+parent_e3(j)^2;
        % THIS IS DIVERGENCE FREE
        parent_rel_div (j)=abs(s11_sg(j)+s22_sg(j)+s33_sg(j))/(abs(s11_sg(j))+abs(s22_sg(j))+abs(s33_sg(j)));
    
        clear V D
    
    
    end
    
    
    
    wd= 'E:\OLD-Data-DISK\Matlab for PTV';
    cd([fold,num2str(WF)]);
    
   save parent_variables_v3 parent X Y Z pTIME...
        parent_totalpix  parent_xpix  parent_ypix  parent_grv   ...
        parent_co parent_Aij parent_Sij parent_rij parent_w1 parent_w2 parent_w3 parent_v1 parent_v2 parent_v3...
        parent_e1 parent_e2 parent_e3 parent_strain parent_enstro...
        parent_mn_shrt_ax    parent_max_shrt_ax   parent_mn_lng_ax   parent_max_lng_ax...
        parent_mn_shrt_ax_mm   parent_max_shrt_ax_mm  parent_mn_lng_ax_mm...
        parent_max_lng_ax_mm  parent_rel_div
    
    
    cd(wd);
end

clear index t
%%%% For children (all the childrens are non broken)
    c_X=[];
    c_Y=[];
    c_Z=[];
    c_T=[];
    
    c_totalpix=[];
    c_xpix=[];
    c_ypix=[];
    c_grv=[];
    c_co=[];
    
    c_mn_shrt_ax=[];
    c_mn_shrt_ax_mm=[];
    c_max_shrt_ax=[];
    c_max_shrt_ax_mm=[];
    c_mn_lng_ax=[];
    c_mn_lng_ax_mm=[];
    c_max_lng_ax=[];
    c_max_lng_ax_mm=[];
    
    c_rel_div=[];
    
    w1_c=[];
    w2_c=[];
    w3_c=[];
    w11_c=[];
    w12_c=[];
    w13_c=[];
    w21_c=[];
    w22_c=[];
    w23_c=[];
    w31_c=[];
    w32_c=[];
    w33_c=[];
    s11_c=[];
    s12_c=[];
    s13_c=[];
    s22_c=[];
    s23_c=[];
    s33_c=[];
    
    for traj_id=1:length(child)
        
        t=time(child(traj_id),:);
        index=find(t>0);
        
        c_T=[c_T;t(index)'];
        cTIME=c_T;
        x=tx(child(traj_id),:);
        c_X=[c_X;x(index)'];%%% This is a bug remover !!
        y=ty(child(traj_id),:);
        c_Y=[c_Y;y(index)'];
        z=tz(child(traj_id),:);
        c_Z=[c_Z;z(index)'];
        
        
        w1_c=[w1_c;w1_rc(child(traj_id),index)'];
        w2_c=[w2_c;w2_rc(child(traj_id),index)'];
        w3_c=[w3_c;w3_rc(child(traj_id),index)'];
        
        w11_c=[w11_c;w11_rc(child(traj_id),index)'];
        w12_c=[w12_c;w12_rc(child(traj_id),index)'];
        w13_c=[w13_c;w13_rc(child(traj_id),index)'];
        w21_c=[w21_c;w21_rc(child(traj_id),index)'];
        w22_c=[w22_c;w22_rc(child(traj_id),index)'];
        w23_c=[w23_c;w23_rc(child(traj_id),index)'];
        w31_c=[w31_c;w31_rc(child(traj_id),index)'];
        w32_c=[w32_c;w32_rc(child(traj_id),index)'];
        w33_c=[w33_c;w33_rc(child(traj_id),index)'];
        
        s11_c=[s11_c;s11_rc(child(traj_id),index)'];
        s12_c=[s12_c;s12_rc(child(traj_id),index)'];
        s13_c=[s13_c;s13_rc(child(traj_id),index)'];
        s22_c=[s22_c;s22_rc(child(traj_id),index)'];
        s23_c=[s23_c;s23_rc(child(traj_id),index)'];
        s33_c=[s33_c;s33_rc(child(traj_id),index)'];
        
        
        ctotalpix=totalpix(child(traj_id),:);
        ctotalpix=ctotalpix(index)';
        ctotalpix=remove_outlayer(ctotalpix,length(ctotalpix));
        c_totalpix=[c_totalpix;ctotalpix];
        
        cxpix=xpix(child(traj_id),:);
        cxpix=cxpix(index)';
        cxpix=remove_outlayer(cxpix,length(cxpix));
        c_xpix=[c_xpix;cxpix];
        
        cypix=ypix(child(traj_id),:);
        cypix=cypix(index)';
        cypix=remove_outlayer(cypix,length(cypix));
        c_ypix=[c_ypix;cypix];
        
        cgrv=grv(child(traj_id),:);
        cgrv=cgrv(index)';
        cgrv=remove_outlayer(cgrv,length(cgrv));
        c_grv=[c_grv;cgrv];
        
        
        cco=co(child(traj_id),:);
        cco=cco(index)';
        cco=remove_outlayer(cco,length(cco));
        c_co=[c_co;cco];
        
        
        cmn_shrt_ax=mn_shrt_ax(child(traj_id),:);
        cmn_shrt_ax=cmn_shrt_ax(index)';
        cmn_shrt_ax=remove_outlayer(cmn_shrt_ax, length(cmn_shrt_ax));
        c_mn_shrt_ax=[c_mn_shrt_ax;cmn_shrt_ax];
        
        cmax_shrt_ax=max_shrt_ax(child(traj_id),:);
        cmax_shrt_ax=cmax_shrt_ax(index)';
        cmax_shrt_ax=remove_outlayer(cmax_shrt_ax,length(cmax_shrt_ax));
        c_max_shrt_ax=[c_max_shrt_ax;cmax_shrt_ax];
        
        cmn_shrt_ax_mm=mn_shrt_ax_mm(child(traj_id),:);
        cmn_shrt_ax_mm=cmn_shrt_ax_mm(index)';
        cmn_shrt_ax_mm=remove_outlayer(cmn_shrt_ax_mm,length(cmn_shrt_ax_mm));
        c_mn_shrt_ax_mm=[c_mn_shrt_ax_mm;cmn_shrt_ax_mm];
        
        cmax_shrt_ax_mm=max_shrt_ax_mm(child(traj_id),:);
        cmax_shrt_ax_mm=cmax_shrt_ax_mm(index)';
        cmax_shrt_ax_mm=remove_outlayer(cmax_shrt_ax_mm,length(cmax_shrt_ax_mm));
        c_max_shrt_ax_mm=[c_max_shrt_ax_mm;cmax_shrt_ax_mm];
        
        cmn_lng_ax=mn_lng_ax(child(traj_id),:);
        cmn_lng_ax= cmn_lng_ax(index)';
        cmn_lng_ax=remove_outlayer(cmn_lng_ax,length(cmn_lng_ax));
        c_mn_lng_ax=[c_mn_lng_ax;cmn_lng_ax];
        
        cmax_lng_ax=max_lng_ax(child(traj_id),:);
        cmax_lng_ax=cmax_lng_ax(index)';
        cmax_lng_ax=remove_outlayer(cmax_lng_ax,length(cmax_lng_ax));
        c_max_lng_ax=[c_max_lng_ax;cmax_lng_ax];
        
        cmn_lng_ax_mm=mn_lng_ax_mm(child(traj_id),:);
        cmn_lng_ax_mm=cmn_lng_ax_mm(index)';
        cmn_lng_ax_mm=remove_outlayer(cmn_lng_ax_mm,length(cmn_lng_ax_mm));
        c_mn_lng_ax_mm=[c_mn_lng_ax_mm;cmn_lng_ax_mm];
        
        cmax_lng_ax_mm=max_lng_ax_mm(child(traj_id),:);
        cmax_lng_ax_mm=cmax_lng_ax_mm(index)';
        cmax_lng_ax_mm=remove_outlayer(cmax_lng_ax_mm,length(cmax_lng_ax_mm));
        c_max_lng_ax_mm=[c_max_lng_ax_mm;cmax_lng_ax_mm];
        
        
    end
    
    %%%% Only smooth
    
    %width=51;%41;
    %order=3;%2;
    
    c_X=smooth(c_X,width,'sgolay',order);
    c_Y=smooth(c_Y,width,'sgolay',order);
    c_Z=smooth(c_Z,width,'sgolay',order);
 
    c_totalpix=smooth(c_totalpix,width,'sgolay',order);
    c_xpix=smooth( c_xpix,width,'sgolay',order);
    c_ypix=smooth( c_ypix,width,'sgolay',order);
    
    c_grv=smooth( c_grv,width,'sgolay',order);
    c_co=smooth( c_co,width,'sgolay',order);
    
    c_mn_shrt_ax=smooth( c_mn_shrt_ax,width,'sgolay',order);
    c_max_shrt_ax=smooth( c_max_shrt_ax,width,'sgolay',order);
    
    c_mn_shrt_ax_mm=smooth( c_mn_shrt_ax_mm,width,'sgolay',order);
    c_max_shrt_ax_mm=smooth( c_max_shrt_ax_mm,width,'sgolay',order);
    
    c_mn_lng_ax=smooth( c_mn_lng_ax,width,'sgolay',order);
    c_max_lng_ax=smooth( c_max_lng_ax,width,'sgolay',order);
    
    c_mn_lng_ax_mm=smooth( c_mn_lng_ax_mm,width,'sgolay',order);
    c_max_lng_ax_mm=smooth( c_max_lng_ax_mm,width,'sgolay',order);
    
    
    c_w1 = smooth(w1_c,width,'sgolay',order);
    c_w2 = smooth(w2_c,width,'sgolay',order);
    c_w3 = smooth(w3_c,width,'sgolay',order);
    
    w11_sg_c=smooth(w11_c,width,'sgolay',order);
    w12_sg_c=smooth(w12_c,width,'sgolay',order);
    w13_sg_c=smooth(w13_c,width,'sgolay',order);
    w21_sg_c=smooth(w21_c,width,'sgolay',order);
    w22_sg_c=smooth(w22_c,width,'sgolay',order);
    w23_sg_c=smooth(w23_c,width,'sgolay',order);
    w31_sg_c=smooth(w31_c,width,'sgolay',order);
    w32_sg_c=smooth(w32_c,width,'sgolay',order);
    w33_sg_c=smooth(w33_c,width,'sgolay',order);
    
    s11_sg_c = smooth(s11_c,width,'sgolay',order);
    s12_sg_c = smooth(s12_c,width,'sgolay',order);
    s13_sg_c = smooth(s13_c,width,'sgolay',order);
    s22_sg_c = smooth(s22_c,width,'sgolay',order);
    s23_sg_c = smooth(s23_c,width,'sgolay',order);
    s33_sg_c = smooth(s33_c,width,'sgolay',order);
    
    
    
    for j=1:length(c_w1)
        
        c_rij(j,1:3,1:3)=[w11_sg_c(j)  w12_sg_c(j) w13_sg_c(j); w21_sg_c(j) w22_sg_c(j) w23_sg_c(j); w31_sg_c(j)  w32_sg_c(j) w33_sg_c(j)];
        c_enstro(j)=c_w1(j)^2+c_w2(j)^2+c_w3(j)^2;
        %A=[s11(i,j) s12(i,j) s13(i,j);s12(i,j) s22(i,j) s23(i,j);s13(i,j) s23(i,j) s33(i,j)];
        %[V,D] = eig(A);
        
        c_Sij(j,1:3,1:3)=[s11_sg_c(j) s12_sg_c(j) s13_sg_c(j);s12_sg_c(j) s22_sg_c(j) s23_sg_c(j);s13_sg_c(j) s23_sg_c(j) s33_sg_c(j)];
        
        c_Aij(j,1:3,1:3)=c_Sij(j,1:3,1:3)+c_rij(j,1:3,1:3);

        
        [V,D] = eig(squeeze(c_Sij(j,1:3,1:3)));
        
        [D,k] = sort(diag(D)); 	% ascending order
        D = diag(D(end:-1:1));  % descending order
        V = V(:,k(end:-1:1)); 	% the same order for eigenvectors
        c_v1(j,1:3)=V(:,1)';
        c_v2(j,1:3)=V(:,2)';
        c_v3(j,1:3)=V(:,3)';
        c_e1(j)=D(1,1);
        c_e2(j)=D(2,2);
        c_e3(j)=D(3,3);
        c_strain(j)=c_e1(j)^2+c_e2(j)^2+c_e3(j)^2;
        % THIS IS DIVERGENCE FREE
        c_rel_div (j)=abs(s11_sg_c(j)+s22_sg_c(j)+s33_sg_c(j))/(abs(s11_sg_c(j))+abs(s22_sg_c(j))+abs(s33_sg_c(j)));
        
        clear V D
    end
    
    
    
    wd= 'E:\OLD-Data-DISK\Matlab for PTV';
    cd([fold,num2str(WF)]);
    
   save children_variables_v3 child c_X c_Y c_Z cTIME...
        c_totalpix c_xpix  c_ypix  c_grv   ...
        c_co c_Aij c_Sij c_rij c_w1 c_w2 c_w3 c_v1 c_v2 c_v3...
        c_e1 c_e2 c_e3 c_strain c_enstro...
        c_mn_shrt_ax    c_max_shrt_ax   c_mn_lng_ax  c_max_lng_ax...
        c_mn_shrt_ax_mm   c_max_shrt_ax_mm  c_mn_lng_ax_mm...
        c_max_lng_ax_mm  c_rel_div
    
    
cd(wd);




%%%%%% Now the non breaking events

 clear index t 


file=dir([fold,num2str(WF),'/nb_variables_*']);
if size(file,1)>0
for num=1:size(file,1)
    
    var=load([fold,num2str(WF),'/nb_variables_',num2str(num)]);
    non_break=var.non_break;
    nb_T=var.nb_T;
    difference=nb_T(2:end)-nb_T(1:end-1); % sometimes a traj begins immediately, but still broken
    index_10=find(difference>1);
    if size(non_break,2)>1 & ~isempty(index_10)
         
        nbX=[];
        nbY=[];
        nbZ=[];
        nbT=[];
        
        nb_totalpix=[];
        nb_xpix=[];
        nb_ypix=[];
        nb_grv=[];
        nb_co=[];
        
        nb_mn_shrt_ax=[];
        nb_mn_shrt_ax_mm=[];
        nb_max_shrt_ax=[];
        nb_max_shrt_ax_mm=[];
        nb_mn_lng_ax=[];
        nb_mn_lng_ax_mm=[];
        nb_max_lng_ax=[];
        nb_max_lng_ax_mm=[];
        
        nb_rel_div=[];
        
        w1_nb=[];
        w2_nb=[];
        w3_nb=[];
        w11_nb=[];
        w12_nb=[];
        w13_nb=[];
        w21_nb=[];
        w22_nb=[];
        w23_nb=[];
        w31_nb=[];
        w32_nb=[];
        w33_nb=[];
        s11_nb=[];
        s12_nb=[];
        s13_nb=[];
        s22_nb=[];
        s23_nb=[];
        s33_nb=[];
        
        for traj_id=1:length(non_break)
            
            t=time(non_break(traj_id),:);
            index=find(t>0);
            
            nbT=[nbT;t(index)'];
            x=tx(non_break(traj_id),:);
            nbX=[nbX;x(index)'];%%% This is a bug remover !!
            y=ty(non_break(traj_id),:);
            nbY=[nbY;y(index)'];
            z=tz(non_break(traj_id),:);
            nbZ=[nbZ;z(index)'];
            
            
            w1_nb=[w1_nb;w1_rc(non_break(traj_id),index)'];
            w2_nb=[w2_nb;w2_rc(non_break(traj_id),index)'];
            w3_nb=[w3_nb;w3_rc(non_break(traj_id),index)'];
            
            w11_nb=[w11_nb;w11_rc(non_break(traj_id),index)'];
            w12_nb=[w12_nb;w12_rc(non_break(traj_id),index)'];
            w13_nb=[w13_nb;w13_rc(non_break(traj_id),index)'];
            w21_nb=[w21_nb;w21_rc(non_break(traj_id),index)'];
            w22_nb=[w22_nb;w22_rc(non_break(traj_id),index)'];
            w23_nb=[w23_nb;w23_rc(non_break(traj_id),index)'];
            w31_nb=[w31_nb;w31_rc(non_break(traj_id),index)'];
            w32_nb=[w32_nb;w32_rc(non_break(traj_id),index)'];
            w33_nb=[w33_nb;w33_rc(non_break(traj_id),index)'];
            
            s11_nb=[s11_nb;s11_rc(non_break(traj_id),index)'];
            s12_nb=[s12_nb;s12_rc(non_break(traj_id),index)'];
            s13_nb=[s13_nb;s13_rc(non_break(traj_id),index)'];
            s22_nb=[s22_nb;s22_rc(non_break(traj_id),index)'];
            s23_nb=[s23_nb;s23_rc(non_break(traj_id),index)'];
            s33_nb=[s33_nb;s33_rc(non_break(traj_id),index)'];
            
            
            nbtotalpix=totalpix(non_break(traj_id),:);
            nbtotalpix=nbtotalpix(index)';
            nbtotalpix=remove_outlayer(nbtotalpix,length(nbtotalpix));
            nb_totalpix=[nb_totalpix;nbtotalpix];
            
            nbxpix=xpix(non_break(traj_id),:);
            nbxpix=nbxpix(index)';
            nbxpix=remove_outlayer(nbxpix,length(nbxpix));
            nb_xpix=[nb_xpix;nbxpix];
            
            nbypix=ypix(non_break(traj_id),:);
            nbypix=nbypix(index)';
            nbypix=remove_outlayer(nbypix,length(nbypix));
            nb_ypix=[nb_ypix;nbypix];
            
            nbgrv=grv(non_break(traj_id),:);
            nbgrv=nbgrv(index)';
            nbgrv=remove_outlayer(nbgrv,length(nbgrv));
            nb_grv=[nb_grv;nbgrv];
            
            
            nbco=co(non_break(traj_id),:);
            nbco=nbco(index)';
            nbco=remove_outlayer(nbco,length(nbco));
            nb_co=[nb_co;nbco];
            
            
            nbmn_shrt_ax=mn_shrt_ax(non_break(traj_id),:);
            nbmn_shrt_ax=nbmn_shrt_ax(index)';
            nbmn_shrt_ax=remove_outlayer(nbmn_shrt_ax, length(nbmn_shrt_ax));
            nb_mn_shrt_ax=[nb_mn_shrt_ax;nbmn_shrt_ax];
            
            nbmax_shrt_ax=max_shrt_ax(non_break(traj_id),:);
            nbmax_shrt_ax=nbmax_shrt_ax(index)';
            nbmax_shrt_ax=remove_outlayer(nbmax_shrt_ax,length(nbmax_shrt_ax));
            nb_max_shrt_ax=[nb_max_shrt_ax;nbmax_shrt_ax];
            
            nbmn_shrt_ax_mm=mn_shrt_ax_mm(non_break(traj_id),:);
            nbmn_shrt_ax_mm=nbmn_shrt_ax_mm(index)';
            nbmn_shrt_ax_mm=remove_outlayer(nbmn_shrt_ax_mm,length(nbmn_shrt_ax_mm));
            nb_mn_shrt_ax_mm=[nb_mn_shrt_ax_mm;nbmn_shrt_ax_mm];
            
            nbmax_shrt_ax_mm=max_shrt_ax_mm(non_break(traj_id),:);
            nbmax_shrt_ax_mm=nbmax_shrt_ax_mm(index)';
            nbmax_shrt_ax_mm=remove_outlayer(nbmax_shrt_ax_mm,length(nbmax_shrt_ax_mm));
            nb_max_shrt_ax_mm=[nb_max_shrt_ax_mm;nbmax_shrt_ax_mm];
            
            nbmn_lng_ax=mn_lng_ax(non_break(traj_id),:);
            nbmn_lng_ax= nbmn_lng_ax(index)';
            nbmn_lng_ax=remove_outlayer(nbmn_lng_ax,length(nbmn_lng_ax));
            nb_mn_lng_ax=[nb_mn_lng_ax;nbmn_lng_ax];
            
            nbmax_lng_ax=max_lng_ax(non_break(traj_id),:);
            nbmax_lng_ax=nbmax_lng_ax(index)';
            nbmax_lng_ax=remove_outlayer(nbmax_lng_ax,length(nbmax_lng_ax));
            nb_max_lng_ax=[nb_max_lng_ax;nbmax_lng_ax];
            
            nbmn_lng_ax_mm=mn_lng_ax_mm(non_break(traj_id),:);
            nbmn_lng_ax_mm=nbmn_lng_ax_mm(index)';
            nbmn_lng_ax_mm=remove_outlayer(nbmn_lng_ax_mm,length(nbmn_lng_ax_mm));
            nb_mn_lng_ax_mm=[nb_mn_lng_ax_mm;nbmn_lng_ax_mm];
            
            nbmax_lng_ax_mm=max_lng_ax_mm(non_break(traj_id),:);
            nbmax_lng_ax_mm=nbmax_lng_ax_mm(index)';
            nbmax_lng_ax_mm=remove_outlayer(nbmax_lng_ax_mm,length(nbmax_lng_ax_mm));
            nb_max_lng_ax_mm=[nb_max_lng_ax_mm;nbmax_lng_ax_mm];
            
            
        end
        
        %%%% Fill and smooth
        
        %width=51;%41;
        %order=3;%2;
        
        
        nbT_fill=fill_time(nbT);
        nbTIME=nbT_fill;%smooth(nbT_fill,width,'sgolay',order);
        nbX=fill_vector(nbT,nbX);
        nbnX=smooth(nbX,width,'sgolay',order);
        nbY=fill_vector(nbT,nbY);
        nbY=smooth(nbY,width,'sgolay',order);
        nbZ=fill_vector(nbT,nbZ);
        nbZ=smooth(nbZ,width,'sgolay',order);
        
        
        
        nb_totalpix=fill_vector(nbT,nb_totalpix);
        nb_totalpix=smooth(nb_totalpix,width,'sgolay',order);
        
        nb_xpix=fill_vector(nbT, nb_xpix);
        nb_xpix=smooth( nb_xpix,width,'sgolay',order);
        
        nb_ypix=fill_vector(nbT, nb_ypix);
        nb_ypix=smooth( nb_ypix,width,'sgolay',order);
        
        nb_grv=fill_vector(nbT, nb_grv);
        nb_grv=smooth( nb_grv,width,'sgolay',order);
        
        nb_co=fill_vector(nbT, nb_co);
        nb_co=smooth( nb_co,width,'sgolay',order);
        
        nb_mn_shrt_ax=fill_vector(nbT, nb_mn_shrt_ax);
        nb_mn_shrt_ax=smooth( nb_mn_shrt_ax,width,'sgolay',order);
        
        nb_max_shrt_ax=fill_vector(nbT, nb_max_shrt_ax);
        nb_max_shrt_ax=smooth( nb_max_shrt_ax,width,'sgolay',order);
        
        nb_mn_shrt_ax_mm=fill_vector(nbT, nb_mn_shrt_ax_mm);
        nb_mn_shrt_ax_mm=smooth( nb_mn_shrt_ax_mm,width,'sgolay',order);
        
        nb_max_shrt_ax_mm=fill_vector(nbT, nb_max_shrt_ax_mm);
        nb_max_shrt_ax_mm=smooth( nb_max_shrt_ax_mm,width,'sgolay',order);
        
        nb_mn_lng_ax=fill_vector(nbT, nb_mn_lng_ax);
        nb_mn_lng_ax=smooth( nb_mn_lng_ax,width,'sgolay',order);
        
        nb_max_lng_ax=fill_vector(nbT, nb_max_lng_ax);
        nb_max_lng_ax=smooth( nb_max_lng_ax,width,'sgolay',order);
        
        nb_mn_lng_ax_mm=fill_vector(nbT, nb_mn_lng_ax_mm);
        nb_mn_lng_ax_mm=smooth( nb_mn_lng_ax_mm,width,'sgolay',order);
        
        nb_max_lng_ax_mm=fill_vector(nbT, nb_max_lng_ax_mm);
        nb_max_lng_ax_mm=smooth( nb_max_lng_ax_mm,width,'sgolay',order);
        
        
        w1_nb=fill_vector(nbT, w1_nb);
        w2_nb=fill_vector(nbT, w2_nb);
        w3_nb=fill_vector(nbT, w3_nb);
        
        
        w11_nb=fill_vector(nbT, w11_nb);
        w12_nb=fill_vector(nbT, w12_nb);
        w13_nb=fill_vector(nbT, w13_nb);
        w21_nb=fill_vector(nbT, w21_nb);
        w22_nb=fill_vector(nbT, w22_nb);
        w23_nb=fill_vector(nbT, w23_nb);
        w31_nb=fill_vector(nbT, w31_nb);
        w32_nb=fill_vector(nbT, w32_nb);
        w33_nb=fill_vector(nbT, w33_nb);
        
        s11_nb=fill_vector(nbT, s11_nb);
        s12_nb=fill_vector(nbT, s12_nb);
        s13_nb=fill_vector(nbT, s13_nb);
        s22_nb=fill_vector(nbT, s22_nb);
        s23_nb=fill_vector(nbT, s23_nb);
        s33_nb=fill_vector(nbT, s33_nb);
        
        
        
        
        nb_w1 = smooth(w1_nb,width,'sgolay',order);
        nb_w2 = smooth(w2_nb,width,'sgolay',order);
        nb_w3 = smooth(w3_nb,width,'sgolay',order);
        
        w11_sg_nb=smooth(w11_nb,width,'sgolay',order);
        w12_sg_nb=smooth(w12_nb,width,'sgolay',order);
        w13_sg_nb=smooth(w13_nb,width,'sgolay',order);
        w21_sg_nb=smooth(w21_nb,width,'sgolay',order);
        w22_sg_nb=smooth(w22_nb,width,'sgolay',order);
        w23_sg_nb=smooth(w23_nb,width,'sgolay',order);
        w31_sg_nb=smooth(w31_nb,width,'sgolay',order);
        w32_sg_nb=smooth(w32_nb,width,'sgolay',order);
        w33_sg_nb=smooth(w33_nb,width,'sgolay',order);
        
        s11_sg_nb = smooth(s11_nb,width,'sgolay',order);
        s12_sg_nb = smooth(s12_nb,width,'sgolay',order);
        s13_sg_nb = smooth(s13_nb,width,'sgolay',order);
        s22_sg_nb = smooth(s22_nb,width,'sgolay',order);
        s23_sg_nb = smooth(s23_nb,width,'sgolay',order);
        s33_sg_nb = smooth(s33_nb,width,'sgolay',order);
        
        
        
        for j=1:length(nb_w1)
            
            nb_rij(j,1:3,1:3)=[w11_sg_nb(j)  w12_sg_nb(j) w13_sg_nb(j); w21_sg_nb(j) w22_sg_nb(j) w23_sg_nb(j); w31_sg_nb(j)  w32_sg_nb(j) w33_sg_nb(j)];
            nb_enstro(j)=nb_w1(j)^2+nb_w2(j)^2+nb_w3(j)^2;
            %A=[s11(i,j) s12(i,j) s13(i,j);s12(i,j) s22(i,j) s23(i,j);s13(i,j) s23(i,j) s33(i,j)];
            %[V,D] = eig(A);
            
            nb_Sij(j,1:3,1:3)=[s11_sg_nb(j) s12_sg_nb(j) s13_sg_nb(j);s12_sg_nb(j) s22_sg_nb(j) s23_sg_nb(j);s13_sg_nb(j) s23_sg_nb(j) s33_sg_nb(j)];
            
            nb_Aij(j,1:3,1:3)=nb_Sij(j,1:3,1:3)+nb_rij(j,1:3,1:3);

            [V,D] = eig(squeeze(nb_Sij(j,1:3,1:3)));
            
            [D,k] = sort(diag(D)); 	% ascending order
            D = diag(D(end:-1:1));  % descending order
            V = V(:,k(end:-1:1)); 	% the same order for eigenvectors
            nb_v1(j,1:3)=V(:,1)';
            nb_v2(j,1:3)=V(:,2)';
            nb_v3(j,1:3)=V(:,3)';
            nb_e1(j)=D(1,1);
            nb_e2(j)=D(2,2);
            nb_e3(j)=D(3,3);
            nb_strain(j)=nb_e1(j)^2+nb_e2(j)^2+nb_e3(j)^2;
            % THIS IS DIVERGENCE FREE
            nb_rel_div (j)=abs(s11_sg_nb(j)+s22_sg_nb(j)+s33_sg_nb(j))/(abs(s11_sg_nb(j))+abs(s22_sg_nb(j))+abs(s33_sg_nb(j)));
            
            clear V D
        end
        
        
        wd= 'E:\OLD-Data-DISK\Matlab for PTV';
        cd([fold,num2str(WF)]);
        
        save (['non_brk_',num2str(num)],'non_break', 'nbX' ,'nbY', 'nbZ' ,'nbTIME',...
            'nb_totalpix' , 'nb_xpix',  'nb_ypix' , 'nb_grv',   ...
            'nb_co','nb_Aij','nb_Sij', 'nb_rij' ,'nb_w1' ,'nb_w2', 'nb_w3', 'nb_v1' ,'nb_v2', 'nb_v3',...
            'nb_e1', 'nb_e2' ,'nb_e3' ,'nb_strain', 'nb_enstro',...
            'nb_mn_shrt_ax' ,   'nb_max_shrt_ax',   'nb_mn_lng_ax',   'nb_max_lng_ax',...
            'nb_mn_shrt_ax_mm',   'nb_max_shrt_ax_mm' , 'nb_mn_lng_ax_mm',...
            'nb_max_lng_ax_mm'  ,'nb_rel_div')
        
        cd(wd);
        
        %end 
        
    else
        
        nbX=[];
        nbY=[];
        nbZ=[];
        nbT=[];
        
        nb_totalpix=[];
        nb_xpix=[];
        nb_ypix=[];
        nb_grv=[];
        nb_co=[];
        
        nb_mn_shrt_ax=[];
        nb_mn_shrt_ax_mm=[];
        nb_max_shrt_ax=[];
        nb_max_shrt_ax_mm=[];
        nb_mn_lng_ax=[];
        nb_mn_lng_ax_mm=[];
        nb_max_lng_ax=[];
        nb_max_lng_ax_mm=[];
        
        nb_rel_div=[];
        
        w1_nb=[];
        w2_nb=[];
        w3_nb=[];
        w11_nb=[];
        w12_nb=[];
        w13_nb=[];
        w21_nb=[];
        w22_nb=[];
        w23_nb=[];
        w31_nb=[];
        w32_nb=[];
        w33_nb=[];
        s11_nb=[];
        s12_nb=[];
        s13_nb=[];
        s22_nb=[];
        s23_nb=[];
        s33_nb=[];
        
        for traj_id=1:length(non_break)
            
            t=time(non_break(traj_id),:);
            index=find(t>0);
            
            nbT=[nbT;t(index)'];
            nbTIME=nbT; 
            x=tx(non_break(traj_id),:);
            nbX=[nbX;x(index)'];%%% This is a bug remover !!
            y=ty(non_break(traj_id),:);
            nbY=[nbY;y(index)'];
            z=tz(non_break(traj_id),:);
            nbZ=[nbZ;z(index)'];
            
            
            w1_nb=[w1_nb;w1_rc(non_break(traj_id),index)'];
            w2_nb=[w2_nb;w2_rc(non_break(traj_id),index)'];
            w3_nb=[w3_nb;w3_rc(non_break(traj_id),index)'];
            
            w11_nb=[w11_nb;w11_rc(non_break(traj_id),index)'];
            w12_nb=[w12_nb;w12_rc(non_break(traj_id),index)'];
            w13_nb=[w13_nb;w13_rc(non_break(traj_id),index)'];
            w21_nb=[w21_nb;w21_rc(non_break(traj_id),index)'];
            w22_nb=[w22_nb;w22_rc(non_break(traj_id),index)'];
            w23_nb=[w23_nb;w23_rc(non_break(traj_id),index)'];
            w31_nb=[w31_nb;w31_rc(non_break(traj_id),index)'];
            w32_nb=[w32_nb;w32_rc(non_break(traj_id),index)'];
            w33_nb=[w33_nb;w33_rc(non_break(traj_id),index)'];
            
            s11_nb=[s11_nb;s11_rc(non_break(traj_id),index)'];
            s12_nb=[s12_nb;s12_rc(non_break(traj_id),index)'];
            s13_nb=[s13_nb;s13_rc(non_break(traj_id),index)'];
            s22_nb=[s22_nb;s22_rc(non_break(traj_id),index)'];
            s23_nb=[s23_nb;s23_rc(non_break(traj_id),index)'];
            s33_nb=[s33_nb;s33_rc(non_break(traj_id),index)'];
            
            
            nbtotalpix=totalpix(non_break(traj_id),:);
            nbtotalpix=nbtotalpix(index)';
            nbtotalpix=remove_outlayer(nbtotalpix,length(nbtotalpix));
            nb_totalpix=[nb_totalpix;nbtotalpix];
            
            nbxpix=xpix(non_break(traj_id),:);
            nbxpix=nbxpix(index)';
            nbxpix=remove_outlayer(nbxpix,length(nbxpix));
            nb_xpix=[nb_xpix;nbxpix];
            
            nbypix=ypix(non_break(traj_id),:);
            nbypix=nbypix(index)';
            nbypix=remove_outlayer(nbypix,length(nbypix));
            nb_ypix=[nb_ypix;nbypix];
            
            nbgrv=grv(non_break(traj_id),:);
            nbgrv=nbgrv(index)';
            nbgrv=remove_outlayer(nbgrv,length(nbgrv));
            nb_grv=[nb_grv;nbgrv];
            
            
            nbco=co(non_break(traj_id),:);
            nbco=nbco(index)';
            nbco=remove_outlayer(nbco,length(nbco));
            nb_co=[nb_co;nbco];
            
            
            nbmn_shrt_ax=mn_shrt_ax(non_break(traj_id),:);
            nbmn_shrt_ax=nbmn_shrt_ax(index)';
            nbmn_shrt_ax=remove_outlayer(nbmn_shrt_ax, length(nbmn_shrt_ax));
            nb_mn_shrt_ax=[nb_mn_shrt_ax;nbmn_shrt_ax];
            
            nbmax_shrt_ax=max_shrt_ax(non_break(traj_id),:);
            nbmax_shrt_ax=nbmax_shrt_ax(index)';
            nbmax_shrt_ax=remove_outlayer(nbmax_shrt_ax,length(nbmax_shrt_ax));
            nb_max_shrt_ax=[nb_max_shrt_ax;nbmax_shrt_ax];
            
            nbmn_shrt_ax_mm=mn_shrt_ax_mm(non_break(traj_id),:);
            nbmn_shrt_ax_mm=nbmn_shrt_ax_mm(index)';
            nbmn_shrt_ax_mm=remove_outlayer(nbmn_shrt_ax_mm,length(nbmn_shrt_ax_mm));
            nb_mn_shrt_ax_mm=[nb_mn_shrt_ax_mm;nbmn_shrt_ax_mm];
            
            nbmax_shrt_ax_mm=max_shrt_ax_mm(non_break(traj_id),:);
            nbmax_shrt_ax_mm=nbmax_shrt_ax_mm(index)';
            nbmax_shrt_ax_mm=remove_outlayer(nbmax_shrt_ax_mm,length(nbmax_shrt_ax_mm));
            nb_max_shrt_ax_mm=[nb_max_shrt_ax_mm;nbmax_shrt_ax_mm];
            
            nbmn_lng_ax=mn_lng_ax(non_break(traj_id),:);
            nbmn_lng_ax= nbmn_lng_ax(index)';
            nbmn_lng_ax=remove_outlayer(nbmn_lng_ax,length(nbmn_lng_ax));
            nb_mn_lng_ax=[nb_mn_lng_ax;nbmn_lng_ax];
            
            nbmax_lng_ax=max_lng_ax(non_break(traj_id),:);
            nbmax_lng_ax=nbmax_lng_ax(index)';
            nbmax_lng_ax=remove_outlayer(nbmax_lng_ax,length(nbmax_lng_ax));
            nb_max_lng_ax=[nb_max_lng_ax;nbmax_lng_ax];
            
            nbmn_lng_ax_mm=mn_lng_ax_mm(non_break(traj_id),:);
            nbmn_lng_ax_mm=nbmn_lng_ax_mm(index)';
            nbmn_lng_ax_mm=remove_outlayer(nbmn_lng_ax_mm,length(nbmn_lng_ax_mm));
            nbarent_mn_lng_ax_mm=[nb_mn_lng_ax_mm;nbmn_lng_ax_mm];
            
            nbmax_lng_ax_mm=max_lng_ax_mm(non_break(traj_id),:);
            nbmax_lng_ax_mm=nbmax_lng_ax_mm(index)';
            nbmax_lng_ax_mm=remove_outlayer(nbmax_lng_ax_mm,length(nbmax_lng_ax_mm));
            nb_max_lng_ax_mm=[nb_max_lng_ax_mm;nbmax_lng_ax_mm];
            
            
        end
        
        %%%% Only smooth
        
        %width=51;%41;
        %order=3;%2;
        
        
        
        nbnX=smooth(nbX,width,'sgolay',order);
        nbY=smooth(nbY,width,'sgolay',order);
        nbZ=smooth(nbZ,width,'sgolay',order);
   
        nb_totalpix=smooth(nb_totalpix,width,'sgolay',order);
        nb_xpix=smooth( nb_xpix,width,'sgolay',order);
        nb_ypix=smooth( nb_ypix,width,'sgolay',order);
        
        nb_grv=smooth( nb_grv,width,'sgolay',order);
        nb_co=smooth( nb_co,width,'sgolay',order);
        
        nb_mn_shrt_ax=smooth( nb_mn_shrt_ax,width,'sgolay',order);
        nb_max_shrt_ax=smooth( nb_max_shrt_ax,width,'sgolay',order);
        
        nb_mn_shrt_ax_mm=smooth( nb_mn_shrt_ax_mm,width,'sgolay',order);
        nb_max_shrt_ax_mm=smooth( nb_max_shrt_ax_mm,width,'sgolay',order);
        
        nb_mn_lng_ax=smooth( nb_mn_lng_ax,width,'sgolay',order);
        nb_max_lng_ax=smooth( nb_max_lng_ax,width,'sgolay',order);
        
        nb_mn_lng_ax_mm=smooth( nb_mn_lng_ax_mm,width,'sgolay',order);
        nb_max_lng_ax_mm=smooth( nb_max_lng_ax_mm,width,'sgolay',order);
        
        nb_w1 = smooth(w1_nb,width,'sgolay',order);
        nb_w2 = smooth(w2_nb,width,'sgolay',order);
        nb_w3 = smooth(w3_nb,width,'sgolay',order);
        
        w11_sg_nb=smooth(w11_nb,width,'sgolay',order);
        w12_sg_nb=smooth(w12_nb,width,'sgolay',order);
        w13_sg_nb=smooth(w13_nb,width,'sgolay',order);
        w21_sg_nb=smooth(w21_nb,width,'sgolay',order);
        w22_sg_nb=smooth(w22_nb,width,'sgolay',order);
        w23_sg_nb=smooth(w23_nb,width,'sgolay',order);
        w31_sg_nb=smooth(w31_nb,width,'sgolay',order);
        w32_sg_nb=smooth(w32_nb,width,'sgolay',order);
        w33_sg_nb=smooth(w33_nb,width,'sgolay',order);
        
        s11_sg_nb = smooth(s11_nb,width,'sgolay',order);
        s12_sg_nb = smooth(s12_nb,width,'sgolay',order);
        s13_sg_nb = smooth(s13_nb,width,'sgolay',order);
        s22_sg_nb = smooth(s22_nb,width,'sgolay',order);
        s23_sg_nb = smooth(s23_nb,width,'sgolay',order);
        s33_sg_nb = smooth(s33_nb,width,'sgolay',order);
        
        for j=1:length(nb_w1)
            
            nb_rij(j,1:3,1:3)=[w11_sg_nb(j)  w12_sg_nb(j) w13_sg_nb(j); w21_sg_nb(j) w22_sg_nb(j) w23_sg_nb(j); w31_sg_nb(j)  w32_sg_nb(j) w33_sg_nb(j)];
            nb_enstro(j)=nb_w1(j)^2+nb_w2(j)^2+nb_w3(j)^2;
            %A=[s11(i,j) s12(i,j) s13(i,j);s12(i,j) s22(i,j) s23(i,j);s13(i,j) s23(i,j) s33(i,j)];
            %[V,D] = eig(A);
            
            nb_Sij(j,1:3,1:3)=[s11_sg_nb(j) s12_sg_nb(j) s13_sg_nb(j);s12_sg_nb(j) s22_sg_nb(j) s23_sg_nb(j);s13_sg_nb(j) s23_sg_nb(j) s33_sg_nb(j)];
            nb_Aij(j,1:3,1:3)=nb_Sij(j,1:3,1:3)+nb_rij(j,1:3,1:3);

            
            [V,D] = eig(squeeze(nb_Sij(j,1:3,1:3)));
            
            [D,k] = sort(diag(D)); 	% ascending order
            D = diag(D(end:-1:1));  % descending order
            V = V(:,k(end:-1:1)); 	% the same order for eigenvectors
            nb_v1(j,1:3)=V(:,1)';
            nb_v2(j,1:3)=V(:,2)';
            nb_v3(j,1:3)=V(:,3)';
            nb_e1(j)=D(1,1);
            nb_e2(j)=D(2,2);
            nb_e3(j)=D(3,3);
            nb_strain(j)=nb_e1(j)^2+nb_e2(j)^2+nb_e3(j)^2;
            % THIS IS DIVERGENCE FREE
            nb_rel_div (j)=abs(s11_sg_nb(j)+s22_sg_nb(j)+s33_sg_nb(j))/(abs(s11_sg_nb(j))+abs(s22_sg_nb(j))+abs(s33_sg_nb(j)));
            
            clear V D
        end
        
        
        wd= 'E:\OLD-Data-DISK\Matlab for PTV';
        cd([fold,num2str(WF)]);
        
        save (['non_brk_',num2str(num)],'non_break', 'nbX' ,'nbY', 'nbZ' ,'nbTIME',...
            'nb_totalpix' , 'nb_xpix',  'nb_ypix' , 'nb_grv',   ...
            'nb_co' ,'nb_Aij','nb_Sij', 'nb_rij' ,'nb_w1' ,'nb_w2', 'nb_w3', 'nb_v1' ,'nb_v2', 'nb_v3',...
            'nb_e1', 'nb_e2' ,'nb_e3' ,'nb_strain', 'nb_enstro',...
            'nb_mn_shrt_ax' ,   'nb_max_shrt_ax',   'nb_mn_lng_ax',   'nb_max_lng_ax',...
            'nb_mn_shrt_ax_mm',   'nb_max_shrt_ax_mm' , 'nb_mn_lng_ax_mm',...
            'nb_max_lng_ax_mm'  ,'nb_rel_div')
        
        cd(wd);
    end
    
    
end
end

clearvars -except wf max_length order width already_there count long_traj_count
end




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



% 
% %%%%%%%%%%%%%%%%%%%%%%%%% Figure section %%%%%%%%%%%%%%%%%%%%%%% 
% 
% sep=zeros(size(child,2),500,7);
% a=1;    
% parent_time=T(T>0);
% X=X(abs(X)>0);Y=Y(abs(Y)>0);Z=Z(abs(Z)>0);
% for j=1:size(child,2) % j is/are the trajectory other than the parent one
%     begi_fr=ctime{j}(1);
%     if parent_time(end)>ctime{j}(end)
%         end_fr=ctime{j}(end);
%     else
%         end_fr=parent_time(end);
%     end
%     frame_length=end_fr-begi_fr+1;
%     cc=0;
%     for gg=1:frame_length;
%         %e=(begi_fr+gg-a)-parent_time(1)+1;% This was for ideal case when
%                                            % parent traj was not fragmented anywhere
%         e=find(ctime{j}(gg)==T);%This for when the parent line is fragmented either before or after breakage!
%         if ~isempty(e)% This is if the parent line is fragmented after breakage frame
%             cc=cc+1;
%         sep(j,cc,1:7)=[X(e) cx{j}(cc) Y(e) cy{j}(cc) Z(e) cz{j}(cc)...
%                        sqrt((X(e)-cx{j}(cc)).^2+(Y(e)-cy{j}(cc)).^2+(Z(e)-cz{j}(cc)).^2)];
%         end
%         %hold on
%         %plot3([X(e) cx{j}(gg)], [Y(e) cy{j}(gg)] ,[Z(e) cz{j}(gg)],'b')
%     end
% end
% 
% %save sep
% 
% 
% 
% %figure;
% %threshold=(mn_shrt_ax_mm(par,1)/2)*1e-3;
% threshold=(parent_max_lng_ax_mm(1)/2)*1e-3;
% for hh=1:size(child,2)
% 
%     sep_dist=sep(hh,:,7);
%     sep_dist=sep_dist(sep_dist>0);
%     if length(sep_dist)>35
%         r=sep_dist(1:35);
%         t=ctime{hh}(1:35);
%     else
%         r=sep_dist;
%         t=ctime{hh}(1:length(r));
%     end
%     
%     p=polyfit(t,r,1);
%     break_frame_X1(hh)=round((threshold-p(2))/p(1));
%    
%     if  break_frame_X1(hh)>ctime{j}(1)  
%     break_frame_X1(hh)=ctime{j}(1);
%     end
%     
% %      hold on
% %      plot(t,r)
% %      clear sep_dist r t
% end
% % xlabel('Frame number')
% % ylabel('r (m)')
% 
% 
% 
% %%% according to balastic regime.......
% threshold=(parent_max_lng_ax_mm(1)/2)*1e-3;
% epsi=9e-5;
% K=11/3;
% C2=2.13;
% parent_time=T(T>0);
% X=X(abs(X)>0);Y=Y(abs(Y)>0);Z=Z(abs(Z)>0);
% for j=1:size(child,2) % j is/are the trajectory other than the parent one
%     
%     %begi_fr=ctime{j}(1);
%     %e=begi_fr-parent_time(1)+1;
%     e=find(ctime{j}(1)==T);
%     if isempty(e)
%         diference=abs(T-ctime{j}(1));
%         a=find(diference==min(diference));
%         e=a;
%     end
%         
%     del_r(j)=sqrt((X(e)-cx{j}(1)).^2+(Y(e)-cy{j}(1)).^2+(Z(e)-cz{j}(1)).^2);
%     tt(j)=ctime{j}(1)-round(sqrt((del_r(j)-threshold)^2/(K*C2*(epsi*threshold)^(2/3)))*250);
%     if tt(j)>0
%         index=find(T==tt(j));
%         if~isempty(index)%% This is when the formula does not work !!!
%                          %%% i.e., it goes back beyond the parent traj
%                          %%% begining time(Example WF 10)
%             break_frame_X2(j)=T(T==tt(j));
%         
%      else
%             break_frame_X2(j)=ctime{j}(1);
%         end
%     end
%     %hold on
%     %plot3([X(e) cx{j}(1)], [Y(e) cy{j}(1)] ,[Z(e) cz{j}(1)],'r')
% end
% 
% 
% 
% % or take as it is for breakage frame
% for kk=1:size(child,2)
% break_frame_ptv(kk)=time(child(kk),1);
% end
% 
% save break_frame break_frame_ptv break_frame_X1 break_frame_X2
% 
% 
% 
% cd(wd)
% 
% for ee=1:size(child,2)
% figure;
% ptv_fr=break_frame_ptv(ee);
% %ptv_fr=break_frame_X2(ee);
% 
% subplot(2,2,1)
% imshow([fold2,'image_agg\','filt_Cam1.',num2str(ptv_fr)])
% title(['breakage at frame ',num2str(ptv_fr)])
% subplot(2,2,2)
% imshow([fold2,'image_agg\','filt_Cam2.',num2str(ptv_fr)])
% title(['breakage at frame ',num2str(ptv_fr)])
% subplot(2,2,4)
% imshow([fold2,'image_agg\','filt_Cam3.',num2str(ptv_fr)])
% title(['breakage at frame ',num2str(ptv_fr)])
% subplot(2,2,3)
% imshow([fold2,'image_agg\','filt_Cam4.',num2str(ptv_fr)])
% title(['breakage at frame ',num2str(ptv_fr)])
% end
% 
% 
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 





% 
% %%% This is Beat's figure part------------------
% %----------------------------------------------------
% a=146;b=-7; %looking down along omega 
% %a=147;b=16; %visible separation
% 
% figure;hold on;
% ind_a=find(xa<1e9);
% ind_b=find(xa(ind_a)<0.0125);
% ind_c=find(za(ind_a(ind_b))>0.002);
% ind=ind_a(ind_b(ind_c));
% scatter3(xa(ind),ya(ind),za(ind),30,'b','filled')
% xlabel('x')
% ylabel('y')
% zlabel('z')
% view(a,b)
% axis equal;
% box on;
% 
% for i=1:tid
%     count=0;
%     for j=1:traj_len(i)
%         if 2*strain(i,j)<1000 & enstro(i,j)<1000 & strain(i,j)>0
%             count=count+1;
%             tx(count)=agg(traj(i,j,1),traj(i,j,2),1);
%             ty(count)=agg(traj(i,j,1),traj(i,j,2),2);
%             tz(count)=agg(traj(i,j,1),traj(i,j,2),3);
%             ve1(count,1:3)=v1(i,j,1:3);
%             ve2(count,1:3)=v2(i,j,1:3);
%             ve3(count,1:3)=v3(i,j,1:3);
%             ei1(count)=e1(i,j);
%             ei2(count)=e2(i,j);
%             ei3(count)=e3(i,j);
%             wo1(count)=w1(i,j);
%             wo2(count)=w2(i,j);
%             wo3(count)=w3(i,j);
%             st(count)=strain(i,j);
%             en(count)=enstro(i,j);
%         end
%     end
%     figure(10+offset);hold on;
%     scatter3(tx(1:count),ty(1:count),st(1:count),15,'r','filled')
%     xlabel('x')
%     ylabel('y')
%     zlabel('s^2')
%     view(a,b)
%     box on;
%     figure(11+offset);hold on;
%     scatter3(tx(1:count),ty(1:count),en(1:count),15,'r','filled')
%     xlabel('x')
%     ylabel('y')
%     zlabel('\omega^2')
%     view(a,b)
%     box on;
%     figure(12+offset);hold on;
%     quiver3(tx(1:count),ty(1:count),tz(1:count),ve1(1:count,1),ve1(1:count,2),ve1(1:count,3),1,'r')
%     quiver3(tx(1:count),ty(1:count),tz(1:count),-ve1(1:count,1),-ve1(1:count,2),-ve1(1:count,3),1,'r')
%     quiver3(tx(1:count),ty(1:count),tz(1:count),ve2(1:count,1),ve2(1:count,2),ve2(1:count,3),1,'g')
%     quiver3(tx(1:count),ty(1:count),tz(1:count),-ve2(1:count,1),-ve2(1:count,2),-ve2(1:count,3),1,'g')
%     xlabel('x')
%     ylabel('y')
%     zlabel('z')
%     title('\lambda_1 (red), \lambda_2 (green)')
%     view(a,b)
%     axis equal;
%     box on;
%     figure(13+offset);hold on;
%     quiver3(tx(1:count),ty(1:count),tz(1:count),ve3(1:count,1),ve3(1:count,2),ve3(1:count,3),1,'b')
%     quiver3(tx(1:count),ty(1:count),tz(1:count),-ve3(1:count,1),-ve3(1:count,2),-ve3(1:count,3),1,'b')
%     xlabel('x')
%     ylabel('y')
%     zlabel('z')
%     title('\lambda_3')
%     view(a,b)
%     axis equal;
%     box on;
%     figure(14+offset);hold on;
%     quiver3(tx(1:count),ty(1:count),tz(1:count),wo1(1:count),wo2(1:count),wo3(1:count),2,'m')
%     xlabel('x')
%     ylabel('y')
%     zlabel('z')
%     title('\omega')
%     view(a,b)
%     axis equal;
%     box on;
%     figure(15+offset);hold on;
%     scatter3(tx(1:count),ty(1:count),ei1(1:count),15,'k','filled')
%     xlabel('x')
%     ylabel('y')
%     zlabel('z')
%     title('\Lambda_1')
%     view(a,b)
%     box on;
% end
% 
% 
% 
% 
% 
% 
% 
