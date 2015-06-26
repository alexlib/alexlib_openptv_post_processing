clear all;
close all;
first=100;
last=410;
radius_sphere=0.006;
tol=2e-3;%tol=0.5e-3

offset=10;

long_traj_count=0;
max_length=0;

agg=zeros(last-first+1,100,8);
xu=zeros(last-first+1,2000,7);


%read and draw aggregate
name_root=['D:\Turbulent_data_backup\flow_induced_aggregates_breakage_at_100rpm\WF6\res_agg/'];
for i=first:last

    name=[name_root,'rt_is_v2.',num2Str(i)];
    fid = fopen(name, 'r');
    num_points = fscanf(fid, '%i', [1 1]);    % It has two rows now.
    tmp = fscanf(fid, '%i %f %f %f %i %i %i %i %f %f %f %f ', [12 num_points]);
    A=tmp';
    A(:,2:4)=A(:,2:4)*0.001;
    agg(i-first+1,1:num_points,1:7)=A(:,[2 3 4 9 10 11 12]);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fclose(fid);
end

xa=[];
ya=[];
za=[];
ta=[];
for i=1:last-first+1
    ind=find(abs(agg(i,:,1))>0 & agg(i,:,2)>0.001);
    xa=[xa;agg(i,ind,1)'];
    ya=[ya;agg(i,ind,2)'];
    za=[za;agg(i,ind,3)'];
    t=repmat(i,length(ind),1);
    ta=[ta;t];
end
figure;hold on;
scatter3(xa,ya,za,30,'b','filled')

%track aggregates using nearest neighbour, 'tol' tolerance
tid=0;
for i=1:last-first
    [2 i]
    ind=find(abs(agg(i,:,1))>0 & agg(i,:,2)>0.001);
    ind2=find(agg(i,ind,8)==0);
    ind=ind(ind2);
    num_points=length(ind);
    for j=1:num_points
        tid=tid+1;%trajectory id
        pid=1;
        traj(tid,pid,1)=i;%frame id
        traj(tid,pid,2)=ind(j);%pos in frame
        agg(traj(tid,pid,1),traj(tid,pid,2),8)=1;%occupied as of now
        %find next ?
        end_of_traj=0;
        while end_of_traj==0
            pos_agg=[agg(traj(tid,pid,1),traj(tid,pid,2),1) agg(traj(tid,pid,1),traj(tid,pid,2),2) agg(traj(tid,pid,1),traj(tid,pid,2),3)];
            ind2=find(abs(agg(traj(tid,pid,1)+1,:,1))>0 & abs(agg(traj(tid,pid,1)+1,:,2))>0.001);
            pos=repmat(pos_agg,length(ind2),1);
            if length(ind2)>1
                delta=squeeze(agg(traj(tid,pid,1)+1,ind2,1:3))-pos;
                dist=sum(delta.^2,2).^0.5;
            else
                delta=squeeze(agg(traj(tid,pid,1)+1,ind2,1:3))-pos';
                dist=sum(delta.^2,1).^0.5;
            end
            ind3=find(dist<tol);
            if ~isempty(ind3)
                ind2=ind2(ind3);
            else
                ind2=[];
            end
            if ~isempty(ind2) && traj(tid,pid,1)<last-first
                pid=pid+1;
                traj(tid,pid,1)=traj(tid,pid-1,1)+1;
                traj(tid,pid,2)=ind2(1);
                agg(traj(tid,pid,1),traj(tid,pid,2),8)=1;%occupied as of now
            else
                end_of_traj=1;
                if pid>10 %pid=3
                    long_traj_count=long_traj_count+1;
                    traj_len(tid)=pid;
                    if pid>max_length
                        max_length=pid;
                    end
                    %clear tx ty tz
                    for k=1:pid
                        tx(tid,k)=agg(traj(tid,k,1),traj(tid,k,2),1);%%% tx(k)
                        ty(tid,k)=agg(traj(tid,k,1),traj(tid,k,2),2);
                        tz(tid,k)=agg(traj(tid,k,1),traj(tid,k,2),3);
                        ttotalpix(tid,k)=agg(traj(tid,k,1),traj(tid,k,2),4);
                        txpix(tid,k)=agg(traj(tid,k,1),traj(tid,k,2),5);
                        typix(tid,k)=agg(traj(tid,k,1),traj(tid,k,2),6);
                        tgrv(tid,k)=agg(traj(tid,k,1),traj(tid,k,2),7);
                    end
                    scatter3(tx,ty,tz,'r');
                else
                    tid=tid-1;
                end
            end
        end
    end
end

figure(400);
hold on;
x1=tx(1,:);
x1=x1(x1~=0);

y1=ty(1,:);
y1=y1(y1~=0);

z1=tz(1,:);
z1=z1(z1~=0);

scatter3(x1,y1,z1,'r')

x2=tx(2,:);
x2=x2(x2~=0);

y2=ty(2,:);
y2=y2(y2~=0);

z2=tz(2,:);
z2=z2(z2~=0);

scatter3(x2,y2,z2,'g')

x3=tx(3,:);
x3=x3(x3~=0);

y3=ty(3,:);
y3=y3(y3~=0);

z3=tz(3,:);
z3=z3(z3~=0);

scatter3(x3,y3,z3,'b')
hold off

traj1_totalpix=ttotalpix(1,:);
traj1_totalpix=traj1_totalpix(traj1_totalpix~=0);
traj2_totalpix=ttotalpix(2,:);
traj2_totalpix=traj2_totalpix(traj2_totalpix~=0);
traj3_totalpix=ttotalpix(3,:);
traj3_totalpix=traj3_totalpix(traj3_totalpix~=0);

traj1=traj(1,:,1);
traj1=traj1(traj1>0);
traj2=traj(2,:,1);
traj2=traj2(traj2>0);
traj3=traj(3,:,1);
traj3=traj3(traj3>0);
figure(500);
hold on
scatter(traj1,traj1_totalpix,'r')
scatter(traj2,traj2_totalpix,'g')
scatter(traj3,traj3_totalpix,'b')
hold off

long_traj_count=long_traj_count

%read flow tracer field
%pos_nozzle=[1.6795e-002  1.671e-002 -5.3560e-003];
name_root=['G:\flow_induced_aggregates_breakage_at_100rpm\WF6\res/'];
for i=first:last
    if mod(i-first,100)==0
        [3 i]
    end
    name=[name_root,'xuap.',num2Str(i)];
    clear A;
    A=load(name);
    si=size(A);
    num_points=si(1,1);
    xu(i-first+1,1:num_points,1:6)=A(:,6:11);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    xu(i-first+1,1:num_points,7)=A(:,15);
end

x=zeros(tid,max_length,100);
y=zeros(tid,max_length,100);
z=zeros(tid,max_length,100);
Aij=zeros(tid,max_length,3,3);
co=zeros(tid,max_length,1);

%get Aij along all aggregate trajectories
for i=1:tid
    for j=1:traj_len(i)
        %find swarm
        ind=find(xu(traj(i,j,1),:,7)==1);
        pos_agg=[agg(traj(i,j,1),traj(i,j,2),1) agg(traj(i,j,1),traj(i,j,2),2) agg(traj(i,j,1),traj(i,j,2),3)];
        pos=repmat(pos_agg,length(ind),1);
        delta=squeeze(xu(traj(i,j,1),ind,1:3))-pos;
        dist=sum(delta.^2,2).^0.5;
        ind2=find(dist<=radius_sphere);
        ind=ind(ind2);
        co(i,j)=length(ind);% this tells me that how many points there are around each point on a trajectory

        if length(ind)>=4
            %get Aij around this point
            Gij=get_G_ij_linear(xu(traj(i,j,1),ind,1),xu(traj(i,j,1),ind,2),xu(traj(i,j,1),ind,3),xu(traj(i,j,1),ind,4),xu(traj(i,j,1),ind,5),xu(traj(i,j,1),ind,6));
            Aij(i,j,1:3,1:3)=Gij(1:3,1:3);
        else
            Aij(i,j,1:3,1:3)=zeros(3,3);
        end
    end
    %compute omega, sij and eigen-stuff
    for j=1:traj_len(i)
        w1(i,j)=Aij(i,j,3,2)-Aij(i,j,2,3);%G(3,2)-G(2,3);%
        w2(i,j)=Aij(i,j,1,3)-Aij(i,j,3,1);%G(1,3)-G(3,1);%
        w3(i,j)=Aij(i,j,2,1)-Aij(i,j,1,2);%G(2,1)-G(1,2);%
        enstro(i,j)=w1(i,j)^2+w2(i,j)^2+w3(i,j)^2;
        s11(i,j)=0.5*(Aij(i,j,1,1)+Aij(i,j,1,1));%(G(1,1)+G(1,1));%
        s12(i,j)=0.5*(Aij(i,j,1,2)+Aij(i,j,2,1));%(G(1,2)+G(2,1));%
        s13(i,j)=0.5*(Aij(i,j,1,3)+Aij(i,j,3,1));%(G(1,3)+G(3,1));%
        s22(i,j)=0.5*(Aij(i,j,2,2)+Aij(i,j,2,2));%(G(2,2)+G(2,2));%
        s23(i,j)=0.5*(Aij(i,j,2,3)+Aij(i,j,3,2));%(G(2,3)+G(3,2));%
        s33(i,j)=0.5*(Aij(i,j,3,3)+Aij(i,j,3,3));%(G(3,3)+G(3,3));%

        A=[s11(i,j) s12(i,j) s13(i,j);s12(i,j) s22(i,j) s23(i,j);s13(i,j) s23(i,j) s33(i,j)];
        [V,D] = eig(A);
        [D,k] = sort(diag(D)); 	% ascending order
        D = diag(D(end:-1:1)); 			% descending order
        V = V(:,k(end:-1:1)); 	% the same order for eigenvectors
        v1(i,j,1:3)=V(:,1)';
        v2(i,j,1:3)=V(:,2)';
        v3(i,j,1:3)=V(:,3)';
        e1(i,j)=D(1,1);
        e2(i,j)=D(2,2);
        e3(i,j)=D(3,3);
        strain(i,j)=e1(i,j)^2+e2(i,j)^2+e3(i,j)^2;
    end
    %remove crap
    w1(i,:)=remove_outlayer(w1(i,:),traj_len(i));
    w2(i,:)=remove_outlayer(w2(i,:),traj_len(i));
    w3(i,:)=remove_outlayer(w3(i,:),traj_len(i));
    s11(i,:)=remove_outlayer(s11(i,:),traj_len(i));
    s12(i,:)=remove_outlayer(s12(i,:),traj_len(i));
    s13(i,:)=remove_outlayer(s13(i,:),traj_len(i));
    s22(i,:)=remove_outlayer(s22(i,:),traj_len(i));
    s23(i,:)=remove_outlayer(s23(i,:),traj_len(i));
    s33(i,:)=remove_outlayer(s33(i,:),traj_len(i));
    
    %low pass filter
    width=41;
    order=2;
    w1(i,:) = smooth(squeeze(w1(i,:)),width,'sgolay',order);
    w2(i,:) = smooth(squeeze(w2(i,:)),width,'sgolay',order);
    w3(i,:) = smooth(squeeze(w3(i,:)),width,'sgolay',order);
    s11(i,:) = smooth(squeeze(s11(i,:)),width,'sgolay',order);
    s12(i,:) = smooth(squeeze(s12(i,:)),width,'sgolay',order);
    s13(i,:) = smooth(squeeze(s13(i,:)),width,'sgolay',order);
    s22(i,:) = smooth(squeeze(s22(i,:)),width,'sgolay',order);
    s23(i,:) = smooth(squeeze(s23(i,:)),width,'sgolay',order);
    s33(i,:) = smooth(squeeze(s33(i,:)),width,'sgolay',order);
    for j=1:traj_len(i)
        enstro(i,j)=w1(i,j)^2+w2(i,j)^2+w3(i,j)^2;
        A=[s11(i,j) s12(i,j) s13(i,j);s12(i,j) s22(i,j) s23(i,j);s13(i,j) s23(i,j) s33(i,j)];
        [V,D] = eig(A);
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
    end
end

save flow_induced_aggregates_breakage_at_100rpm_WF6 xa ya za ta co Aij w1 w2 w3 v1 v2 v3 e1 e2 e3 strain enstro traj tx ty tz ttotalpix txpix typix tgrv



% Time series of Strain
strain_on_traj1=strain(1,:);
strain_on_traj1=strain_on_traj1(strain_on_traj1>0);
traj1=traj(1,:,1);
traj1=traj1(traj1>0);

strain_on_traj2=strain(2,:);
strain_on_traj2=strain_on_traj2(strain_on_traj2>0);
traj2=traj(2,:,1);
traj2=traj2(traj2>0);

strain_on_traj3=strain(3,:);
strain_on_traj3=strain_on_traj3(strain_on_traj3>0);
traj3=traj(3,:,1);
traj3=traj3(traj3>0);

tau_eta=sqrt(1/(2*mean([strain_on_traj1 strain_on_traj2 strain_on_traj3])));

figure(1000);
hold on;
scatter(traj1,strain_on_traj1*tau_eta^2,'r');
scatter(traj2,strain_on_traj2*tau_eta^2,'g');
scatter(traj3,strain_on_traj3*tau_eta^2,'b');
hold off

xlabel('\tau_\eta')
ylabel('2S^2_{ij}\cdot\tau_\eta^2')


% Time series of Enstrophy

enstro_on_traj1=enstro(1,:);
enstro_on_traj1=enstro_on_traj1(enstro_on_traj1>0);
traj1=traj(1,:,1);
traj1=traj1(traj1>0);

enstro_on_traj2=enstro(2,:);
enstro_on_traj2=enstro_on_traj2(enstro_on_traj2>0);
traj2=traj(2,:,1);
traj2=traj2(traj2>0);

enstro_on_traj3=enstro(3,:);
enstro_on_traj3=enstro_on_traj3(enstro_on_traj3>0);
traj3=traj(3,:,1);
traj3=traj3(traj3>0);

figure(2000);
hold on;
scatter(traj1,enstro_on_traj1,'r');
scatter(traj2,enstro_on_traj2,'g');
scatter(traj3,enstro_on_traj3,'b');
hold off

xlabel('\tau_\eta')
ylabel('\omega^2')


% Angle between eigen vector and trajectory

traj_3=[x3(1) y3(1) z3(1)];

lambda1x=v1(3,:,1);
lambda1x=lambda1x(lambda1x~=0);
lambda1y=v1(3,:,2);
lambda1y=lambda1y(lambda1y~=0);
lambda1z=v1(3,:,3);
lambda1z=lambda3z(lambda1z~=0);

lambda2x=v2(3,:,1);
lambda2x=lambda2x(lambda2x~=0);
lambda2y=v2(3,:,2);
lambda2y=lambda2y(lambda2y~=0);
lambda2z=v2(3,:,3);
lambda2z=lambda2z(lambda2z~=0);

lambda3x=v3(3,:,1);
lambda3x=lambda3x(lambda3x~=0);
lambda3y=v3(3,:,2);
lambda3y=lambda3y(lambda3y~=0);
lambda3z=v3(3,:,3);
lambda3z=lambda3z(lambda3z~=0);

lambda1=[lambda1x(1) lambda1y(1) lambda1z(1) ];
lambda2=[lambda2x(1) lambda2y(1) lambda2z(1) ];
lambda3=[lambda3x(1) lambda3y(1) lambda3z(1) ];


cos_traj3_lambda1=(traj_3(1).*lambda1(1)+traj_3(2).*lambda1(2)+traj_3(3).*lambda1(3))...
    ./(sqrt(traj_3(1)^2+traj_3(2)^2+traj_3(3)^3).*sqrt(lambda1(1)^2+lambda1(2)^2+lambda1(3)^2));


cos_traj3_lambda2=(traj_3(1).*lambda2(1)+traj_3(2).*lambda2(2)+traj_3(3).*lambda2(3))...
    ./(sqrt(traj_3(1)^2+traj_3(2)^2+traj_3(3)^3).*sqrt(lambda2(1)^2+lambda2(2)^2+lambda2(3)^2));

cos_traj3_lambda3=(traj_3(1).*lambda3(1)+traj_3(2).*lambda3(2)+traj_3(3).*lambda3(3))...
    ./(sqrt(traj_3(1)^2+traj_3(2)^2+traj_3(3)^3).*sqrt(lambda3(1)^2+lambda3(2)^2+lambda3(3)^2));






%Separation vector
all_points=[xa ya za];
rr=[(all_points(2:end,1)-all_points(1:end-1,1))  (all_points(2:end,2)-all_points(1:end-1,2)) (all_points(2:end,3)-all_points(1:end-1,3))];





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
    quiver3(tx(1:count),ty(1:count),tz(1:count),ve1(1:count,1)',ve1(1:count,2)',ve1(1:count,3)',1,'r')
    quiver3(tx(1:count),ty(1:count),tz(1:count),-ve1(1:count,1)',-ve1(1:count,2)',-ve1(1:count,3)',1,'r')
    quiver3(tx(1:count),ty(1:count),tz(1:count),ve2(1:count,1)',ve2(1:count,2)',ve2(1:count,3)',1,'g')
    quiver3(tx(1:count),ty(1:count),tz(1:count),-ve2(1:count,1)',-ve2(1:count,2)',-ve2(1:count,3)',1,'g')
    xlabel('x')
    ylabel('y')
    zlabel('z')
    title('\lambda_1 (red), \lambda_2 (green)')
    view(a,b)
    axis equal;
    box on;
    figure(13+offset);hold on;
    quiver3(tx(1:count),ty(1:count),tz(1:count),ve3(1:count,1)',ve3(1:count,2)',ve3(1:count,3)',1,'b')
    quiver3(tx(1:count),ty(1:count),tz(1:count),-ve3(1:count,1)',-ve3(1:count,2)',-ve3(1:count,3)',1,'b')
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





