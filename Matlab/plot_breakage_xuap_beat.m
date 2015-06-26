first=300;
last=400;
radius_sphere=0.007;
tol=1.0e-3;

offset=10;

long_traj_count=0;
max_length=10;

agg=zeros(last-first+1,100,4);
xu=zeros(last-first+1,2000,7);

%read and draw aggregate
name_root=['G:\flow_induced_aggregates_breakage_at_100rpm\WF6\res_agg\'];
for i=first:last
    if mod(i-first,100)==0
        [1 i]
    end
    name=[name_root,'rt_is.',num2Str(i)];
    fid = fopen(name, 'r');
    num_points = fscanf(fid, '%i', [1 1]);    % It has two rows now.
    tmp = fscanf(fid, '%i %f %f %f %i %i %i %i', [8 num_points]);
    A=tmp';
    A(:,2:4)=A(:,2:4)*0.001;
    agg(i-first+1,1:num_points,1:3)=A(:,2:4);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
% figure;%hold on;
% scatter3(xa,ya,za,30,'b','filled')


figure;hold on;
%track aggregates using nearest neighbour, 'tol' tolerance
tid=0;
for i=1:last-first
    [2 i]
    ind=find(abs(agg(i,:,1))>0 & agg(i,:,2)>0.001);
    ind2=find(agg(i,ind,4)==0);
    ind=ind(ind2);
    num_points=length(ind);
    for j=1:num_points
        tid=tid+1;
        pid=1;
        traj(tid,pid,1)=i;%frame id
        traj(tid,pid,2)=ind(j);%pos in frame
        agg(traj(tid,pid,1),traj(tid,pid,2),4)=1;%occupied as of now
        %find next ?
        end_of_traj=0;
        while end_of_traj==0
            pos_agg=[agg(traj(tid,pid,1),traj(tid,pid,2),1) agg(traj(tid,pid,1),traj(tid,pid,2),2) agg(traj(tid,pid,1),traj(tid,pid,2),3)];
            ind2=find(abs(agg(traj(tid,pid,1)+1,:,1))>0 & abs(agg(traj(tid,pid,1)+1,:,2))>0.001 & agg(traj(tid,pid,1)+1,:,4)<1);
            pos=repmat(pos_agg,length(ind2),1);
            if length(ind2>0)
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
                traj(tid,pid,2)=ind2(1);
                agg(traj(tid,pid,1),traj(tid,pid,2),4)=1;%occupied as of now
            else
                end_of_traj=1;
                if pid>20
                    long_traj_count=long_traj_count+1;
                    traj_len(tid)=pid;
                    if pid>max_length
                        max_length=pid;
                    end
                    clear tx ty tz
                    for k=1:pid
                        tx(k)=agg(traj(tid,k,1),traj(tid,k,2),1);
                        ty(k)=agg(traj(tid,k,1),traj(tid,k,2),2);
                        tz(k)=agg(traj(tid,k,1),traj(tid,k,2),3);
                    end
                    hold on
                    scatter3(tx,ty,tz,'k.');
                else
                    tid=tid-1;
                end
            end
        end
    end
end
long_traj_count=long_traj_count

%read flow tracer field
%pos_nozzle=[1.6795e-002  1.671e-002 -5.3560e-003];
name_root=['\\Ifu-gwh-disk\hydromechanik3\Beat_Colloid/WF1/res/'];
for i=first:last
    if mod(i-first,100)==0
        [3 i]
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

%get Aij along all aggregate trajectories
for i=1:tid
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

save breakage_case_50rpm_Sept_2011 xa ya za ta co Aij w1 w2 w3 v1 v2 v3 e1 e2 e3 strain enstro

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




