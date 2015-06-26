

clear all
clear all
clear all
%close all

%WF=[1 2 6 10 11 12 13 15 16 18 21 22 23 25 26 27 28 29 30 31 32 33 34 35 36 37 38 41 42 45 46 50 51 52 53 55 56 58 59 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 79 80 81 82 83 84 86 87 88 89 90 92 93 95 97 98  99 100 101 102 103 104 105 106 107 108 109 110 111 112 113];

WF=[2  ];
%fold='C:\ETH\Debashish\data';
fold='D:\Turbulent_data_backup\flow_induced_aggregates_breakage_at_100rpm\WF';


for j=1:1%length(WF)
    
    % par_var=[fold,'/parent_variables'];
    % child_var=[fold,'/children_variables'];
    % break_var=[fold,'/break_frame'];
    
    par_var=[fold,num2str(WF(j)),'/parent_variables'];
    child_var=[fold,num2str(WF(j)),'/children_variables'];
    break_var=[fold,num2str(WF(j)),'/break_frame'];
    
    
    parent=load(par_var);
    children=load(child_var);
    break_frame=load(break_var);
    
    break_fr_ptv=break_frame.break_frame_ptv;
    break_fr_X1=break_frame.break_frame_X1;
    break_fr_X2=break_frame.break_frame_X2;
    
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
    
    figure(1);
    plot3(x_,y_,z_)
    title('Parent trajectory')
    
    
    c_time=1;
    dt=1/450;
    B(c_time,1:3,1:3)=[1 0 0;0 1 0;0 0 1];
    B(c_time,1:3,1:3)=B(c_time,1:3,1:3)*1e-3;
    
    int_t0_tb=break_fr_ptv(1)-time(1)+1;
    int_window=10*35;% 2*tau_B
    frame_break=break_fr_ptv(1)-time(1)+1;
    frame_0=frame_break-int_window+1;
    
    if int_window>int_t0_tb
        int_window=int_t0_tb;
        frame_0=1;
    end
    % Integration over time
    for frame=frame_0:frame_break
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
    
    cX1=children.cx{1}(1); %%% Children coordinate
    cY1=children.cy{1}(1);
    cZ1=children.cz{1}(1);
    
    cX2=parent.X(frame_break);
    cY2=parent.Y(frame_break);
    cZ2=parent.Z(frame_break);
    
    r=[cX1-cX2;cY1-cY2;cZ1-cZ2];
    r=r/sum(r.^2)^0.5; %%% Separation @ breakage between parent and fragment
    
    a=W1(end,:);
    cos_W1_r(j)=dot(a',r);%%% Angle between separation vetor and W1
    cos_lambda1_r(j)=dot(lam_1(end,:)',r);%%% Angle between separation vetor and lambda_1
    cos_W1_l1(j)=dot(a',lam_1(end,:)');%%% Angle between W1 and lambda_1
    
    b=W2(end,:);
    cos_W2_r(j)=dot(b',r);
    cos_lambda2_r(j)=dot(lam_2(end,:)',r);
    cos_W2_l2(j)=dot(b',lam_2(end,:)');
    
    c=W3(end,:);
    cos_W3_r(j)=dot(c',r);
    cos_lambda3_r(j)=dot(lam_3(end,:)',r);
    cos_W3_l3(j)=dot(c',lam_3(end,:)');
    
    w1_stat(j)=w1(end);
    w2_stat(j)=w2(end);
    w3_stat(j)=w3(end);
end


for kk=1:10:size(W1,1)
    
    R11=W1(kk,1);%l1x;
    R12=W1(kk,2);%l1y;
    R13=W1(kk,3);%l1z;
    
    R21=W2(kk,1);%l2x;
    R22=W2(kk,2);%l2y;
    R23=W2(kk,3);%l2z;
    
    R31=W3(kk,1);%l3x;
    R32=W3(kk,2);%l3y;
    R33=W3(kk,3);%l3z;
    
    % Rotation Matrix
    R=[R11 R12 R13; R21 R22 R23;R31 R32 R33]';
    % mesh gridpoints of elipsoid
    [xe,ye,ze]=ellipsoid(0,0,0,w1(kk)/2,w2(kk)/2,w3(kk)/2);
    %%%% Markus-----translation; keep ellipsoids apart from each other
    %[xe,ye,ze]=ellipsoid(x_(kk)*100,y_(kk)/100,z_(kk),w1(kk)/2,w2(kk)/2,w3(kk)/2);
    
    % rotate coordinates
    for i=1:size(xe,1)
        for j=1:size(xe,1)
            dummy = R*[xe(i,j) ye(i,j) ze(i,j)]';
            x_r(i,j) = dummy(1)+x_(kk);
            y_r(i,j) = dummy(2)+y_(kk);
            z_r(i,j) = dummy(3)+z_(kk);
        end
    end
    
    %figure(100);
    
    hold on
    box on
    surf(x_r,y_r,z_r);
    scatter3(x_(kk),y_(kk),z_(kk),100,'filled')
    plot3(x_(kk),y_(kk),z_(kk),'r')
    shading flat
    alpha(0.2)
    %     tit=['frame=',num2str(kk+frame_0-1)];
    %     title(tit)
    xlabel('\bf x')
    ylabel('\bf y')
    zlabel('\bf z')
    %hold on;
    %quiver3([0 0 0]',[0 0 0]',[0 0 0]',V(:,1)*w1(kk)/2,V(:,2)*w2(kk)/2,V(:,3)*w3(kk)/2,0)
    %     line([-V(1,1)*w1(kk)/2 V(1,1)*w1(kk)/2],[-V(2,1)*w1(kk)/2 V(2,1)*w1(kk)/2],[-V(3,1)*w1(kk)/2 V(3,1)*w1(kk)/2],'Color','r')
    %     line([-V(1,2)*w2(kk)/2 V(1,2)*w2(kk)/2],[-V(2,2)*w2(kk)/2 V(2,2)*w2(kk)/2],[-V(3,2)*w2(kk)/2 V(3,2)*w2(kk)/2],'Color','g')
    %     line([-V(1,3)*w3(kk)/2 V(1,3)*w3(kk)/2],[-V(2,3)*w3(kk)/2 V(2,3)*w3(kk)/2],[-V(3,3)*w3(kk)/2 V(3,3)*w3(kk)/2],'Color','b')
    %quiver3([0 0 0]',[0 0 0]',[0 0 0]',V(:,1),V(:,2),V(:,3),0)
    %axis equal
    %[V D.^0.5]
end




%
% figure('Name','Eigen values','NumberTitle','on');
% plot(log(w1),'r')
% hold on
% plot(log(w2),'g')
% plot(log(w3),'b')
% xlabel('Frame')
% ylabel('\bf ln (l_{i})')
% legend('l1','l2','l3')
















