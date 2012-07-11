




function store=inter_v2(x1_0,y1_0,z1_0);




dl=0.0007;
frame_rate=430;
%frame_rate=5000;
delta_t=1/frame_rate;

load data_rsl_p7mm
load velocity_component

%load data_rsl_p7mm_5kHz
%load velocity_component_5kHz



[ux,uy,uz] = gradient(m_u,dl);
[vx,vy,vz] = gradient(m_v,dl);
[wx,wy,wz] = gradient(m_w,dl);







%--------------------Authetic-------------------
% Frame 19730(potential suspect of the breakage); x_t=587; y_t= 409;
%              x_rt=8.831; y_rt=31.519;z_rt=11.263;
% x1_0=8.831/1000;
% y1_0=31.519/1000;
% z1_0=11.263/1000;
%--------------------------------------------

% x1_0=8.831/1000;
% y1_0=29.519/1000;
% z1_0=11.263/1000;

%--------------------------------------------
time_step=1;
store=repmat(nan,[8,time_step]);

%[x1_0,y1_0,z1_0]=search(n,xpix,ypix);

% [x1_0,y1_0,z1_0]=distance(varargin);
for time=1:time_step;

    if time>1
        x1_0=x1_1(time-1);
        y1_0=y1_1(time-1);
        z1_0=z1_1(time-1);
    end

    ux1=interp3(m_x,m_y,m_z,ux,x1_0,y1_0,z1_0);
    uy1=interp3(m_x,m_y,m_z,uy,x1_0,y1_0,z1_0);
    uz1=interp3(m_x,m_y,m_z,uz,x1_0,y1_0,z1_0);

    vx1=interp3(m_x,m_y,m_z,vx,x1_0,y1_0,z1_0);
    vy1=interp3(m_x,m_y,m_z,vy,x1_0,y1_0,z1_0);
    vz1=interp3(m_x,m_y,m_z,vz,x1_0,y1_0,z1_0);

    wx1=interp3(m_x,m_y,m_z,wx,x1_0,y1_0,z1_0);
    wy1=interp3(m_x,m_y,m_z,wy,x1_0,y1_0,z1_0);
    wz1=interp3(m_x,m_y,m_z,wz,x1_0,y1_0,z1_0);

    s_11=0.5*(ux1+ux1);
    s_12=0.5*(uy1+vx1);
    s_13=0.5*(uz1+wx1);
    s_22=0.5*(vy1+vy1);
    s_23=0.5*(vz1+wy1);
    s_33=0.5*(wz1+wz1);

    strain=s_11.*s_11+s_22.*s_22+s_33.*s_33+2*s_12.*s_12+s_13.*s_13+s_23.*s_23;

    %scatter3(x1_1(time),y1_1(time),z1_1(time),log10(strain(time)))

    u1=interp3(m_x,m_y,m_z,m_u,x1_0,y1_0,z1_0);
    v1=interp3(m_x,m_y,m_z,m_v,x1_0,y1_0,z1_0);
    w1=interp3(m_x,m_y,m_z,m_w,x1_0,y1_0,z1_0);
    
    vel=(u1.^2+v1.^2+w1.^2).^0.5;
   
    % Corrected pixel; without streak
    
%     total_pix=correct(total_pix,v1);
%     xpix=correct(xpix,u1);
%     ypix=correct(ypix,v1);
    
    
    %%----additional parameter----
    
%     total_pix_mm=total_pix/40;
%     xpix_mm=xpix/40;
%     ypix_mm=ypix/40;
%     
%     ecentricity=sqrt(abs(1-(xpix_mm/ypix_mm)^2));
%     area=3.14*0.25*xpix_mm*ypix_mm;
%     circumference=3.14*(xpix_mm/2+ypix_mm/2)*(1+((3*((xpix_mm-ypix_mm)/(xpix_mm+ypix_mm))^2)/(10+sqrt(4-3*((xpix_mm-ypix_mm)/(xpix_mm+ypix_mm))^2))));
%     
    
    %%%-----------------------------------------------
    
    x1_1(time)=x1_0+u1*delta_t;
    y1_1(time)=y1_0+v1*delta_t;
    z1_1(time)=z1_0+w1*delta_t;



    strain_mat=[s_11 s_12 s_13; s_12 s_22 s_23; s_13 s_23 s_33];
    if norm(strain_mat)>0
        [e_vector,e_value]=eig(strain_mat);

        [e_value,col_num]=sort(diag(e_value));
        e_value=diag(e_value(end:-1:1));
        e_vector=e_vector(:,col_num(end:-1:1));


        e_val_1=e_value(1,1);
        e_val_2=e_value(2,2);
        e_val_3=e_value(3,3);

        e_vec_1x=e_vector(1,1);
        e_vec_1y=e_vector(2,1);
        e_vec_1z=e_vector(3,1);


        e_vec_2x=e_vector(1,2);
        e_vec_2y=e_vector(2,2);
        e_vec_2z=e_vector(3,2);

        e_vec_3x=e_vector(1,3);
        e_vec_3y=e_vector(2,3);
        e_vec_3z=e_vector(3,3);

    else
        e_val_1=0;
        e_val_2=0;
        e_val_3=0;

        e_vec_1x=0;
        e_vec_1y=0;
        e_vec_1z=0;


        e_vec_2x=0;
        e_vec_2y=0;
        e_vec_2z=0;

        e_vec_3x=0;
        e_vec_3y=0;
        e_vec_3z=0;
    end

    


%     if vel~=0
%         abscostheta=abs((e_vec_1x.*u1+e_vec_1y.*v1+e_vec_1z.*w1)...
%             ./(e_vec_1x.^2+e_vec_1y.^2+e_vec_1z.^2).^0.5...
%             ./(u1.^2+v1.^2+w1.^2).^0.5);
% 
%         theta_radian=(acos(abscostheta));
%         theta_degree=rad2deg(theta_radian);
% 
%     else
%         vel=-999;
%         theta_degree=-999;
%         %error('Velocity at that point is zero and can not be a denominator !!!')
% 
%     end
store(:,time)=[x1_0 y1_0 z1_0 strain e_val_1 e_val_2 e_val_3   vel];%theta_degree
 
end
%theta_degree(zero_vel)=NaN;

%theta_degree

 %look_up=[x1_1; y1_1; z1_1; strain; e_val_1;e_val_2;e_val_3;theta_degree ];
%  fprintf (' x\t       y\t        z\t     Strain\t    lambda_1\t lambda_2\t lambda_3\t  angle\t    totalpix\t    xpix\t    ypix\t     \n')
%  fprintf ('%2.5f\t %2.5f\t %2.5f\t %2.5f\t %2.5f\t %2.5f\t %2.5f\t %2.5f\t %5.2f\t  %7.2f\t  %7.2f\t \n',store)

% figure;
% scatter3(x1_1,y1_1,z1_1)
% hold on;
% draw_channel(1)
%
% figure;
% sus_frame=imread('E:\PTV\Working_folder\WF_1\Breakage_event_stat\img\cam1.19740.tif');
% imshow(sus_frame)
% 


% colorbar
%
%  xlim([-0.003 0.022 ])
%  ylim([0 0.045])
%  zlim([0.00 0.012])

% figure;
% plot(y1_1,log10(strain))
% xlabel('Stream wise distance (m)');
% ylabel('log10_{strain}')
% grid on
