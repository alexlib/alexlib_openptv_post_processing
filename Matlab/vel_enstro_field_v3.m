%function res=vel_enstro_field_v2(first,last,minx,maxx,miny,maxy,minz,maxz,dl,data_rsl_p7mm)

 clear all
 clc
 first=1;
 last=10000;%10918;
 minx=0.01;%-0.005;
 maxx=0.055;%0.025;
 miny=0;%-0.010;
 maxy=0.055;%0.025;
 minz=0;%0.005;
 maxz=0.020;%0.020;
 dl=0.5/1000;
 savee=0;



cond_c=zeros(length(minx:dl:maxx),length(miny:dl:maxy),length(minz:dl:maxz));
cond_u=zeros(length(minx:dl:maxx),length(miny:dl:maxy),length(minz:dl:maxz));
cond_v=zeros(length(minx:dl:maxx),length(miny:dl:maxy),length(minz:dl:maxz));
cond_w=zeros(length(minx:dl:maxx),length(miny:dl:maxy),length(minz:dl:maxz));


%if savee ==0
    for ii=first:1:last
        ii
        %f=load(['E:\WF_turbulent\exp5\res\xuap.',num2Str(ii)]);
        f=load(['E:\mean_flow_analysis_in_turbulent_box\rpm_1000\res\xuap.',num2Str(ii)]);
        s=size(f);
        if s(1,1)>0
            x=f(:,6);
            y=f(:,7);
            z=f(:,8);
            u=f(:,9);
            v=f(:,10);
            w=f(:,11);
            q=f(:,15);
            vel=(u.^2+v.^2+w.^2).^0.5;

            ind=find(q>0 & x>minx & x<maxx  &  y>miny & y<maxy  & z> minz & z<maxz & vel<1);

            x=x(ind);
            y=y(ind);
            z=z(ind);
            u=u(ind);
            v=v(ind);
            w=w(ind);


            for n=1:length(x)
                ind_x=round((x(n)-minx)/dl);
                ind_y=round((y(n)-miny)/dl);
                ind_z=round((z(n)-minz)/dl);
                if ind_x >0 & ind_y >0 & ind_z>0

                cond_c(ind_x,ind_y,ind_z)=cond_c(ind_x,ind_y,ind_z)+1;
                cond_u(ind_x,ind_y,ind_z)=cond_u(ind_x,ind_y,ind_z)+u(n);
                cond_v(ind_x,ind_y,ind_z)=cond_v(ind_x,ind_y,ind_z)+v(n);
                cond_w(ind_x,ind_y,ind_z)=cond_w(ind_x,ind_y,ind_z)+w(n);
                
                end

            end


        end
    end
    %save data_rsl_p5mm_exp5 cond_c cond_u cond_v cond_w
    %save rsl_2mm_rpm1000 cond_c cond_u cond_v cond_w 
%else
    %load data_rsl_p1mm_exp5;
    %load rsl_2mm_rpm1000
%end

%%----------------------Averaging stuff-------------------%%%5
i_c=0;
for i=minx:dl:maxx
    i_c=i_c+1;
    j_c=0;
    for j=miny:dl:maxy
        j_c=j_c+1;
        k_c=0;
        for k=minz:dl:maxz
            k_c=k_c+1;
            if  cond_c(i_c,j_c,k_c)> 5
                cond_u(i_c,j_c,k_c)=cond_u(i_c,j_c,k_c)/cond_c(i_c,j_c,k_c);
                cond_v(i_c,j_c,k_c)=cond_v(i_c,j_c,k_c)/cond_c(i_c,j_c,k_c);
                cond_w(i_c,j_c,k_c)=cond_w(i_c,j_c,k_c)/cond_c(i_c,j_c,k_c);
            else
                cond_u(i_c,j_c,k_c)=0;
                cond_v(i_c,j_c,k_c)=0;
                cond_w(i_c,j_c,k_c)=0;
            end

            % Boundary condition %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             if j<0.0245 % chamber
%                 if k< 0.000 | k>0.012 | i<-0.003 | i>0.022
% 
%                     cond_u(i_c,j_c,k_c)=0;
%                     cond_v(i_c,j_c,k_c)=0;
%                     cond_w(i_c,j_c,k_c)=0;
%                  end
% 
%             else % channel
% 
%                 if k<0.009 | k>0.012  | i<0.0075 | i>0.0105
%                     cond_u(i_c,j_c,k_c)=0;
%                     cond_v(i_c,j_c,k_c)=0;
%                     cond_w(i_c,j_c,k_c)=0;
% 
%                 end
%             end
        end
    end
end

%%-------------Mesh grid------------%%
i_c=0;
for i=minx:dl:maxx
    i_c=i_c+1;
    j_c=0;
    for j=miny:dl:maxy
        j_c=j_c+1;
        k_c=0;
        for k=minz:dl:maxz
            k_c=k_c+1;
            m_x(j_c,i_c,k_c)=i;%-minx;
            m_y(j_c,i_c,k_c)=j;%-miny;
            m_z(j_c,i_c,k_c)=k;%-minz; 
            m_c(j_c,i_c,k_c)=cond_c(i_c,j_c,k_c);
            m_u(j_c,i_c,k_c)=cond_u(i_c,j_c,k_c);
            m_v(j_c,i_c,k_c)=cond_v(i_c,j_c,k_c);
            m_w(j_c,i_c,k_c)=cond_w(i_c,j_c,k_c); 
            m_vel(j_c,i_c,k_c)=(m_u(j_c,i_c,k_c)^2+m_v(j_c,i_c,k_c)^2+m_w(j_c,i_c,k_c)^2)^0.5;
        end
    end
end

%save velocity_component_tb_p5mm_exp5.mat  m_c m_x m_y m_z m_u m_v m_w m_vel
if savee==0;
save rpm1000_gridp5mm.mat  cond_c cond_u cond_v cond_w m_c m_x m_y m_z m_u m_v m_w m_vel
else
    load rpm1000_gridp5mm
end



% % for s=1:5
% %     m_u   = smooth3(m_u,'gaussian',3);
% %     m_v   = smooth3(m_v,'gaussian',3);
% %     m_w   = smooth3(m_w ,'gaussian',3);
% %     m_vel = smooth3(m_vel,'gaussian',3);
% % 
% % end
% 
% ij=find(m_vel==0);
% [ux,uy,uz] = gradient(m_u,dl);
% [vx,vy,vz] = gradient(m_v,dl);
% [wx,wy,wz] = gradient(m_w,dl);
% 
% 
% 
% 
% %% ---------------------------------Eigen value------------------------%%
% 
% i_c=0;
% for i=minx:dl:maxx
%     i_c=i_c+1;
%     j_c=0;
%     for j=miny:dl:maxy
%         j_c=j_c+1;
%         k_c=0;
%         for k=minz:dl:maxz
%             k_c=k_c+1;
%             
%             s11=1/2*(ux(j_c,i_c,k_c)+ux(j_c,i_c,k_c));
%             s12=1/2*(uy(j_c,i_c,k_c)+vx(j_c,i_c,k_c));
%             s13=1/2*(uz(j_c,i_c,k_c)+wx(j_c,i_c,k_c));
%             s22=1/2*(vy(j_c,i_c,k_c)+vy(j_c,i_c,k_c));
%             s23=1/2*(vz(j_c,i_c,k_c)+wy(j_c,i_c,k_c));
%             s33=1/2*(wz(j_c,i_c,k_c)+wz(j_c,i_c,k_c));
%             A=[s11 s12 s13;s12 s22 s23;s13 s23 s33];
%             if norm(A)>0
%                 [V,D] = eig(A);
%                 [D,k] = sort(diag(D));     % ascending order
%                 D = (diag(D(end:-1:1)));   % descending order and again makes it a diagonal matrix
%                 V = V(:,k(end:-1:1));      % the same order(descending)for eigenvectors
%                 
%                 m_eig_value1(j_c,i_c,k_c)=D(1,1);
%                 m_eig_value2(j_c,i_c,k_c)=D(2,2);
%                 m_eig_value3(j_c,i_c,k_c)=D(3,3);
%                 
%                 m_eig_vec1x(j_c,i_c,k_c)=V(1,1);
%                 m_eig_vec1y(j_c,i_c,k_c)=V(2,1);
%                 m_eig_vec1z(j_c,i_c,k_c)=V(3,1);
%                 
%                 m_eig_vec2x(j_c,i_c,k_c)=V(1,2);
%                 m_eig_vec2y(j_c,i_c,k_c)=V(2,2);
%                 m_eig_vec2z(j_c,i_c,k_c)=V(3,2);
%                 
%                 m_eig_vec3x(j_c,i_c,k_c)=V(1,3);
%                 m_eig_vec3y(j_c,i_c,k_c)=V(2,3);
%                 m_eig_vec3z(j_c,i_c,k_c)=V(3,3);
%                 
%                               
%             else
%                 m_eig_value1(j_c,i_c,k_c)=0;
%                 m_eig_value2(j_c,i_c,k_c)=0;
%                 m_eig_value3(j_c,i_c,k_c)=0;
%                 m_eig_vec1x(j_c,i_c,k_c)=0;
%                 m_eig_vec1y(j_c,i_c,k_c)=0;
%                 m_eig_vec1z(j_c,i_c,k_c)=0;
%                 m_eig_vec2x(j_c,i_c,k_c)=0;
%                 m_eig_vec2y(j_c,i_c,k_c)=0;
%                 m_eig_vec2z(j_c,i_c,k_c)=0;
%                 m_eig_vec3x(j_c,i_c,k_c)=0;
%                 m_eig_vec3y(j_c,i_c,k_c)=0;
%                 m_eig_vec3z(j_c,i_c,k_c)=0;
%                 
%             end
%         end
%     end
% end
% 
% 
% for s=1:5
%     m_eig_value1  = smooth3(m_eig_value1,'gaussian',3);
%     m_eig_value2  = smooth3(m_eig_value2,'gaussian',3);
%     m_eig_value3  = smooth3(m_eig_value3,'gaussian',3);
%     m_eig_vec1x = smooth3(m_eig_vec1x,'gaussian',3);
%     m_eig_vec1y = smooth3(m_eig_vec1y,'gaussian',3);
%     m_eig_vec1z = smooth3(m_eig_vec1z,'gaussian',3);
%     m_eig_vec2x = smooth3(m_eig_vec2x,'gaussian',3);
%     m_eig_vec2y = smooth3(m_eig_vec2y,'gaussian',3);
%     m_eig_vec2z = smooth3(m_eig_vec2z,'gaussian',3);
%     m_eig_vec3x = smooth3(m_eig_vec3x,'gaussian',3);
%     m_eig_vec3y = smooth3(m_eig_vec3y,'gaussian',3);
%     m_eig_vec3z = smooth3(m_eig_vec3z,'gaussian',3);
% 
% end
% 
% 
%     
% abscostheta=abs((m_eig_vec3x.*m_u+m_eig_vec3y.*m_v+m_eig_vec3z.*m_w)...
%     ./(m_eig_vec3x.^2+m_eig_vec3y.^2+m_eig_vec3z.^2).^0.5...
%     ./(m_u.^2+m_v.^2+m_w.^2).^0.5);
% theta_radian=(acos(abscostheta));
% theta_degree=rad2deg(theta_radian);
% 
% theta_degree(ij)=NaN;
% save theta_degree
% 
% % eig_vec1=[m_eig_vec1x(j_c,i_c,k_c);m_eig_vec1y(j_c,i_c,k_c);m_eig_vec1z(j_c,i_c,k_c)].* m_eig_value1(j_c,i_c,k_c);
% % eig_vec2=[m_eig_vec2x(j_c,i_c,k_c);m_eig_vec2y(j_c,i_c,k_c);m_eig_vec2z(j_c,i_c,k_c)].* m_eig_value2(j_c,i_c,k_c);
% % eig_vec3=[m_eig_vec3x(j_c,i_c,k_c);m_eig_vec3y(j_c,i_c,k_c);m_eig_vec3z(j_c,i_c,k_c)].* m_eig_value3(j_c,i_c,k_c);
% % 
% % vel=[m_u(j_c,i_c,k_c);m_v(j_c,i_c,k_c);m_w(j_c,i_c,k_c)];
% % mag_vel=sqrt((m_u(j_c,i_c,k_c))^2+(m_v(j_c,i_c,k_c))^2+(m_w(j_c,i_c,k_c)^2));
% 
% 
% % mag_eig_vec1=sqrt((m_eig_vec1x(j_c,i_c,k_c))^2+(m_eig_vec1y(j_c,i_c,k_c))^2+(m_eig_vec1z(j_c,i_c,k_c))^2);
% % mag_eig_vec2=sqrt((m_eig_vec2x(j_c,i_c,k_c))^2+(m_eig_vec2y(j_c,i_c,k_c))^2+(m_eig_vec2z(j_c,i_c,k_c))^2);
% % mag_eig_vec3=sqrt((m_eig_vec3x(j_c,i_c,k_c))^2+(m_eig_vec3y(j_c,i_c,k_c))^2+(m_eig_vec3z(j_c,i_c,k_c))^2);
% % 
% % 
% % costheta=dot(vel,eig_vec1)/(mag_vel*mag_eig_vec1);
% % 
% % theta_radian=(acos(costheta));
% % theta_degree=rad2deg(theta_radian)
% 
% %--------------------------------------------------------------------------
% 
% %     %%%%%%%%LAMBDA 2%%%%%%%%%%%%%%%%%%%%%%%%%
% %     i_c=0;
% %     for i=minx:dl:maxx
% %         i_c=i_c+1;
% %         j_c=0;
% %         for j=miny:dl:maxy
% %             j_c=j_c+1;
% %             k_c=0;
% %             for k=minz:dl:maxz
% %                 k_c=k_c+1;
% %                 uM=[ux(j_c,i_c,k_c) uy(j_c,i_c,k_c) uz(j_c,i_c,k_c);vx(j_c,i_c,k_c) vy(j_c,i_c,k_c) vz(j_c,i_c,k_c);wx(j_c,i_c,k_c) wy(j_c,i_c,k_c) wz(j_c,i_c,k_c)];
% %                 for i=1:3
% %                     for j=1:3
% %                         s(i,j)=0.5*(uM(i,j)+uM(j,i));
% %                         o(i,j)=0.5*(uM(i,j)-uM(j,i));
% %                     end
% %                 end
% %                 M=zeros(3,3);
% %                 for i=1:3
% %                     for j=1:3
% %                         for k=1:3
% %                             M(i,j)=M(i,j)+o(i,k)*o(k,j)+s(i,k)*s(k,j);
% %                         end
% %                     end
% %                 end
% %                 [V,D] = eig(M);
% %                 [D,k] = sort(diag(D)); 	% ascending order
% %                 D = diag(D(end:-1:1)); 	% descending order
% %                 V = V(:,k(end:-1:1)); 	% the same order for eigenvectors
% %                 vort(j_c,i_c,k_c)=D(2,2);
% %             end
% %         end
% %     end
% %     vort = smooth3(vort,'box',3);
% %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %%%%%%%%enstro&strain %%%%%%%%%%%%%%%%%%%%%%%%%
% s_11=0.5*(ux+ux);
% s_12=0.5*(uy+vx);
% s_13=0.5*(uz+wx);
% s_22=0.5*(vy+vy);
% s_23=0.5*(vz+wy);
% s_33=0.5*(wz+wz);
% w1=wy-vz;
% w2=uz-wx;
% w3=vx-uy;
% strain=s_11.*s_11+s_22.*s_22+s_33.*s_33+2*(s_12.*s_12+s_13.*s_13+s_23.*s_23);
% enstro=w1.*w1+w2.*w2+w3.*w3;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %%% Reynolds number%%%%
% 
% %% Flux field%%%%%%
% %close all;
% 
% %.............Positive value of velocity gradient for plotting....Debashish......%
% value=1e-8; % this extremely low value is aasigned for the zeros and negative values as they are hard to treat by log 
% positive_vz=vz;
% location=find( positive_vz<=0);
% positive_vz(location)=value; %fills up the negatives and zeros with negligible figure
% %........................................................................
% % Problem:
% %This works fine for vz(=dv/dz) but for vx(=dv/dx) it shows a picture where
% %half of the graphics are filled with almost no value(blue colored zone) and comparison can not be made with numerical... 
% %.......................................................................
% %........................................................................
% 
% 
% % The following code does not replace all the negatives and zeros rather it
% % takes the absolute values of them and replaces only those whose value are
% % smaller than 1e-7 by 1e-7. This comply with the numericals as numericals
% % also account for abs value.
% %.......................Beat.................................................
% positive_vx=abs(vx);
% location=find(abs(vx)<1e-7);
% positive_vx(location)=1e-7;
% 
% %................................................................%
% 
% parameter1=m_v;
% parameter2=log10(positive_vz);
% parameter3=log10(abs(positive_vx));
% 
% %--------------------------------Scaling for quiver stuff-----------------%
% fact=0.0003;
% m_eig_vec1x = fact*m_eig_vec1x;
% m_eig_vec1y = fact*m_eig_vec1y;
% m_eig_vec1z = fact*m_eig_vec1z;
% m_eig_vec2x = fact*m_eig_vec2x;
% m_eig_vec2y = fact*m_eig_vec2y;
% m_eig_vec2z = fact*m_eig_vec2z;
% m_eig_vec3x = fact*m_eig_vec3x;
% m_eig_vec3y = fact*m_eig_vec3y;
% m_eig_vec3z = fact*m_eig_vec3z;
% %-------------------------------------------------------------------------%
% 
% %save Analysis_parameter.mat m_x m_y m_z m_eig_vec1x m_eig_vec1y m_eig_vec1z m_eig_vec2x m_eig_vec2y m_eig_vec2z m_eig_vec3x m_eig_vec3y m_eig_vec3z m_eig_value1 m_eig_value2 m_eig_value3
% 
% %----------------------------------------------------
% 
% %----Angle between eigen vector and velocity gradient vector ......................%
% 
% 
% 
% %----------------------------------------------------%
% 
% figure;
% hold on;
% box on;
% % ---------------------------------------Quiver----------------------------%


   %quiver3(m_x(1:5:end,1:5:end,1:5:end),m_y(1:5:end,1:5:end,1:5:end),m_z(1:5:end,1:5:end,1:5:end),m_u(1:5:end,1:5:end,1:5:end),m_v(1:5:end,1:5:end,1:5:end),m_w(1:5:end,1:5:end,1:5:end),2);
      quiverc3D(m_x(1:2:end,1:2:end,1:2:end),m_y(1:2:end,1:2:end,1:2:end),m_z(1:2:end,1:2:end,1:2:end),m_u(1:2:end,1:2:end,1:2:end),m_v(1:2:end,1:2:end,1:2:end),m_w(1:2:end,1:2:end,1:2:end),5);  

      axis equal
%   %quiver3(m_x(10:2:end,1:end,1:end),m_y(10:2:end,1:end,1:end),m_z(10:2:end,1:end,1:end),m_u(10:2:end,1:end,1:end),m_v(10:2:end,1:end,1:end),m_w(10:2:end,1:end,1:end),5);
%   %draw_channel(1)
%   %geometry(1)
%   
%   
% % quiver3(m_x,m_y,m_z,m_eig_vec1x.*m_eig_value1,m_eig_vec1y.*m_eig_value1,m_eig_vec1z.*m_eig_value1,0,'r');
% % quiver3(m_x,m_y,m_z,-m_eig_vec1x.*m_eig_value1,-m_eig_vec1y.*m_eig_value1,-m_eig_vec1z.*m_eig_value1,0,'r'); 
% 
% % quiver3(m_x,m_y,m_z,m_eig_vec2x.*m_eig_value2,m_eig_vec2y.*m_eig_value2,m_eig_vec2z.*m_eig_value2,0,'g');
% % quiver3(m_x,m_y,m_z,-m_eig_vec2x.*m_eig_value2,-m_eig_vec2y.*m_eig_value2,-m_eig_vec2z.*m_eig_value2,0,'g');
% 
% % quiver3(m_x,m_y,m_z,m_eig_vec3x.*m_eig_value3,m_eig_vec3y.*m_eig_value3,m_eig_vec3z.*m_eig_value3,0,'b');
% % quiver3(m_x,m_y,m_z,-m_eig_vec3x.*m_eig_value3,-m_eig_vec3y.*m_eig_value3,-m_eig_vec3z.*m_eig_value3,0,'b');
% 
% 
% 
% 
% 
% %------------------- Analysis from different "Angle of view" ---------------%
% % -----------------It is saved in "Eigen_frame_analysis.m" file--------------------%
% 
% %------------------------------------- Slice -----------------------------%
% %xslice =[.006:.001:.012];
% 
% % xslice =[.0075:.0005:.0105];
% % yslice = [.02:.005:0.04]; 
% % zslice =[.001:.002:0.009]; 
% 
% % slice(m_x,m_y,m_z,m_eig_value1,xslice,yslice,zslice)
% % slice(m_x,m_y,m_z,m_eig_value2,xslice,yslice,zslice)
% % slice(m_x,m_y,m_z,m_eig_value3,xslice,yslice,zslice)
% % 
% % slice(m_x,m_y,m_z,m_eig_value1.^2+m_eig_value2.^2+m_eig_value3.^2,xslice,yslice,zslice)
% 
% 
% 
% 
% % h=slice(m_x(30:end-5,12:end-12,1:end),m_y(30:end-5,12:end-12,1:end),m_z(30:end-5,12:end-12,1:end),m_vel(30:end-5,12:end-12,1:end),xslice,yslice,zslice);
% % %h=slice(m_x,m_y,m_z,strain,xslice,yslice,zslice);
% % set(h,'FaceColor','interp',...
% % 	'EdgeColor','none',...
% % 	'DiffuseStrength',.8)
% %  alpha(0.5)
% %  colorbar
% %  
%  
%  
%  
%  
% %  lightangle(-45,45)
% % colormap (jet(24))
% % set(gcf,'Renderer','zbuffer')
% 
% %  colormap (flipud(jet(24)))
% % colorbar('horiz')
% 
% 
%  
% 
%  
%  
% strain(find(strain==0))=1E-1;
% %slice(m_x(20:end-20,:,:),m_y(20:end-20,:,:),m_z(20:end-20,:,:),log10((strain(20:end-20,:,:))),xslice,yslice,zslice)
% 
% %-------------------------------------------------------------------------%
% 
% 
% 
% % camlight 
% % lighting gouraud
% % camlight headlight
% %-------------------------------------------------------------------------%
% 
% 
% %       [sx sy sz] = meshgrid([0.15 0.25],0:0.02:0.2,0:0.003:0.03);
% %       [sx sy sz] = meshgrid(0.05,0:0.01:0.1,0:0.003:0.03);
% %       daspect([1 1 1])
% %       h=streamline(m_x,m_y,m_z,m_u,m_v,m_w,sx,sy,sz);
% %       set(h,'Color','blue')
% % 
% %     verts = stream3(m_x,m_y,m_z,m_u,m_v,m_w,sx,sy,sz);
% %     cav = curl(m_x,m_y,m_z,m_u,m_v,m_w);
% %     %[CURLX, CURLY, CURLZ]= curl(m_x,m_y,m_z,m_u,m_v,m_w);
% %     spd = sqrt(m_u.^2 + m_v.^2 + m_w.^2).*.1;
% %     streamribbon(verts,m_x,m_y,m_z,cav,spd);%
% %     axis tight
% %     shading interp
% %     view(3)
% %     camlight; %lighting gouraud
% %     lighting phong
% 
% %CURLX-w1
% 
% 
% 
% %threshold=60;
% figure;
% threshold=log10(2000);
% 
% 
% hpatch = patch(isosurface(m_x(20:end-20,:,:),m_y(20:end-20,:,:),m_z(20:end-20,:,:),log10((strain(20:end-20,:,:))),threshold));
% isonormals(m_x(20:end-20,:,:),m_y(20:end-20,:,:),m_z(20:end-20,:,:),strain(20:end-20,:,:),hpatch)
% set(hpatch,'FaceColor','red','EdgeColor','none')
% hpatch = patch(isocaps(m_x(20:end-20,:,:),m_y(20:end-20,:,:),m_z(20:end-20,:,:),strain(20:end-20,:,:),threshold));
% isocaps(m_x(20:end-20,:,:),m_y(20:end-20,:,:),m_z(20:end-20,:,:),strain(20:end-20,:,:),hpatch)
% set(hpatch,'FaceColor','interp','EdgeColor','none')
% 
% 
% 
% % max_value=max(max(parameter2))
% 
% 
% 
% 
% 
% %camlight left; camlight; lighting phong
% set(gcf,'Renderer','zbuffer');
% %set(gca, 'CLim', [-0.04, 0.24])
% 
% %daspect([1 1 1])
% %view([13,18])
% view([-24,37])
% xlabel('x(mm)','FontName','Times New Roman','Fontsize',14)
% ylabel('y(mm)','FontName','Times New Roman','Fontsize',14)
% zlabel('z(mm)','FontName','Times New Roman','Fontsize',14)
% %     tit=['time: ',num2Str(floor(10*loop/20+0.5)/10),' sec'];
% %     title(tit);
% %axis([0.00075 0.0105 0.018 0.03 0.00 0.012])
% axis equal
% %colorbar
% 
% 
% %SaveSingleFigure(loop);

% 
% 
% 
% % 
% % [dummyx dummyy]=meshgrid(1:10,1:10);
% % figure;surf(dummyx,dummyy,dummyx.*dummyy)
% % 
% % 
% %figure;contourf(squeeze(m_x(34,:,:)),squeeze(m_z(34,:,:)),squeeze(log10(abs(vz(34,:,:)))),-3:.5:3)
