% function res=cosines_a_l_i(first,last,tail,stepSize);
%pdfsPointOne(0,5990,5,5)

% cosines_a_l_i(10000,15990,3,1) % from script.m
first = 1;
last = 2990; 
tail = 3;
stepSize = 1;


%% Assign the working data
curdir = cd;
hfiles{1} = [curdir(1),':\Matlab\Alex\HDF\0104_water.hdf'];
% or
hfiles{2} = [curdir(1),':\Matlab\Alex\HDF\0104_hompol.hdf'];

counter=0;
pos_diva_pos_ua_counter=0;
pos_diva_neg_ua_counter=0;
neg_diva_weak_counter=0;
neg_diva_strong_counter=0;

if 1<2
    statistic=zeros(200000,30);
    %%%%1017288
    for i=first:stepSize:last
         f1 = hdfread(hfiles{1},int2str(i),...
            'Fields','x,y,z,u,v,w,ax,ay,az,w1,w2,w3,s11,s12,s13,s22,s23,s33,ww1,ww2,ww3,wws,sss,R,Q,diss,div,trajnum,t,alx,aly,alz',...
            'FirstRecord',1);
        s = length(f1{1});
        
        if s(1,1)>0
            
            
            [...
                xr,yr,zr,...                % 1-3
                u,v,w,...                   % 4-6
                ax,ay,az,...                % 7-9
                w1,w2,w3,...                % 10-12
                s11,s12,s13,s22,s23,s33,... % 13-18
                ww1,ww2,ww3,...             % 19-21
                wws,sss,R,Q,diss,...        % 22,23,24,25,26
                reldiv,trajnum,...             % 27,28
                age,...                       % 29 (!)
                alx,aly,alz] = ...          % 30-32
                deal(f1{1:32});
            
%             xr=f1(:,1);
%             yr=f1(:,2);
%             zr=f1(:,3);
%             u=f1(:,4);
%             v=f1(:,5);
%             w=f1(:,6);
            absu=(u.^2+v.^2+w.^2).^0.5;
%             ax=f1(:,19);
%             ay=f1(:,20);
%             az=f1(:,21);
            absa=(ax.^2+ay.^2+az.^2).^0.5;
%             alx=f1(:,16);
%             aly=f1(:,17);
%             alz=f1(:,18);
            absal=(alx.^2+aly.^2+alz.^2).^0.5;
            
%             ux=f1(:,7);
%             uy=f1(:,8);
%             uz=f1(:,9);
%             vx=f1(:,10);
%             vy=f1(:,11);
%             vz=f1(:,12);
%             wx=f1(:,13);
%             wy=f1(:,14);
%             wz=f1(:,15); 
            ux  =   s11;
            uy  =   s12 - 0.5 * w3;
            uz  =   s13 + 0.5 * w2;
            vx  =   s12 + 0.5 * w3;
            vy  =   s22;
            vz  =   s23 - 0.5 * w1;
            wx  =   s13 - 0.5 * w2;
            wy  =   s23 + 0.5 * w1;
            wz  =   s33;
            
%             
%             w1=wy-vz;
%             w2=uz-wx;
%             w3=vx-uy;
            enstrophy=w1.^2+w2.^2+w3.^2;
            absw=enstrophy.^0.5;
%             s11=0.5*(ux+ux);
%             s12=0.5*(uy+vx);
%             s13=0.5*(uz+wx);
%             s22=0.5*(vy+vy);
%             s23=0.5*(vz+wy);
%             s33=0.5*(wz+wz);
            
            wS1=w1.*s11+w2.*s12+w3.*s13;
            wS2=w1.*s12+w2.*s22+w3.*s23;
            wS3=w1.*s13+w2.*s23+w3.*s33;
            absws=(wS1.^2+wS2.^2+wS3.^2).^0.5;
            wws=w1.*wS1+w2.*wS2+w3.*wS3;

            
            acx=u.*ux+v.*uy+w.*uz;%%%%%%%%%%%%%%%%%%%%
            acy=u.*vx+v.*vy+w.*vz;%%%%%%%%%%%%%%%%%%%%
            acz=u.*wx+v.*wy+w.*wz;%%%%%%%%%%%%%%%%%%%%
            absac=(acx.^2+acy.^2+acz.^2).^0.5;%%%%%%%%%%%%%%%%%%%%%%%%%
            aLamx=w2.*w-w3.*v;%%%%%%%%%%%%%%%%%%%%%%%
            aLamy=w3.*u-w1.*w;%%%%%%%%%%%%%%%%%%%%%%%
            aLamz=w1.*v-w2.*u;%%%%%%%%%%%%%%%%%%%%%%%%%
            absaLam=(aLamx.^2+aLamy.^2+aLamz.^2).^0.5;%%%%%%%%%%%%%%%%%
%             aBx=0.5*f2(:,45);%acx-aLamx;%%%%%%%%%%%%%%%%%%%%%
%             aBy=0.5*f2(:,46);%acy-aLamy;%%%%%%%%%%%%%%%%%%%%%%%
%             aBz=0.5*f2(:,47);%acz-aLamz;%%%%%%%%%%%%%%%%%%%%%%%%%
%             absaB=(aBx.^2+aBy.^2+aBz.^2).^0.5;%%%%%%%%%%%%%%%%%%%
            unx=u./absu;
            uny=v./absu;
            unz=w./absu;
            aun=ax.*unx+ay.*uny+az.*unz;
            aPx=aun.*unx;%%%%%%%%%%%%%%%%%%%%%%%%
            aPy=aun.*uny;%%%%%%%%%%%%%%%%%%%%%%%%%
            aPz=aun.*unz;%%%%%%%%%%%%%%%%%%%%%%%%%%
            absaP=(aPx.^2+aPy.^2+aPz.^2).^0.5;%%%%%%%%%%%%%%%%%%%%%
            aOx=ax-aPx;%%%%%%%%%%%%%%%%%%%%
            aOy=ay-aPy;%%%%%%%%%%%%%%%%%%%
            aOz=az-aPz;%%%%%%%%%%%%%%%%%%%%
            absaO=(aOx.^2+aOy.^2+aOz.^2).^0.5;%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             cux=unx.*f2(:,36)+uny.*f2(:,37)+unz.*f2(:,38);%%%(u.*ux+v.*uy+w.*uz)./(absu.^2);
%             cuy=unx.*f2(:,39)+uny.*f2(:,40)+unz.*f2(:,41);%%%(u.*vx+v.*vy+w.*vz)./(absu.^2);
%             cuz=unx.*f2(:,42)+uny.*f2(:,43)+unz.*f2(:,44);%%%(u.*wx+v.*wy+w.*wz)./(absu.^2);
%             curvGrad=(cux.^2+cuy.^2+cuz.^2).^0.5;
            curvAcc=(((v.*az-w.*ay).^2+(w.*ax-u.*az).^2+(u.*ay-v.*ax).^2).^0.5)./(absu.^3);
            
            % reldiv=abs(ux+vy+wz)./(abs(ux)+abs(vy)+abs(wz));
            % age=t; %f1(:,32);
            
%             axx=f2(:,22);
%             axy=f2(:,23);
%             axz=f2(:,24);
%             ayx=f2(:,25);
%             ayy=f2(:,26);
%             ayz=f2(:,27);
%             azx=f2(:,28);
%             azy=f2(:,29);
%             azz=f2(:,30); 
            
%             aw1=azy-ayz;
%             aw2=axz-azx;
%             aw3=ayx-axy;
%             aenstrophy=aw1.^2+aw2.^2+aw3.^2;
%             absaw=aenstrophy.^0.5;
%             as11=0.5*(axx+axx);
%             as12=0.5*(axy+ayx);
%             as13=0.5*(axz+azx);
%             as22=0.5*(ayy+ayy);
%             as23=0.5*(ayz+azy);
%             as33=0.5*(azz+azz);
            
%             diva=as11+as22+as33;
          
            
            %       so we have: 
            %       absu,     1
            %       absa,     2 
            %       absal,    3 
            %       absac,    4 
            %       absaLam,  5 
            %       absaB,    6 
            %       absaP,    7 
            %       absaO,    8
            
            ende = 0;
            leng = 1;
            for iii=1:s(1,1)-1
                start=ende+1;
                ende=start;
                leng=1;
                searching=1;
                while ende<s(1,1)-1 & searching==1
                    if age(ende+1)==age(ende)+1
                        ende=ende+1;
                        leng=leng+1;
                    else
                        searching=0;
                    end
                end
                if ende<s(1,1)
                    for ii=start:ende
                        if age(ii)>tail & age(ii)<leng-tail+1 & reldiv(ii)<0.1
                            counter=counter+1;
                            
                            A=[s11(ii) s12(ii) s13(ii);s12(ii) s22(ii) s23(ii);s13(ii) s23(ii) s33(ii)];
                            [V,D] = eig(A);
                            [D,k] = sort(diag(D)); 	% ascending order
                            D = diag(D(end:-1:1)); 			% descending order		
                            V = V(:,k(end:-1:1)); 	% the same order for eigenvectors
                            strain=D(1,1)^2+D(2,2)^2+D(3,3)^2;
                            Q=0.25*(enstrophy(ii)-2*strain);
                       
                            statistic(counter,1)=(acx(ii)*V(1,1)+acy(ii)*V(2,1)+acz(ii)*V(3,1))/(absac(ii));
                            statistic(counter,2)=(acx(ii)*V(1,2)+acy(ii)*V(2,2)+acz(ii)*V(3,2))/(absac(ii));
                            statistic(counter,3)=(acx(ii)*V(1,3)+acy(ii)*V(2,3)+acz(ii)*V(3,3))/(absac(ii));
                            
                            statistic(counter,16)=(alx(ii)*ax(ii)+aly(ii)*ay(ii)+alz(ii)*az(ii))/(absal(ii)*absa(ii));
                            statistic(counter,17)=(acx(ii)*ax(ii)+acy(ii)*ay(ii)+acz(ii)*az(ii))/(absac(ii)*absa(ii));
                            statistic(counter,18)=(alx(ii)*acx(ii)+aly(ii)*acy(ii)+alz(ii)*acz(ii))/(absal(ii)*absac(ii));
                            
                            %if absaP(ii)/absaO(ii)>0.8 %high parallel less
                            %alingment
                            %if Q>0
                            %if strain>15 %more aligned but not strong
                            %if u(ii)*ax(ii)+v(ii)*ay(ii)+w(ii)*az(ii)>1e-4
                            %%NICE!
                            %if (u(ii)*ax(ii)+v(ii)*ay(ii)+w(ii)*az(ii))/(absa(ii)*absu(ii))>0.7
                            
                            if Q<-1 %pos diva
                                if (u(ii)*ax(ii)+v(ii)*ay(ii)+w(ii)*az(ii))/(absa(ii)*absu(ii))>0.7 %pos ua
                                    pos_diva_pos_ua_counter=pos_diva_pos_ua_counter+1;
                                    statistic(pos_diva_pos_ua_counter,4)=(acx(ii)*V(1,1)+acy(ii)*V(2,1)+acz(ii)*V(3,1))/(absac(ii));
                                    statistic(pos_diva_pos_ua_counter,5)=(acx(ii)*V(1,2)+acy(ii)*V(2,2)+acz(ii)*V(3,2))/(absac(ii));
                                    statistic(pos_diva_pos_ua_counter,6)=(acx(ii)*V(1,3)+acy(ii)*V(2,3)+acz(ii)*V(3,3))/(absac(ii));
                                    statistic(pos_diva_pos_ua_counter,19)=(alx(ii)*ax(ii)+aly(ii)*ay(ii)+alz(ii)*az(ii))/(absal(ii)*absa(ii));
                                    statistic(pos_diva_pos_ua_counter,20)=(acx(ii)*ax(ii)+acy(ii)*ay(ii)+acz(ii)*az(ii))/(absac(ii)*absa(ii));
                                    statistic(pos_diva_pos_ua_counter,21)=(alx(ii)*acx(ii)+aly(ii)*acy(ii)+alz(ii)*acz(ii))/(absal(ii)*absac(ii));
                                end
                                if (u(ii)*ax(ii)+v(ii)*ay(ii)+w(ii)*az(ii))/(absa(ii)*absu(ii))<-0.7 %neg ua
                                    pos_diva_neg_ua_counter=pos_diva_neg_ua_counter+1;
                                    statistic(pos_diva_neg_ua_counter,7)=(acx(ii)*V(1,1)+acy(ii)*V(2,1)+acz(ii)*V(3,1))/(absac(ii));
                                    statistic(pos_diva_neg_ua_counter,8)=(acx(ii)*V(1,2)+acy(ii)*V(2,2)+acz(ii)*V(3,2))/(absac(ii));
                                    statistic(pos_diva_neg_ua_counter,9)=(acx(ii)*V(1,3)+acy(ii)*V(2,3)+acz(ii)*V(3,3))/(absac(ii));
                                    statistic(pos_diva_neg_ua_counter,22)=(alx(ii)*ax(ii)+aly(ii)*ay(ii)+alz(ii)*az(ii))/(absal(ii)*absa(ii));
                                    statistic(pos_diva_neg_ua_counter,23)=(acx(ii)*ax(ii)+acy(ii)*ay(ii)+acz(ii)*az(ii))/(absac(ii)*absa(ii));
                                    statistic(pos_diva_neg_ua_counter,24)=(alx(ii)*acx(ii)+aly(ii)*acy(ii)+alz(ii)*acz(ii))/(absal(ii)*absac(ii));
                                end 
                            end
                            %if absaP(ii)/absaO(ii)<0.2 %low parrallel more
                            %alignment
                            %if Q<0
                            %if strain<3 %less aligned but not strong
                            %if u(ii)*ax(ii)+v(ii)*ay(ii)+w(ii)*az(ii)<-1e-4
                            %%NICE!
                            %if (u(ii)*ax(ii)+v(ii)*ay(ii)+w(ii)*az(ii))/(absa(ii)*absu(ii))<-0.7
                            if Q>1.5  
                                if Q<5 %neg diva weak
                                    neg_diva_weak_counter=neg_diva_weak_counter+1;
                                    statistic(neg_diva_weak_counter,10)=(acx(ii)*V(1,1)+acy(ii)*V(2,1)+acz(ii)*V(3,1))/(absac(ii));
                                    statistic(neg_diva_weak_counter,11)=(acx(ii)*V(1,2)+acy(ii)*V(2,2)+acz(ii)*V(3,2))/(absac(ii));
                                    statistic(neg_diva_weak_counter,12)=(acx(ii)*V(1,3)+acy(ii)*V(2,3)+acz(ii)*V(3,3))/(absac(ii));
                                    statistic(neg_diva_weak_counter,25)=(alx(ii)*ax(ii)+aly(ii)*ay(ii)+alz(ii)*az(ii))/(absal(ii)*absa(ii));
                                    statistic(neg_diva_weak_counter,26)=(acx(ii)*ax(ii)+acy(ii)*ay(ii)+acz(ii)*az(ii))/(absac(ii)*absa(ii));
                                    statistic(neg_diva_weak_counter,27)=(alx(ii)*acx(ii)+aly(ii)*acy(ii)+alz(ii)*acz(ii))/(absal(ii)*absac(ii));
                                end
                                if Q>5 %neg diva strong
                                    neg_diva_strong_counter=neg_diva_strong_counter+1;
                                    statistic(neg_diva_strong_counter,13)=(acx(ii)*V(1,1)+acy(ii)*V(2,1)+acz(ii)*V(3,1))/(absac(ii));
                                    statistic(neg_diva_strong_counter,14)=(acx(ii)*V(1,2)+acy(ii)*V(2,2)+acz(ii)*V(3,2))/(absac(ii));
                                    statistic(neg_diva_strong_counter,15)=(acx(ii)*V(1,3)+acy(ii)*V(2,3)+acz(ii)*V(3,3))/(absac(ii));
                                    statistic(neg_diva_strong_counter,28)=(alx(ii)*ax(ii)+aly(ii)*ay(ii)+alz(ii)*az(ii))/(absal(ii)*absa(ii));
                                    statistic(neg_diva_strong_counter,29)=(acx(ii)*ax(ii)+acy(ii)*ay(ii)+acz(ii)*az(ii))/(absac(ii)*absa(ii));
                                    statistic(neg_diva_strong_counter,30)=(alx(ii)*acx(ii)+aly(ii)*acy(ii)+alz(ii)*acz(ii))/(absal(ii)*absac(ii));
                                end
                            end                           
                        end
                    end
                end
            end
        end
    end
    
    %save resultOne statistic;
else
    %load resultOne;
end


figure;

subplot(1,2,1);
title('a_c')
hold on;
grid on;
box on;
[N,X] = nhist(statistic(1:counter,1),100); %%%cos(a,\lambda_1)
plot(X,N,'r');
[N,X] = nhist(statistic(1:counter,2),100); %%%cos(a,\lambda_2)
plot(X,N,'g');
[N,X] = nhist(statistic(1:counter,3),100); %%%cos(a,\lambda_3)
plot(X,N,'b');
xlabel('cos(a_c,\lambda_i)')
ylabel('pdf');

subplot(1,2,2);
hold on;
grid on;
box on;
[N,X] = nhist(statistic(1:counter,16),100); %%%cos(a,a_l)
plot(X,N,'r');
[N,X] = nhist(statistic(1:counter,17),100); %%%cos(a,a_c)
plot(X,N,'g');
[N,X] = nhist(statistic(1:counter,18),100); %%%cos(a_l,a_c)
plot(X,N,'b');
xlabel('r,g,b, cos(a,a_l), cos(a,a_c), cos(a_l,a_c)')
ylabel('pdf');


figure;

subplot(2,2,4);
hold on;
grid on;
box on;
[N,X] = nhist(statistic(1:pos_diva_pos_ua_counter,4),50); %%%cos(a,\lambda_1)
plot(X,N,'r');
[N,X] = nhist(statistic(1:pos_diva_pos_ua_counter,5),50); %%%cos(a,\lambda_2)
plot(X,N,'g');
[N,X] = nhist(statistic(1:pos_diva_pos_ua_counter,6),50); %%%cos(a,\lambda_3)
plot(X,N,'b');

subplot(2,2,2);
hold on;
grid on;
box on;
[N,X] = nhist(statistic(1:pos_diva_neg_ua_counter,7),50); %%%cos(a,\lambda_1)
plot(X,N,'r');
[N,X] = nhist(statistic(1:pos_diva_neg_ua_counter,8),50); %%%cos(a,\lambda_2)
plot(X,N,'g');
[N,X] = nhist(statistic(1:pos_diva_neg_ua_counter,9),50); %%%cos(a,\lambda_3)
plot(X,N,'b');

subplot(2,2,1);
title('a_c')
hold on;
grid on;
box on;
[N,X] = nhist(statistic(1:neg_diva_weak_counter,10),50); %%%cos(a,\lambda_1)
plot(X,N,'r');
[N,X] = nhist(statistic(1:neg_diva_weak_counter,11),50); %%%cos(a,\lambda_2)
plot(X,N,'g');
[N,X] = nhist(statistic(1:neg_diva_weak_counter,12),50); %%%cos(a,\lambda_3)
plot(X,N,'b');

subplot(2,2,3);
hold on;
grid on;
box on;
[N,X] = nhist(statistic(1:neg_diva_strong_counter,13),50); %%%cos(a,\lambda_1)
plot(X,N,'r');
[N,X] = nhist(statistic(1:neg_diva_strong_counter,14),50); %%%cos(a,\lambda_2)
plot(X,N,'g');
[N,X] = nhist(statistic(1:neg_diva_strong_counter,15),50); %%%cos(a,\lambda_3)
plot(X,N,'b');



figure;
subplot(2,2,4);
hold on;
grid on;
box on;
[N,X] = nhist(statistic(1:pos_diva_pos_ua_counter,19),50); %%%cos(a,\lambda_1)
plot(X,N,'r');
[N,X] = nhist(statistic(1:pos_diva_pos_ua_counter,20),50); %%%cos(a,\lambda_2)
plot(X,N,'g');
[N,X] = nhist(statistic(1:pos_diva_pos_ua_counter,21),50); %%%cos(a,\lambda_3)
plot(X,N,'b');

subplot(2,2,2);
hold on;
grid on;
box on;
[N,X] = nhist(statistic(1:pos_diva_neg_ua_counter,22),50); %%%cos(a,\lambda_1)
plot(X,N,'r');
[N,X] = nhist(statistic(1:pos_diva_neg_ua_counter,23),50); %%%cos(a,\lambda_2)
plot(X,N,'g');
[N,X] = nhist(statistic(1:pos_diva_neg_ua_counter,24),50); %%%cos(a,\lambda_3)
plot(X,N,'b');

subplot(2,2,1);
hold on;
grid on;
box on;
[N,X] = nhist(statistic(1:neg_diva_weak_counter,25),50); %%%cos(a,\lambda_1)
plot(X,N,'r');
[N,X] = nhist(statistic(1:neg_diva_weak_counter,26),50); %%%cos(a,\lambda_2)
plot(X,N,'g');
[N,X] = nhist(statistic(1:neg_diva_weak_counter,27),50); %%%cos(a,\lambda_3)
plot(X,N,'b');

subplot(2,2,3);
hold on;
grid on;
box on;
[N,X] = nhist(statistic(1:neg_diva_strong_counter,28),50); %%%cos(a,\lambda_1)
plot(X,N,'r');
[N,X] = nhist(statistic(1:neg_diva_strong_counter,29),50); %%%cos(a,\lambda_2)
plot(X,N,'g');
[N,X] = nhist(statistic(1:neg_diva_strong_counter,30),50); %%%cos(a,\lambda_3)
plot(X,N,'b');



counter
pos_diva_pos_ua_counter
pos_diva_neg_ua_counter
neg_diva_weak_counter
neg_diva_strong_counter
   
