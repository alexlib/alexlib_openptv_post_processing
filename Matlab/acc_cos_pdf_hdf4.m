function res=acc_cos_pdf_hdf4(first,last,tail,stepSize,threshold);
%acc_cos_pdf(10000,19900,0,20,1.37)
% modified to read from HDF files instead of trajAcc.*
% Nov. 21, 2007

% from readme.txt
%{
low= both strain & enstro < average (2/3 events, 1/3 of a^2)
strain= strain>enstro & strain>average (1/3 events, 1/3 of a^2)
enstro= enstro>strain & enstro>average (1/3 events, 1/3 of a^2)


the color code of all figures in this folder is as follows:

if only black, red, green, blue:
all,low, strain, enstro

if black, red,green, blue,cyan:
a,a_{\perp},a_{||},a_l,a_c

if black, blue,cyan:
then like above, but a_{\perp},a_{||} are meaningless
%}

if nargin == 0
    first = 1;
    last = 2990;
    tail = 3;
    stepSize = 1;
    threshold = 1.37;
end

counter=0;
count_low=0;
count_enstro=0;
count_strain=0;

hist_absa=zeros(700000,4);
hist_absao=zeros(700000,4);
hist_absap=zeros(700000,4);
hist_absal=zeros(700000,4);
hist_absac=zeros(700000,4);
hist_curv=zeros(700000,4);

cos_alac=zeros(700000,4);
cos_apac=zeros(700000,4);
ratio_apac=zeros(700000,4);

cos_a_w=zeros(700000,4);
cos_ao_w=zeros(700000,4);
cos_ap_w=zeros(700000,4);
cos_ac_w=zeros(700000,4);
cos_al_w=zeros(700000,4);

cos_a_l1=zeros(700000,4);
cos_ao_l1=zeros(700000,4);
cos_ap_l1=zeros(700000,4);
cos_ac_l1=zeros(700000,4);
cos_al_l1=zeros(700000,4);

cos_a_l2=zeros(700000,4);
cos_ao_l2=zeros(700000,4);
cos_ap_l2=zeros(700000,4);
cos_ac_l2=zeros(700000,4);
cos_al_l2=zeros(700000,4);

cos_a_l3=zeros(700000,4);
cos_ao_l3=zeros(700000,4);
cos_ap_l3=zeros(700000,4);
cos_ac_l3=zeros(700000,4);
cos_al_l3=zeros(700000,4);

cos_a_u=zeros(700000,4);
cos_ao_u=zeros(700000,4);
cos_ap_u=zeros(700000,4);
cos_ac_u=zeros(700000,4);
cos_al_u=zeros(700000,4);

mean_s=3.4;
allpoints=0;




%% Assign the working data
curdir = cd;
hfiles{1} = [curdir(1),':\Matlab\Alex\HDF\0104_water.hdf'];
% or
hfiles{2} = [curdir(1),':\Matlab\Alex\HDF\0104_hompol.hdf'];

for i=first:stepSize:last
    %     clear f;
    %     i
    %     name=['D:/data/Riso_micro_PTV/trajPoint.',num2Str(i)];
    %     f=load(name);
    %     s=size(f);
    %     if s(1,1)>0
    %         xr=f(:,1);
    %         yr=f(:,2);
    %         zr=f(:,3);
    %         u=f(:,4);
    %         v=f(:,5);
    %         w=f(:,6);
    %         absu=(u.^2+v.^2+w.^2).^0.5;
    %         ax=f(:,7);%%%%%%%%%%%
    %         ay=f(:,8);%%%%%%%%%%%
    %         az=f(:,9);%%%%%%%%%%%
    %         absa=(ax.^2+ay.^2+az.^2).^0.5;%%%%%%%%%%%%%
    %         alx=f(:,19);%%%%%%%%%%%
    %         aly=f(:,20);%%%%%%%%%%%
    %         alz=f(:,21);%%%%%%%%%%%
    %         absal=(alx.^2+aly.^2+alz.^2).^0.5;%%%%%%%%%%%%%%%
    %         w1=f(:,10);
    %         w2=f(:,11);
    %         w3=f(:,12);
    %         absw=(w1.^2+w2.^2+w3.^2).^0.5;
    %         enstro=w1.^2+w2.^2+w3.^2;
    %         s11=f(:,13);
    %         s12=f(:,14);
    %         s13=f(:,15);
    %         s22=f(:,16);
    %         s23=f(:,17);
    %         s33=f(:,18);
    %         ux=s11;
    %         uy=s12-0.5*w3;
    %         uz=s13+0.5*w2;
    %         vx=s12+0.5*w3;
    %         vy=s22;
    %         vz=s23-0.5*w1;
    %         wx=s13-0.5*w2;
    %         wy=s23+0.5*w1;
    %         wz=s33;
    %         acx=u.*ux+v.*uy+w.*uz;%%%%%%%%%%%%%%%%%%%%
    %         acy=u.*vx+v.*vy+w.*vz;%%%%%%%%%%%%%%%%%%%%
    %         acz=u.*wx+v.*wy+w.*wz;%%%%%%%%%%%%%%%%%%%%
    %         absac=(acx.^2+acy.^2+acz.^2).^0.5;%%%%%%%%%%%%%%%%%%%%%%%%%
    %         aLamx=w2.*w-w3.*v;%%%%%%%%%%%%%%%%%%%%%%%
    %         aLamy=w3.*u-w1.*w;%%%%%%%%%%%%%%%%%%%%%%%
    %         aLamz=w1.*v-w2.*u;%%%%%%%%%%%%%%%%%%%%%%%%%
    %         absaLam=(aLamx.^2+aLamy.^2+aLamz.^2).^0.5;%%%%%%%%%%%%%%%%%
    %         aBx=acx-aLamx;%%%%%%%%%%%%%%%%%%%%%
    %         aBy=acy-aLamy;%%%%%%%%%%%%%%%%%%%%%%%
    %         aBz=acz-aLamz;%%%%%%%%%%%%%%%%%%%%%%%%%
    %         absaB=(aBx.^2+aBy.^2+aBz.^2).^0.5;%%%%%%%%%%%%%%%%%%%
    %         unx=u./absu;
    %         uny=v./absu;
    %         unz=w./absu;
    %         aun=ax.*unx+ay.*uny+az.*unz;
    %         aPx=aun.*unx;%%%%%%%%%%%%%%%%%%%%%%%%
    %         aPy=aun.*uny;%%%%%%%%%%%%%%%%%%%%%%%%%
    %         aPz=aun.*unz;%%%%%%%%%%%%%%%%%%%%%%%%%%
    %         absaP=(aPx.^2+aPy.^2+aPz.^2).^0.5;%%%%%%%%%%%%%%%%%%%%%
    %         aOx=ax-aPx;%%%%%%%%%%%%%%%%%%%%
    %         aOy=ay-aPy;%%%%%%%%%%%%%%%%%%%
    %         aOz=az-aPz;%%%%%%%%%%%%%%%%%%%%
    %         absaO=(aOx.^2+aOy.^2+aOz.^2).^0.5;%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %         cux=(u.*ux+v.*uy+w.*uz)./(absu.^2);
    %         cuy=(u.*vx+v.*vy+w.*vz)./(absu.^2);
    %         cuz=(u.*wx+v.*wy+w.*wz)./(absu.^2);
    %         curvGrad=(cux.^2+cuy.^2+cuz.^2).^0.5;
    %         curv=(((v.*az-w.*ay).^2+(w.*ax-u.*az).^2+(u.*ay-v.*ax).^2).^0.5)./(absu.^3);
    %
    %         ua=u.*acx+v.*acy+w.*acz;
    %
    %         daxdx=f(:,22);
    %         daxdy=f(:,23);
    %         daxdz=f(:,24);
    %         daydx=f(:,25);
    %         daydy=f(:,26);
    %         daydz=f(:,27);
    %         dazdx=f(:,28);
    %         dazdy=f(:,29);
    %         dazdz=f(:,30);
    %
    %         reldiv=f(:,31);
    %         age=f(:,32);
    %
    %         %       so we have:
    %         %       absu,     1
    %         %       absa,     2
    %         %       absal,    3
    %         %       absac,    4
    %         %       absaLam,  5
    %         %       absaB,    6
    %         %       absaP,    7
    %         %       absaO,    8
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
        enstro=w1.^2+w2.^2+w3.^2;
        absw=enstro.^0.5;
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
        curv=(((v.*az-w.*ay).^2+(w.*ax-u.*az).^2+(u.*ay-v.*ax).^2).^0.5)./(absu.^3);
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

        ende=0;
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
                    allpoints=allpoints+1;
                    if age(ii)>tail & age(ii)<leng-tail+1 & reldiv(ii)<0.2 & absa(ii)^2/8e-4<100 & curv(ii)<3000
                        uM=[ux(ii) uy(ii) uz(ii);vx(ii) vy(ii) vz(ii);wx(ii) wy(ii) wz(ii)];
                        A=[s11(ii) s12(ii) s13(ii);s12(ii) s22(ii) s23(ii);s13(ii) s23(ii) s33(ii)];
                        [V,D] = eig(A);
                        [D,k] = sort(diag(D)); 	% ascending order
                        D = diag(D(end:-1:1)); 			% descending order
                        V = V(:,k(end:-1:1)); 	% the same order for eigenvectors
                        strain=D(1,1)^2+D(2,2)^2+D(3,3)^2;


                        if counter<700000+1 %& ua(ii)<0
                            counter=counter+1;
                            hist_absa(counter,1)=absa(ii);
                            hist_absao(counter,1)=absaO(ii);
                            hist_absap(counter,1)=absaP(ii);
                            hist_absal(counter,1)=absal(ii);
                            hist_absac(counter,1)=absac(ii);
                            hist_curv(counter,1)=curv(ii);

                            cos_a_w(counter,1)=[ax(ii) ay(ii) az(ii)]*[w1(ii) w2(ii) w3(ii)]'/(absa(ii)*absw(ii));
                            cos_ao_w(counter,1)=[aOx(ii) aOy(ii) aOz(ii)]*[w1(ii) w2(ii) w3(ii)]'/(absaO(ii)*absw(ii));
                            cos_ap_w(counter,1)=[aPx(ii) aPy(ii) aPz(ii)]*[w1(ii) w2(ii) w3(ii)]'/(absaP(ii)*absw(ii));
                            cos_ac_w(counter,1)=[acx(ii) acy(ii) acz(ii)]*[w1(ii) w2(ii) w3(ii)]'/(absac(ii)*absw(ii));
                            cos_al_w(counter,1)=[alx(ii) aly(ii) alz(ii)]*[w1(ii) w2(ii) w3(ii)]'/(absal(ii)*absw(ii));

                            cos_a_l1(counter,1)=[ax(ii) ay(ii) az(ii)]*V(:,1)/(absa(ii)*(sum(V(:,1).^2))^0.5);
                            cos_ao_l1(counter,1)=[aOx(ii) aOy(ii) aOz(ii)]*V(:,1)/(absaO(ii)*(sum(V(:,1).^2))^0.5);
                            cos_ap_l1(counter,1)=[aPx(ii) aPy(ii) aPz(ii)]*V(:,1)/(absaP(ii)*(sum(V(:,1).^2))^0.5);
                            cos_ac_l1(counter,1)=[acx(ii) acy(ii) acz(ii)]*V(:,1)/(absac(ii)*(sum(V(:,1).^2))^0.5);
                            cos_al_l1(counter,1)=[alx(ii) aly(ii) alz(ii)]*V(:,1)/(absal(ii)*(sum(V(:,1).^2))^0.5);

                            cos_a_l2(counter,1)=[ax(ii) ay(ii) az(ii)]*V(:,2)/(absa(ii)*(sum(V(:,2).^2))^0.5);
                            cos_ao_l2(counter,1)=[aOx(ii) aOy(ii) aOz(ii)]*V(:,2)/(absaO(ii)*(sum(V(:,2).^2))^0.5);
                            cos_ap_l2(counter,1)=[aPx(ii) aPy(ii) aPz(ii)]*V(:,2)/(absaP(ii)*(sum(V(:,2).^2))^0.5);
                            cos_ac_l2(counter,1)=[acx(ii) acy(ii) acz(ii)]*V(:,2)/(absac(ii)*(sum(V(:,2).^2))^0.5);
                            cos_al_l2(counter,1)=[alx(ii) aly(ii) alz(ii)]*V(:,2)/(absal(ii)*(sum(V(:,2).^2))^0.5);

                            cos_a_l3(counter,1)=[ax(ii) ay(ii) az(ii)]*V(:,3)/(absa(ii)*(sum(V(:,3).^2))^0.5);
                            cos_ao_l3(counter,1)=[aOx(ii) aOy(ii) aOz(ii)]*V(:,3)/(absaO(ii)*(sum(V(:,3).^2))^0.5);
                            cos_ap_l3(counter,1)=[aPx(ii) aPy(ii) aPz(ii)]*V(:,3)/(absaP(ii)*(sum(V(:,3).^2))^0.5);
                            cos_ac_l3(counter,1)=[acx(ii) acy(ii) acz(ii)]*V(:,3)/(absac(ii)*(sum(V(:,3).^2))^0.5);
                            cos_al_l3(counter,1)=[alx(ii) aly(ii) alz(ii)]*V(:,3)/(absal(ii)*(sum(V(:,3).^2))^0.5);

                            cos_a_u(counter,1)=[ax(ii) ay(ii) az(ii)]*[u(ii) v(ii) w(ii)]'/(absa(ii)*absu(ii));
                            cos_ao_u(counter,1)=[aOx(ii) aOy(ii) aOz(ii)]*[u(ii) v(ii) w(ii)]'/(absaO(ii)*absu(ii));
                            cos_ap_u(counter,1)=[aPx(ii) aPy(ii) aPz(ii)]*[u(ii) v(ii) w(ii)]'/(absaP(ii)*absu(ii));
                            cos_ac_u(counter,1)=[acx(ii) acy(ii) acz(ii)]*[u(ii) v(ii) w(ii)]'/(absac(ii)*absu(ii));
                            cos_al_u(counter,1)=[alx(ii) aly(ii) alz(ii)]*[u(ii) v(ii) w(ii)]'/(absal(ii)*absu(ii));

                            cos_alac(counter,1)=[alx(ii) aly(ii) alz(ii)]*[acx(ii) acy(ii) acz(ii)]'/(absal(ii)*absac(ii));
                            cos_apac(counter,1)=[aPx(ii) aPy(ii) aPz(ii)]*[acx(ii) acy(ii) acz(ii)]'/(absaP(ii)*absac(ii));
                            ratio_apac(counter,1)=absaP(ii)/absac(ii);
                        end

                        if enstro(ii)<10^threshold & 2*strain<10^threshold %& ua(ii)<0
                            if count_low<700000+1
                                count_low=count_low+1;
                                hist_absa(count_low,2)=absa(ii);
                                hist_absao(count_low,2)=absaO(ii);
                                hist_absap(count_low,2)=absaP(ii);
                                hist_absal(count_low,2)=absal(ii);
                                hist_absac(count_low,2)=absac(ii);
                                hist_curv(count_low,2)=curv(ii);

                                cos_a_w(count_low,2)=[ax(ii) ay(ii) az(ii)]*[w1(ii) w2(ii) w3(ii)]'/(absa(ii)*absw(ii));
                                cos_ao_w(count_low,2)=[aOx(ii) aOy(ii) aOz(ii)]*[w1(ii) w2(ii) w3(ii)]'/(absaO(ii)*absw(ii));
                                cos_ap_w(count_low,2)=[aPx(ii) aPy(ii) aPz(ii)]*[w1(ii) w2(ii) w3(ii)]'/(absaP(ii)*absw(ii));
                                cos_ac_w(count_low,2)=[acx(ii) acy(ii) acz(ii)]*[w1(ii) w2(ii) w3(ii)]'/(absac(ii)*absw(ii));
                                cos_al_w(count_low,2)=[alx(ii) aly(ii) alz(ii)]*[w1(ii) w2(ii) w3(ii)]'/(absal(ii)*absw(ii));

                                cos_a_l1(count_low,2)=[ax(ii) ay(ii) az(ii)]*V(:,1)/(absa(ii)*(sum(V(:,1).^2))^0.5);
                                cos_ao_l1(count_low,2)=[aOx(ii) aOy(ii) aOz(ii)]*V(:,1)/(absaO(ii)*(sum(V(:,1).^2))^0.5);
                                cos_ap_l1(count_low,2)=[aPx(ii) aPy(ii) aPz(ii)]*V(:,1)/(absaP(ii)*(sum(V(:,1).^2))^0.5);
                                cos_ac_l1(count_low,2)=[acx(ii) acy(ii) acz(ii)]*V(:,1)/(absac(ii)*(sum(V(:,1).^2))^0.5);
                                cos_al_l1(count_low,2)=[alx(ii) aly(ii) alz(ii)]*V(:,1)/(absal(ii)*(sum(V(:,1).^2))^0.5);

                                cos_a_l2(count_low,2)=[ax(ii) ay(ii) az(ii)]*V(:,2)/(absa(ii)*(sum(V(:,2).^2))^0.5);
                                cos_ao_l2(count_low,2)=[aOx(ii) aOy(ii) aOz(ii)]*V(:,2)/(absaO(ii)*(sum(V(:,2).^2))^0.5);
                                cos_ap_l2(count_low,2)=[aPx(ii) aPy(ii) aPz(ii)]*V(:,2)/(absaP(ii)*(sum(V(:,2).^2))^0.5);
                                cos_ac_l2(count_low,2)=[acx(ii) acy(ii) acz(ii)]*V(:,2)/(absac(ii)*(sum(V(:,2).^2))^0.5);
                                cos_al_l2(count_low,2)=[alx(ii) aly(ii) alz(ii)]*V(:,2)/(absal(ii)*(sum(V(:,2).^2))^0.5);

                                cos_a_l3(count_low,2)=[ax(ii) ay(ii) az(ii)]*V(:,3)/(absa(ii)*(sum(V(:,3).^2))^0.5);
                                cos_ao_l3(count_low,2)=[aOx(ii) aOy(ii) aOz(ii)]*V(:,3)/(absaO(ii)*(sum(V(:,3).^2))^0.5);
                                cos_ap_l3(count_low,2)=[aPx(ii) aPy(ii) aPz(ii)]*V(:,3)/(absaP(ii)*(sum(V(:,3).^2))^0.5);
                                cos_ac_l3(count_low,2)=[acx(ii) acy(ii) acz(ii)]*V(:,3)/(absac(ii)*(sum(V(:,3).^2))^0.5);
                                cos_al_l3(count_low,2)=[alx(ii) aly(ii) alz(ii)]*V(:,3)/(absal(ii)*(sum(V(:,3).^2))^0.5);

                                cos_a_u(count_low,2)=[ax(ii) ay(ii) az(ii)]*[u(ii) v(ii) w(ii)]'/(absa(ii)*absu(ii));
                                cos_ao_u(count_low,2)=[aOx(ii) aOy(ii) aOz(ii)]*[u(ii) v(ii) w(ii)]'/(absaO(ii)*absu(ii));
                                cos_ap_u(count_low,2)=[aPx(ii) aPy(ii) aPz(ii)]*[u(ii) v(ii) w(ii)]'/(absaP(ii)*absu(ii));
                                cos_ac_u(count_low,2)=[acx(ii) acy(ii) acz(ii)]*[u(ii) v(ii) w(ii)]'/(absac(ii)*absu(ii));
                                cos_al_u(count_low,2)=[alx(ii) aly(ii) alz(ii)]*[u(ii) v(ii) w(ii)]'/(absal(ii)*absu(ii));

                                cos_alac(count_low,2)=[alx(ii) aly(ii) alz(ii)]*[acx(ii) acy(ii) acz(ii)]'/(absal(ii)*absac(ii));
                                cos_apac(count_low,2)=[aPx(ii) aPy(ii) aPz(ii)]*[acx(ii) acy(ii) acz(ii)]'/(absaP(ii)*absac(ii));
                                ratio_apac(count_low,2)=absaP(ii)/absac(ii);
                            end
                        end
                        if 2*strain<enstro(ii) & enstro(ii)>10^threshold %& ua(ii)<0
                            if count_enstro<700000+1
                                count_enstro=count_enstro+1;
                                hist_absa(count_enstro,3)=absa(ii);
                                hist_absao(count_enstro,3)=absaO(ii);
                                hist_absap(count_enstro,3)=absaP(ii);
                                hist_absal(count_enstro,3)=absal(ii);
                                hist_absac(count_enstro,3)=absac(ii);
                                hist_curv(count_enstro,3)=curv(ii);

                                cos_a_w(count_enstro,3)=[ax(ii) ay(ii) az(ii)]*[w1(ii) w2(ii) w3(ii)]'/(absa(ii)*absw(ii));
                                cos_ao_w(count_enstro,3)=[aOx(ii) aOy(ii) aOz(ii)]*[w1(ii) w2(ii) w3(ii)]'/(absaO(ii)*absw(ii));
                                cos_ap_w(count_enstro,3)=[aPx(ii) aPy(ii) aPz(ii)]*[w1(ii) w2(ii) w3(ii)]'/(absaP(ii)*absw(ii));
                                cos_ac_w(count_enstro,3)=[acx(ii) acy(ii) acz(ii)]*[w1(ii) w2(ii) w3(ii)]'/(absac(ii)*absw(ii));
                                cos_al_w(count_enstro,3)=[alx(ii) aly(ii) alz(ii)]*[w1(ii) w2(ii) w3(ii)]'/(absal(ii)*absw(ii));

                                cos_a_l1(count_enstro,3)=[ax(ii) ay(ii) az(ii)]*V(:,1)/(absa(ii)*(sum(V(:,1).^2))^0.5);
                                cos_ao_l1(count_enstro,3)=[aOx(ii) aOy(ii) aOz(ii)]*V(:,1)/(absaO(ii)*(sum(V(:,1).^2))^0.5);
                                cos_ap_l1(count_enstro,3)=[aPx(ii) aPy(ii) aPz(ii)]*V(:,1)/(absaP(ii)*(sum(V(:,1).^2))^0.5);
                                cos_ac_l1(count_enstro,3)=[acx(ii) acy(ii) acz(ii)]*V(:,1)/(absac(ii)*(sum(V(:,1).^2))^0.5);
                                cos_al_l1(count_enstro,3)=[alx(ii) aly(ii) alz(ii)]*V(:,1)/(absal(ii)*(sum(V(:,1).^2))^0.5);

                                cos_a_l2(count_enstro,3)=[ax(ii) ay(ii) az(ii)]*V(:,2)/(absa(ii)*(sum(V(:,2).^2))^0.5);
                                cos_ao_l2(count_enstro,3)=[aOx(ii) aOy(ii) aOz(ii)]*V(:,2)/(absaO(ii)*(sum(V(:,2).^2))^0.5);
                                cos_ap_l2(count_enstro,3)=[aPx(ii) aPy(ii) aPz(ii)]*V(:,2)/(absaP(ii)*(sum(V(:,2).^2))^0.5);
                                cos_ac_l2(count_enstro,3)=[acx(ii) acy(ii) acz(ii)]*V(:,2)/(absac(ii)*(sum(V(:,2).^2))^0.5);
                                cos_al_l2(count_enstro,3)=[alx(ii) aly(ii) alz(ii)]*V(:,2)/(absal(ii)*(sum(V(:,2).^2))^0.5);

                                cos_a_l3(count_enstro,3)=[ax(ii) ay(ii) az(ii)]*V(:,3)/(absa(ii)*(sum(V(:,3).^2))^0.5);
                                cos_ao_l3(count_enstro,3)=[aOx(ii) aOy(ii) aOz(ii)]*V(:,3)/(absaO(ii)*(sum(V(:,3).^2))^0.5);
                                cos_ap_l3(count_enstro,3)=[aPx(ii) aPy(ii) aPz(ii)]*V(:,3)/(absaP(ii)*(sum(V(:,3).^2))^0.5);
                                cos_ac_l3(count_enstro,3)=[acx(ii) acy(ii) acz(ii)]*V(:,3)/(absac(ii)*(sum(V(:,3).^2))^0.5);
                                cos_al_l3(count_enstro,3)=[alx(ii) aly(ii) alz(ii)]*V(:,3)/(absal(ii)*(sum(V(:,3).^2))^0.5);

                                cos_a_u(count_enstro,3)=[ax(ii) ay(ii) az(ii)]*[u(ii) v(ii) w(ii)]'/(absa(ii)*absu(ii));
                                cos_ao_u(count_enstro,3)=[aOx(ii) aOy(ii) aOz(ii)]*[u(ii) v(ii) w(ii)]'/(absaO(ii)*absu(ii));
                                cos_ap_u(count_enstro,3)=[aPx(ii) aPy(ii) aPz(ii)]*[u(ii) v(ii) w(ii)]'/(absaP(ii)*absu(ii));
                                cos_ac_u(count_enstro,3)=[acx(ii) acy(ii) acz(ii)]*[u(ii) v(ii) w(ii)]'/(absac(ii)*absu(ii));
                                cos_al_u(count_enstro,3)=[alx(ii) aly(ii) alz(ii)]*[u(ii) v(ii) w(ii)]'/(absal(ii)*absu(ii));

                                cos_alac(count_enstro,3)=[alx(ii) aly(ii) alz(ii)]*[acx(ii) acy(ii) acz(ii)]'/(absal(ii)*absac(ii));
                                cos_apac(count_enstro,3)=[aPx(ii) aPy(ii) aPz(ii)]*[acx(ii) acy(ii) acz(ii)]'/(absaP(ii)*absac(ii));
                                ratio_apac(count_enstro,3)=absaP(ii)/absac(ii);
                            end
                        end
                        if 2*strain>enstro(ii) & 2*strain>10^threshold %& ua(ii)<0
                            if count_strain<700000+1
                                count_strain=count_strain+1;
                                hist_absa(count_strain,4)=absa(ii);
                                hist_absao(count_strain,4)=absaO(ii);
                                hist_absap(count_strain,4)=absaP(ii);
                                hist_absal(count_strain,4)=absal(ii);
                                hist_absac(count_strain,4)=absac(ii);
                                hist_curv(count_strain,4)=curv(ii);

                                cos_a_w(count_strain,4)=[ax(ii) ay(ii) az(ii)]*[w1(ii) w2(ii) w3(ii)]'/(absa(ii)*absw(ii));
                                cos_ao_w(count_strain,4)=[aOx(ii) aOy(ii) aOz(ii)]*[w1(ii) w2(ii) w3(ii)]'/(absaO(ii)*absw(ii));
                                cos_ap_w(count_strain,4)=[aPx(ii) aPy(ii) aPz(ii)]*[w1(ii) w2(ii) w3(ii)]'/(absaP(ii)*absw(ii));
                                cos_ac_w(count_strain,4)=[acx(ii) acy(ii) acz(ii)]*[w1(ii) w2(ii) w3(ii)]'/(absac(ii)*absw(ii));
                                cos_al_w(count_strain,4)=[alx(ii) aly(ii) alz(ii)]*[w1(ii) w2(ii) w3(ii)]'/(absal(ii)*absw(ii));

                                cos_a_l1(count_strain,4)=[ax(ii) ay(ii) az(ii)]*V(:,1)/(absa(ii)*(sum(V(:,1).^2))^0.5);
                                cos_ao_l1(count_strain,4)=[aOx(ii) aOy(ii) aOz(ii)]*V(:,1)/(absaO(ii)*(sum(V(:,1).^2))^0.5);
                                cos_ap_l1(count_strain,4)=[aPx(ii) aPy(ii) aPz(ii)]*V(:,1)/(absaP(ii)*(sum(V(:,1).^2))^0.5);
                                cos_ac_l1(count_strain,4)=[acx(ii) acy(ii) acz(ii)]*V(:,1)/(absac(ii)*(sum(V(:,1).^2))^0.5);
                                cos_al_l1(count_strain,4)=[alx(ii) aly(ii) alz(ii)]*V(:,1)/(absal(ii)*(sum(V(:,1).^2))^0.5);

                                cos_a_l2(count_strain,4)=[ax(ii) ay(ii) az(ii)]*V(:,2)/(absa(ii)*(sum(V(:,2).^2))^0.5);
                                cos_ao_l2(count_strain,4)=[aOx(ii) aOy(ii) aOz(ii)]*V(:,2)/(absaO(ii)*(sum(V(:,2).^2))^0.5);
                                cos_ap_l2(count_strain,4)=[aPx(ii) aPy(ii) aPz(ii)]*V(:,2)/(absaP(ii)*(sum(V(:,2).^2))^0.5);
                                cos_ac_l2(count_strain,4)=[acx(ii) acy(ii) acz(ii)]*V(:,2)/(absac(ii)*(sum(V(:,2).^2))^0.5);
                                cos_al_l2(count_strain,4)=[alx(ii) aly(ii) alz(ii)]*V(:,2)/(absal(ii)*(sum(V(:,2).^2))^0.5);

                                cos_a_l3(count_strain,4)=[ax(ii) ay(ii) az(ii)]*V(:,3)/(absa(ii)*(sum(V(:,3).^2))^0.5);
                                cos_ao_l3(count_strain,4)=[aOx(ii) aOy(ii) aOz(ii)]*V(:,3)/(absaO(ii)*(sum(V(:,3).^2))^0.5);
                                cos_ap_l3(count_strain,4)=[aPx(ii) aPy(ii) aPz(ii)]*V(:,3)/(absaP(ii)*(sum(V(:,3).^2))^0.5);
                                cos_ac_l3(count_strain,4)=[acx(ii) acy(ii) acz(ii)]*V(:,3)/(absac(ii)*(sum(V(:,3).^2))^0.5);
                                cos_al_l3(count_strain,4)=[alx(ii) aly(ii) alz(ii)]*V(:,3)/(absal(ii)*(sum(V(:,3).^2))^0.5);

                                cos_a_u(count_strain,4)=[ax(ii) ay(ii) az(ii)]*[u(ii) v(ii) w(ii)]'/(absa(ii)*absu(ii));
                                cos_ao_u(count_strain,4)=[aOx(ii) aOy(ii) aOz(ii)]*[u(ii) v(ii) w(ii)]'/(absaO(ii)*absu(ii));
                                cos_ap_u(count_strain,4)=[aPx(ii) aPy(ii) aPz(ii)]*[u(ii) v(ii) w(ii)]'/(absaP(ii)*absu(ii));
                                cos_ac_u(count_strain,4)=[acx(ii) acy(ii) acz(ii)]*[u(ii) v(ii) w(ii)]'/(absac(ii)*absu(ii));
                                cos_al_u(count_strain,4)=[alx(ii) aly(ii) alz(ii)]*[u(ii) v(ii) w(ii)]'/(absal(ii)*absu(ii));

                                cos_alac(count_strain,4)=[alx(ii) aly(ii) alz(ii)]*[acx(ii) acy(ii) acz(ii)]'/(absal(ii)*absac(ii));
                                cos_apac(count_strain,4)=[aPx(ii) aPy(ii) aPz(ii)]*[acx(ii) acy(ii) acz(ii)]'/(absaP(ii)*absac(ii));
                                ratio_apac(count_strain,4)=absaP(ii)/absac(ii);
                            end
                        end
                    end
                end
            end
        end
    end
end


figure;
hold on;
grid on;
box on;
[n,x]=nhist(hist_absa(1:counter,1));
plot(x,n,'k');
[n,x]=nhist(hist_absa(1:count_low,2));
plot(x,n,'r');
[n,x]=nhist(hist_absa(1:count_enstro,3));
plot(x,n,'g');
[n,x]=nhist(hist_absa(1:count_strain,4));
plot(x,n,'b');
xlabel('|a|','FontSize',12,'FontName','Times New Roman')
ylabel('pdf','FontSize',12,'FontName','Times New Roman')
set(gca,'XScale','log')
set(gca,'YScale','log')

figure;
hold on;
grid on;
box on;
[n,x]=nhist(hist_absao(1:counter,1),1000);
plot(x,n,'k');
[n,x]=nhist(hist_absao(1:count_low,2),1000);
plot(x,n,'r');
[n,x]=nhist(hist_absao(1:count_enstro,3),1000);
plot(x,n,'g');
[n,x]=nhist(hist_absao(1:count_strain,4),1000);
plot(x,n,'b');
xlabel('|a_{\perp}|','FontSize',12,'FontName','Times New Roman')
ylabel('pdf','FontSize',12,'FontName','Times New Roman')
set(gca,'XScale','log')
set(gca,'YScale','log')

figure;
hold on;
grid on;
box on;
[n,x]=nhist(hist_absap(1:counter,1),1000);
plot(x,n,'k');
[n,x]=nhist(hist_absap(1:count_low,2),1000);
plot(x,n,'r');
[n,x]=nhist(hist_absap(1:count_enstro,3),1000);
plot(x,n,'g');
[n,x]=nhist(hist_absap(1:count_strain,4),1000);
plot(x,n,'b');
xlabel('|a_{||}|','FontSize',12,'FontName','Times New Roman')
ylabel('pdf','FontSize',12,'FontName','Times New Roman')
set(gca,'XScale','log')
set(gca,'YScale','log')

figure;
hold on;
grid on;
box on;
[n,x]=nhist(hist_absal(1:counter,1),1000);
plot(x,n,'k');
[n,x]=nhist(hist_absal(1:count_low,2),1000);
plot(x,n,'r');
[n,x]=nhist(hist_absal(1:count_enstro,3),1000);
plot(x,n,'g');
[n,x]=nhist(hist_absal(1:count_strain,4),1000);
plot(x,n,'b');
xlabel('|a_l|','FontSize',12,'FontName','Times New Roman')
ylabel('pdf','FontSize',12,'FontName','Times New Roman')
set(gca,'XScale','log')
set(gca,'YScale','log')

figure;
hold on;
grid on;
box on;
[n,x]=nhist(hist_absac(1:counter,1),1000);
plot(x,n,'k');
[n,x]=nhist(hist_absac(1:count_low,2),1000);
plot(x,n,'r');
[n,x]=nhist(hist_absac(1:count_enstro,3),1000);
plot(x,n,'g');
[n,x]=nhist(hist_absac(1:count_strain,4),1000);
plot(x,n,'b');
xlabel('|a_c|','FontSize',12,'FontName','Times New Roman')
ylabel('pdf','FontSize',12,'FontName','Times New Roman')
set(gca,'XScale','log')
set(gca,'YScale','log')

figure;
hold on;
grid on;
box on;
[n,x]=nhist(hist_curv(1:counter,1),1000);
plot(x,n,'k');
[n,x]=nhist(hist_curv(1:count_low,2),1000);
plot(x,n,'r');
[n,x]=nhist(hist_curv(1:count_enstro,3),1000);
plot(x,n,'g');
[n,x]=nhist(hist_curv(1:count_strain,4),1000);
plot(x,n,'b');
xlabel('curvature','FontSize',12,'FontName','Times New Roman')
ylabel('pdf','FontSize',12,'FontName','Times New Roman')
set(gca,'XScale','log')
set(gca,'YScale','log')



figure;
hold on;
grid on;
box on;
[n,x]=nhist(cos_a_w(1:counter,1),50);
plot(x,n,'k');
[n,x]=nhist(cos_ao_w(1:counter,1),50);
plot(x,n,'r');
[n,x]=nhist(cos_ap_w(1:counter,1),50);
plot(x,n,'g');
[n,x]=nhist(cos_al_w(1:counter,1),50);
plot(x,n,'b');
[n,x]=nhist(cos_ac_w(1:counter,1),50);
plot(x,n,'c');
xlabel('cos(a,\omega) all','FontSize',12,'FontName','Times New Roman')
ylabel('pdf','FontSize',12,'FontName','Times New Roman')
axis([-1 1 0 2])

figure;
hold on;
grid on;
box on;
[n,x]=nhist(cos_a_w(1:count_low,2),50);
plot(x,n,'k');
[n,x]=nhist(cos_ao_w(1:count_low,2),50);
plot(x,n,'r');
[n,x]=nhist(cos_ap_w(1:count_low,2),50);
plot(x,n,'g');
[n,x]=nhist(cos_al_w(1:count_low,2),50);
plot(x,n,'b');
[n,x]=nhist(cos_ac_w(1:count_low,2),50);
plot(x,n,'c');
xlabel('cos(a,\omega) low','FontSize',12,'FontName','Times New Roman')
ylabel('pdf','FontSize',12,'FontName','Times New Roman')
axis([-1 1 0 2])

figure;
hold on;
grid on;
box on;
[n,x]=nhist(cos_a_w(1:count_enstro,3),50);
plot(x,n,'k');
[n,x]=nhist(cos_ao_w(1:count_enstro,3),50);
plot(x,n,'r');
[n,x]=nhist(cos_ap_w(1:count_enstro,3),50);
plot(x,n,'g');
[n,x]=nhist(cos_al_w(1:count_enstro,3),50);
plot(x,n,'b');
[n,x]=nhist(cos_ac_w(1:count_enstro,3),50);
plot(x,n,'c');
xlabel('cos(a,\omega) enstro','FontSize',12,'FontName','Times New Roman')
ylabel('pdf','FontSize',12,'FontName','Times New Roman')
axis([-1 1 0 2])

figure;
hold on;
grid on;
box on;
[n,x]=nhist(cos_a_w(1:count_strain,4),50);
plot(x,n,'k');
[n,x]=nhist(cos_ao_w(1:count_strain,4),50);
plot(x,n,'r');
[n,x]=nhist(cos_ap_w(1:count_strain,4),50);
plot(x,n,'g');
[n,x]=nhist(cos_al_w(1:count_strain,4),50);
plot(x,n,'b');
[n,x]=nhist(cos_ac_w(1:count_strain,4),50);
plot(x,n,'c');
xlabel('cos(a,\omega) strain','FontSize',12,'FontName','Times New Roman')
ylabel('pdf','FontSize',12,'FontName','Times New Roman')
axis([-1 1 0 2])

figure;
hold on;
grid on;
box on;
[n,x]=nhist(cos_a_u(1:counter,1),50);
plot(x,n,'k');
[n,x]=nhist(cos_al_u(1:counter,1),50);
plot(x,n,'b');
[n,x]=nhist(cos_ac_u(1:counter,1),50);
plot(x,n,'c');
xlabel('cos(u,a) all','FontSize',12,'FontName','Times New Roman')
ylabel('pdf','FontSize',12,'FontName','Times New Roman')
axis([-1 1 0 2])

figure;
hold on;
grid on;
box on;
[n,x]=nhist(cos_a_u(1:count_low,2),50);
plot(x,n,'k');
[n,x]=nhist(cos_al_u(1:count_low,2),50);
plot(x,n,'b');
[n,x]=nhist(cos_ac_u(1:count_low,2),50);
plot(x,n,'c');
xlabel('cos(u,a) low','FontSize',12,'FontName','Times New Roman')
ylabel('pdf','FontSize',12,'FontName','Times New Roman')
axis([-1 1 0 2])

figure;
hold on;
grid on;
box on;
[n,x]=nhist(cos_a_u(1:count_enstro,3),50);
plot(x,n,'k');
[n,x]=nhist(cos_al_u(1:count_enstro,3),50);
plot(x,n,'b');
[n,x]=nhist(cos_ac_u(1:count_enstro,3),50);
plot(x,n,'c');
xlabel('cos(u,a) enstro','FontSize',12,'FontName','Times New Roman')
ylabel('pdf','FontSize',12,'FontName','Times New Roman')
axis([-1 1 0 2])

figure;
hold on;
grid on;
box on;
[n,x]=nhist(cos_a_u(1:count_strain,4),50);
plot(x,n,'k');
[n,x]=nhist(cos_al_u(1:count_strain,4),50);
plot(x,n,'b');
[n,x]=nhist(cos_ac_u(1:count_strain,4),50);
plot(x,n,'c');
xlabel('cos(u,a) strain','FontSize',12,'FontName','Times New Roman')
ylabel('pdf','FontSize',12,'FontName','Times New Roman')
axis([-1 1 0 2])

figure;
subplot(1,3,1);
hold on;
grid on;
box on;
[n,x]=nhist(cos_a_l1(1:counter,1),50);
plot(x,n,'k');
[n,x]=nhist(cos_ao_l1(1:counter,1),50);
plot(x,n,'r');
[n,x]=nhist(cos_ap_l1(1:counter,1),50);
plot(x,n,'g');
[n,x]=nhist(cos_al_l1(1:counter,1),50);
plot(x,n,'b');
[n,x]=nhist(cos_ac_l1(1:counter,1),50);
plot(x,n,'c');
xlabel('cos(a,\lambda_1) all','FontSize',12,'FontName','Times New Roman')
ylabel('pdf','FontSize',12,'FontName','Times New Roman')
axis([-1 1 0 2])

subplot(1,3,2);
hold on;
grid on;
box on;
[n,x]=nhist(cos_a_l2(1:counter,1),50);
plot(x,n,'k');
[n,x]=nhist(cos_ao_l2(1:counter,1),50);
plot(x,n,'r');
[n,x]=nhist(cos_ap_l2(1:counter,1),50);
plot(x,n,'g');
[n,x]=nhist(cos_al_l2(1:counter,1),50);
plot(x,n,'b');
[n,x]=nhist(cos_ac_l2(1:counter,1),50);
plot(x,n,'c');
xlabel('cos(a,\lambda_2) all','FontSize',12,'FontName','Times New Roman')
ylabel('pdf','FontSize',12,'FontName','Times New Roman')
axis([-1 1 0 2])

subplot(1,3,3);
hold on;
grid on;
box on;
[n,x]=nhist(cos_a_l3(1:counter,1),50);
plot(x,n,'k');
[n,x]=nhist(cos_ao_l3(1:counter,1),50);
plot(x,n,'r');
[n,x]=nhist(cos_ap_l3(1:counter,1),50);
plot(x,n,'g');
[n,x]=nhist(cos_al_l3(1:counter,1),50);
plot(x,n,'b');
[n,x]=nhist(cos_ac_l3(1:counter,1),50);
plot(x,n,'c');
xlabel('cos(a,\lambda_3) all','FontSize',12,'FontName','Times New Roman')
ylabel('pdf','FontSize',12,'FontName','Times New Roman')
axis([-1 1 0 2])

figure;
subplot(1,3,1);
hold on;
grid on;
box on;
[n,x]=nhist(cos_a_l1(1:count_low,2),50);
plot(x,n,'k');
[n,x]=nhist(cos_ao_l1(1:count_low,2),50);
plot(x,n,'r');
[n,x]=nhist(cos_ap_l1(1:count_low,2),50);
plot(x,n,'g');
[n,x]=nhist(cos_al_l1(1:count_low,2),50);
plot(x,n,'b');
[n,x]=nhist(cos_ac_l1(1:count_low,2),50);
plot(x,n,'c');
xlabel('cos(a,\lambda_1) low','FontSize',12,'FontName','Times New Roman')
ylabel('pdf','FontSize',12,'FontName','Times New Roman')
axis([-1 1 0 2])

subplot(1,3,2);
hold on;
grid on;
box on;
[n,x]=nhist(cos_a_l2(1:count_low,2),50);
plot(x,n,'k');
[n,x]=nhist(cos_ao_l2(1:count_low,2),50);
plot(x,n,'r');
[n,x]=nhist(cos_ap_l2(1:count_low,2),50);
plot(x,n,'g');
[n,x]=nhist(cos_al_l2(1:count_low,2),50);
plot(x,n,'b');
[n,x]=nhist(cos_ac_l2(1:count_low,2),50);
plot(x,n,'c');
xlabel('cos(a,\lambda_2) low','FontSize',12,'FontName','Times New Roman')
ylabel('pdf','FontSize',12,'FontName','Times New Roman')
axis([-1 1 0 2])

subplot(1,3,3);
hold on;
grid on;
box on;
[n,x]=nhist(cos_a_l3(1:count_low,2),50);
plot(x,n,'k');
[n,x]=nhist(cos_ao_l3(1:count_low,2),50);
plot(x,n,'r');
[n,x]=nhist(cos_ap_l3(1:count_low,2),50);
plot(x,n,'g');
[n,x]=nhist(cos_al_l3(1:count_low,2),50);
plot(x,n,'b');
[n,x]=nhist(cos_ac_l3(1:count_low,2),50);
plot(x,n,'c');
xlabel('cos(a,\lambda_3) low','FontSize',12,'FontName','Times New Roman')
ylabel('pdf','FontSize',12,'FontName','Times New Roman')
axis([-1 1 0 2])

figure;
subplot(1,3,1);
hold on;
grid on;
box on;
[n,x]=nhist(cos_a_l1(1:count_enstro,3),50);
plot(x,n,'k');
[n,x]=nhist(cos_ao_l1(1:count_enstro,3),50);
plot(x,n,'r');
[n,x]=nhist(cos_ap_l1(1:count_enstro,3),50);
plot(x,n,'g');
[n,x]=nhist(cos_al_l1(1:count_enstro,3),50);
plot(x,n,'b');
[n,x]=nhist(cos_ac_l1(1:count_enstro,3),50);
plot(x,n,'c');
xlabel('cos(a,\lambda_1) enstro','FontSize',12,'FontName','Times New Roman')
ylabel('pdf','FontSize',12,'FontName','Times New Roman')
axis([-1 1 0 2])

subplot(1,3,2);
hold on;
grid on;
box on;
[n,x]=nhist(cos_a_l2(1:count_enstro,3),50);
plot(x,n,'k');
[n,x]=nhist(cos_ao_l2(1:count_enstro,3),50);
plot(x,n,'r');
[n,x]=nhist(cos_ap_l2(1:count_enstro,3),50);
plot(x,n,'g');
[n,x]=nhist(cos_al_l2(1:count_enstro,3),50);
plot(x,n,'b');
[n,x]=nhist(cos_ac_l2(1:count_enstro,3),50);
plot(x,n,'c');
xlabel('cos(a,\lambda_2) enstro','FontSize',12,'FontName','Times New Roman')
ylabel('pdf','FontSize',12,'FontName','Times New Roman')
axis([-1 1 0 2])

subplot(1,3,3);
hold on;
grid on;
box on;
[n,x]=nhist(cos_a_l3(1:count_enstro,3),50);
plot(x,n,'k');
[n,x]=nhist(cos_ao_l3(1:count_enstro,3),50);
plot(x,n,'r');
[n,x]=nhist(cos_ap_l3(1:count_enstro,3),50);
plot(x,n,'g');
[n,x]=nhist(cos_al_l3(1:count_enstro,3),50);
plot(x,n,'b');
[n,x]=nhist(cos_ac_l3(1:count_enstro,3),50);
plot(x,n,'c');
xlabel('cos(a,\lambda_3) enstro','FontSize',12,'FontName','Times New Roman')
ylabel('pdf','FontSize',12,'FontName','Times New Roman')
axis([-1 1 0 2])

figure;
subplot(1,3,1);
hold on;
grid on;
box on;
[n,x]=nhist(cos_a_l1(1:count_strain,4),50);
plot(x,n,'k');
[n,x]=nhist(cos_ao_l1(1:count_strain,4),50);
plot(x,n,'r');
[n,x]=nhist(cos_ap_l1(1:count_strain,4),50);
plot(x,n,'g');
[n,x]=nhist(cos_al_l1(1:count_strain,4),50);
plot(x,n,'b');
[n,x]=nhist(cos_ac_l1(1:count_strain,4),50);
plot(x,n,'c');
xlabel('cos(a,\lambda_1) strain','FontSize',12,'FontName','Times New Roman')
ylabel('pdf','FontSize',12,'FontName','Times New Roman')
axis([-1 1 0 2])

subplot(1,3,2);
hold on;
grid on;
box on;
[n,x]=nhist(cos_a_l2(1:count_strain,4),50);
plot(x,n,'k');
[n,x]=nhist(cos_ao_l2(1:count_strain,4),50);
plot(x,n,'r');
[n,x]=nhist(cos_ap_l2(1:count_strain,4),50);
plot(x,n,'g');
[n,x]=nhist(cos_al_l2(1:count_strain,4),50);
plot(x,n,'b');
[n,x]=nhist(cos_ac_l2(1:count_strain,4),50);
plot(x,n,'c');
xlabel('cos(a,\lambda_2) strain','FontSize',12,'FontName','Times New Roman')
ylabel('pdf','FontSize',12,'FontName','Times New Roman')
axis([-1 1 0 2])

subplot(1,3,3);
hold on;
grid on;
box on;
[n,x]=nhist(cos_a_l3(1:count_strain,4),50);
plot(x,n,'k');
[n,x]=nhist(cos_ao_l3(1:count_strain,4),50);
plot(x,n,'r');
[n,x]=nhist(cos_ap_l3(1:count_strain,4),50);
plot(x,n,'g');
[n,x]=nhist(cos_al_l3(1:count_strain,4),50);
plot(x,n,'b');
[n,x]=nhist(cos_ac_l3(1:count_strain,4),50);
plot(x,n,'c');
xlabel('cos(a,\lambda_3) strain','FontSize',12,'FontName','Times New Roman')
ylabel('pdf','FontSize',12,'FontName','Times New Roman')
axis([-1 1 0 2])


figure;
hold on;
grid on;
box on;
[n,x]=nhist(cos_alac(1:counter,1),50);
plot(x,n,'k');
[n,x]=nhist(cos_alac(1:count_low,2),50);
plot(x,n,'r');
[n,x]=nhist(cos_alac(1:count_enstro,3),50);
plot(x,n,'g');
[n,x]=nhist(cos_alac(1:count_strain,4),50);
plot(x,n,'b');
xlabel('cos(a_l,a_c) ','FontSize',12,'FontName','Times New Roman')
ylabel('pdf','FontSize',12,'FontName','Times New Roman')
axis([-1 1 0 2])

figure;
hold on;
grid on;
box on;
[n,x]=nhist(cos_apac(1:counter,1),50);
plot(x,n,'k');
[n,x]=nhist(cos_apac(1:count_low,2),50);
plot(x,n,'r');
[n,x]=nhist(cos_apac(1:count_enstro,3),50);
plot(x,n,'g');
[n,x]=nhist(cos_apac(1:count_strain,4),50);
plot(x,n,'b');
xlabel('cos(a_{||},a_c) ','FontSize',12,'FontName','Times New Roman')
ylabel('pdf','FontSize',12,'FontName','Times New Roman')
axis([-1 1 0 2])


figure;
hold on;
grid on;
box on;
[n,x]=nhist(ratio_apac(1:counter,1),50);
plot(x,n,'k');
[n,x]=nhist(ratio_apac(1:count_low,2),50);
plot(x,n,'r');
[n,x]=nhist(ratio_apac(1:count_enstro,3),50);
plot(x,n,'g');
[n,x]=nhist(ratio_apac(1:count_strain,4),50);
plot(x,n,'b');
xlabel('|a_{||}|/|a_c|','FontSize',12,'FontName','Times New Roman')
ylabel('pdf','FontSize',12,'FontName','Times New Roman')
set(gca,'XScale','log')
set(gca,'YScale','log')