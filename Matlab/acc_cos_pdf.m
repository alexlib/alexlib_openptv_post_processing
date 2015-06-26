function res=acc_cos_pdf(first,last,tail,stepSize,threshold,Galilean_vel);
%acc_cos_pdf(10000,19900,0,20,1.37,0);0.015*10^6)

%Hi Alex, I outcomment everything that is not needed for the figures
%threshold is chosen in my case that 10^1.37`=23.4~<w^2>=<2s^2>

av_curv=67.6413;
av_curv_ac=133.4948;
av_curv_high=2.0039e-011;
av_curv_ac_high=5.7013e-005;

counter=0;
count_low=0;
count_enstro=0;
count_strain=0;

hist_absa=zeros(650000,4);
hist_absao=zeros(650000,4);
hist_absap=zeros(650000,4);
hist_absal=zeros(650000,4);
hist_absac=zeros(650000,4);
% hist_absus=zeros(650000,4);
% hist_absuo=zeros(650000,4);
hist_curv=zeros(650000,4);
hist_curvac=zeros(650000,4);

% cos_alac=zeros(650000,4);
% cos_usac=zeros(650000,4);
% cos_uoac=zeros(650000,4);
% cos_usuo=zeros(650000,4);
% cos_apac=zeros(650000,4);
% cos_apus=zeros(650000,4);
% cos_aouo=zeros(650000,4);

cos_a_w=zeros(650000,4);
cos_ao_w=zeros(650000,4);
cos_ap_w=zeros(650000,4);
cos_ac_w=zeros(650000,4);
% cos_al_w=zeros(650000,4);
% cos_us_w=zeros(650000,4);

cos_a_u=zeros(650000,4);
cos_ac_u=zeros(650000,4);

cos_a_l1=zeros(650000,4);
% cos_ao_l1=zeros(650000,4);
% cos_ap_l1=zeros(650000,4);
cos_ac_l1=zeros(650000,4);
% cos_al_l1=zeros(650000,4);
% cos_us_l1=zeros(650000,4);
% cos_uo_l1=zeros(650000,4);

cos_a_l2=zeros(650000,4);
% cos_ao_l2=zeros(650000,4);
% cos_ap_l2=zeros(650000,4);
cos_ac_l2=zeros(650000,4);
% cos_al_l2=zeros(650000,4);
% cos_us_l2=zeros(650000,4);
% cos_uo_l2=zeros(650000,4);

cos_a_l3=zeros(650000,4);
% cos_ao_l3=zeros(650000,4);
% cos_ap_l3=zeros(650000,4);
cos_ac_l3=zeros(650000,4);
% cos_al_l3=zeros(650000,4);
% cos_us_l3=zeros(650000,4);
% cos_uo_l3=zeros(650000,4);

% cos_a_u=zeros(650000,4);
% cos_ao_u=zeros(650000,4);
% cos_ap_u=zeros(650000,4);
% cos_ac_u=zeros(650000,4);
% cos_al_u=zeros(650000,4);
% cos_us_u=zeros(650000,4);

mean_s=3.4;
allpoints=0;


for i=first:stepSize:last
    clear f;
    i
    name=['D:/data/Riso_micro_PTV/trajPoint.',num2Str(i)];
    f=load(name);
    s=size(f);
    if s(1,1)>0
        xr=f(:,1);
        yr=f(:,2);
        zr=f(:,3);
        u=f(:,4)+Galilean_vel;
        v=f(:,5)+Galilean_vel;
        w=f(:,6)+Galilean_vel;
        absu=(u.^2+v.^2+w.^2).^0.5;
        ax=f(:,7);%%%%%%%%%%%
        ay=f(:,8);%%%%%%%%%%%
        az=f(:,9);%%%%%%%%%%%
        absa=(ax.^2+ay.^2+az.^2).^0.5;%%%%%%%%%%%%%
        alx=f(:,19);%%%%%%%%%%%
        aly=f(:,20);%%%%%%%%%%%
        alz=f(:,21);%%%%%%%%%%%
        absal=(alx.^2+aly.^2+alz.^2).^0.5;%%%%%%%%%%%%%%%
        w1=f(:,10);
        w2=f(:,11);
        w3=f(:,12);
        absw=(w1.^2+w2.^2+w3.^2).^0.5;
        enstro=w1.^2+w2.^2+w3.^2;
        s11=f(:,13);
        s12=f(:,14);
        s13=f(:,15);
        s22=f(:,16);
        s23=f(:,17);
        s33=f(:,18);
        usx=u.*s11+v.*s12+w.*s13;
        usy=u.*s12+v.*s22+w.*s23;
        usz=u.*s13+v.*s23+w.*s33;
        absus=(usx.^2+usy.^2+usz.^2).^0.5;
        ux=s11;
        uy=s12-0.5*w3;
        uz=s13+0.5*w2;
        vx=s12+0.5*w3;
        vy=s22;
        vz=s23-0.5*w1;
        wx=s13-0.5*w2;
        wy=s23+0.5*w1;
        wz=s33;
%         o11=ux-s11;
%         o12=uy-s12;
%         o13=uz-s13;
%         o21=vx-s12;
%         o22=vy-s22;
%         o23=vz-s23;
%         o31=wx-s13;
%         o32=wy-s23;
%         o33=wz-s33;
%         uox=u.*o11+v.*o12+w.*o13;
%         uoy=u.*o21+v.*o22+w.*o23;
%         uoz=u.*o31+v.*o32+w.*o33;
%         absuo=(uox.^2+uoy.^2+uoz.^2).^0.5;
        acx=u.*ux+v.*uy+w.*uz;%%%%%%%%%%%%%%%%%%%%
        acy=u.*vx+v.*vy+w.*vz;%%%%%%%%%%%%%%%%%%%%
        acz=u.*wx+v.*wy+w.*wz;%%%%%%%%%%%%%%%%%%%%
        absac=(acx.^2+acy.^2+acz.^2).^0.5;%%%%%%%%%%%%%%%%%%%%%%%%%
%         aLamx=w2.*w-w3.*v;%%%%%%%%%%%%%%%%%%%%%%%
%         aLamy=w3.*u-w1.*w;%%%%%%%%%%%%%%%%%%%%%%%
%         aLamz=w1.*v-w2.*u;%%%%%%%%%%%%%%%%%%%%%%%%%
%         absaLam=(aLamx.^2+aLamy.^2+aLamz.^2).^0.5;%%%%%%%%%%%%%%%%%
%         aBx=acx-aLamx;%%%%%%%%%%%%%%%%%%%%%
%         aBy=acy-aLamy;%%%%%%%%%%%%%%%%%%%%%%%
%         aBz=acz-aLamz;%%%%%%%%%%%%%%%%%%%%%%%%%
%         absaB=(aBx.^2+aBy.^2+aBz.^2).^0.5;%%%%%%%%%%%%%%%%%%%
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
        cux=(u.*ux+v.*uy+w.*uz)./(absu.^2);
        cuy=(u.*vx+v.*vy+w.*vz)./(absu.^2);
        cuz=(u.*wx+v.*wy+w.*wz)./(absu.^2);
        curvGrad=(cux.^2+cuy.^2+cuz.^2).^0.5;
        curv=(((v.*az-w.*ay).^2+(w.*ax-u.*az).^2+(u.*ay-v.*ax).^2).^0.5)./(absu.^3);
        curvac=(((v.*acz-w.*acy).^2+(w.*acx-u.*acz).^2+(u.*acy-v.*acx).^2).^0.5)./(absu.^3);
        
%         ua=u.*acx+v.*acy+w.*acz;
        
%         daxdx=f(:,22);
%         daxdy=f(:,23);
%         daxdz=f(:,24);
%         daydx=f(:,25);
%         daydy=f(:,26);
%         daydz=f(:,27);
%         dazdx=f(:,28);
%         dazdy=f(:,29);
%         dazdz=f(:,30);
        
        reldiv=f(:,31);
        age=f(:,32);

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
                    if age(ii)>tail & age(ii)<leng-tail+1 & reldiv(ii)<0.2 & absa(ii)^2/8e-4<100 & curv(ii)<3000 %& absu(ii)>0.015*1.5%%
                        uM=[ux(ii) uy(ii) uz(ii);vx(ii) vy(ii) vz(ii);wx(ii) wy(ii) wz(ii)];
                        A=[s11(ii) s12(ii) s13(ii);s12(ii) s22(ii) s23(ii);s13(ii) s23(ii) s33(ii)];
                        [V,D] = eig(A);
                        [D,k] = sort(diag(D)); 	% ascending order
                        D = diag(D(end:-1:1)); 			% descending order
                        V = V(:,k(end:-1:1)); 	% the same order for eigenvectors
                        strain=D(1,1)^2+D(2,2)^2+D(3,3)^2;
                        
                        ua  = [u(ii) v(ii) w(ii)]*[ax(ii) ay(ii) az(ii)]';
                        
                        if counter<650000+1 %& ua<-0.5e-3
                            counter=counter+1;
                            hist_absa(counter,1)=absa(ii);
                            hist_absao(counter,1)=absaO(ii);
                            hist_absap(counter,1)=absaP(ii);
                            hist_absal(counter,1)=absal(ii);
                            hist_absac(counter,1)=absac(ii);
                            hist_curv(counter,1)=curv(ii);
                            hist_curvac(counter,1)=curvac(ii);
%                             hist_absus(counter,1)=absus(ii);
%                             hist_absuo(counter,1)=absuo(ii);
                            
                            cos_a_w(counter,1)=[ax(ii) ay(ii) az(ii)]*[w1(ii) w2(ii) w3(ii)]'/(absa(ii)*absw(ii));
                            cos_ao_w(counter,1)=[aOx(ii) aOy(ii) aOz(ii)]*[w1(ii) w2(ii) w3(ii)]'/(absaO(ii)*absw(ii));
                            cos_ap_w(counter,1)=[aPx(ii) aPy(ii) aPz(ii)]*[w1(ii) w2(ii) w3(ii)]'/(absaP(ii)*absw(ii));
                            cos_ac_w(counter,1)=[acx(ii) acy(ii) acz(ii)]*[w1(ii) w2(ii) w3(ii)]'/(absac(ii)*absw(ii));
%                             cos_al_w(counter,1)=[alx(ii) aly(ii) alz(ii)]*[w1(ii) w2(ii) w3(ii)]'/(absal(ii)*absw(ii));
%                             cos_us_w(counter,1)=[usx(ii) usy(ii) usz(ii)]*[w1(ii) w2(ii) w3(ii)]'/(absus(ii)*absw(ii));
%                             cos_uo_w(counter,1)=[uox(ii) uoy(ii) uoz(ii)]*[w1(ii) w2(ii) w3(ii)]'/(absuo(ii)*absw(ii));
                            cos_a_u(counter,1)=[ax(ii) ay(ii) az(ii)]*[u(ii) v(ii) w(ii)]'/(absa(ii)*absu(ii));
                            cos_ac_u(counter,1)=[acx(ii) acy(ii) acz(ii)]*[u(ii) v(ii) w(ii)]'/(absac(ii)*absu(ii));
                            
                            cos_a_l1(counter,1)=[ax(ii) ay(ii) az(ii)]*V(:,1)/(absa(ii)*(sum(V(:,1).^2))^0.5);
%                             cos_ao_l1(counter,1)=[aOx(ii) aOy(ii) aOz(ii)]*V(:,1)/(absaO(ii)*(sum(V(:,1).^2))^0.5);
%                             cos_ap_l1(counter,1)=[aPx(ii) aPy(ii) aPz(ii)]*V(:,1)/(absaP(ii)*(sum(V(:,1).^2))^0.5);
                            cos_ac_l1(counter,1)=[acx(ii) acy(ii) acz(ii)]*V(:,1)/(absac(ii)*(sum(V(:,1).^2))^0.5);
%                             cos_al_l1(counter,1)=[alx(ii) aly(ii) alz(ii)]*V(:,1)/(absal(ii)*(sum(V(:,1).^2))^0.5);
%                             cos_us_l1(counter,1)=[usx(ii) usy(ii) usz(ii)]*V(:,1)/(absus(ii)*(sum(V(:,1).^2))^0.5);
%                             cos_uo_l1(counter,1)=[uox(ii) uoy(ii) uoz(ii)]*V(:,1)/(absuo(ii)*(sum(V(:,1).^2))^0.5);
                            
                            cos_a_l2(counter,1)=[ax(ii) ay(ii) az(ii)]*V(:,2)/(absa(ii)*(sum(V(:,2).^2))^0.5);
%                             cos_ao_l2(counter,1)=[aOx(ii) aOy(ii) aOz(ii)]*V(:,2)/(absaO(ii)*(sum(V(:,2).^2))^0.5);
%                             cos_ap_l2(counter,1)=[aPx(ii) aPy(ii) aPz(ii)]*V(:,2)/(absaP(ii)*(sum(V(:,2).^2))^0.5);
                            cos_ac_l2(counter,1)=[acx(ii) acy(ii) acz(ii)]*V(:,2)/(absac(ii)*(sum(V(:,2).^2))^0.5);
%                             cos_al_l2(counter,1)=[alx(ii) aly(ii) alz(ii)]*V(:,2)/(absal(ii)*(sum(V(:,2).^2))^0.5);
%                             cos_us_l2(counter,1)=[usx(ii) usy(ii) usz(ii)]*V(:,2)/(absus(ii)*(sum(V(:,2).^2))^0.5);
%                             cos_uo_l2(counter,1)=[uox(ii) uoy(ii) uoz(ii)]*V(:,2)/(absuo(ii)*(sum(V(:,2).^2))^0.5);
                            
                            cos_a_l3(counter,1)=[ax(ii) ay(ii) az(ii)]*V(:,3)/(absa(ii)*(sum(V(:,3).^2))^0.5);
%                             cos_ao_l3(counter,1)=[aOx(ii) aOy(ii) aOz(ii)]*V(:,3)/(absaO(ii)*(sum(V(:,3).^2))^0.5);
%                             cos_ap_l3(counter,1)=[aPx(ii) aPy(ii) aPz(ii)]*V(:,3)/(absaP(ii)*(sum(V(:,3).^2))^0.5);
                            cos_ac_l3(counter,1)=[acx(ii) acy(ii) acz(ii)]*V(:,3)/(absac(ii)*(sum(V(:,3).^2))^0.5);
%                             cos_al_l3(counter,1)=[alx(ii) aly(ii) alz(ii)]*V(:,3)/(absal(ii)*(sum(V(:,3).^2))^0.5);
%                             cos_us_l3(counter,1)=[usx(ii) usy(ii) usz(ii)]*V(:,3)/(absus(ii)*(sum(V(:,3).^2))^0.5);
%                             cos_uo_l3(counter,1)=[uox(ii) uoy(ii) uoz(ii)]*V(:,3)/(absuo(ii)*(sum(V(:,3).^2))^0.5);
                            
%                             cos_a_u(counter,1)=[ax(ii) ay(ii) az(ii)]*[u(ii) v(ii) w(ii)]'/(absa(ii)*absu(ii));
%                             cos_ao_u(counter,1)=[aOx(ii) aOy(ii) aOz(ii)]*[u(ii) v(ii) w(ii)]'/(absaO(ii)*absu(ii));
%                             cos_ap_u(counter,1)=[aPx(ii) aPy(ii) aPz(ii)]*[u(ii) v(ii) w(ii)]'/(absaP(ii)*absu(ii));
%                             cos_ac_u(counter,1)=[acx(ii) acy(ii) acz(ii)]*[u(ii) v(ii) w(ii)]'/(absac(ii)*absu(ii));
%                             cos_al_u(counter,1)=[alx(ii) aly(ii) alz(ii)]*[u(ii) v(ii) w(ii)]'/(absal(ii)*absu(ii));
%                             cos_us_u(counter,1)=[usx(ii) usy(ii) usz(ii)]*[u(ii) v(ii) w(ii)]'/(absus(ii)*absu(ii));
%                             cos_uo_u(counter,1)=[uox(ii) uoy(ii) uoz(ii)]*[u(ii) v(ii) w(ii)]'/(absuo(ii)*absu(ii));
%                             
%                             cos_alac(counter,1)=[alx(ii) aly(ii) alz(ii)]*[acx(ii) acy(ii) acz(ii)]'/(absal(ii)*absac(ii));
%                             cos_usac(counter,1)=[usx(ii) usy(ii) usz(ii)]*[acx(ii) acy(ii) acz(ii)]'/(absus(ii)*absac(ii));
%                             cos_uoac(counter,1)=[uox(ii) uoy(ii) uoz(ii)]*[acx(ii) acy(ii) acz(ii)]'/(absuo(ii)*absac(ii));
%                             cos_usuo(counter,1)=[usx(ii) usy(ii) usz(ii)]*[uox(ii) uoy(ii) uoz(ii)]'/(absus(ii)*absuo(ii));
%                             cos_apac(counter,1)=[aPx(ii) aPy(ii) aPz(ii)]*[acx(ii) acy(ii) acz(ii)]'/(absaP(ii)*absac(ii));
%                             cos_apus(counter,1)=[aPx(ii) aPy(ii) aPz(ii)]*[usx(ii) usy(ii) usz(ii)]'/(absaP(ii)*absus(ii));
%                             cos_aouo(counter,1)=[aOx(ii) aOy(ii) aOz(ii)]*[uox(ii) uoy(ii) uoz(ii)]'/(absaO(ii)*absuo(ii));
                        end
                        
                        if enstro(ii)<10^threshold & 2*strain<10^threshold %& ua(ii)<0
                            if count_low<650000+1
                                count_low=count_low+1;
                                hist_absa(count_low,2)=absa(ii);
%                                 hist_absao(count_low,2)=absaO(ii);
%                                 hist_absap(count_low,2)=absaP(ii);
%                                 hist_absal(count_low,2)=absal(ii);
%                                 hist_absac(count_low,2)=absac(ii);
                                hist_curv(count_low,2)=curv(ii);
                                hist_curvac(count_low,2)=curvac(ii);
%                                 hist_absus(count_low,2)=absus(ii);
%                                 hist_absuo(count_low,2)=absuo(ii);

%                                 cos_a_w(count_low,2)=[ax(ii) ay(ii) az(ii)]*[w1(ii) w2(ii) w3(ii)]'/(absa(ii)*absw(ii));
%                                 cos_ao_w(count_low,2)=[aOx(ii) aOy(ii) aOz(ii)]*[w1(ii) w2(ii) w3(ii)]'/(absaO(ii)*absw(ii));
%                                 cos_ap_w(count_low,2)=[aPx(ii) aPy(ii) aPz(ii)]*[w1(ii) w2(ii) w3(ii)]'/(absaP(ii)*absw(ii));
%                                 cos_ac_w(count_low,2)=[acx(ii) acy(ii) acz(ii)]*[w1(ii) w2(ii) w3(ii)]'/(absac(ii)*absw(ii));
%                                 cos_al_w(count_low,2)=[alx(ii) aly(ii) alz(ii)]*[w1(ii) w2(ii) w3(ii)]'/(absal(ii)*absw(ii));
%                                 cos_us_w(count_low,2)=[usx(ii) usy(ii) usz(ii)]*[w1(ii) w2(ii) w3(ii)]'/(absus(ii)*absw(ii));
%                                 cos_uo_w(count_low,2)=[uox(ii) uoy(ii) uoz(ii)]*[w1(ii) w2(ii) w3(ii)]'/(absuo(ii)*absw(ii));
%                                 
%                                 cos_a_l1(count_low,2)=[ax(ii) ay(ii) az(ii)]*V(:,1)/(absa(ii)*(sum(V(:,1).^2))^0.5);
%                                 cos_ao_l1(count_low,2)=[aOx(ii) aOy(ii) aOz(ii)]*V(:,1)/(absaO(ii)*(sum(V(:,1).^2))^0.5);
%                                 cos_ap_l1(count_low,2)=[aPx(ii) aPy(ii) aPz(ii)]*V(:,1)/(absaP(ii)*(sum(V(:,1).^2))^0.5);
%                                 cos_ac_l1(count_low,2)=[acx(ii) acy(ii) acz(ii)]*V(:,1)/(absac(ii)*(sum(V(:,1).^2))^0.5);
%                                 cos_al_l1(count_low,2)=[alx(ii) aly(ii) alz(ii)]*V(:,1)/(absal(ii)*(sum(V(:,1).^2))^0.5);
%                                 cos_us_l1(count_low,2)=[usx(ii) usy(ii) usz(ii)]*V(:,1)/(absus(ii)*(sum(V(:,1).^2))^0.5);
%                                 cos_uo_l1(count_low,2)=[uox(ii) uoy(ii) uoz(ii)]*V(:,1)/(absuo(ii)*(sum(V(:,1).^2))^0.5);
% 
%                                 cos_a_l2(count_low,2)=[ax(ii) ay(ii) az(ii)]*V(:,2)/(absa(ii)*(sum(V(:,2).^2))^0.5);
%                                 cos_ao_l2(count_low,2)=[aOx(ii) aOy(ii) aOz(ii)]*V(:,2)/(absaO(ii)*(sum(V(:,2).^2))^0.5);
%                                 cos_ap_l2(count_low,2)=[aPx(ii) aPy(ii) aPz(ii)]*V(:,2)/(absaP(ii)*(sum(V(:,2).^2))^0.5);
%                                 cos_ac_l2(count_low,2)=[acx(ii) acy(ii) acz(ii)]*V(:,2)/(absac(ii)*(sum(V(:,2).^2))^0.5);
%                                 cos_al_l2(count_low,2)=[alx(ii) aly(ii) alz(ii)]*V(:,2)/(absal(ii)*(sum(V(:,2).^2))^0.5);
%                                 cos_us_l2(count_low,2)=[usx(ii) usy(ii) usz(ii)]*V(:,2)/(absus(ii)*(sum(V(:,2).^2))^0.5);
%                                 cos_uo_l2(count_low,2)=[uox(ii) uoy(ii) uoz(ii)]*V(:,2)/(absuo(ii)*(sum(V(:,2).^2))^0.5);
% 
%                                 cos_a_l3(count_low,2)=[ax(ii) ay(ii) az(ii)]*V(:,3)/(absa(ii)*(sum(V(:,3).^2))^0.5);
%                                 cos_ao_l3(count_low,2)=[aOx(ii) aOy(ii) aOz(ii)]*V(:,3)/(absaO(ii)*(sum(V(:,3).^2))^0.5);
%                                 cos_ap_l3(count_low,2)=[aPx(ii) aPy(ii) aPz(ii)]*V(:,3)/(absaP(ii)*(sum(V(:,3).^2))^0.5);
%                                 cos_ac_l3(count_low,2)=[acx(ii) acy(ii) acz(ii)]*V(:,3)/(absac(ii)*(sum(V(:,3).^2))^0.5);
%                                 cos_al_l3(count_low,2)=[alx(ii) aly(ii) alz(ii)]*V(:,3)/(absal(ii)*(sum(V(:,3).^2))^0.5);
%                                 cos_us_l3(count_low,2)=[usx(ii) usy(ii) usz(ii)]*V(:,3)/(absus(ii)*(sum(V(:,3).^2))^0.5);
%                                 cos_uo_l3(count_low,2)=[uox(ii) uoy(ii) uoz(ii)]*V(:,3)/(absuo(ii)*(sum(V(:,3).^2))^0.5);
%                                 
%                                 cos_a_u(count_low,2)=[ax(ii) ay(ii) az(ii)]*[u(ii) v(ii) w(ii)]'/(absa(ii)*absu(ii));
%                                 cos_ao_u(count_low,2)=[aOx(ii) aOy(ii) aOz(ii)]*[u(ii) v(ii) w(ii)]'/(absaO(ii)*absu(ii));
%                                 cos_ap_u(count_low,2)=[aPx(ii) aPy(ii) aPz(ii)]*[u(ii) v(ii) w(ii)]'/(absaP(ii)*absu(ii));
%                                 cos_ac_u(count_low,2)=[acx(ii) acy(ii) acz(ii)]*[u(ii) v(ii) w(ii)]'/(absac(ii)*absu(ii));
%                                 cos_al_u(count_low,2)=[alx(ii) aly(ii) alz(ii)]*[u(ii) v(ii) w(ii)]'/(absal(ii)*absu(ii));
%                                 cos_us_u(count_low,2)=[usx(ii) usy(ii) usz(ii)]*[u(ii) v(ii) w(ii)]'/(absus(ii)*absu(ii));
%                                 cos_uo_u(count_low,2)=[uox(ii) uoy(ii) uoz(ii)]*[u(ii) v(ii) w(ii)]'/(absuo(ii)*absu(ii));
%                                 
%                                 cos_alac(count_low,2)=[alx(ii) aly(ii) alz(ii)]*[acx(ii) acy(ii) acz(ii)]'/(absal(ii)*absac(ii));
%                                 cos_usac(count_low,2)=[usx(ii) usy(ii) usz(ii)]*[acx(ii) acy(ii) acz(ii)]'/(absus(ii)*absac(ii));
%                                 cos_uoac(count_low,2)=[uox(ii) uoy(ii) uoz(ii)]*[acx(ii) acy(ii) acz(ii)]'/(absuo(ii)*absac(ii));
%                                 cos_usuo(count_low,2)=[usx(ii) usy(ii) usz(ii)]*[uox(ii) uoy(ii) uoz(ii)]'/(absus(ii)*absuo(ii));
%                                 cos_apac(count_low,2)=[aPx(ii) aPy(ii) aPz(ii)]*[acx(ii) acy(ii) acz(ii)]'/(absaP(ii)*absac(ii));
%                                 cos_apus(count_low,2)=[aPx(ii) aPy(ii) aPz(ii)]*[usx(ii) usy(ii) usz(ii)]'/(absaP(ii)*absus(ii));
%                                 cos_aouo(count_low,2)=[aOx(ii) aOy(ii) aOz(ii)]*[uox(ii) uoy(ii) uoz(ii)]'/(absaO(ii)*absuo(ii));
                            end
                        end
                        if 2*strain<enstro(ii) & enstro(ii)>10^threshold %& ua(ii)<0
                            if count_enstro<650000+1
                                count_enstro=count_enstro+1;
                                hist_absa(count_enstro,3)=absa(ii);
%                                 hist_absao(count_enstro,3)=absaO(ii);
%                                 hist_absap(count_enstro,3)=absaP(ii);
%                                 hist_absal(count_enstro,3)=absal(ii);
%                                 hist_absac(count_enstro,3)=absac(ii);
                                hist_curv(count_enstro,3)=curv(ii);
                                hist_curvac(count_enstro,3)=curvac(ii);
%                                 hist_absus(count_enstro,3)=absus(ii);
%                                 hist_absuo(count_enstro,3)=absuo(ii);
                                
                                cos_a_w(count_enstro,3)=[ax(ii) ay(ii) az(ii)]*[w1(ii) w2(ii) w3(ii)]'/(absa(ii)*absw(ii));
                                cos_ao_w(count_enstro,3)=[aOx(ii) aOy(ii) aOz(ii)]*[w1(ii) w2(ii) w3(ii)]'/(absaO(ii)*absw(ii));
                                cos_ap_w(count_enstro,3)=[aPx(ii) aPy(ii) aPz(ii)]*[w1(ii) w2(ii) w3(ii)]'/(absaP(ii)*absw(ii));
                                cos_ac_w(count_enstro,3)=[acx(ii) acy(ii) acz(ii)]*[w1(ii) w2(ii) w3(ii)]'/(absac(ii)*absw(ii));
%                                 cos_al_w(count_enstro,3)=[alx(ii) aly(ii) alz(ii)]*[w1(ii) w2(ii) w3(ii)]'/(absal(ii)*absw(ii));
%                                 cos_us_w(count_enstro,3)=[usx(ii) usy(ii) usz(ii)]*[w1(ii) w2(ii) w3(ii)]'/(absus(ii)*absw(ii));
%                                 cos_uo_w(count_enstro,3)=[uox(ii) uoy(ii) uoz(ii)]*[w1(ii) w2(ii) w3(ii)]'/(absuo(ii)*absw(ii));
                                cos_a_u(count_enstro,3)=[ax(ii) ay(ii) az(ii)]*[u(ii) v(ii) w(ii)]'/(absa(ii)*absu(ii));
                                cos_ac_u(count_enstro,3)=[acx(ii) acy(ii) acz(ii)]*[u(ii) v(ii) w(ii)]'/(absac(ii)*absu(ii));
                                
                                cos_a_l1(count_enstro,3)=[ax(ii) ay(ii) az(ii)]*V(:,1)/(absa(ii)*(sum(V(:,1).^2))^0.5);
%                                 cos_ao_l1(count_enstro,3)=[aOx(ii) aOy(ii) aOz(ii)]*V(:,1)/(absaO(ii)*(sum(V(:,1).^2))^0.5);
%                                 cos_ap_l1(count_enstro,3)=[aPx(ii) aPy(ii) aPz(ii)]*V(:,1)/(absaP(ii)*(sum(V(:,1).^2))^0.5);
                                cos_ac_l1(count_enstro,3)=[acx(ii) acy(ii) acz(ii)]*V(:,1)/(absac(ii)*(sum(V(:,1).^2))^0.5);
%                                 cos_al_l1(count_enstro,3)=[alx(ii) aly(ii) alz(ii)]*V(:,1)/(absal(ii)*(sum(V(:,1).^2))^0.5);
%                                 cos_us_l1(count_enstro,3)=[usx(ii) usy(ii) usz(ii)]*V(:,1)/(absus(ii)*(sum(V(:,1).^2))^0.5);
%                                 cos_uo_l1(count_enstro,3)=[uox(ii) uoy(ii) uoz(ii)]*V(:,1)/(absuo(ii)*(sum(V(:,1).^2))^0.5);

                                cos_a_l2(count_enstro,3)=[ax(ii) ay(ii) az(ii)]*V(:,2)/(absa(ii)*(sum(V(:,2).^2))^0.5);
%                                 cos_ao_l2(count_enstro,3)=[aOx(ii) aOy(ii) aOz(ii)]*V(:,2)/(absaO(ii)*(sum(V(:,2).^2))^0.5);
%                                 cos_ap_l2(count_enstro,3)=[aPx(ii) aPy(ii) aPz(ii)]*V(:,2)/(absaP(ii)*(sum(V(:,2).^2))^0.5);
                                cos_ac_l2(count_enstro,3)=[acx(ii) acy(ii) acz(ii)]*V(:,2)/(absac(ii)*(sum(V(:,2).^2))^0.5);
%                                 cos_al_l2(count_enstro,3)=[alx(ii) aly(ii) alz(ii)]*V(:,2)/(absal(ii)*(sum(V(:,2).^2))^0.5);
%                                 cos_us_l2(count_enstro,3)=[usx(ii) usy(ii) usz(ii)]*V(:,2)/(absus(ii)*(sum(V(:,2).^2))^0.5);
%                                 cos_uo_l2(count_enstro,3)=[uox(ii) uoy(ii) uoz(ii)]*V(:,2)/(absuo(ii)*(sum(V(:,2).^2))^0.5);

                                cos_a_l3(count_enstro,3)=[ax(ii) ay(ii) az(ii)]*V(:,3)/(absa(ii)*(sum(V(:,3).^2))^0.5);
%                                 cos_ao_l3(count_enstro,3)=[aOx(ii) aOy(ii) aOz(ii)]*V(:,3)/(absaO(ii)*(sum(V(:,3).^2))^0.5);
%                                 cos_ap_l3(count_enstro,3)=[aPx(ii) aPy(ii) aPz(ii)]*V(:,3)/(absaP(ii)*(sum(V(:,3).^2))^0.5);
                                cos_ac_l3(count_enstro,3)=[acx(ii) acy(ii) acz(ii)]*V(:,3)/(absac(ii)*(sum(V(:,3).^2))^0.5);
%                                 cos_al_l3(count_enstro,3)=[alx(ii) aly(ii) alz(ii)]*V(:,3)/(absal(ii)*(sum(V(:,3).^2))^0.5);
%                                 cos_us_l3(count_enstro,3)=[usx(ii) usy(ii) usz(ii)]*V(:,3)/(absus(ii)*(sum(V(:,3).^2))^0.5);
%                                 cos_uo_l3(count_enstro,3)=[uox(ii) uoy(ii) uoz(ii)]*V(:,3)/(absuo(ii)*(sum(V(:,3).^2))^0.5);
                                
%                                 cos_a_u(count_enstro,3)=[ax(ii) ay(ii) az(ii)]*[u(ii) v(ii) w(ii)]'/(absa(ii)*absu(ii));
%                                 cos_ao_u(count_enstro,3)=[aOx(ii) aOy(ii) aOz(ii)]*[u(ii) v(ii) w(ii)]'/(absaO(ii)*absu(ii));
%                                 cos_ap_u(count_enstro,3)=[aPx(ii) aPy(ii) aPz(ii)]*[u(ii) v(ii) w(ii)]'/(absaP(ii)*absu(ii));
%                                 cos_ac_u(count_enstro,3)=[acx(ii) acy(ii) acz(ii)]*[u(ii) v(ii) w(ii)]'/(absac(ii)*absu(ii));
%                                 cos_al_u(count_enstro,3)=[alx(ii) aly(ii) alz(ii)]*[u(ii) v(ii) w(ii)]'/(absal(ii)*absu(ii));
%                                 cos_us_u(count_enstro,3)=[usx(ii) usy(ii) usz(ii)]*[u(ii) v(ii) w(ii)]'/(absus(ii)*absu(ii));
%                                 cos_uo_u(count_enstro,3)=[uox(ii) uoy(ii) uoz(ii)]*[u(ii) v(ii) w(ii)]'/(absuo(ii)*absu(ii));
%                                 
%                                 cos_alac(count_enstro,3)=[alx(ii) aly(ii) alz(ii)]*[acx(ii) acy(ii) acz(ii)]'/(absal(ii)*absac(ii));
%                                 cos_usac(count_enstro,3)=[usx(ii) usy(ii) usz(ii)]*[acx(ii) acy(ii) acz(ii)]'/(absus(ii)*absac(ii));
%                                 cos_uoac(count_enstro,3)=[uox(ii) uoy(ii) uoz(ii)]*[acx(ii) acy(ii) acz(ii)]'/(absuo(ii)*absac(ii));
%                                 cos_usuo(count_enstro,3)=[usx(ii) usy(ii) usz(ii)]*[uox(ii) uoy(ii) uoz(ii)]'/(absus(ii)*absuo(ii));
%                                 cos_apac(count_enstro,3)=[aPx(ii) aPy(ii) aPz(ii)]*[acx(ii) acy(ii) acz(ii)]'/(absaP(ii)*absac(ii));
%                                 cos_apus(count_enstro,3)=[aPx(ii) aPy(ii) aPz(ii)]*[usx(ii) usy(ii) usz(ii)]'/(absaP(ii)*absus(ii));
%                                 cos_aouo(count_enstro,3)=[aOx(ii) aOy(ii) aOz(ii)]*[uox(ii) uoy(ii) uoz(ii)]'/(absaO(ii)*absuo(ii));
                            end
                        end
                        %if 2*strain>enstro(ii) & 2*strain>10^threshold %& ua(ii)<0
                        if 2*strain>10^threshold & enstro(ii)>10^threshold %& ua(ii)<0
                            if count_strain<650000+1
                                count_strain=count_strain+1;
                                hist_absa(count_strain,4)=absa(ii);
%                                 hist_absao(count_strain,4)=absaO(ii);
%                                 hist_absap(count_strain,4)=absaP(ii);
%                                 hist_absal(count_strain,4)=absal(ii);
%                                 hist_absac(count_strain,4)=absac(ii);
                                hist_curv(count_strain,4)=curv(ii);
                                hist_curvac(count_strain,4)=curvac(ii);
%                                 hist_absus(count_strain,4)=absus(ii);
%                                 hist_absuo(count_strain,4)=absuo(ii);
                                
                                cos_a_w(count_strain,4)=[ax(ii) ay(ii) az(ii)]*[w1(ii) w2(ii) w3(ii)]'/(absa(ii)*absw(ii));
                                cos_ao_w(count_strain,4)=[aOx(ii) aOy(ii) aOz(ii)]*[w1(ii) w2(ii) w3(ii)]'/(absaO(ii)*absw(ii));
                                cos_ap_w(count_strain,4)=[aPx(ii) aPy(ii) aPz(ii)]*[w1(ii) w2(ii) w3(ii)]'/(absaP(ii)*absw(ii));
                                cos_ac_w(count_strain,4)=[acx(ii) acy(ii) acz(ii)]*[w1(ii) w2(ii) w3(ii)]'/(absac(ii)*absw(ii));
%                                 cos_al_w(count_strain,4)=[alx(ii) aly(ii) alz(ii)]*[w1(ii) w2(ii) w3(ii)]'/(absal(ii)*absw(ii));
%                                 cos_us_w(count_strain,4)=[usx(ii) usy(ii) usz(ii)]*[w1(ii) w2(ii) w3(ii)]'/(absus(ii)*absw(ii));
%                                 cos_uo_w(count_strain,4)=[uox(ii) uoy(ii) uoz(ii)]*[w1(ii) w2(ii) w3(ii)]'/(absuo(ii)*absw(ii));
                                cos_a_u(count_strain,4)=[ax(ii) ay(ii) az(ii)]*[u(ii) v(ii) w(ii)]'/(absa(ii)*absu(ii));
                                cos_ac_u(count_strain,4)=[acx(ii) acy(ii) acz(ii)]*[u(ii) v(ii) w(ii)]'/(absac(ii)*absu(ii));                                

                                cos_a_l1(count_strain,4)=[ax(ii) ay(ii) az(ii)]*V(:,1)/(absa(ii)*(sum(V(:,1).^2))^0.5);
%                                 cos_ao_l1(count_strain,4)=[aOx(ii) aOy(ii) aOz(ii)]*V(:,1)/(absaO(ii)*(sum(V(:,1).^2))^0.5);
%                                 cos_ap_l1(count_strain,4)=[aPx(ii) aPy(ii) aPz(ii)]*V(:,1)/(absaP(ii)*(sum(V(:,1).^2))^0.5);
                                cos_ac_l1(count_strain,4)=[acx(ii) acy(ii) acz(ii)]*V(:,1)/(absac(ii)*(sum(V(:,1).^2))^0.5);
%                                 cos_al_l1(count_strain,4)=[alx(ii) aly(ii) alz(ii)]*V(:,1)/(absal(ii)*(sum(V(:,1).^2))^0.5);
%                                 cos_us_l1(count_strain,4)=[usx(ii) usy(ii) usz(ii)]*V(:,1)/(absus(ii)*(sum(V(:,1).^2))^0.5);
%                                 cos_uo_l1(count_strain,4)=[uox(ii) uoy(ii) uoz(ii)]*V(:,1)/(absuo(ii)*(sum(V(:,1).^2))^0.5);
% 
                                cos_a_l2(count_strain,4)=[ax(ii) ay(ii) az(ii)]*V(:,2)/(absa(ii)*(sum(V(:,2).^2))^0.5);
%                                 cos_ao_l2(count_strain,4)=[aOx(ii) aOy(ii) aOz(ii)]*V(:,2)/(absaO(ii)*(sum(V(:,2).^2))^0.5);
%                                 cos_ap_l2(count_strain,4)=[aPx(ii) aPy(ii) aPz(ii)]*V(:,2)/(absaP(ii)*(sum(V(:,2).^2))^0.5);
                                cos_ac_l2(count_strain,4)=[acx(ii) acy(ii) acz(ii)]*V(:,2)/(absac(ii)*(sum(V(:,2).^2))^0.5);
%                                 cos_al_l2(count_strain,4)=[alx(ii) aly(ii) alz(ii)]*V(:,2)/(absal(ii)*(sum(V(:,2).^2))^0.5);
%                                 cos_us_l2(count_strain,4)=[usx(ii) usy(ii) usz(ii)]*V(:,2)/(absus(ii)*(sum(V(:,2).^2))^0.5);
%                                 cos_uo_l2(count_strain,4)=[uox(ii) uoy(ii) uoz(ii)]*V(:,2)/(absuo(ii)*(sum(V(:,2).^2))^0.5);
% 
                                cos_a_l3(count_strain,4)=[ax(ii) ay(ii) az(ii)]*V(:,3)/(absa(ii)*(sum(V(:,3).^2))^0.5);
%                                 cos_ao_l3(count_strain,4)=[aOx(ii) aOy(ii) aOz(ii)]*V(:,3)/(absaO(ii)*(sum(V(:,3).^2))^0.5);
%                                 cos_ap_l3(count_strain,4)=[aPx(ii) aPy(ii) aPz(ii)]*V(:,3)/(absaP(ii)*(sum(V(:,3).^2))^0.5);
                                cos_ac_l3(count_strain,4)=[acx(ii) acy(ii) acz(ii)]*V(:,3)/(absac(ii)*(sum(V(:,3).^2))^0.5);
%                                 cos_al_l3(count_strain,4)=[alx(ii) aly(ii) alz(ii)]*V(:,3)/(absal(ii)*(sum(V(:,3).^2))^0.5);
%                                 cos_us_l3(count_strain,4)=[usx(ii) usy(ii) usz(ii)]*V(:,3)/(absus(ii)*(sum(V(:,3).^2))^0.5);
%                                 cos_uo_l3(count_strain,4)=[uox(ii) uoy(ii) uoz(ii)]*V(:,3)/(absuo(ii)*(sum(V(:,3).^2))^0.5);
%                                 
%                                 cos_a_u(count_strain,4)=[ax(ii) ay(ii) az(ii)]*[u(ii) v(ii) w(ii)]'/(absa(ii)*absu(ii));
%                                 cos_ao_u(count_strain,4)=[aOx(ii) aOy(ii) aOz(ii)]*[u(ii) v(ii) w(ii)]'/(absaO(ii)*absu(ii));
%                                 cos_ap_u(count_strain,4)=[aPx(ii) aPy(ii) aPz(ii)]*[u(ii) v(ii) w(ii)]'/(absaP(ii)*absu(ii));
%                                 cos_ac_u(count_strain,4)=[acx(ii) acy(ii) acz(ii)]*[u(ii) v(ii) w(ii)]'/(absac(ii)*absu(ii));
%                                 cos_al_u(count_strain,4)=[alx(ii) aly(ii) alz(ii)]*[u(ii) v(ii) w(ii)]'/(absal(ii)*absu(ii));
%                                 cos_us_u(count_strain,4)=[usx(ii) usy(ii) usz(ii)]*[u(ii) v(ii) w(ii)]'/(absus(ii)*absu(ii));
%                                 cos_uo_u(count_strain,4)=[uox(ii) uoy(ii) uoz(ii)]*[u(ii) v(ii) w(ii)]'/(absuo(ii)*absu(ii));
%                                 
%                                 cos_alac(count_strain,4)=[alx(ii) aly(ii) alz(ii)]*[acx(ii) acy(ii) acz(ii)]'/(absal(ii)*absac(ii));
%                                 cos_usac(count_strain,4)=[usx(ii) usy(ii) usz(ii)]*[acx(ii) acy(ii) acz(ii)]'/(absus(ii)*absac(ii));
%                                 cos_uoac(count_strain,4)=[uox(ii) uoy(ii) uoz(ii)]*[acx(ii) acy(ii) acz(ii)]'/(absuo(ii)*absac(ii));
%                                 cos_usuo(count_strain,4)=[usx(ii) usy(ii) usz(ii)]*[uox(ii) uoy(ii) uoz(ii)]'/(absus(ii)*absuo(ii));
%                                 cos_apac(count_strain,4)=[aPx(ii) aPy(ii) aPz(ii)]*[acx(ii) acy(ii) acz(ii)]'/(absaP(ii)*absac(ii));
%                                 cos_apus(count_strain,4)=[aPx(ii) aPy(ii) aPz(ii)]*[usx(ii) usy(ii) usz(ii)]'/(absaP(ii)*absus(ii));
%                                 cos_aouo(count_strain,4)=[aOx(ii) aOy(ii) aOz(ii)]*[uox(ii) uoy(ii) uoz(ii)]'/(absaO(ii)*absuo(ii));
                            end
                        end
                    end
                end
            end
        end
    end
end

%figure 4a
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
xlabel('|a|,low,enstro,strain','FontSize',12,'FontName','Times New Roman')
ylabel('pdf','FontSize',12,'FontName','Times New Roman')
set(gca,'XScale','log')
set(gca,'YScale','log')

%figure 5a
figure;
hold on;
grid on;
box on;
[n,x]=nhist(hist_absa(1:counter,1));
plot(x,n,'k');
[n,x]=nhist(hist_absap(1:counter,1),1000);
plot(x,n,'g');
[n,x]=nhist(hist_absao(1:counter,1),1000);
plot(x,n,'r');
[n,x]=nhist(hist_absal(1:counter,1),1000);
plot(x,n,'c');
[n,x]=nhist(hist_absac(1:counter,1),1000);
plot(x,n,'b');
xlabel('|a|,||,\perp,a_l,a_c','FontSize',12,'FontName','Times New Roman')
ylabel('pdf','FontSize',12,'FontName','Times New Roman')
set(gca,'XScale','log')
set(gca,'YScale','log')

% figure;
% hold on;
% grid on;
% box on;
% [n,x]=nhist(hist_absao(1:counter,1),1000);
% plot(x,n,'k');
% [n,x]=nhist(hist_absao(1:count_low,2),1000);
% plot(x,n,'r');
% [n,x]=nhist(hist_absao(1:count_enstro,3),1000);
% plot(x,n,'g');
% [n,x]=nhist(hist_absao(1:count_strain,4),1000);
% plot(x,n,'b');
% xlabel('|a_{\perp}|','FontSize',12,'FontName','Times New Roman')
% ylabel('pdf','FontSize',12,'FontName','Times New Roman')
% set(gca,'XScale','log')
% set(gca,'YScale','log')

% figure;
% hold on;
% grid on;
% box on;
% [n,x]=nhist(hist_absap(1:counter,1),1000);
% plot(x,n,'k');
% [n,x]=nhist(hist_absap(1:count_low,2),1000);
% plot(x,n,'r');
% [n,x]=nhist(hist_absap(1:count_enstro,3),1000);
% plot(x,n,'g');
% [n,x]=nhist(hist_absap(1:count_strain,4),1000);
% plot(x,n,'b');
% xlabel('|a_{||}|','FontSize',12,'FontName','Times New Roman')
% ylabel('pdf','FontSize',12,'FontName','Times New Roman')
% set(gca,'XScale','log')
% set(gca,'YScale','log')

% figure;
% hold on;
% grid on;
% box on;
% [n,x]=nhist(hist_absal(1:counter,1),1000);
% plot(x,n,'k');
% [n,x]=nhist(hist_absal(1:count_low,2),1000);
% plot(x,n,'r');
% [n,x]=nhist(hist_absal(1:count_enstro,3),1000);
% plot(x,n,'g');
% [n,x]=nhist(hist_absal(1:count_strain,4),1000);
% plot(x,n,'b');
% xlabel('|a_l|','FontSize',12,'FontName','Times New Roman')
% ylabel('pdf','FontSize',12,'FontName','Times New Roman')
% set(gca,'XScale','log')
% set(gca,'YScale','log')


%figure 5b
figure;
hold on;
grid on;
box on;
[n,x]=nhist(hist_curv(1:counter,1),1000);
ind=find(n.*x==max(n.*x));
val=x(ind)
[n,x]=nhist(hist_curv(1:counter,1)/val,1000);
plot(x,n,'k');
[n,x]=nhist(hist_curv(1:count_low,2)/val,1000);
plot(x,n,'r');
[n,x]=nhist(hist_curv(1:count_enstro,3)/val,1000);
plot(x,n,'g');
[n,x]=nhist(hist_curv(1:count_strain,4)/val,1000);
plot(x,n,'b');
xlabel('curvature','FontSize',12,'FontName','Times New Roman')
ylabel('pdf','FontSize',12,'FontName','Times New Roman')
set(gca,'XScale','log')
set(gca,'YScale','log')

%figure 5c
figure;
hold on;
grid on;
box on;
[n,x]=nhist(hist_curvac(1:counter,1),1000);
ind=find(n.*x==max(n.*x));
val=x(ind)
[n,x]=nhist(hist_curvac(1:counter,1)/val,1000);
plot(x,n,'k');
[n,x]=nhist(hist_curvac(1:count_low,2)/val,1000);
plot(x,n,'r');
[n,x]=nhist(hist_curvac(1:count_enstro,3)/val,1000);
plot(x,n,'g');
[n,x]=nhist(hist_curvac(1:count_strain,4)/val,1000);
plot(x,n,'b');
xlabel('curvature a_c','FontSize',12,'FontName','Times New Roman')
ylabel('pdf','FontSize',12,'FontName','Times New Roman')
set(gca,'XScale','log')
set(gca,'YScale','log')

% figure;
% hold on;
% grid on;
% box on;
% [n,x]=nhist(hist_absus(1:counter,1),1000);
% plot(x,n,'k');
% [n,x]=nhist(hist_absus(1:count_low,2),1000);
% plot(x,n,'r');
% [n,x]=nhist(hist_absus(1:count_enstro,3),1000);
% plot(x,n,'g');
% [n,x]=nhist(hist_absus(1:count_strain,4),1000);
% plot(x,n,'b');
% xlabel('|us|','FontSize',12,'FontName','Times New Roman')
% ylabel('pdf','FontSize',12,'FontName','Times New Roman')
% set(gca,'XScale','log')
% set(gca,'YScale','log')

% figure;
% hold on;
% grid on;
% box on;
% [n,x]=nhist(hist_absuo(1:counter,1),1000);
% plot(x,n,'k');
% [n,x]=nhist(hist_absuo(1:count_low,2),1000);
% plot(x,n,'r');
% [n,x]=nhist(hist_absuo(1:count_enstro,3),1000);
% plot(x,n,'g');
% [n,x]=nhist(hist_absuo(1:count_strain,4),1000);
% plot(x,n,'b');
% xlabel('|uo|','FontSize',12,'FontName','Times New Roman')
% ylabel('pdf','FontSize',12,'FontName','Times New Roman')
% set(gca,'XScale','log')
% set(gca,'YScale','log')

%figure 6d?
figure
hold on;
grid on;
box on;
[n,x]=nhist(cos_a_u(1:counter,1),50);
plot(x,n,'k');
[n,x]=nhist(cos_a_u(1:count_enstro,3),50);
plot(x,n,'k--');
[n,x]=nhist(cos_a_u(1:count_strain,4),50);
plot(x,n,'k:');
[n,x]=nhist(cos_ac_u(1:counter,1),50);
plot(x,n,'g');
[n,x]=nhist(cos_ac_u(1:count_enstro,3),50);
plot(x,n,'g--');
[n,x]=nhist(cos_ac_u(1:count_strain,4),50);
plot(x,n,'g:');

xlabel('cos(a,u), cos(a_c,u)','FontSize',12,'FontName','Times New Roman')
ylabel('pdf','FontSize',12,'FontName','Times New Roman')
axis([-1 1 0 1])

%figure 7a
figure;
subplot(1,2,1)
hold on;
grid on;
box on;
[n,x]=nhist(cos_a_w(1:counter,1),50);
plot(x,n,'k');
[n,x]=nhist(cos_ao_w(1:counter,1),50);
plot(x,n,'r');
[n,x]=nhist(cos_ap_w(1:counter,1),50);
plot(x,n,'c');
[n,x]=nhist(cos_ac_w(1:counter,1),50);
plot(x,n,'g');
[n,x]=nhist(cos_a_w(1:count_enstro,3),50);
plot(x,n,'k--');
[n,x]=nhist(cos_ao_w(1:count_enstro,3),50);
plot(x,n,'r--');
[n,x]=nhist(cos_ap_w(1:count_enstro,3),50);
plot(x,n,'c--');
[n,x]=nhist(cos_ac_w(1:count_enstro,3),50);
plot(x,n,'g--');
xlabel('cos(a,\omega)','FontSize',12,'FontName','Times New Roman')
ylabel('pdf','FontSize',12,'FontName','Times New Roman')
axis([-1 1 0 2])

%figure 7b
subplot(1,2,2)
hold on;
grid on;
box on;
[n,x]=nhist(cos_a_w(1:counter,1),50);
plot(x,n,'k');
[n,x]=nhist(cos_ao_w(1:counter,1),50);
plot(x,n,'r');
[n,x]=nhist(cos_ap_w(1:counter,1),50);
plot(x,n,'c');
[n,x]=nhist(cos_ac_w(1:counter,1),50);
plot(x,n,'g');
[n,x]=nhist(cos_a_w(1:count_strain,4),50);
plot(x,n,'k--');
[n,x]=nhist(cos_ao_w(1:count_strain,4),50);
plot(x,n,'r--');
[n,x]=nhist(cos_ap_w(1:count_strain,4),50);
plot(x,n,'c--');
[n,x]=nhist(cos_ac_w(1:count_strain,4),50);
plot(x,n,'g--');
xlabel('cos(a,\omega)','FontSize',12,'FontName','Times New Roman')
ylabel('pdf','FontSize',12,'FontName','Times New Roman')
axis([-1 1 0 2])

% figure;
% hold on;
% grid on;
% box on;
% [n,x]=nhist(cos_a_w(1:count_low,2),50);
% plot(x,n,'k');
% [n,x]=nhist(cos_ao_w(1:count_low,2),50);
% plot(x,n,'r');
% [n,x]=nhist(cos_ap_w(1:count_low,2),50);
% plot(x,n,'g');
% [n,x]=nhist(cos_al_w(1:count_low,2),50);
% plot(x,n,'b');
% [n,x]=nhist(cos_ac_w(1:count_low,2),50);
% plot(x,n,'c');
% [n,x]=nhist(cos_us_w(1:count_low,2),50);
% plot(x,n,'m');
% xlabel('cos(a,\omega) low','FontSize',12,'FontName','Times New Roman')
% ylabel('pdf','FontSize',12,'FontName','Times New Roman')
% axis([-1 1 0 2])

% figure;
% hold on;
% grid on;
% box on;
% [n,x]=nhist(cos_a_w(1:count_enstro,3),50);
% plot(x,n,'k');
% [n,x]=nhist(cos_ao_w(1:count_enstro,3),50);
% plot(x,n,'r');
% [n,x]=nhist(cos_ap_w(1:count_enstro,3),50);
% plot(x,n,'g');
% [n,x]=nhist(cos_al_w(1:count_enstro,3),50);
% plot(x,n,'b');
% [n,x]=nhist(cos_ac_w(1:count_enstro,3),50);
% plot(x,n,'c');
% [n,x]=nhist(cos_us_w(1:count_enstro,3),50);
% plot(x,n,'m');
% xlabel('cos(a,\omega) enstro','FontSize',12,'FontName','Times New Roman')
% ylabel('pdf','FontSize',12,'FontName','Times New Roman')
% axis([-1 1 0 2])

% figure;
% hold on;
% grid on;
% box on;
% [n,x]=nhist(cos_a_w(1:count_strain,4),50);
% plot(x,n,'k');
% [n,x]=nhist(cos_ao_w(1:count_strain,4),50);
% plot(x,n,'r');
% [n,x]=nhist(cos_ap_w(1:count_strain,4),50);
% plot(x,n,'g');
% [n,x]=nhist(cos_al_w(1:count_strain,4),50);
% plot(x,n,'b');
% [n,x]=nhist(cos_ac_w(1:count_strain,4),50);
% plot(x,n,'c');
% [n,x]=nhist(cos_us_w(1:count_strain,4),50);
% plot(x,n,'m');
% xlabel('cos(a,\omega) strain','FontSize',12,'FontName','Times New Roman')
% ylabel('pdf','FontSize',12,'FontName','Times New Roman')
% axis([-1 1 0 2])

% figure;
% hold on;
% grid on;
% box on;
% [n,x]=nhist(cos_a_u(1:counter,1),50);
% plot(x,n,'k');
% [n,x]=nhist(cos_al_u(1:counter,1),50);
% plot(x,n,'b');
% [n,x]=nhist(cos_ac_u(1:counter,1),50);
% plot(x,n,'c');
% [n,x]=nhist(cos_us_u(1:counter,1),50);
% plot(x,n,'m');
% xlabel('cos(u,a) all','FontSize',12,'FontName','Times New Roman')
% ylabel('pdf','FontSize',12,'FontName','Times New Roman')
% axis([-1 1 0 2])
% 
% figure;
% hold on;
% grid on;
% box on;
% [n,x]=nhist(cos_a_u(1:count_low,2),50);
% plot(x,n,'k');
% [n,x]=nhist(cos_al_u(1:count_low,2),50);
% plot(x,n,'b');
% [n,x]=nhist(cos_ac_u(1:count_low,2),50);
% plot(x,n,'c');
% [n,x]=nhist(cos_us_u(1:count_low,2),50);
% plot(x,n,'m');
% xlabel('cos(u,a) low','FontSize',12,'FontName','Times New Roman')
% ylabel('pdf','FontSize',12,'FontName','Times New Roman')
% axis([-1 1 0 2])
% 
% figure;
% hold on;
% grid on;
% box on;
% [n,x]=nhist(cos_a_u(1:count_enstro,3),50);
% plot(x,n,'k');
% [n,x]=nhist(cos_al_u(1:count_enstro,3),50);
% plot(x,n,'b');
% [n,x]=nhist(cos_ac_u(1:count_enstro,3),50);
% plot(x,n,'c');
% [n,x]=nhist(cos_us_u(1:count_enstro,3),50);
% plot(x,n,'m');
% xlabel('cos(u,a) enstro','FontSize',12,'FontName','Times New Roman')
% ylabel('pdf','FontSize',12,'FontName','Times New Roman')
% axis([-1 1 0 2])

% figure;
% hold on;
% grid on;
% box on;
% [n,x]=nhist(cos_a_u(1:count_strain,4),50);
% plot(x,n,'k');
% [n,x]=nhist(cos_al_u(1:count_strain,4),50);
% plot(x,n,'b');
% [n,x]=nhist(cos_ac_u(1:count_strain,4),50);
% plot(x,n,'c');
% [n,x]=nhist(cos_us_u(1:count_strain,4),50);
% plot(x,n,'m');
% xlabel('cos(u,a) strain','FontSize',12,'FontName','Times New Roman')
% ylabel('pdf','FontSize',12,'FontName','Times New Roman')
% axis([-1 1 0 2])

%figure 8
figure;
subplot(1,3,1);
hold on;
grid on;
box on;
[n,x]=nhist(cos_a_l1(1:counter,1),50);
plot(x,n,'k');
[n,x]=nhist(cos_ac_l1(1:counter,1),50);
plot(x,n,'g');
[n,x]=nhist(cos_a_l1(1:count_enstro,3),50);
plot(x,n,'k--');
[n,x]=nhist(cos_ac_l1(1:count_enstro,3),50);
plot(x,n,'g--');
[n,x]=nhist(cos_a_l1(1:count_strain,4),50);
plot(x,n,'k:');
[n,x]=nhist(cos_ac_l1(1:count_strain,4),50);
plot(x,n,'g:');
xlabel('cos(a,\lambda_1) all','FontSize',12,'FontName','Times New Roman')
ylabel('pdf','FontSize',12,'FontName','Times New Roman')
axis([-1 1 0 2])
subplot(1,3,2);
hold on;
grid on;
box on;
[n,x]=nhist(cos_a_l2(1:counter,1),50);
plot(x,n,'k');
[n,x]=nhist(cos_ac_l2(1:counter,1),50);
plot(x,n,'g');
[n,x]=nhist(cos_a_l2(1:count_enstro,3),50);
plot(x,n,'k--');
[n,x]=nhist(cos_ac_l2(1:count_enstro,3),50);
plot(x,n,'g--');
[n,x]=nhist(cos_a_l2(1:count_strain,4),50);
plot(x,n,'k:');
[n,x]=nhist(cos_ac_l2(1:count_strain,4),50);
plot(x,n,'g:');
xlabel('cos(a,\lambda_2) all','FontSize',12,'FontName','Times New Roman')
axis([-1 1 0 2])
subplot(1,3,3);
hold on;
grid on;
box on;
[n,x]=nhist(cos_a_l3(1:counter,1),50);
plot(x,n,'k');
[n,x]=nhist(cos_ac_l3(1:counter,1),50);
plot(x,n,'g');
[n,x]=nhist(cos_a_l3(1:count_enstro,3),50);
plot(x,n,'k--');
[n,x]=nhist(cos_ac_l3(1:count_enstro,3),50);
plot(x,n,'g--');
[n,x]=nhist(cos_a_l3(1:count_strain,4),50);
plot(x,n,'k:');
[n,x]=nhist(cos_ac_l3(1:count_strain,4),50);
plot(x,n,'g:');
xlabel('cos(a,\lambda_3) all','FontSize',12,'FontName','Times New Roman')
axis([-1 1 0 2])


% figure;
% subplot(1,3,1);
% hold on;
% grid on;
% box on;
% [n,x]=nhist(cos_a_l1(1:counter,1),50);
% plot(x,n,'k');
% [n,x]=nhist(cos_ao_l1(1:counter,1),50);
% plot(x,n,'r');
% [n,x]=nhist(cos_ap_l1(1:counter,1),50);
% plot(x,n,'g');
% [n,x]=nhist(cos_al_l1(1:counter,1),50);
% plot(x,n,'b');
% [n,x]=nhist(cos_ac_l1(1:counter,1),50);
% plot(x,n,'c');
% [n,x]=nhist(cos_us_l1(1:counter,1),50);
% plot(x,n,'m');
% [n,x]=nhist(cos_uo_l1(1:counter,1),50);
% plot(x,n,'y');
% xlabel('cos(a,\lambda_1) all','FontSize',12,'FontName','Times New Roman')
% ylabel('pdf','FontSize',12,'FontName','Times New Roman')
% axis([-1 1 0 2])
% 
% subplot(1,3,2);
% hold on;
% grid on;
% box on;
% [n,x]=nhist(cos_a_l2(1:counter,1),50);
% plot(x,n,'k');
% [n,x]=nhist(cos_ao_l2(1:counter,1),50);
% plot(x,n,'r');
% [n,x]=nhist(cos_ap_l2(1:counter,1),50);
% plot(x,n,'g');
% [n,x]=nhist(cos_al_l2(1:counter,1),50);
% plot(x,n,'b');
% [n,x]=nhist(cos_ac_l2(1:counter,1),50);
% plot(x,n,'c');
% [n,x]=nhist(cos_us_l2(1:counter,1),50);
% plot(x,n,'m');
% [n,x]=nhist(cos_uo_l2(1:counter,1),50);
% plot(x,n,'y');
% xlabel('cos(a,\lambda_2) all','FontSize',12,'FontName','Times New Roman')
% ylabel('pdf','FontSize',12,'FontName','Times New Roman')
% axis([-1 1 0 2])
% 
% subplot(1,3,3);
% hold on;
% grid on;
% box on;
% [n,x]=nhist(cos_a_l3(1:counter,1),50);
% plot(x,n,'k');
% [n,x]=nhist(cos_ao_l3(1:counter,1),50);
% plot(x,n,'r');
% [n,x]=nhist(cos_ap_l3(1:counter,1),50);
% plot(x,n,'g');
% [n,x]=nhist(cos_al_l3(1:counter,1),50);
% plot(x,n,'b');
% [n,x]=nhist(cos_ac_l3(1:counter,1),50);
% plot(x,n,'c');
% [n,x]=nhist(cos_us_l3(1:counter,1),50);
% plot(x,n,'m');
% [n,x]=nhist(cos_uo_l3(1:counter,1),50);
% plot(x,n,'y');
% xlabel('cos(a,\lambda_3) all','FontSize',12,'FontName','Times New Roman')
% ylabel('pdf','FontSize',12,'FontName','Times New Roman')
% axis([-1 1 0 2])
% 
% figure;
% subplot(1,3,1);
% hold on;
% grid on;
% box on;
% [n,x]=nhist(cos_a_l1(1:count_low,2),50);
% plot(x,n,'k');
% [n,x]=nhist(cos_ao_l1(1:count_low,2),50);
% plot(x,n,'r');
% [n,x]=nhist(cos_ap_l1(1:count_low,2),50);
% plot(x,n,'g');
% [n,x]=nhist(cos_al_l1(1:count_low,2),50);
% plot(x,n,'b');
% [n,x]=nhist(cos_ac_l1(1:count_low,2),50);
% plot(x,n,'c');
% [n,x]=nhist(cos_us_l1(1:count_low,2),50);
% plot(x,n,'m');
% [n,x]=nhist(cos_uo_l1(1:count_low,2),50);
% plot(x,n,'y');
% xlabel('cos(a,\lambda_1) low','FontSize',12,'FontName','Times New Roman')
% ylabel('pdf','FontSize',12,'FontName','Times New Roman')
% axis([-1 1 0 2])
% 
% subplot(1,3,2);
% hold on;
% grid on;
% box on;
% [n,x]=nhist(cos_a_l2(1:count_low,2),50);
% plot(x,n,'k');
% [n,x]=nhist(cos_ao_l2(1:count_low,2),50);
% plot(x,n,'r');
% [n,x]=nhist(cos_ap_l2(1:count_low,2),50);
% plot(x,n,'g');
% [n,x]=nhist(cos_al_l2(1:count_low,2),50);
% plot(x,n,'b');
% [n,x]=nhist(cos_ac_l2(1:count_low,2),50);
% plot(x,n,'c');
% [n,x]=nhist(cos_us_l2(1:count_low,2),50);
% plot(x,n,'m');
% [n,x]=nhist(cos_uo_l2(1:count_low,2),50);
% plot(x,n,'y');
% xlabel('cos(a,\lambda_2) low','FontSize',12,'FontName','Times New Roman')
% ylabel('pdf','FontSize',12,'FontName','Times New Roman')
% axis([-1 1 0 2])
% 
% subplot(1,3,3);
% hold on;
% grid on;
% box on;
% [n,x]=nhist(cos_a_l3(1:count_low,2),50);
% plot(x,n,'k');
% [n,x]=nhist(cos_ao_l3(1:count_low,2),50);
% plot(x,n,'r');
% [n,x]=nhist(cos_ap_l3(1:count_low,2),50);
% plot(x,n,'g');
% [n,x]=nhist(cos_al_l3(1:count_low,2),50);
% plot(x,n,'b');
% [n,x]=nhist(cos_ac_l3(1:count_low,2),50);
% plot(x,n,'c');
% [n,x]=nhist(cos_us_l3(1:count_low,2),50);
% plot(x,n,'m');
% [n,x]=nhist(cos_uo_l3(1:count_low,2),50);
% plot(x,n,'y');
% xlabel('cos(a,\lambda_3) low','FontSize',12,'FontName','Times New Roman')
% ylabel('pdf','FontSize',12,'FontName','Times New Roman')
% axis([-1 1 0 2])

% figure;
% subplot(1,3,1);
% hold on;
% grid on;
% box on;
% [n,x]=nhist(cos_a_l1(1:count_enstro,3),50);
% plot(x,n,'k');
% [n,x]=nhist(cos_ao_l1(1:count_enstro,3),50);
% plot(x,n,'r');
% [n,x]=nhist(cos_ap_l1(1:count_enstro,3),50);
% plot(x,n,'g');
% [n,x]=nhist(cos_al_l1(1:count_enstro,3),50);
% plot(x,n,'b');
% [n,x]=nhist(cos_ac_l1(1:count_enstro,3),50);
% plot(x,n,'c');
% [n,x]=nhist(cos_us_l1(1:count_enstro,3),50);
% plot(x,n,'m');
% [n,x]=nhist(cos_uo_l1(1:count_enstro,3),50);
% plot(x,n,'y');
% xlabel('cos(a,\lambda_1) enstro','FontSize',12,'FontName','Times New Roman')
% ylabel('pdf','FontSize',12,'FontName','Times New Roman')
% axis([-1 1 0 2])
% 
% subplot(1,3,2);
% hold on;
% grid on;
% box on;
% [n,x]=nhist(cos_a_l2(1:count_enstro,3),50);
% plot(x,n,'k');
% [n,x]=nhist(cos_ao_l2(1:count_enstro,3),50);
% plot(x,n,'r');
% [n,x]=nhist(cos_ap_l2(1:count_enstro,3),50);
% plot(x,n,'g');
% [n,x]=nhist(cos_al_l2(1:count_enstro,3),50);
% plot(x,n,'b');
% [n,x]=nhist(cos_ac_l2(1:count_enstro,3),50);
% plot(x,n,'c');
% [n,x]=nhist(cos_us_l2(1:count_enstro,3),50);
% plot(x,n,'m');
% [n,x]=nhist(cos_uo_l2(1:count_enstro,3),50);
% plot(x,n,'y');
% xlabel('cos(a,\lambda_2) enstro','FontSize',12,'FontName','Times New Roman')
% ylabel('pdf','FontSize',12,'FontName','Times New Roman')
% axis([-1 1 0 2])
% 
% subplot(1,3,3);
% hold on;
% grid on;
% box on;
% [n,x]=nhist(cos_a_l3(1:count_enstro,3),50);
% plot(x,n,'k');
% [n,x]=nhist(cos_ao_l3(1:count_enstro,3),50);
% plot(x,n,'r');
% [n,x]=nhist(cos_ap_l3(1:count_enstro,3),50);
% plot(x,n,'g');
% [n,x]=nhist(cos_al_l3(1:count_enstro,3),50);
% plot(x,n,'b');
% [n,x]=nhist(cos_ac_l3(1:count_enstro,3),50);
% plot(x,n,'c');
% [n,x]=nhist(cos_us_l3(1:count_enstro,3),50);
% plot(x,n,'m');
% [n,x]=nhist(cos_uo_l3(1:count_enstro,3),50);
% plot(x,n,'y');
% xlabel('cos(a,\lambda_3) enstro','FontSize',12,'FontName','Times New Roman')
% ylabel('pdf','FontSize',12,'FontName','Times New Roman')
% axis([-1 1 0 2])

% figure;
% subplot(1,3,1);
% hold on;
% grid on;
% box on;
% [n,x]=nhist(cos_a_l1(1:count_strain,4),50);
% plot(x,n,'k');
% [n,x]=nhist(cos_ao_l1(1:count_strain,4),50);
% plot(x,n,'r');
% [n,x]=nhist(cos_ap_l1(1:count_strain,4),50);
% plot(x,n,'g');
% [n,x]=nhist(cos_al_l1(1:count_strain,4),50);
% plot(x,n,'b');
% [n,x]=nhist(cos_ac_l1(1:count_strain,4),50);
% plot(x,n,'c');
% [n,x]=nhist(cos_us_l1(1:count_strain,4),50);
% plot(x,n,'m');
% [n,x]=nhist(cos_uo_l1(1:count_strain,4),50);
% plot(x,n,'y');
% xlabel('cos(a,\lambda_1) strain','FontSize',12,'FontName','Times New Roman')
% ylabel('pdf','FontSize',12,'FontName','Times New Roman')
% axis([-1 1 0 2])

% subplot(1,3,2);
% hold on;
% grid on;
% box on;
% [n,x]=nhist(cos_a_l2(1:count_strain,4),50);
% plot(x,n,'k');
% [n,x]=nhist(cos_ao_l2(1:count_strain,4),50);
% plot(x,n,'r');
% [n,x]=nhist(cos_ap_l2(1:count_strain,4),50);
% plot(x,n,'g');
% [n,x]=nhist(cos_al_l2(1:count_strain,4),50);
% plot(x,n,'b');
% [n,x]=nhist(cos_ac_l2(1:count_strain,4),50);
% plot(x,n,'c');
% [n,x]=nhist(cos_us_l2(1:count_strain,4),50);
% plot(x,n,'m');
% [n,x]=nhist(cos_uo_l2(1:count_strain,4),50);
% plot(x,n,'y');
% xlabel('cos(a,\lambda_2) strain','FontSize',12,'FontName','Times New Roman')
% ylabel('pdf','FontSize',12,'FontName','Times New Roman')
% axis([-1 1 0 2])

% subplot(1,3,3);
% hold on;
% grid on;
% box on;
% [n,x]=nhist(cos_a_l3(1:count_strain,4),50);
% plot(x,n,'k');
% [n,x]=nhist(cos_ao_l3(1:count_strain,4),50);
% plot(x,n,'r');
% [n,x]=nhist(cos_ap_l3(1:count_strain,4),50);
% plot(x,n,'g');
% [n,x]=nhist(cos_al_l3(1:count_strain,4),50);
% plot(x,n,'b');
% [n,x]=nhist(cos_ac_l3(1:count_strain,4),50);
% plot(x,n,'c');
% [n,x]=nhist(cos_us_l3(1:count_strain,4),50);
% plot(x,n,'m');
% [n,x]=nhist(cos_uo_l3(1:count_strain,4),50);
% plot(x,n,'y');
% xlabel('cos(a,\lambda_3) strain','FontSize',12,'FontName','Times New Roman')
% ylabel('pdf','FontSize',12,'FontName','Times New Roman')
% axis([-1 1 0 2])


% figure;
% hold on;
% grid on;
% box on;
% [n,x]=nhist(cos_alac(1:counter,1),50);
% plot(x,n,'k');
% [n,x]=nhist(cos_alac(1:count_low,2),50);
% plot(x,n,'r');
% [n,x]=nhist(cos_alac(1:count_enstro,3),50);
% plot(x,n,'g');
% [n,x]=nhist(cos_alac(1:count_strain,4),50);
% plot(x,n,'b');
% xlabel('cos(a_l,a_c) ','FontSize',12,'FontName','Times New Roman')
% ylabel('pdf','FontSize',12,'FontName','Times New Roman')
% axis([-1 1 0 10])

% figure;
% hold on;
% grid on;
% box on;
% [n,x]=nhist(cos_usac(1:counter,1),50);
% plot(x,n,'k');
% [n,x]=nhist(cos_usac(1:count_low,2),50);
% plot(x,n,'r');
% [n,x]=nhist(cos_usac(1:count_enstro,3),50);
% plot(x,n,'g');
% [n,x]=nhist(cos_usac(1:count_strain,4),50);
% plot(x,n,'b');
% xlabel('cos(us,a_c) ','FontSize',12,'FontName','Times New Roman')
% ylabel('pdf','FontSize',12,'FontName','Times New Roman')
% axis([-1 1 0 8])
% 
% figure;
% hold on;
% grid on;
% box on;
% [n,x]=nhist(cos_uoac(1:counter,1),50);
% plot(x,n,'k');
% [n,x]=nhist(cos_uoac(1:count_low,2),50);
% plot(x,n,'r');
% [n,x]=nhist(cos_uoac(1:count_enstro,3),50);
% plot(x,n,'g');
% [n,x]=nhist(cos_uoac(1:count_strain,4),50);
% plot(x,n,'b');
% xlabel('cos(uo,a_c) ','FontSize',12,'FontName','Times New Roman')
% ylabel('pdf','FontSize',12,'FontName','Times New Roman')
% axis([-1 1 0 8])

% figure;
% hold on;
% grid on;
% box on;
% [n,x]=nhist(cos_usuo(1:counter,1),50);
% plot(x,n,'k');
% [n,x]=nhist(cos_usuo(1:count_low,2),50);
% plot(x,n,'r');
% [n,x]=nhist(cos_usuo(1:count_enstro,3),50);
% plot(x,n,'g');
% [n,x]=nhist(cos_usuo(1:count_strain,4),50);
% plot(x,n,'b');
% xlabel('cos(us,uo) ','FontSize',12,'FontName','Times New Roman')
% ylabel('pdf','FontSize',12,'FontName','Times New Roman')
% axis([-1 1 0 2])
% 
% figure;
% hold on;
% grid on;
% box on;
% [n,x]=nhist(cos_apac(1:counter,1),50);
% plot(x,n,'k');
% [n,x]=nhist(cos_apac(1:count_low,2),50);
% plot(x,n,'r');
% [n,x]=nhist(cos_apac(1:count_enstro,3),50);
% plot(x,n,'g');
% [n,x]=nhist(cos_apac(1:count_strain,4),50);
% plot(x,n,'b');
% xlabel('cos(a_{||},a_c) ','FontSize',12,'FontName','Times New Roman')
% ylabel('pdf','FontSize',12,'FontName','Times New Roman')
% axis([-1 1 0 2])

% figure;
% hold on;
% grid on;
% box on;
% [n,x]=nhist(cos_apus(1:counter,1),50);
% plot(x,n,'k');
% [n,x]=nhist(cos_apus(1:count_low,2),50);
% plot(x,n,'r');
% [n,x]=nhist(cos_apus(1:count_enstro,3),50);
% plot(x,n,'g');
% [n,x]=nhist(cos_apus(1:count_strain,4),50);
% plot(x,n,'b');
% xlabel('cos(a_{||},us) ','FontSize',12,'FontName','Times New Roman')
% ylabel('pdf','FontSize',12,'FontName','Times New Roman')
% axis([-1 1 0 2])
% 
% figure;
% hold on;
% grid on;
% box on;
% [n,x]=nhist(cos_aouo(1:counter,1),50);
% plot(x,n,'k');
% [n,x]=nhist(cos_aouo(1:count_low,2),50);
% plot(x,n,'r');
% [n,x]=nhist(cos_aouo(1:count_enstro,3),50);
% plot(x,n,'g');
% [n,x]=nhist(cos_aouo(1:count_strain,4),50);
% plot(x,n,'b');
% xlabel('cos(a_{\perp},uo) ','FontSize',12,'FontName','Times New Roman')
% ylabel('pdf','FontSize',12,'FontName','Times New Roman')
% axis([-1 1 0 2])


