function res=RQ_grid(first,last,tail,stepSize);
%RQ_grid(1,6000,0,200)

counter=0;
jRQ=zeros(81,41);
st=zeros(2000000,1);
mean_s=0.69;
%mean_s=0.4;

allpoints=0;


for i=first:stepSize:last
    clear f;
    if mod(i-first,10)==0
    i 
    counter/allpoints
    end
    name=['X:/rotating/PTV08/080811/run1/coarse_grid_raw.',num2Str(i)];
    f=load(name);
    s=size(f);
    if s(1,1)>0
        xr=f(:,1);
        yr=f(:,2);
        zr=f(:,3);
        u=f(:,4);
        v=f(:,5);
        w=f(:,6);
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
        s11=f(:,13);
        s12=f(:,14);
        s13=f(:,15);
        s22=f(:,16);
        s23=f(:,17);
        s33=f(:,18);
        ux=s11;
        uy=s12-0.5*w3;
        uz=s13+0.5*w2;
        vx=s12+0.5*w3;
        vy=s22;
        vz=s23-0.5*w1;
        wx=s13-0.5*w2;
        wy=s23+0.5*w1;
        wz=s33;
        acx=u.*ux+v.*uy+w.*uz;%%%%%%%%%%%%%%%%%%%%
        acy=u.*vx+v.*vy+w.*vz;%%%%%%%%%%%%%%%%%%%%
        acz=u.*wx+v.*wy+w.*wz;%%%%%%%%%%%%%%%%%%%%
        absac=(acx.^2+acy.^2+acz.^2).^0.5;%%%%%%%%%%%%%%%%%%%%%%%%%
        aLamx=w2.*w-w3.*v;%%%%%%%%%%%%%%%%%%%%%%%
        aLamy=w3.*u-w1.*w;%%%%%%%%%%%%%%%%%%%%%%%
        aLamz=w1.*v-w2.*u;%%%%%%%%%%%%%%%%%%%%%%%%%
        absaLam=(aLamx.^2+aLamy.^2+aLamz.^2).^0.5;%%%%%%%%%%%%%%%%%
        aBx=acx-aLamx;%%%%%%%%%%%%%%%%%%%%%
        aBy=acy-aLamy;%%%%%%%%%%%%%%%%%%%%%%%
        aBz=acz-aLamz;%%%%%%%%%%%%%%%%%%%%%%%%%
        absaB=(aBx.^2+aBy.^2+aBz.^2).^0.5;%%%%%%%%%%%%%%%%%%%
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
        
        daxdx=f(:,22);
        daxdy=f(:,23);
        daxdz=f(:,24);
        daydx=f(:,25);
        daydy=f(:,26);
        daydz=f(:,27);
        dazdx=f(:,28);
        dazdy=f(:,29);
        dazdz=f(:,30);
        
        reldiv=f(:,31);
       
        %rel_error=(1/3)*(abs(ax-alx-acx)./(abs(ax)+abs(alx)+abs(acx))+abs(ay-aly-acy)./(abs(ay)+abs(aly)+abs(acy))+abs(az-alz-acz)./(abs(az)+abs(alz)+abs(acz)));
        
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
        for ii=1:s(1,1)-1

            allpoints=allpoints+1;
            if isnan(u(ii))==0 & yr(ii)>0.03 %& rel_error(ii)<0.8 %& reldiv(ii)<1 %absa(ii)^2/8e-4<100 %& curv(ii)<3000
                counter=counter+1;

                uM=[ux(ii) uy(ii) uz(ii);vx(ii) vy(ii) vz(ii);wx(ii) wy(ii) wz(ii)];
                A=[s11(ii) s12(ii) s13(ii);s12(ii) s22(ii) s23(ii);s13(ii) s23(ii) s33(ii)];
                [V,D] = eig(A);
                [D,k] = sort(diag(D)); 	% ascending order
                D = diag(D(end:-1:1)); 	% descending order
                V = V(:,k(end:-1:1)); 	% the same order for eigenvectors
                strain=D(1,1)^2+D(2,2)^2+D(3,3)^2;
                if counter<2000000 & strain<100
                    st(counter)=strain^0.5;
                end

                Q=-(1/2)*(trace(uM^2))/mean_s^2;
                R=-(1/3)*(trace(uM^3))/mean_s^3;

                indJ=floor(R*10+0.5)+20;%+/-100
                indI=floor(Q*10+0.5)+40;%+  100
                if indI>0 &indJ>0 & indI<82 & indJ<42 & strain<100
                    jRQ(indI,indJ)=jRQ(indI,indJ)+1;
                end
            end
        end
    end
end

success=counter/allpoints
counter=counter
if counter>2000000
    counter=2000000;
end

x=0:40;
x=(x-20)/10;
y=0:80;
y=(y-40)/10;
figure
hold on;
contourf(x,y,log(jRQ+1)/log(10),20);
h=['R'];
xlabel(h,'FontSize',12,'FontName','Times New Roman')
ylabel('Q','FontSize',12,'FontName','Times New Roman');
title('log_{10}[p(R,Q)]','FontSize',12,'FontName','Times New Roman')
colorbar
set(gca,'FontSize',12,'FontName','Times New Roman')
axis([-2 2 -4 4])

figure
hold on;
[n,x]=nhist(st(1:counter),1000);
plot(x,n.*x);
h=['s'];
xlabel(h,'FontSize',12,'FontName','Times New Roman')
ylabel('pdf','FontSize',12,'FontName','Times New Roman');
set(gca,'XScale','log')
set(gca,'YScale','log')
mean_s=mean(st(1:counter))







