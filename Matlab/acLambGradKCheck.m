function res=acLambGradKCheck(first,last,tail,stepSize);

ac=zeros(1500000,3);
LambGrad=zeros(1500000,3);
tot=0;

for i=first:stepSize:last
   clear f;
   i
   %name='../JFMruns/40points4mmNoKriging/trajPoint.';
   %name='F:/mai10/trajPoint.';
   %name='F:/26April2004_water/fromAlex/trajAcc.'; %processed from Alex,linear 2,4mm
   name='F:/26April2004_water/fromAlex/trajAcc.';
   le=length(name);
   ext=int2str(i);
   nd=length(ext);
   for j=1:nd
      name(le+j)=ext(j);
   end
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
      ux=f(:,7);
      uy=f(:,8);
      uz=f(:,9);
      vx=f(:,10);
      vy=f(:,11);
      vz=f(:,12);
      wx=f(:,13);
      wy=f(:,14);
      wz=f(:,15);
      w1=wy-vz;
      w2=uz-wx;
      w3=vx-uy;
      enstrophy=w1.^2+w2.^2+w3.^2;
      s11=0.5*(ux+ux);
      s12=0.5*(uy+vx);
      s13=0.5*(uz+wx);
      s22=0.5*(vy+vy);
      s23=0.5*(vz+wy);
      s33=0.5*(wz+wz);
      
      wS1=w1.*s11+w2.*s12+w3.*s13;
      wS2=w1.*s12+w2.*s22+w3.*s23;
      wS3=w1.*s13+w2.*s23+w3.*s33;
      wws=w1.*wS1+w2.*wS2+w3.*wS3;
      
      alx=f(:,16);
      aly=f(:,17);
      alz=f(:,18);
      absal=(alx.^2+aly.^2+alz.^2).^0.5;%%%
      
      acx=u.*ux+v.*uy+w.*uz;
      acy=u.*vx+v.*vy+w.*vz;
      acz=u.*wx+v.*wy+w.*wz;
      absac=(acx.^2+acy.^2+acz.^2).^0.5;  
      
      ax=f(:,19);
      ay=f(:,20);
      az=f(:,21);
%       ax=alx+acx;
%       ay=aly+acy;
%       az=alz+acz;
      absa=(ax.^2+ay.^2+az.^2).^0.5;
      
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
      
      reldiv=abs(ux+vy+wz)./(abs(ux)+abs(vy)+abs(wz));
      age=f(:,32);
      
      nx=f(:,33);
      ny=f(:,34);
      nz=f(:,35);      
      dnxdx=f(:,36);
      dnxdy=f(:,37);
      dnxdz=f(:,38);
      dnydx=f(:,39);
      dnydy=f(:,40);
      dnydz=f(:,41);
      dnzdx=f(:,42);
      dnzdy=f(:,43);
      dnzdz=f(:,44);
      
      dkdx=f(:,45);
      dkdy=f(:,46);
      dkdz=f(:,47);
      
      cnx=nx.*dnxdx+ny.*dnxdy+nz.*dnxdz;
      cny=nx.*dnydx+ny.*dnydy+nz.*dnydz;
      cnz=nx.*dnzdx+ny.*dnzdy+nz.*dnzdz;
      curvGrad=(cnx.^2+cny.^2+cnz.^2).^0.5;
      
      curv=(((v.*az-w.*ay).^2+(w.*ax-u.*az).^2+(u.*ay-v.*ax).^2).^0.5)./(absu.^3); 
      curvO=absaO./absu.^2; 
      
      reldiv=abs(ux+vy+wz)./(abs(ux)+abs(vy)+abs(wz));
      age=f(:,32);
      
      ende=0;
      for iii=1:s(1,1)-1
         start=ende+1;
         ende=start;
         length=1;
         searching=1;
         while ende<s(1,1)-1 & searching==1
            if age(ende+1)==age(ende)+1
               ende=ende+1;
               length=length+1;
            else
               searching=0;
            end
         end
         counter=0;
         oSquare=0;
         if ende<s(1,1)
            for ii=start:ende
               if age(ii)>tail & age(ii)<length-tail+1 & reldiv(ii)<0.1
                   tot=tot+1;
                   LaGr=cross([w1(ii) w2(ii) w3(ii)]',[u(ii) v(ii) w(ii)]')+0.5*[dkdx(ii) dkdy(ii) dkdz(ii)]';
                   LambGrad(tot,1)=LaGr(1);
                   ac(tot,1)=acx(ii);
                   LambGrad(tot,2)=LaGr(2);
                   ac(tot,2)=acy(ii);
                   LambGrad(tot,3)=LaGr(3);
                   ac(tot,3)=acz(ii);
               end
            end
         end
      end
      tot
   end
end
tot


limit = 0.04;
ind = LambGrad(1:tot,1) > -limit & LambGrad(1:tot,1) < limit & ac(1:tot,1) > -limit & ac(1:tot,1) < limit;
hf = figure
hold on
[n,x,bins] = histmulti5([LambGrad(ind,1),ac(ind,1)]);
contourf(x(:,1),x(:,2),n',20)
grid on
set(gca,'xlim',[-limit limit], 'ylim',[-limit limit]);
xlabel('\omega x u+1/2\nablau^2, x');
ylabel('a_c, x');

ind = LambGrad(1:tot,2) > -limit & LambGrad(1:tot,2) < limit & ac(1:tot,2) > -limit & ac(1:tot,2) < limit;
hf = figure
hold on
[n,x,bins] = histmulti5([LambGrad(ind,2),ac(ind,2)]);
contourf(x(:,1),x(:,2),n',20)
grid on
set(gca,'xlim',[-limit limit], 'ylim',[-limit limit]);
xlabel('\omega x u+1/2\nablau^2, y');
ylabel('a_c, y');

ind = LambGrad(1:tot,3) > -limit & LambGrad(1:tot,3) < limit & ac(1:tot,3) > -limit & ac(1:tot,3) < limit;
hf = figure
hold on
[n,x,bins] = histmulti5([LambGrad(ind,3),ac(ind,3)]);
contourf(x(:,1),x(:,2),n',20)
grid on
set(gca,'xlim',[-limit limit], 'ylim',[-limit limit]);
xlabel('\omega x u+1/2\nablau^2, z');
ylabel('a_c, z');

   
