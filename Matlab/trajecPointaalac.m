function res=trajecPointaalac(first,last);

cosineAlAc=zeros(41,1);
jsxLi=zeros(101,101);
jsyLi=zeros(101,101);
jszLi=zeros(101,101);
good=0;

for i=first:last
   clear f;
   i
   name='trajPoint.';
   le=length(name);
   ext=int2str(i);
   nd=length(ext);
   for j=1:nd
      name(le+j)=ext(j);
   end
   f=load(name);
   s=size(f);
   
   xr=f(:,1);
   yr=f(:,2);
   zr=f(:,3);
   u=f(:,4);
   v=f(:,5);
   w=f(:,6);
   ax=f(:,7);
   ay=f(:,8);
   az=f(:,9);
   absa=(ax.^2+ay.^2+az.^2).^0.5;
   alx=f(:,30);%%%
   aly=f(:,31);%%%
   alz=f(:,32);%%%
   absal=(alx.^2+aly.^2+alz.^2).^0.5;%%%
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
   acx=u.*ux+v.*uy+w.*uz;
   acy=u.*vx+v.*vy+w.*vz;
   acz=u.*wx+v.*wy+w.*wz;
   absac=(acx.^2+acy.^2+acz.^2).^0.5;  
   %alx=ax-acx;%%%
   %aly=ay-acy;%%%
   %alz=az-acz;%%%
   %absal=(alx.^2+aly.^2+alz.^2).^0.5;%%%
   reldiv=f(:,27);
   
   for ii=1:s(1,1)
      diffx=ax(ii)-alx(ii)-acx(ii);
      diffy=ay(ii)-aly(ii)-acy(ii);
      diffz=az(ii)-alz(ii)-acz(ii);
      absDiff=(diffx^2+diffy^2+diffz^2)^0.5;
      if reldiv(ii)<0.1 %& absDiff<0.01 
         index=floor(((alx(ii)*acx(ii)+aly(ii)*acy(ii)+alz(ii)*acz(ii))/(absal(ii)*absac(ii))+1)/0.05)+1;
         if index>0 &index<42
            cosineAlAc(index)=cosineAlAc(index)+1;
         end   
         good=good+1;
         indJ=floor((ax(ii)+0.01)*5000)+1;
         indI=floor((alx(ii)+acx(ii)+0.01)*5000)+1;
         if indI>0 &indJ>0 & indI<102 & indJ<102
            jsxLi(indI,indJ)=jsxLi(indI,indJ)+1;
         end
         indJ=floor((ay(ii)+0.01)*5000)+1;
         indI=floor((aly(ii)+acy(ii)+0.01)*5000)+1;
         if indI>0 &indJ>0 & indI<102 & indJ<102
            jsyLi(indI,indJ)=jsyLi(indI,indJ)+1;
         end
         indJ=floor((az(ii)+0.01)*5000)+1;
         indI=floor((alz(ii)+acz(ii)+0.01)*5000)+1;
         if indI>0 &indJ>0 & indI<102 & indJ<102
            jszLi(indI,indJ)=jszLi(indI,indJ)+1;
         end
      end      
   end
   good
end
good

save cosineAlAc
save jsxLi
save jsyLi
save jszLi

figure
s=size(cosineAlAc)
fact=(s(1,1)/2)/sum(cosineAlAc)
x=-1:0.05:1;
plot(x,cosineAlAc*fact,'k');
xlabel('cos(a_l,a_c)','FontSize',22);
ylabel('pdf','FontSize',22);
figure
x=-0.01:0.0002:0.01;
contourf(x,x,jsxLi,10);
xlabel('a_{x}, rel. div<0.1','FontSize',22);
ylabel('a_{l,x}+a_{c,x}','FontSize',22);
axis equal
axis([-0.01 0.01 -0.01 0.01])
colorbar;
figure
contourf(x,x,jsyLi,10);
xlabel('a_{y}, rel. div<0.1','FontSize',22);
ylabel('a_{l,y}+a_{c,y}','FontSize',22);
axis equal
axis([-0.01 0.01 -0.01 0.01])
colorbar;
figure
contourf(x,x,jszLi,10);
xlabel('a_{z}, rel. div<0.1','FontSize',22);
ylabel('a_{l,z}+a_{c,z}','FontSize',22);
axis equal
axis([-0.01 0.01 -0.01 0.01])
colorbar;


   
   
