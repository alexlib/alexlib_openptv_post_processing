function res=trajecPointDivCheck(first,last);

vywz=zeros(101,101);
uxwz=zeros(101,101);
uxvy=zeros(101,101);
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
   
   s11=f(:,13);
   s12=f(:,14);
   s13=f(:,15);
   s22=f(:,16);
   s23=f(:,17);
   s33=f(:,18);
   reldiv=f(:,27);
   age=f(:,29);

   for ii=1:s(1,1)
      if age(ii)>-1  
         indJ=floor((s22(ii)+s33(ii)+2)*25)+1;
         indI=floor((-s11(ii)+2)*25)+1;
         if indI>0 &indJ>0 & indI<102 & indJ<102
            vywz(indI,indJ)=vywz(indI,indJ)+1;
         end
         indJ=floor((s11(ii)+s33(ii)+2)*25)+1;
         indI=floor((-s22(ii)+2)*25)+1;
         if indI>0 &indJ>0 & indI<102 & indJ<102
            uxwz(indI,indJ)=uxwz(indI,indJ)+1;
         end
         indJ=floor((s11(ii)+s22(ii)+2)*25)+1;
         indI=floor((-s33(ii)+2)*25)+1;
         if indI>0 &indJ>0 & indI<102 & indJ<102
            uxvy(indI,indJ)=uxvy(indI,indJ)+1;
         end
      end      
   end
end

save vywz
save uxwz
save uxvy

figure
x=-2:0.04:2;
contourf(x,x,vywz,10);
xlabel('v_y+w_z','FontSize',22);
ylabel('-u_x','FontSize',22);
axis equal
axis([-2 2 -2 2])
colorbar;
figure
contourf(x,x,uxwz,10);
xlabel('u_x+w_z','FontSize',22);
ylabel('-v_y','FontSize',22);
axis equal
axis([-2 2 -2 2])
colorbar;
figure
contourf(x,x,uxvy,10);
xlabel('u_x+v_y','FontSize',22);
ylabel('-w_z','FontSize',22);
axis equal
axis([-2 2 -2 2])
colorbar;


   
   
