function res=trajecPointCoswl_i(first,last);

if 1>2
   counter=0;
   eiCounter=0;
   counter1=0;
   counter2=0;
   counter3=0;
   counter4=0;
   coswl1=zeros(41,1);
   coswl2=zeros(41,1);
   coswl3=zeros(41,1);
   ei=zeros(2000000,3);
   
   good=0;
   
   for i=first:1:last
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
      si=size(f);
      
      xr=f(:,1);
      yr=f(:,2);
      zr=f(:,3);
      u=f(:,4);
      v=f(:,5);
      w=f(:,6);
      w1=f(:,10);
      w2=f(:,11);
      w3=f(:,12);
      
      s11=f(:,13);
      s12=f(:,14);
      s13=f(:,15);
      s22=f(:,16);
      s23=f(:,17);
      s33=f(:,18);
      
      absw=(w1.*w1+w2.*w2+w3.*w3).^0.5;
      divref=f(:,27);
      
      for ii=1:si(1,1)
         if divref(ii)<0.1
            counter=counter+1;
            good=good+1;
            
            A=[s11(ii) s12(ii) s13(ii);s12(ii) s22(ii) s23(ii);s13(ii) s23(ii) s33(ii)];
            
            [V,D]=eig(A);
            
            if D(1,1)<D(2,2)
               hV=V(:,1);
               hL=D(1,1);
               V(:,1)=V(:,2);
               D(1,1)=D(2,2);
               V(:,2)=hV;
               D(2,2)=hL;
            end
            if D(2,2)<D(3,3)
               hV=V(:,2);
               hL=D(2,2);
               V(:,2)=V(:,3);
               D(2,2)=D(3,3);
               V(:,3)=hV;
               D(3,3)=hL;
            end
            if D(1,1)<D(2,2)
               hV=V(:,1);
               hL=D(1,1);
               V(:,1)=V(:,2);
               D(1,1)=D(2,2);
               V(:,2)=hV;
               D(2,2)=hL;
            end
            
            if abs(D(1,1))<20 & abs(D(2,2))<20 & abs(D(3,3))<20
               eiCounter=eiCounter+1;
               ei(eiCounter,1)=D(1,1);
               ei(eiCounter,2)=D(2,2);
               ei(eiCounter,3)=D(3,3);
            end
            
            hiwl1=(w1(ii)*V(1,1)+w2(ii)*V(2,1)+w3(ii)*V(3,1))/absw(ii);
            hiwl2=(w1(ii)*V(1,2)+w2(ii)*V(2,2)+w3(ii)*V(3,2))/absw(ii);
            hiwl3=(w1(ii)*V(1,3)+w2(ii)*V(2,3)+w3(ii)*V(3,3))/absw(ii);
            
            index=floor((hiwl1+1)/0.045);
            if index>0 & index<42
               coswl1(index)=coswl1(index)+1;
            end
            
            index=floor((hiwl2+1)/0.045);
            if index>0 & index<42
               coswl2(index)=coswl2(index)+1;
            end
            
            index=floor((hiwl3+1)/0.045);
            if index>0 & index<42
               coswl3(index)=coswl3(index)+1;
            end
            
         end
      end
      good=good;
   end
   
   save coswl1
   save coswl2
   save coswl3
else
   load coswl1
   load coswl2
   load coswl3
end

figure
s=size(coswl1)
fact=(s(1,1)/2)/sum(coswl1)
x=-1:0.05:1;
plot(x,coswl1*fact,'r');
hold on;
s=size(coswl2)
fact=(s(1,1)/2)/sum(coswl2)
x=-1:0.05:1;
plot(x,coswl2*fact,'g');
s=size(coswl3)
fact=(s(1,1)/2)/sum(coswl3)
x=-1:0.05:1;
plot(x,coswl3*fact,'b');


[n,xout]=hist(ei(1:eiCounter,1),1000);;
figure;
plot(xout,20*n./sum(n),'r');
%xlim([-5 5]);
title('lambda 1,2,3');
hold on;
[n,xout]=hist(ei(1:eiCounter,2),100);;
plot(xout,20*n./sum(n),'g');
[n,xout]=hist(ei(1:eiCounter,3),100);;
plot(xout,20*n./sum(n),'b');%


%%%%%%%%5Arkady figure







