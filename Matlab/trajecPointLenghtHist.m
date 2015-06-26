function res=trajecPointLenghtHist(first,last);
   
ageHistmai10=zeros(100,1);
divHistmai10=zeros(100,1);
enstro1Histmai10=zeros(100,1);
strain1Histmai10=zeros(100,1);
enstro2Histmai10=zeros(100,1);
strain2Histmai10=zeros(100,1);
eKin1Histmai10=zeros(100,1);
eKin2Histmai10=zeros(100,1);
asq1Histmai10=zeros(100,1);
asq2Histmai10=zeros(100,1);
pointsForStatmai10=zeros(100,1);
totalPointsmai10=zeros(100,1);

   
if 1<2
   
   
   counter=0;
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
      u=f(:,4);
      v=f(:,5);
      w=f(:,6);
      ax=f(:,7);
      ay=f(:,8);
      az=f(:,9);
      asq=ax.^2+ay.^2+az.^2;
      eKin=0.5*(u.^2+v.^2+w.^2);
      w1=f(:,10);
      w2=f(:,11);
      w3=f(:,12);
      enstro=w1.^2+w2.^2+w3.^2;
      strain=f(:,26)/2.4;
      divref=f(:,27);
      age=f(:,29);
      
      for ii=1:s(1,1)-1
         if age(ii)==0
            counter=0;
            counter2=0;
            div=0;
            stra1=0;
            enst1=0;
            stra2=0;
            enst2=0;
            eKin1=0;
            eKin2=0;
            asq1=0;
            asq2=0;
         end
         counter=counter+1;
         div=div+divref(ii);
         enst1=enst1+enstro(ii);
         stra1=stra1+strain(ii);
         eKin1=eKin1+eKin(ii);
         asq1=asq1+asq(ii);
         if divref(ii)<0.1 & strain(ii)<400 & enstro(ii)<400 & eKin(ii)<0.005
            enst2=enst2+enstro(ii);
            stra2=stra2+strain(ii);
            eKin2=eKin2+eKin(ii);
            asq2=asq2+asq(ii);
            counter2=counter2+1;
         end   
         if age(ii+1)==0
            div=div/counter;
            enst1=enst1/counter;
            stra1=stra1/counter;
            eKin1=eKin1/counter;
            asq1=asq1/counter;
            if counter2>0
               enst2=enst2/counter2;
               stra2=stra2/counter2;
               eKin2=eKin2/counter2;
               asq2=asq2/counter2;
            else
               enst2=0;
               stra2=0;
               eKin2=0;
               asq2=0;
            end
            if age(ii)<100
               ageHistmai10(age(ii)+1,1)=ageHistmai10(age(ii)+1,1)+1;
               divHistmai10(age(ii)+1,1)=(divHistmai10(age(ii)+1,1)*(ageHistmai10(age(ii)+1,1)-1)+div)/ageHistmai10(age(ii)+1,1);
               enstro1Histmai10(age(ii)+1,1)=(enstro1Histmai10(age(ii)+1,1)*(ageHistmai10(age(ii)+1,1)-1)+enst1)/ageHistmai10(age(ii)+1,1);
               strain1Histmai10(age(ii)+1,1)=(strain1Histmai10(age(ii)+1,1)*(ageHistmai10(age(ii)+1,1)-1)+stra1)/ageHistmai10(age(ii)+1,1);
               enstro2Histmai10(age(ii)+1,1)=(enstro2Histmai10(age(ii)+1,1)*(ageHistmai10(age(ii)+1,1)-1)+enst2)/ageHistmai10(age(ii)+1,1);
               strain2Histmai10(age(ii)+1,1)=(strain2Histmai10(age(ii)+1,1)*(ageHistmai10(age(ii)+1,1)-1)+stra2)/ageHistmai10(age(ii)+1,1);
               eKin1Histmai10(age(ii)+1,1)=(eKin1Histmai10(age(ii)+1,1)*(ageHistmai10(age(ii)+1,1)-1)+eKin1)/ageHistmai10(age(ii)+1,1);
               eKin2Histmai10(age(ii)+1,1)=(eKin2Histmai10(age(ii)+1,1)*(ageHistmai10(age(ii)+1,1)-1)+eKin2)/ageHistmai10(age(ii)+1,1);
               asq1Histmai10(age(ii)+1,1)=(asq1Histmai10(age(ii)+1,1)*(ageHistmai10(age(ii)+1,1)-1)+asq1)/ageHistmai10(age(ii)+1,1);
               asq2Histmai10(age(ii)+1,1)=(asq2Histmai10(age(ii)+1,1)*(ageHistmai10(age(ii)+1,1)-1)+asq2)/ageHistmai10(age(ii)+1,1);
            end
         end
      end
   end
   
   save ageHistmai10
   save divHistmai10
   save enstro1Histmai10
   save strain1Histmai10
   save enstro2Histmai10
   save strain2Histmai10
   save eKin1Histmai10
   save eKin2Histmai10
   save asq1Histmai10
   save asq2Histmai10
end
load ageHistmai10
load divHistmai10
load enstro1Histmai10
load strain1Histmai10
load enstro2Histmai10
load strain2Histmai10
load eKin1Histmai10
load eKin2Histmai10
load asq1Histmai10
load asq2Histmai10

for i=1:100
   pointsForStatmai10(i,1)=ageHistmai10(i,1)*i;
   for j=i:100
      totalPointsmai10(i,1)=totalPointsmai10(i,1)+pointsForStatmai10(i,1);
   end
end
for l=1:0
   for i=3:99
      dummy1(i)=(1/3)*(enstro1Histmai10(i-1,1)+enstro1Histmai10(i,1)+enstro1Histmai10(i+1,1));
      dummy2(i)=(1/3)*(enstro2Histmai10(i-1,1)+enstro2Histmai10(i,1)+enstro2Histmai10(i+1,1));
      dummy3(i)=(1/3)*(strain1Histmai10(i-1,1)+strain1Histmai10(i,1)+strain1Histmai10(i+1,1));
      dummy4(i)=(1/3)*(strain2Histmai10(i-1,1)+strain2Histmai10(i,1)+strain2Histmai10(i+1,1));   
   end
   for i=3:99
      enstro1Histmai10(i,1)=dummy1(i);
      enstro2Histmai10(i,1)=dummy2(i);
      strain1Histmai10(i,1)=dummy3(i);
      strain2Histmai10(i,1)=dummy4(i);
   end
end

x=0:0.0714:7.1;
semilogy(x,ageHistmai10,'-rv')
hold on;
semilogy(x,pointsForStatmai10,'-b*');
semilogy(x,totalPointsmai10,'-ks');

h = legend('#trajectories','#points','total points',0);
xlabel('\tau_{\eta}','FontSize',22);
ylabel('#','FontSize',22);
grid on;
figure
plot(x,divHistmai10)
xlabel('\tau_{\eta}','FontSize',22);
ylabel('<rel. div>','FontSize',22);
grid on;


figure
plot(x,enstro1Histmai10,':ko')
hold on;
plot(x,strain1Histmai10,':k')
plot(x,enstro2Histmai10,'-ko','LineWidth',2)
plot(x,strain2Histmai10,'-k','LineWidth',2)
h = legend('<\omega^2>, enstrophy','<s^2>, strain','<\omega^2>, enstrophy, rel.div<0.1','<s^2>, strain,rel.div<0.1',0);
xlabel('\tau_{\eta}','FontSize',22);
ylabel('\omega^2, s^2','FontSize',22);
grid on;

figure
plot(x,eKin2Histmai10,'-ko','LineWidth',2)
hold on;
plot(x,eKin1Histmai10,':k')
h = legend('e_{kin}, rel.div<0.1','e_{kin}',0);
xlabel('\tau_{\eta}','FontSize',22);
ylabel('e_{kin}','FontSize',22);
grid on;

figure
plot(x,asq2Histmai10,'-ko','LineWidth',2)
hold on;
plot(x,asq1Histmai10,':k')
h = legend('a^2, rel.div<0.1','a^2',0);
xlabel('\tau_{\eta}','FontSize',22);
ylabel('a^2','FontSize',22);
grid on;



