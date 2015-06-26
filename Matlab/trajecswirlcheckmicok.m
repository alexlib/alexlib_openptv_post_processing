function res=trajecswirlcheckmicok(first,last,len,minT,maxT);
numT=0;

% ------compute RQ diagram according to Chong anmd Perry for P=0
lenn =100;
Rbound=zeros(2*lenn,1);
Qbound=zeros(2*lenn,1);
zero=zeros(2*lenn,1);
begg=5;
Qbound(1)=-begg;
for i=2:lenn
    Qbound(i)= Qbound(1)+begg*i/lenn;
end
for i=lenn:2*lenn-1
    Qbound(i)= Qbound(i-lenn+1);
end

for i=1:lenn
    Rbound(i)=-sqrt(-4/27*Qbound(i)^3);
end
for i=lenn:2*lenn-1
    Rbound(i)=sqrt(-4/27*Qbound(i)^3);
end
%------------------------------- end RQ diagram
for n=first:last
   n;
  name='D:\michele\trajPoints_a\trajPoint.';
   le=length(name);
   ext=int2str(n);
   nd=length(ext);
   for j=1:nd
      name(le+j)=ext(j);
   end
   f=load(name);
   si=size(f);
   %------si is the (varying from file to file) length 
   %------of each trajpoint file  (rem: swi(1,1) is the number of rows 
   %------while swi(1,2) is the number of columns (fixed = 32)
   
   swi=zeros(si(1,1),1);
   b=zeros(si(1,1),1);
   swiswi=zeros(si(1,1),1);
   bb=zeros(si(1,1),1);
   sig=zeros(si(1,1),1);
   Rc=zeros(si(1,1),1);
   Qc=zeros(si(1,1),1);
     
   Rbeg=zeros(si(1,1),1);
   Qbeg=zeros(si(1,1),1);
   swibeg=zeros(si(1,1),1);
   xbeg=zeros(si(1,1),1);
   ybeg=zeros(si(1,1),1);
   zbeg=zeros(si(1,1),1);
   bbeg=zeros(si(1,1),1);
   sigbeg=zeros(si(1,1),1);
   zbeg(:)=-70;
      
    
   x=1000*f(:,1);
   y=1000*f(:,2);
   z=1000*f(:,3);
   
   xx=1000*f(:,1);
   yy=1000*f(:,2);
   zz=1000*f(:,3);
   
   w1=f(:,10);
   w2=f(:,11);
   w3=f(:,12);
   
   s11=f(:,13);
   s12=f(:,14);
   s13=f(:,15);
   s22=f(:,16);
   s23=f(:,17);
   s33=f(:,18);
   
   R=f(:,24);
   Q=f(:,25);
   RR=f(:,24);
   QQ=f(:,25);
   wws=f(:,22);
   age=f(:,29); 
   div=f(:,27);
      
   
   
   begi=0;
   ende=0;
   while begi<si(1,1)-2   %stop to look for beginning traj. 2 rows
                          %before the end of fole of trajpoint...
      begi=ende+1;
      ende=begi;
      fou=0;    %start a traj: for the whole traj fou = 0
      
      while fou==0 & ende<si(1,1)-1
            if age(ende+1)>age(ende)
%***********computation velocity gradients matrix A (chong and perry)
                ii=ende;
                begi;
                ende;
                a11=s11(ii);
                a12=s12(ii)-w3(ii)/2.;
                a13=s13(ii)+w2(ii)/2.;
   
                a21=s12(ii)+w3(ii)/2.;
			    a22=s22(ii);
			    a23=s23(ii)-w1(ii)/2.;
		        
                a31=s13(ii)-w2(ii)/2.;
		        a32=s23(ii)+w1(ii)/2.;
			    a33=s33(ii);
			    
                A=zeros(3,3);
                A=[a11 a12 a13
                a21 a22 a23
                a31 a32 a33];
			
                P=-trace(A);
                Qcalc=0.5*(P^2.-trace(A^2));    
                Rcalc=-det(A);
                lam=eig(A);
                
                sigg=real(lam(2))  ;
                sig(ii)=sigg;
                Rc(ii)=Rcalc;
                Qc(ii)=Qcalc;
                
                swirl=abs(imag(lam(2)));
                swi(ii)=swirl;
        
                swiswi(ii)=swi(ii);
    
    %check: if the eigenvalues are all real it does not make sense,
    % so let's run the check only for points where 2 eigv are c.c.
    
                if swirl~=0
                    if real(lam(1))==real(lam(2))      
                      bbb=real(lam(3));
                      B=[sigg swirl 0
                        -swirl sigg 0
                        0 0 bbb];
                    else
                       bbb=real(lam(1));
                       B=[bbb 0 0
                         0 sigg swirl 
                         0 -swirl sigg ];    
                    end   
                 
                     b(ii)=bbb;
                     bb(ii)=bbb;
                     lam;           % eig of A
                     bbb;          
                     lam2=eig(B);   % eig of B
                 end
    %  fine controllo
    
                  if 1>2      
                  if div(begi)>0.1
                     w1(begi)=0;
                     w2(begi)=0;
                     w3(begi)=0;
                     wws(begi)=0;
                     
                     %%%%%%%%%%%%%
                        x(begi)=0;
                        y(begi)=0;
                        z(begi)=-70;
                     R(begi)=0;
                     Q(begi)=0;         
                     swi(begi)=0;
                     b(begi)=0;
                     sig(begi)=0;
                     Rc(begi)=0;
                     Qc(begi)=0;  
                     
                     %%%%%%%%%%%%%
                 end
                 end
                                 
                    
                    if 3>2                        
                        if swi(ii)==0;                            
                             x(ii)=0;
                             y(ii)=0;
                             z(ii)=-70;
                             Rc(ii)=0;
                             Qc(ii)=0; 
                             R(ii)=0;
                             Q(ii)=0;
                             sig(ii)=0;
                             b(ii)=0;
                         end
                     end
                                                                  
                 
                if ii>begi+1
                    % scelgo l'inizio di ogni porzione di traiettoria in cui 
                    % gli autovalori sono complessi
                    if 3>2
                        if swi(ii-1)==0
                            if x(ii)~=0          
                             xbeg(ii)=x(ii);
                             ybeg(ii)=y(ii);
                             zbeg(ii)=z(ii);
                             Rbeg(ii)=R(ii);
                             Qbeg(ii)=Q(ii);
                             swibeg(ii)=swi(ii);            
                             bbeg(ii)=b(ii);
                             sigbeg(ii)=sig(ii);
                            end                  
                        end
                    end
                end
               
               
                ende;
                if 1>2
                if div(ende)>0.1    
                      x(ende)=0;
                      y(ende)=0;
                      z(ende)=-70;
                      R(ende)=0;
                      Q(ende)=0;
                      swi(ende)=0;
                      b(ende)=0;
                      sig(ende)=0;             
                end
                end
                ende=ende+1;    
             % I'm still within the same traj
             
            else       % if age(ende+1)>age(ende) 
                fou=1;      %uscito dalla traiettoria                   
            end        % if age(ende+1)>age(ende) 
      end   %while fou==0 & ende<si(1,1)-1
      
      if ende-begi>len
         numT=numT+1;
         if numT>minT-1 & numT<maxT+1
            if min(z(begi:ende))>-73
                                
         
                figure(1)
                plot3(x(begi:ende),y(begi:ende),z(begi:ende),'rd');
                hold on
                 plot3(xx(begi:ende),yy(begi:ende),zz(begi:ende),'LineWidth',0.5);
                hold on
        %         plot3(xx(begi:ende),yy(begi:ende),zz(begi:ende),'bd');
        %        hold on
               plot3(xbeg(begi:ende),ybeg(begi:ende),zbeg(begi:ende),'bd'); 
                hold on
                plot3(x(begi),y(begi),z(begi),'kx');
               
                xlabel('x (mm)','FontSize',22)
                ylabel('y (mm)','FontSize',22)
                zlabel('z (mm)','FontSize',22)                                                 
                
            end
            hold on; 
            box on;
            grid on;
                                                           
         end
      end            
      
       %%%%%%%%%%%%%%%%%%
         if 1>2
         if ende-begi>len
         numT=numT+1;
            if numT>minT-1 & numT<maxT+1
              if min(z(begi:ende))>-73
                                
              figure(2)           
              %plot(swiswi(begi:ende),bb(begi:ende),'LineWidth',0.5); 
              %hold on
              plot(swi(begi:ende),b(begi:ende),'bd'); 
              hold on
              plot(swi(begi),b(begi),'yx');
              hold on            
              plot(swibeg(begi:ende),bbeg(begi:ende),'y*');   
              xlabel('lam ','FontSize',22)
              ylabel('b ','FontSize',22)
                
              end
              hold on; 
              box on;
              grid on;                                                           
            end
        end       
        end
      %*******************
      if ende-begi>len
         numT=numT+1;
         if numT>minT-1 & numT<maxT+1
            if min(z(begi:ende))>-73
                    
               figure(3)
               plot3(QQ(begi:ende),RR(begi:ende),swi(begi:ende),'yd'); 
               hold on
               plot3(Qc(begi:ende),Rc(begi:ende),swi(begi:ende),'kx'); 
               hold on
               plot3(Q(begi:ende),R(begi:ende),swi(begi:ende),'bd'); 
               hold on
               plot3(Q(begi),R(begi),swi(begi),'kx');
               hold on            
               plot3(Qbeg(begi:ende),Rbeg(begi:ende),swibeg(begi:ende),'k+');   
               hold on  
               plot3(Rbound(:),Qbound(:),zero(:),'rd');
               xlabel('R ','FontSize',22)
               ylabel('Q ','FontSize',22)
               zlabel('lam','FontSize',22)
              

            end
            hold on;
            box on;
            grid on;                                    
            
         end
      end
      %*******************

      if ende-begi>len
         numT=numT+1;
         if numT>minT-1 & numT<maxT+1
            if min(z(begi:ende))>-73
                     
               figure(4)
%              plot3(QQ(begi:ende),RR(begi:ende),swi(begi:ende),'LineWidth',0.5); 
%              hold on
               plot3(b(begi:ende),sig(begi:ende),swi(begi:ende),'kd'); 
               hold on
               plot3(bb(begi),sig(begi),swiswi(begi),'kx');
%              hold on
               plot3(bbeg(begi:ende),sigbeg(begi:ende),swibeg(begi:ende),'yd');   
               hold on 
%              plot3(Qbeg(begi:ende),Rbeg(begi:ende),swibeg(begi:ende),'kd');   
               
               xlabel('b ','FontSize',22)
               ylabel('sig ','FontSize',22)
               zlabel('lam','FontSize',22)
            end
            hold on;
            box on;
            grid on;                                    
           
         end
      end

      %*******************          
                  
      
   end  %while begi<si(1,1)-2   keep looking at the same frame 
   maxAge=max(age);
 
end  % for n=first:last  repeat for another frame

numT=numT;
     