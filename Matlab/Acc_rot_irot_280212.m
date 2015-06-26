function  Acc_rot_irot_280212(first,last,tail,stepSize)
% version 28.02.12
% Alex, updated for polymers data
% reprocessed 

if ~nargin
    first = 100001; 
    last  = 103000;
    tail = 1;
    stepSize = 1;
    
end

counter=0;
range=50;
statistic=zeros(1200000,3);

if 1<2
    %%%%1017288
    for i=first:stepSize:last
        clear f;
        fprintf(1,'%d\n',i);
        %name='../../Mai10/trajPoint.';
        %name='D:/PTVAlexPolymer75rpm/trajPoint.';
        %name='D:/PTVAlexFeb17/trajPoint.';
        % name='C:\acceleration\april26_alex/trajAcc.';
        name ='/Users/alex/Documents/PTV/0104_water/n_trajPoint.';
        
        
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
            ax=f(:,19);
            ay=f(:,20);
            az=f(:,21);
            absa=(ax.^2+ay.^2+az.^2).^0.5;
            alx=f(:,16);
            aly=f(:,17);
            alz=f(:,18);
            absal=(alx.^2+aly.^2+alz.^2).^0.5;
            
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
            absw=enstrophy.^0.5;
            s11=0.5*(ux+ux);
            s12=0.5*(uy+vx);
            s13=0.5*(uz+wx);
            s22=0.5*(vy+vy);
            s23=0.5*(vz+wy);
            s33=0.5*(wz+wz);
            
            wS1=w1.*s11+w2.*s12+w3.*s13;
            wS2=w1.*s12+w2.*s22+w3.*s23;
            wS3=w1.*s13+w2.*s23+w3.*s33;
            absws=(wS1.^2+wS2.^2+wS3.^2).^0.5;
            wws=w1.*wS1+w2.*wS2+w3.*wS3;
            leng=1;
            
            acx=u.*ux+v.*uy+w.*uz;%%%%%%%%%%%%%%%%%%%%
            acy=u.*vx+v.*vy+w.*vz;%%%%%%%%%%%%%%%%%%%%
            acz=u.*wx+v.*wy+w.*wz;%%%%%%%%%%%%%%%%%%%%
            absac=(acx.^2+acy.^2+acz.^2).^0.5;%%%%%%%%%%%%%%%%%%%%%%%%%
            aLamx=w2.*w-w3.*v;%%%%%%%%%%%%%%%%%%%%%%%
            aLamy=w3.*u-w1.*w;%%%%%%%%%%%%%%%%%%%%%%%
            aLamz=w1.*v-w2.*u;%%%%%%%%%%%%%%%%%%%%%%%%%
            absaLam=(aLamx.^2+aLamy.^2+aLamz.^2).^0.5;%%%%%%%%%%%%%%%%%
%             aBx=0.5*f(:,45);%acx-aLamx;%%%%%%%%%%%%%%%%%%%%%
%             aBy=0.5*f(:,46);%acy-aLamy;%%%%%%%%%%%%%%%%%%%%%%%
%             aBz=0.5*f(:,47);%acz-aLamz;%%%%%%%%%%%%%%%%%%%%%%%%%


            aBx=acx-aLamx;%
        aBy=acy-aLamy;%
        aBz=acz-aLamz;
        
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
            
            axx=f(:,22);
            axy=f(:,23);
            axz=f(:,24);
            ayx=f(:,25);
            ayy=f(:,26);
            ayz=f(:,27);
            azx=f(:,28);
            azy=f(:,29);
            azz=f(:,30); 
            
            aw1=azy-ayz;
            aw2=axz-azx;
            aw3=ayx-axy;
            aenstrophy=aw1.^2+aw2.^2+aw3.^2;
            absaw=aenstrophy.^0.5;
            as11=0.5*(axx+axx);
            as12=0.5*(axy+ayx);
            as13=0.5*(axz+azx);
            as22=0.5*(ayy+ayy);
            as23=0.5*(ayz+azy);
            as33=0.5*(azz+azz);
            
            diva=as11+as22+as33;
            
            reldiv=abs(ux+vy+wz)./(abs(ux)+abs(vy)+abs(wz));
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
                while ende<s(1,1)-1 && searching==1
                    if age(ende+1)==age(ende)+1
                        ende=ende+1;
                        leng=leng+1;
                    else
                        searching=0;
                    end
                end
                if ende<s(1,1)
                    for ii=start:ende
                        if age(ii)>tail && age(ii)<leng-tail+1 && reldiv(ii)<0.1 
                            
                            
                            A=[s11(ii) s12(ii) s13(ii);s12(ii) s22(ii) s23(ii);s13(ii) s23(ii) s33(ii)];
                            [V,D] = eig(A);
                            [D,k] = sort(diag(D)); 	% ascending order
                            D = diag(D(end:-1:1)); 			% descending order		
                            V = V(:,k(end:-1:1)); 	% the same order for eigenvectors
                            strain=D(1,1)^2+D(2,2)^2+D(3,3)^2;
                            Q=0.25*(enstrophy(ii)-2*strain);
                            
                            aA=[as11(ii) as12(ii) as13(ii);...
                                as12(ii) as22(ii) as23(ii);...
                                as13(ii) as23(ii) as33(ii)];
                            [aV,aD] = eig(aA);
                            [aD,k] = sort(diag(aD)); 	% ascending order
                            aD = diag(aD(end:-1:1)); 			% descending order		
                            aV = aV(:,k(end:-1:1)); 	% the same order for eigenvectors
                            astrain=aD(1,1)^2+aD(2,2)^2+aD(3,3)^2;
                            aQ=0.25*(aenstrophy(ii)-2*astrain);
                            
                            % if aenstrophy(ii)<0.1*10e4 && 2*astrain<0.1*10e4
                                counter=counter+1;
                                statistic(counter,1)=astrain;
                                statistic(counter,2)=aenstrophy(ii);
                                statistic(counter,3)=aQ;
                            % end
                        end
                    end
                end
            end
        end
    end
    
    %save asaas;
else
    %load resultOne;
end

figure;
hold on;
grid on;
box on;
[N,X] = nhist(statistic(1:counter,1),10000); %%%
plot(X,N,'r');
[N,X] = nhist(statistic(1:counter,2),10000); %%%
plot(X,N,'b');
title('red=a_{strain}, blue=a_{enstrophy}')

figure;
hold on;
grid on;
box on;
[N,X] = nhist(statistic(1:counter,3),1000); %%%
plot(X,N,'r');
xlabel('Q_a')
ylabel('pdf')
title('red=a_{strain}, blue=a_{enstrophy}')