function res=trajPoint(first,last,minLength)


numTraj=0;
length=zeros(1000000,1);
for i=first:last

    if mod(i,50)==0
        i
    end
    %name=['X:/rotating/PTV08/080725/run1/g_trajPoint.',num2Str(i)];
    %name=['X:/rotating/PTV08/080811/run1/g_trajPoint.',num2Str(i)];
    name=['X:/rotating/PTV08/080827/run2/g_trajPoint.',num2Str(i)];
    f=load(name);
    s=size(f);
    if s(1,1)>0
        x=f(:,1);
        y=f(:,2);
        z=f(:,3);
        ax=f(:,7);
        ay=f(:,8);
        az=f(:,9);
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
        enstro=w1.^2+w2.^2+w3.^2;
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
        age=f(:,32);

        ende=0;
        for ii=1:s(1,1)
            begi=ende+1;
            ende=begi;
            found=0;
            while found==0 && ende<s(1,1)-1
                if age(ende+1)>age(ende)
                    ende=ende+1;
                else
                    found=1;
                    numTraj=numTraj+1;
                    if numTraj==1
                        %figure(1);hold on;box on;grid on;
                        figure(2);hold on;box on;grid on;
                        %figure(3);hold on;box on;grid on;
                    end
                    %plot
                    length(numTraj)=ende-begi+1;
                    if length(numTraj)>minLength-1
                        %figure(1)
                        %plot3(x(begi:ende),y(begi:ende),z(begi:ende),'b');
                        if mean(absa(begi:ende)<0.02) & mean(enstro(begi:ende)<1)
                            %figure(2)
                            plot3(x(begi:ende),z(begi:ende),y(begi:ende),'b');
                        else
                            %figure(3)
                            plot3(x(begi:ende),z(begi:ende),y(begi:ende),'r');
                        end
                    end
                end
            end
        end
    end
end
xlabel('x','FontSize',12,'FontName','Times New Roman');
ylabel('z','FontSize',12,'FontName','Times New Roman');
zlabel('y','FontSize',12,'FontName','Times New Roman');

figure;
hold on;
box on;
[nout,xout]=nhist(length(1:numTraj),50);
plot(xout,nout);
h=['#'];
xlabel(h,'FontSize',12,'FontName','Times New Roman')
ylabel('pdf','FontSize',12,'FontName','Times New Roman');
