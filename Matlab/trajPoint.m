%function res=trajPoint(first,last,minLength)

numTraj=0;
leng=zeros(1000000,1);
first=1;
last=800;
minLength=5;

for i=first:last
    if mod(i,5)==0
        i
    end

    %name=['E:\PTV\Working_folder\WF_1\res_full\trajPoint.',num2Str(i)];
    name=['E:\PTV\Working_folder\Exp21b_1211201\res\trajPoint.',num2Str(i)];
    f=load(name);
    s=size(f);
    clear x y z;

    if s(1,1)>0
        x=f(:,1);
        y=f(:,2);
        z=f(:,3);

        u=f(:,4);
        v=f(:,5);
        w=f(:,6);

        age=f(:,10);


        ende=0;

        for ii=1:s(1,1)
            begi=ende+1;
            ende=begi; % propagating variable
            found=0;
            while found==0 && ende<s(1,1)-1
                if age(ende+1)>age(ende)
                    ende=ende+1;
                else
                    found=1;
                    numTraj=numTraj+1;
                    if numTraj==1
                        figure;hold on;box on;grid on;

                    end
                    %plot
                    leng(numTraj)=ende-begi+1;
                    %[s numTraj leng(numTraj)]

                    if leng(numTraj)>minLength-1
                        %x_i=x(begi+2);
                        %x_e=x(ende);
                        y_i=y(begi+2);
                        y_e=y(ende);
                        %if y_i<.0225 %&& y_e>.0265
                        %if  y_i<0.028 & x_i>9e-3 & x_i<9.5e-31 & y_e<0.028 & x_e>9e-3 & x_e<9.5e-31  % y_i<.022 && y_e>.0265
                        %if y_i>.0225 & y_i<0.03
                            pp=plot3(x(begi+2:ende),y(begi+2:ende),z(begi+2:ende));
                            set(pp,'linewidth',0.1)
                            scatter3(x(begi+2),y(begi+2),z(begi+2),10,'g','filled')
                            scatter3(x(ende),y(ende),z(ende),10,'r','filled')
                            %drawnow
                        %else
                            %if y_i>.03
                                pp=plot3(x(begi+2:ende),y(begi+2:ende),z(begi+2:ende),'k');
                                set(pp,'linewidth',2)
                                scatter3(x(begi+2),y(begi+2),z(begi+2),50,'g','filled')
                                scatter3(x(ende),y(ende),z(ende),10,'r','filled')
                            %end
                        %end
                        %end
                    end
                end
            end
        end
    end


    % set(gca,'xlim',[0 0.018],'ylim',[ 0.015 0.045],'view',[20 60])
    %
    % h=gcf;
    % B(i)=getframe(h);
    %
end

% figure;
% movie(B)
%draw_channel(1)

% cL=[1 0 1];
%
%
%
% %left_face
% x5=[0 0 0 0];
% y5=[0 0.0245 0.0245 0];
% z5=[0 0 0.012 0.012];
% patch(x5,y5,z5,cL)
%
%
% %Right face
% x6=[0.018 0.018 0.018 0.018];
% y6=[0 0.0245 0.0245 0];
% z6=[0 0 0.012 0.012];
% patch(x6,y6,z6,cL)
%
%
%
% % Bottom
% x7=[0 0 0.018 0.018];
% y7=[0 0.0245 0.0245 0];
% z7=[0 0 0 0];
% patch(x7,y7,z7,cL)
%
%
%
%
% % Front face
% x9=[0 0 0.018 0.018];
% y9=[0 0 0 0];
% z9=[0 0.012 0.012 0];
% patch(x9,y9,z9,cL)
%
%
% % Left face channel
% x10=[0.0075 0.0075 0.0075 0.0075 0.0105 0.0105];
% y10=[0.0245 0.0245 0.045 0.045     0.045 0.0245];
% z10=[0.009 0.012 0.012 0.009      0.009 0.009];
% patch(x10,y10,z10,cL)
% %plot3([0.0075 0.0075 ],[0.045 0.0245 ],[0.009 0.009 ],'k')
%
% % Edge
% x11=[0 0 0.0075 0.0075 0.0105 0.0105 0.018 0.018 ];
% y11=[0.0245 0.0245 0.0245 0.0245 0.02450 0.0245 0.0245 0.0245 ];
% z11=[0 0.012 0.012 0.009 0.009 0.012 0.012 0];
% patch(x11,y11,z11,cL)
%
%
% % Right face channel
% x12=[ 0.0105 0.0105 0.0105 0.0105];
% y12=[ 0.045 0.045 0.0245 0.0245];
% z12=[0.009 0.012 0.012 0.009];
% patch(x12,y12,z12,cL)
% alpha(0.05)
%
% xlabel('x','FontSize',12,'FontName','Times New Roman');
% ylabel('y','FontSize',12,'FontName','Times New Roman');
% zlabel('z','FontSize',12,'FontName','Times New Roman');
%tit=['time: ',num2Str(floor(10*first/20+0.5)/10),' sec'];
%title(tit);
%daspect([1 1 0.3])
% axis([-0.0003 0.22 0.00 0.445 0 0.018])

figure;
hold on;
box on;
[nout,xout]=nhist(leng(1:numTraj),50);
plot(xout,nout);
figure;
bar(xout,nout)
h=['#'];
xlabel(h,'FontSize',12,'FontName','Times New Roman')
ylabel('pdf','FontSize',12,'FontName','Times New Roman')


xlim([-0.003 0.022 ])
ylim([0 0.045])
zlim([0.00 0.012])


view([50,60])
lighting phong
material shiny

% h = camlight('left');
% for i = 1:20;
% 	camorbit(10,0)
% 	camlight(h,'left')
% 	drawnow;
% end

