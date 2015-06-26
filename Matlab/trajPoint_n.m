%function res=trajPoint_n(first,last,minLength)
numTraj=0;
leng=zeros(1000000,1);

first=19700;
last=19800;
minLength=20;

for i=first:last
    if mod(i,5)==0 
        i          
    end
 
    name=['E:\PTV\Working_folder\WF_1\res\trajPoint.',num2Str(i)];
    f=load(name);
    s=size(f);
    clear x y z;
    name1=['E:\PTV\Working_folder\WF_1\3d2d_info_diagnosis\rt_is2dTrajPoint_2_trajPoint3d2d_',num2Str(i),'.mat'];
    f1=load(name1);
    s1=size(f1);
    
   
if s(1,1)>0 
        x=f(:,1);
        y=f(:,2);
        z=f(:,3);
        
        age=f(:,10);
        
        n_cam1=f1.a1_t;
        

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
                        figure;hold on;box on;grid on;
                    end
                    %plot
                    leng(numTraj)=ende-begi+1;
                    %[s numTraj leng(numTraj)]
                    if leng(numTraj)>minLength-1
                        
                            
                        scatter3(x(begi+2:ende),y(begi+2:ende),z(begi+2:ende),n_cam1(begi+2:ende),'b');
                        scatter3(x(begi+2),y(begi+2),z(begi+2),n_cam1(begi+2),'g','filled')
                        scatter3(x(ende),y(ende),z(ende),n_cam1(ende),'r','filled')
                        
                    end
                end
            end
        end
    end
end


xlabel('x','FontSize',12,'FontName','Times New Roman');
ylabel('y','FontSize',12,'FontName','Times New Roman');
zlabel('z','FontSize',12,'FontName','Times New Roman');
%tit=['time: ',num2Str(floor(10*first/20+0.5)/10),' sec'];
%title(tit);
%daspect([1 1 0.3])
% axis([-0.0003 0.22 0.00 0.445 0 0.018])

% figure;
% hold on;
% box on;
% [nout,xout]=nhist(leng(1:numTraj),50);
% plot(xout,nout);
% h=['#'];
% xlabel(h,'FontSize',12,'FontName','Times New Roman')
% ylabel('pdf','FontSize',12,'FontName','Times New Roman')


% xlim([-0.003 0.022 ])
% ylim([0 0.045])
% zlim([0.00 0.012])
