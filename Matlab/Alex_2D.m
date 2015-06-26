% 2D function
% Initialization

directoryName = 'Alex/res_06_17_08_1'; %update the file name we run
d=0.0128; % lid velocity mm/s
scene='num 1'; % number of the scene
dx =-0.045:0.005:0.025;
dy=-0.045:0.005:0.04;

% Read the xuap files.
tmp = readXUAPFiles(directoryName);

%creating array that contain data on the middle section of x and y of the
%cavity the length of the section taken is 8 mm

clear x y 
% vel = struct('u',zeros(10,3),'v',zeros(10,3));
for i=1:length(tmp)
    ind =find(abs(tmp(i).uf) > 10e-5 );
    vel(i).u=zeros(11,3);
    vel(i).v=zeros(11,3);
    for j=1:length(ind)
        k=ind(j);
        if (tmp(i).yf(k)>=0.016 && tmp(i).yf(k)<=0.024)
            if (tmp(i).xf(k)>=-0.02 && tmp(i).xf(k)<=0.05)
                x=(tmp(i).xf(k)+0.02)/0.007;
                vel(i).v(floor(x)+1,1)=tmp(i).vf(k)+vel(i).v(floor(x)+1,1);
                vel(i).v(floor(x)+1,3)=vel(i).v(floor(x)+1,3)+1;
            end
        end
        if (tmp(i).xf(k)>=0.011 && tmp(i).xf(k)<=0.019)
            if (tmp(i).yf(k)>=-0.02 && tmp(i).yf(k)<=0.055)
                y=(tmp(i).yf(k)+0.02)/0.0075;
                vel(i).u(floor(y)+1,1)=tmp(i).uf(k)+vel(i).u(floor(y)+1,1);
                vel(i).u(floor(y)+1,3)=vel(i).u(floor(y)+1,3)+1;
            end
        end
    end
end
% dividing by number of particles on the same place at the same time
for i=1:length(vel)
    for j=1:10
        if ~isnan(vel(i).u(j,1))
            vel(i).u(j,1)=vel(i).u(j,1)/vel(i).u(j,3);
        end
        if ~isnan(vel(i).v(j,1))
            vel(i).v(j,1)=vel(i).v(j,1)/vel(i).v(j,3);
        end
    end
end
%averaging velocity calculations

u_av=zeros(10,3);
v_av=zeros(10,3);

for i=1:length(vel)
    for j=1:10
        if ~isnan(vel(i).u(j,1))
            u_av(j,1)=u_av(j,1)+vel(i).u(j,1);
            u_av(j,3)=u_av(j,3)+1;
        end
        if ~isnan(vel(i).v(j,1))
            v_av(j,1)=v_av(j,1)+vel(i).v(j,1);
            v_av(j,3)=v_av(j,3)+1;
        end
    end
end
% averaging
for i=1:10
    u_av(i,1)=u_av(i,1)/u_av(i,3);
    v_av(i,1)=v_av(i,1)/v_av(i,3);
end

% normalization
u_av1=[u_av(:,1)/d];
v_av1=[v_av(:,1)/d];
u_av(:,3)=0;
v_av(:,3)=0;

%adding the walls and lid
u_av1=[0;u_av1;0.5];
v_av1=[0;v_av1;0];
%putting the data in the middle of the graph 
u_av2=u_av1+0.5;
v_av2=v_av1+0.5;

%plot data

h=figure;
plot(u_av2,0:1/(length(u_av2)-1):1,'rs','MarkerFaceColor','g')
hold on
plot(0:1/(length(v_av2)-1):1,v_av2,'rs','MarkerFaceColor','b')
title( [scene,'   av'] )
fileName = [directoryName(6:end),'_av.fig'];
saveas(h,fileName,'fig')
close(h)
clear h;
%rms calculations

for i=1:length(vel)
    for j=1:10
        if ~isnan(vel(i).u(j,1))
            u_av(j,2)=((vel(i).u(j,1)/d)-u_av(j,1))^2+u_av(j,2);
            u_av(j,3)=u_av(j,3)+1;
        end
        if ~isnan(vel(i).v(j,1))
            v_av(j,2)=((vel(i).v(j,1)/d)-v_av(j,1))^2+v_av(j,2);
            v_av(j,3)=v_av(j,3)+1;
        end
    end
end
% normalization
for i=1:10
    u_av(i,2)=u_av(i,2)/u_av(i,3);
    v_av(i,2)=v_av(i,2)/v_av(i,3);
end
%adding the walls and lid
u_rms1=[0;u_av(:,2);0];
v_rms1=[0;v_av(:,2);0];
% scaling the rms data
u_rms1=u_rms1.^0.5;
v_rms1=v_rms1.^0.5;
u_rms1=u_rms1*1;
v_rms1=v_rms1*1;
% bringing the data into the middle
u_rms2=u_rms1+0.5;
v_rms2=v_rms1+0.5;

%plot data

h=figure;
plot(u_rms2,0:1/(length(u_rms2)-1):1,'rs','MarkerFaceColor','g');
hold on
plot(0:1/(length(v_rms2)-1):1,v_rms2,'rs','MarkerFaceColor','b');
title( [scene,'   rms'] )
fileName = [directoryName(6:end),'_rms.fig'];
saveas(h,fileName,'fig')
close(h)
clear h;

% % plot of the cavity
% x=[];
% y=[];
% for i=1:length(tmp)  %the first file is zero file so we start at the second one
%     ind= find(abs(tmp(i).uf) > 10e-5); 
%     x=[x,(tmp(i).xf(ind))']; % joining the x arrays into one in order to craete pdf
%     y=[y,(tmp(i).yf(ind))'];
% end
% plot(x,y,'rs')



