function res=trajectory_from_ptv_is(first,last)
figure;
% first=1;
% lasz=5457;
near_pair=0.001;
frac_of_dist=5;%0.4;%1.2;%0.4
min_traj_length=5;% Beat==5

ptv=zeros(last-first+1,2000,6);
ptv(:,:,1)=-1;
ptv(:,:,2)=-1;
num_traj=0;
traj_bin=[];
num_breakages=0;
breakage=[];
traj_len=[];%% Added by me to see the pdf of trajectory length

name_root=['D:\Flow\WF1\'];
for i=first:last
    
    name=[name_root,'res\ptv_is.',num2Str(i)];
    
    fid = fopen(name, 'r');
    num_points = fscanf(fid, '%i', [1 1]);    % It has two rows now.
    tmp = fscanf(fid, '%i %i %f %f %f  ', [5 num_points]);
    A=tmp';
    A(:,3:5)=A(:,3:5)*0.001;
    ptv(i-first+1,1:num_points,1:5)=A;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fclose(fid);
end

go=0;
for i=1:last-first+1
    [i last-first+1]
    j=1;
    while j<1000
        %find trajectories
        if ptv(i,j,6)==0
            go=1;
            traj=[0 0];
            traj(1,1)=i;
            traj(1,2)=j;
            num_in_traj=1;
        end
        while go==1
            traj=find_next(ptv,traj,num_in_traj);
            old_num_in_traj=num_in_traj;
            num_in_traj=traj(length(traj),1);
            traj=traj(1:length(traj)-1,:);
            if old_num_in_traj==num_in_traj | num_in_traj>last-first-i+1 | traj(1:length(traj),1)>last-first-1
                go=0;
                %plot traj or store it
                if length(traj)>5
                    
                    
                    num_traj=num_traj+1;
                    traj_len=[traj_len;num_in_traj];%% added by me to see the pdf of trajectory length
                     
                    
                    %%%here we can filter the positions
                    ptv=filter_traj(ptv,traj);
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    traj_bin=plot_traj(traj_bin,ptv,traj,first+i-1);
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    ptv=update(ptv,traj);
                    
                    
                end
            end
        end
        j=j+1;
    end
end


%%pdf of trajectory length
figure(100);
[N,X]=nhist(traj_len);
plot(X,N)

%figure(200);box on;hold on;
xlabel('x')
ylabel('y')
zlabel('z')
view([50,60])
box on;grid on;

function traj=find_next(ptv,traj,num_in_traj)

i=traj(num_in_traj,1);
j=traj(num_in_traj,2);
if ptv(i,j,2)>-1
    num_in_traj=num_in_traj+1;
    traj(num_in_traj,1)=i+1;
    traj(num_in_traj,2)=ptv(i,j,2)+1;
end
traj=[traj;[num_in_traj 0]];

function traj_bin=plot_traj(traj_bin,ptv,traj,time)


    
    hold on;

for i=1:length(traj)
    x(i)=ptv(traj(i,1),traj(i,2),3);
    y(i)=ptv(traj(i,1),traj(i,2),4);
    z(i)=ptv(traj(i,1),traj(i,2),5);
    %totpix(i)=ptv(traj(i,1),traj(i,2),6);
    id(i)=i;
    frame_i(i)=traj(i,1);
    j(i)=traj(i,2);
end

%%%% Geometry condition %%%%%%


 %for i=1:length(traj)
     %if y(i)>0.03 | x(i)<9e-3 | x(i)>9.4e-3 | z(i)<0.01
      
      
%         %if y(i)<0.0224 | y(i)>0.0229 | x(i)<9.35e-3 | x(i)>9.7e-3
%         %if y(i)<0.0212 | y(i)>0.0224 | x(i)<2.8e-3 | x(i)>4.2e-3
%         %if y(i)<0.026
%         ok=0;
     %end
 %end

    
        plot3(x,y,z,'b','LineWidth',1);
        %scatter3(x,y,z)%5,totpix)
        
       
     traj_bin=[traj_bin; [frame_i' j' id' x' y' z']];%%%%%%%%%%%%%%%%%%%%%
%     figure(300);hold on;
%     scatter3(time+id,y,totpix)


function ptv=update(ptv,traj)

for i=1:length(traj)
    ptv(traj(i,1),traj(i,2),6)=1;
end

function ptv=filter_traj(ptv,traj)

for i=1:length(traj)
    if i==1
        nx(i)=0.5*ptv(traj(i,1),traj(i,2),3)+0.5*ptv(traj(i+1,1),traj(i+1,2),3);
        ny(i)=0.5*ptv(traj(i,1),traj(i,2),4)+0.5*ptv(traj(i+1,1),traj(i+1,2),4);
        nz(i)=0.5*ptv(traj(i,1),traj(i,2),5)+0.5*ptv(traj(i+1,1),traj(i+1,2),5);
    end
    if i==length(traj)
        nx(i)=0.5*ptv(traj(i,1),traj(i,2),3)+0.5*ptv(traj(i-1,1),traj(i-1,2),3);
        ny(i)=0.5*ptv(traj(i,1),traj(i,2),4)+0.5*ptv(traj(i-1,1),traj(i-1,2),4);
        nz(i)=0.5*ptv(traj(i,1),traj(i,2),5)+0.5*ptv(traj(i-1,1),traj(i-1,2),5);
    end
    if i>1 & i<length(traj)
        order=min(i-1,length(traj)-i);
        if order>4
            order=4;
        end
        switch order
            case 1
                weight=[2 1];
            case 2
                weight=[6 4 1];
            case 3
                weight=[20 15 6 1];
            case 4
                weight=[70 56 28 8 1];
        end
        su=0;
        fx=0;fy=0;fz=0;
        for j=i-order:i+order
            su=su+weight(abs(i-j)+1);
            fx=fx+ptv(traj(j,1),traj(j,2),3)*weight(abs(i-j)+1);
            fy=fy+ptv(traj(j,1),traj(j,2),4)*weight(abs(i-j)+1);
            fz=fz+ptv(traj(j,1),traj(j,2),5)*weight(abs(i-j)+1);
        end
        nx(i)=fx/su;
        ny(i)=fy/su;
        nz(i)=fz/su;
    end
end

for i=1:length(traj)
    ptv(traj(i,1),traj(i,2),3)=nx(i);
    ptv(traj(i,1),traj(i,2),4)=ny(i);
    ptv(traj(i,1),traj(i,2),5)=nz(i);
end






