clear all
clear all
close all

tol=1.00e-3;

long_traj_count=0;
max_length=500;
tau_eta=0.14;
frame_rate=250;
factor=tau_eta*frame_rate;
WF=6;

direc=['D:\Turbulent_data_backup\flow_induced_aggregates_breakage_at_100rpm\WF',num2str(WF),'\res_agg'];
d = [direc,'/rt_is_v2*'];
dd=dir(d);
filename1=dd(1).name;
first = eval(filename1(10:end));

if first==1
    last = size(dd,1);
else
    filename2=dd(end).name;
    last=eval(filename2(10:end));
end

%   first=300;
%   last=400;
agg=repmat(zeros,[last-first+1,200,9]);
load_data=0;
if load_data==0
    name_root=['D:\Turbulent_data_backup\flow_induced_aggregates_breakage_at_100rpm\WF',num2str(WF),'\res_agg'];
    for i=first:last
        
        name=[name_root,'\rt_is_v2.',num2Str(i)];
        fid = fopen(name, 'r');
        num_points = fscanf(fid, '%i', [1 1]);  % It has two rows now
        tmp = fscanf(fid, '%i %f %f %f %i %i %i %i %f %f %f %f ', [12 num_points]);
        A=tmp';
        A(:,2:4)=A(:,2:4)*0.001;
        agg(i-first+1,1:num_points,1:7)=A(:,[2 3 4 9 10 11 12]);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        agg(i-first+1,1:num_points,8)=i;
        fclose(fid);
    end
    
    xa=[];
    ya=[];
    za=[];
    ta=[];
    for i=1:last-first+1
        
        xa=[xa;agg(i,:,1)'];
        ya=[ya;agg(i,:,2)'];
        za=[za;agg(i,:,3)'];
        t=repmat(i,length(agg(i,:,1)),1);
        ta=[ta;t];
        
    end
    figure;
    scatter3(xa,ya,za,30,'b')
    title('Trajectories are not linked in time')
    
    
    tid=0;
    for i=1:last-first
       
        ind=find(abs(agg(i,:,1))>0 & abs(agg(i,:,2))>0);
        %ind=find(agg(i,:,4)>0);
        
        ind2=find(agg(i,ind,9)==0);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ind=ind(ind2);
        num_points=length(ind);
        for j=1:num_points
            tid=tid+1;
            pid=1;
            traj(tid,pid,1)=i;%frame id
            traj(tid,pid,2)=ind(j);%pos in frame
            agg(traj(tid,pid,1),traj(tid,pid,2),9)=1;%occupied as of now
            %find next ?
            end_of_traj=0;
            while end_of_traj==0
                pos_agg=[agg(traj(tid,pid,1),traj(tid,pid,2),1) agg(traj(tid,pid,1),traj(tid,pid,2),2) agg(traj(tid,pid,1),traj(tid,pid,2),3)];
                ind2=find(abs(agg(traj(tid,pid,1)+1,:,1))>0 & abs(agg(traj(tid,pid,1)+1,:,2))>0 & agg(traj(tid,pid,1)+1,:,9)<1);
                %ind2=find(agg(traj(tid,pid,1)+1,:,4)>0  & agg(traj(tid,pid,1)+1,:,9)<1);
                pos=repmat(pos_agg,length(ind2),1);
                if length(ind2>0)%%%%%
                    if length(ind2)>1
                        delta=squeeze(agg(traj(tid,pid,1)+1,ind2,1:3))-pos;
                        dist=sum(delta.^2,2).^0.5;
                    else
                        delta=squeeze(agg(traj(tid,pid,1)+1,ind2,1:3))-pos';
                        dist=sum(delta.^2,1).^0.5;
                    end
                    ind3=find(dist<tol);
                    if length(ind3)>0
                        ind2=ind2(ind3);
                    else
                        ind2=[];
                    end
                else
                    ind2=[];
                end
                
                if length(ind2)>0 & traj(tid,pid,1)<last-first 
                    pid=pid+1;
                    traj(tid,pid,1)=traj(tid,pid-1,1)+1;
                    traj(tid,pid,2)=ind2(1);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    agg(traj(tid,pid,1),traj(tid,pid,2),9)=1;%occupied as of now;
                else
                    end_of_traj=1;
                    
                    if pid>20
                        
                        long_traj_count=long_traj_count+1;
                        
                        traj_len(tid)=pid;
                        
                        if pid>max_length 
                            max_length=pid;
                        end
                        
                        for k=1:pid
                            tx(tid,k)=agg(traj(tid,k,1),traj(tid,k,2),1);
                            ty(tid,k)=agg(traj(tid,k,1),traj(tid,k,2),2);
                            tz(tid,k)=agg(traj(tid,k,1),traj(tid,k,2),3);
                            totalpix(tid,k)=agg(traj(tid,k,1),traj(tid,k,2),4);
                            xpix(tid,k)=agg(traj(tid,k,1),traj(tid,k,2),5);
                            ypix(tid,k)=agg(traj(tid,k,1),traj(tid,k,2),6);
                            grv(tid,k)=agg(traj(tid,k,1),traj(tid,k,2),7);
                            time(tid,k)=agg(traj(tid,k,1),traj(tid,k,2),8);
                     
                        end
                        
                    else
                        
                        tid=tid-1;
                    end
                end
            end
        end
    end
    long_traj_count=long_traj_count
    
    save G:\flow_induced_aggregates_breakage_at_100rpm\WF6\abc xa ya za ta tx ty tz totalpix xpix ypix grv time tol tid traj traj_len
    
else
    load G:\flow_induced_aggregates_breakage_at_100rpm\WF6\abc
end



figure;
for traj_id=1:tid
    t=time(traj_id,:);
    t=t(t>0);
    x=tx(traj_id,:);
    x=x(1:numel(t));%%% This is a bug remover !!
    y=ty(traj_id,:);
    y=y(1:numel(t));
    z=tz(traj_id,:);
    z=z(1:numel(t));
    
    hold on
    scatter3(x,y,z,'k.')
    plot3(x,y,z,'k')
    scatter3(x(1),y(1),z(1),50,'g','filled')
    text(x(1),y(1),z(1),[' T (',num2str(traj_id),')', '   ', 'f',num2str(t(1))]);%num2str(traj(traj_id,1))
    scatter3(x(end),y(end),z(end),50,'r','filled')
    text(x(end),y(end),z(end),[' f',num2str(t(end))]);
    tit=['tol=',num2str(tol*1000),'mm','    ', 'Num\_traj=',num2str(tid)];
    title(tit)
    hold off
end












