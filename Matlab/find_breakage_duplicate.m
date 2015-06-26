function breakage=find_breakage(first,last,plot_mode)
%find_breakage(21030,21055)
first=21030;
last=21055;
plot_mode=0;


near_pair=0.001;
frac_of_dist=0.4;
min_traj_length=5;

ptv=zeros(last-first+1,1000,6);
ptv(:,:,1)=-1;
ptv(:,:,2)=-1;
num_traj=0;
traj_bin=[];
num_breakages=0;
breakage=[];


for i=first:last
    name=['E:\PTV\Working_folder\WF_1\res\ptv_is.',num2Str(i)];
    fid = fopen(name, 'r');
    num_points = fscanf(fid, '%i', [1 1]);    % It has two rows now.
    tmp = fscanf(fid, '%i %i %f %f %f', [5 num_points]);
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
                if length(traj)>2
                    num_traj=num_traj+1;
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    traj_bin=plot_traj(traj_bin,ptv,traj,100,first+i-1,plot_mode);
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    ptv=update(ptv,traj);
                end
            end
        end
        j=j+1;
    end
end
ind_beg=find(traj_bin(:,3)==1);
si=size(ind_beg);
s=si(1,1);
ind_end=ind_beg-1;
ind_end=ind_end(2:s);
ind_end(s)=length(traj_bin(:,3));
si=size(ind_end);
if si(1,2)>1
    ind_end=ind_end';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
traj_ind=[ind_beg ind_end ind_end-ind_beg+1];%[begin end length]
si=size(traj_ind);
le_traj_ind=si(1,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for i=1:le_traj_ind
    if mod(i,50)==0
        [i le_traj_ind]
    end
    %if beg close to something backtrapolate and on long enough trajectory
    if traj_ind(i,3)>=min_traj_length
        ind_beg=traj_ind(i,1);
        ind_i=traj_bin(ind_beg,1);
        ind_j=traj_bin(ind_beg,2);
        xb=ptv(ind_i,ind_j,3);
        yb=ptv(ind_i,ind_j,4);
        zb=ptv(ind_i,ind_j,5);
        %scatter3(xb,yb,zb,'r')
        ind=find(abs(ptv(ind_i,:,3))>0 & abs(ptv(ind_i,:,4))>0 & abs(ptv(ind_i,:,5)>0));
        dx=xb-ptv(ind_i,ind,3);
        dy=yb-ptv(ind_i,ind,4);
        dz=zb-ptv(ind_i,ind,5);
        dist=(dx.^2+dy.^2+dz.^2).^0.5;
        ind_dist=find(dist<near_pair & dist>0);
        near=frac_of_dist*dist(ind_dist);
        if length(ind_dist)>0
            ind_j2=ind(ind_dist);
            %check if ind_i,ind_j2 is element of long enough traj_bin
            for j=1:length(traj_bin)
                if traj_bin(j,1)==ind_i & traj_bin(j,2)==ind_j2
                    ind_beg2=j;
                    for k=1:le_traj_ind
                        if ind_beg2>=traj_ind(k,1) & ind_beg2<=traj_ind(k,2) & traj_ind(k,3)>=min_traj_length
                            i2=k;
                            %scatter3(ptv(ind_i,ind_j ,3),ptv(ind_i,ind_j ,4),ptv(ind_i,ind_j ,5),'c')
                            %scatter3(ptv(ind_i,ind_j2,3),ptv(ind_i,ind_j2,4),ptv(ind_i,ind_j2,5),'m')
                            %if backtrapolate close to something at x
                            ok=1;
                            const=backtrapolate(i ,ind_beg ,traj_bin,traj_ind);
                            if sum(sum(abs(const(:,1:2))))==0
                                ok=0;
                            end
                            [bx by bz]=draw_backtr(const,3,100,'r');
                            const=backtrapolate(i2,ind_beg2,traj_bin,traj_ind);
                            if sum(sum(abs(const(:,1:2))))==0
                                ok=0;
                            end
                            if ok==1
                                [bx2 by2 bz2]=draw_backtr(const,3,100,'k');
                                bxm=0.5*(bx+bx2);
                                bym=0.5*(by+by2);
                                bzm=0.5*(bz+bz2);
                                %scatter3(bxm,bym,bzm,'g');
                                %check if any of the bxm points are close to any
                                found=0;
                                for n=2:6 %loop through bxm points
                                    if found==0 & ind_i-(n-1)>=1
                                        ind=find(abs(ptv(ind_i-(n-1),:,3))>0 & abs(ptv(ind_i-(n-1),:,4))>0 & abs(ptv(ind_i-(n-1),:,5)>0));
                                        dx=bxm(n)-ptv(ind_i-(n-1),ind,3);
                                        dy=bym(n)-ptv(ind_i-(n-1),ind,4);
                                        dz=bzm(n)-ptv(ind_i-(n-1),ind,5);
                                        dist=(dx.^2+dy.^2+dz.^2).^0.5;
                                        ind_dist=find(dist<near & dist>0);
                                        if length(ind_dist)>0
                                            ind_i3=ind_i-(n-1);
                                            ind_j3=ind(ind_dist);
                                            %check if ind_i3,ind_j3 is element of traj_bin
                                            for m=1:length(traj_bin)
                                                if traj_bin(m,1)==ind_i3 & traj_bin(m,2)==ind_j3
                                                    ind_end3=m;
                                                    for o=1:le_traj_ind
                                                        if ind_end3>=traj_ind(o,1) & ind_end3<=traj_ind(o,2) & traj_ind(o,3)>=min_traj_length
                                                            i3=o;
                                                            %extrapolate from x if extrapolate close to midpoint of pair we have breakage
                                                            const=extrapolate(i3,ind_end3,traj_bin,traj_ind);
                                                            [ex ey ez]=draw_extr(const,3,100);
                                                            dx=bxm(1)-ex(n-1);
                                                            dy=bym(1)-ey(n-1);
                                                            dz=bzm(1)-ez(n-1);
                                                            dist=(dx.^2+dy.^2+dz.^2).^0.5;
                                                            yes=0;                                                            
                                                            if dist<near & not(i==i2) & not(i==i3) & not(i2==i3)
                                                                %check if begin on traj 2 is not further away than gap
                                                                time_gap=ind_i-ind_i3;
                                                                if ind_beg2-traj_ind(i2,1)<=time_gap
                                                                    pos_ok=1;
                                                                else
                                                                    pos_ok=0;
                                                                end
                                                                if pos_ok==1
                                                                    yes=1;
                                                                end
                                                            end
                                                            if dist<near & i2==i3
                                                                %if red is on i2 and between blue and cyan
                                                                %aaaattttention_call_beat=1
                                                                %yes=1;
                                                            end%
                                                            
                                                            if yes==1
                                                                %if breakage update list with i of frame and x,y,z, and (im_x,im_y) of end
                                                                %and beg
                                                                %points
                                                                if plot_mode>0
                                                                    figure(100);hold on;
                                                                    scatter3(ptv(ind_i,ind_j ,3),ptv(ind_i,ind_j ,4),ptv(ind_i,ind_j ,5),'c')
                                                                    scatter3(ptv(ind_i,ind_j2,3),ptv(ind_i,ind_j2,4),ptv(ind_i,ind_j2,5),'m')
                                                                    scatter3(ptv(ind_i3,ind_j3,3),ptv(ind_i3,ind_j3,4),ptv(ind_i3,ind_j3,5),'b')
                                                                    scatter3(ex(n-1),ey(n-1),ez(n-1),'r');
                                                                end
                                                                num_breakages=num_breakages+1;
                                                                %get the past, frames i'2, and pos j'2
                                                                past=get_past(traj_ind,traj_bin,i3,ind_end3,first);
                                                                %get the future 1, frames i'2, and pos j'2
                                                                fut1=get_future(traj_ind,traj_bin,i,ind_beg,first);
                                                                %get the future 2, frames i'2, and pos j'2
                                                                fut2=get_future(traj_ind,traj_bin,i2,ind_beg2,first);
                                                                breakage=[breakage;[past' fut1' fut2']];
                                                            end
                                                        end
                                                    end
                                                end
                                            end
                                            found=1;%so it stops if it has found one point that is close
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end


num_breakages=num_breakages

%%%%%%%%%%%%%%%FUNCTIONS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function traj=find_next(ptv,traj,num_in_traj)

i=traj(num_in_traj,1);
j=traj(num_in_traj,2);
if ptv(i,j,2)>-1
    num_in_traj=num_in_traj+1;
    traj(num_in_traj,1)=i+1;
    traj(num_in_traj,2)=ptv(i,j,2)+1;
end
traj=[traj;[num_in_traj 0]];

function traj_bin=plot_traj(traj_bin,ptv,traj,fig_id,time,plot_mode)

if plot_mode==2
    figure(fig_id); hold on;
end
for i=1:length(traj)
    x(i)=ptv(traj(i,1),traj(i,2),3);
    y(i)=ptv(traj(i,1),traj(i,2),4);
    z(i)=ptv(traj(i,1),traj(i,2),5);
    id(i)=i;
    frame_i(i)=traj(i,1);
    j(i)=traj(i,2);
end
ok=1;
for i=1:length(traj)
    %if y(i)>0.03 | x(i)<9e-3 | x(i)>9.4e-3 | z(i)<0.01
    %if y(i)<0.0224 | y(i)>0.0229 | x(i)<9.35e-3 | x(i)>9.7e-3
    if y(i)<0.0212 | y(i)>0.0224 | x(i)<2.8e-3 | x(i)>4.2e-3
        ok=0;
    end
end
if 1<2%ok==1
    if plot_mode==2
        plot3(x,y,z);
    end
    traj_bin=[traj_bin; [frame_i' j' id' x' y' z']];%%%%%%%%%%%%%%%%%%%%%
end


function ptv=update(ptv,traj)

for i=1:length(traj)
    ptv(traj(i,1),traj(i,2),6)=1;
end

function const=extrapolate(index,ind_end,traj_bin,traj_ind)


theor_ind_end=traj_ind(index,2);
len=traj_ind(index,3)-(theor_ind_end-ind_end);
if len>7
    len=7;
end
ind_beg=ind_end-len+1;

for comp=1:3
    count=0;
    for i=ind_beg:ind_end
        count=count+1;
        a(count,1)=(i-ind_end)^2;
        a(count,2)=(i-ind_end)^1;
        a(count,3)=(i-ind_end)^0;
        y(count,1)=traj_bin(i,comp+3);
    end
    c=(a'*a)\(a'*y);
    const(comp,:)=c;
end


function [ex ey ez]=draw_extr(const,n,figid)

count=0;
for i=1:5
    count=count+1;
    x(count)=const(1,1)*i^2+const(1,2)*i^1+const(1,3)*i^0;
    y(count)=const(2,1)*i^2+const(2,2)*i^1+const(2,3)*i^0;
    z(count)=const(3,1)*i^2+const(3,2)*i^1+const(3,3)*i^0;
end

%figure(figid)
%scatter3(x,y,z,'k')
ex=x';
ey=y';
ez=z';

function const=backtrapolate(index,ind_beg,traj_bin,traj_ind)


theor_ind_beg=traj_ind(index,1);
len=traj_ind(index,3)-(ind_beg-theor_ind_beg);
if len>7
    len=7;
end
ind_end=ind_beg+len-1;

for comp=1:3
    count=0;
    for i=ind_beg:ind_end
        count=count+1;
        a(count,1)=(i-ind_beg)^2;
        a(count,2)=(i-ind_beg)^1;
        a(count,3)=(i-ind_beg)^0;
        y(count,1)=traj_bin(i,comp+3);
    end
    if det(a'*a)>0
        c=(a'*a)\(a'*y);
    else
        c=zeros(1,3);
    end
    const(comp,:)=c;
end


function [ex ey ez]=draw_backtr(const,n,figid,col)

count=0;
for i=0:-1:-5
    count=count+1;
    x(count)=const(1,1)*i^2+const(1,2)*i^1+const(1,3)*i^0;
    y(count)=const(2,1)*i^2+const(2,2)*i^1+const(2,3)*i^0;
    z(count)=const(3,1)*i^2+const(3,2)*i^1+const(3,3)*i^0;
end

%figure(figid)
%scatter3(x,y,z,col)
ex=x';
ey=y';
ez=z';

function past=get_past(traj_ind,traj_bin,i3,ind_end3,first)

count=0;
for i=traj_ind(i3,1):ind_end3
    count=count+1;
    past(count,1)=traj_bin(i,1)+first-1;
    past(count,2)=traj_bin(i,2);
end
if count>10
   past(1:10,:)=past(count-9:count,:);
else
   past(10-count+1:10,:)=past(1:count,:);
   past(1:10-count,:)=0;
end
past=past(1:10,:);

function fut=get_future(traj_ind,traj_bin,i,ind_beg,first)

count=0;
for i=ind_beg:traj_ind(i,2)
    count=count+1;
    fut(count,1)=traj_bin(i,1)+first-1;
    fut(count,2)=traj_bin(i,2);
end
if count>10
   %fut(1:10,:)=fut(1:count,:);
else
   fut(count+1:10,:)=0;
end
fut=fut(1:10,:);



