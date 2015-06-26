function res=gluing_ptv_is_v1(first,last,columns,tol,fig_id)
%gluing_ptv_is_v1(1,3000,9,0.3,1)

global traj;
global traj_bin;
global ptv;

global log;
log=[];

global eps;
global gluedim;

global tot_glues;
global suc_jump;
suc_jump=0;
tot_glues=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
max_num_per_frame=100;
min_length=5;
max_jump=30;
dx         = 0.2e-3;
dy         = 0.2e-3;
dz         = 0.4e-3;
d_totalpix = 10;
d_xpx      = 1;
d_ypx      = 1;
d_sumgrv   = 7000;
eps=    [0 0 dx dy dz d_totalpix d_xpx d_ypx d_sumgrv];
gluedim=[    3  4  5  6                              ];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


ptv=zeros(last-first+1,max_num_per_frame,columns+1);
ptv(:,:,1)=-10;
ptv(:,:,2)=-10;
num_traj=0;
traj_bin=[];


name_root=['E:\PTV\Working_folder\Exp25b_09112010\res\'];

for i=first:last
    name=[name_root,'ptv_is_v2.',num2Str(i)];
    fid = fopen(name, 'r');
    num_points = fscanf(fid, '%i', [1 1]);    % It has two rows now.
    tmp = fscanf(fid, '%i %i %f %f %f %f %f %f %f', [columns num_points]);
    A=tmp';
    A(:,3:5)=A(:,3:5)*0.001;
    ptv(i-first+1,1:num_points,1:columns)=A;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fclose(fid);
end

go=0;
for i=1:last-first+1
    [i last-first+1 tot_glues suc_jump]
    j=1;
    while j<max_num_per_frame
        %find trajectories
        go=0;
        if ptv(i,j,columns+1)==0
            go=1;
            num_in_traj=init_traj(i,j);
        end
        while go==1
            old_num_in_traj=num_in_traj;
            num_in_traj=find_next(num_in_traj,columns);
            si=size(traj);
            if num_in_traj>last-first-i+1 | si(1,1)>last-first-1
                go=0;
                %filter, plot, etc....
                if si(1,1)>min_length
                    num_traj=num_traj+1;
                    ptv=filter_traj();
                    plot_traj(fig_id,first+i-1);
                    update(columns);
                end
            elseif old_num_in_traj==num_in_traj
                go=0;
                if si(1,1)>min_length
                    %HERE IS THE ACTUAL CHECK IF IT CAN BE GLUED
                    go=glue_traj(max_jump,i,first,last,columns,tol);
                end
                if go==0 
                    %filter, plot, etc....
                    if si(1,1)>min_length
                        num_traj=num_traj+1;
                        ptv=filter_traj();
                        plot_traj(fig_id,first+i-1);
                        update(columns);
                    end
                end
            end
        end
        j=j+1;
    end
end



%%%%%%%%%%%%%%%FUNCTIONS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function num_in_traj=init_traj(i,j)

global traj;

traj=[0 0];
traj(1,1)=i;
traj(1,2)=j;
num_in_traj=1;

function num_in_traj=find_next(num_in_traj,columns)

global traj;
global ptv;

i=traj(num_in_traj,1);
j=traj(num_in_traj,2);
if ptv(i,j,2)>-1 & ptv(i+1,ptv(i,j,2)+1,columns+1)==0 
    num_in_traj=num_in_traj+1;
    traj(num_in_traj,1)=i+1; 
    traj(num_in_traj,2)=ptv(i,j,2)+1;%% C -> Matlab
end

function go=glue_traj(max_jump,start_frame,first,last,columns,tol)

global traj;
global ptv;

global log;

global eps;
global gluedim;

global tot_glues;
global suc_jump;

for i=1:length(gluedim)
    accuracy(i)=eps(gluedim(i));%[dx dy dz d_totalpix];
end

si=size(traj);
%determine order, which is used to attempt gluing jump
order=floor(si(1,1)/4)-1;
if order<0
    order=0;
end
if order>2
    order=2;
end

%prepare polynomial fit for jumps
en=si(1,1);
be=en-max_jump;
if be<1
    be=1;
end
be=be+start_frame-1;
en=en+start_frame-1;

A=zeros(en-be+1,order+1);
y=zeros(en-be+1,length(gluedim));
count=0;
for frame=be:en
    count=count+1;
    for or=0:order
        A(count,or+1)=frame^or;
    end
    i=traj(frame-start_frame+1,1);
    j=traj(frame-start_frame+1,2);
    for n=1:length(gluedim)%3:6
        y(count,n)=ptv(i,j,gluedim(n));
    end
end
X=(A'*A)\A'*y;

%loop through jump size until max jump is reached
for jump=1:max_jump
   %jump and check in n-dim space whether anything is
   %close enough or not, i.e. x,u?,size?
   frame=en+jump;
   if frame>last-first+1
       go=0;
       break;
   end
   %howmany points in frame n?
   ind=find(ptv(frame,:,1)>-10);
   num_part=length(ind);
   proj=[];
   for n=1:length(gluedim)%3:6
       proj(n)=0;
       for or=0:order
           proj(n)=proj(n)+X(or+1,n)*frame^or;
       end
   end
   proj=repmat(proj,num_part,1);
   acc=repmat(accuracy,num_part,1);
   if num_part>1
       dist=(proj-squeeze(ptv(frame,1:num_part,gluedim)))./acc;
   else %%%fix for unwanted transposed
       dist=(proj-squeeze(ptv(frame,1:num_part,gluedim))')./acc;
   end
   %exclude those canditates, which are already part of another trajectory
   ind_occ=find(ptv(frame,1:num_part,columns+1)==1);
   dist=sum(dist.^2,2).^0.5;
   dist(ind_occ)=1e9;
   mini=min(dist);
   ind=find(dist==mini);
   if length(ind)>0 & mini<tol*sqrt(length(gluedim)) & length(mini) > 0
       %if close enough create new points in gap etc...
       go=1;
       tot_glues=tot_glues+1;
       suc_jump=jump;
       log=[log;[en en+jump jump]];
       %create points in-between
       i=traj(en-start_frame+1,1);
       j=traj(en-start_frame+1,2);
       j_ar(1)=j;
       j_ar(jump+1)=ind;
       p_g(1,1:columns)=squeeze(ptv(i,j,1:columns))';
       p_g(jump+1,1:columns)=squeeze(ptv(frame,ind,1:columns))';
       for fr=en+1:en+jump-1
           %howmany points in frame fr?
           ind=find(ptv(fr,:,1)>-10);
           num_part=length(ind);
           j_ar(fr-en+1)=num_part+1;
           w_e=(fr-en)/jump;
           w_b=1-w_e;
           p_g(fr-en+1,1:columns)=w_b*p_g(1,1:columns)+w_e*p_g(jump+1,1:columns);
       end
       %update ptv array, the links in future and past that is..
       ptv(en,j_ar(1),2)=j_ar(2)-1;
       for fr=2:jump
           ptv(en+fr-1,j_ar(fr),1:columns)=p_g(fr,1:columns);
           if fr>1 & fr<jump+1
               ptv(en+fr-1,j_ar(fr),1)=j_ar(fr-1)-1;
               ptv(en+fr-1,j_ar(fr),2)=j_ar(fr+1)-1;
           end
       end
       ptv(en+jump,j_ar(jump+1),1)=j_ar(jump)-1;
       
       break; %%breaks of for loop through jump
   else
       go=0;
   end
end


function res=plot_traj(fig_id,time)

global traj;
global traj_bin;
global ptv;

figure(fig_id); hold on;
for i=1:length(traj)
    x(i)=ptv(traj(i,1),traj(i,2),3);
    y(i)=ptv(traj(i,1),traj(i,2),4);
    z(i)=ptv(traj(i,1),traj(i,2),5);
    totalpix(i)=ptv(traj(i,1),traj(i,2),6);
    xpx(i)=ptv(traj(i,1),traj(i,2),7);
    ypx(i)=ptv(traj(i,1),traj(i,2),8);
    sumgrv(i)=ptv(traj(i,1),traj(i,2),9);
    id(i)=i;
    frame_i(i)=traj(i,1);
    j(i)=traj(i,2);
end
ok=1;

if 1<2 %min(x)<0.01
    plot3(x,y,totalpix,'b');
    %scatter3(x,y,z,5,totalpix);
    traj_bin=[traj_bin; [frame_i' j' id' x' y' z']];%%%%%%%%%%%%%%%%%%%%%
end

% figure(fig_id+1); hold on;
% if length(traj)>3 & min(x)<0.01
%     scatter3(x,y,z,5,xpx);
% end
% figure(fig_id+2); hold on;
% if length(traj)>3 & min(x)<0.01
%     scatter3(x,y,z,5,ypx);
% end
% figure(fig_id+3); hold on;
% if length(traj)>3 & min(x)<0.01
%     scatter3(x,y,z,5,sumgrv);
% end


function res=update(columns)

global traj;
global ptv;

for i=1:length(traj)
    ptv(traj(i,1),traj(i,2),columns+1)=1;
end

function ptv=filter_traj()

global traj;
global ptv;

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






