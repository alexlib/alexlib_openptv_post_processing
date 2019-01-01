function [ptv, tot_glues, suc_jump, go] = glue_traj(traj, ptv, log, eps, gluedim, tot_glues, suc_jump, max_jump,start_frame,first,last,columns,tol)

% global traj;
% global ptv;

% global log;

% global eps;
% global gluedim;

% global tot_glues;
% global suc_jump;

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
if det(A'*A)>1e-15
    X=(A'*A)\A'*y;
    nogo=0;
else
    nogo=1;
end

%loop through jump size until max jump is reached

if nogo==0
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
            %create points in-betwen
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
            %update ptv array
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
end