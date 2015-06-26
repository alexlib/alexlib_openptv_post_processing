function res=gluing_ptv_is_v2()
%gluing_ptv_is_v2(1)

%%should plot nice figures

global traj;
global traj_bin;
global ptv;

global log;
log=[];

global eps;
global gluedim;

global tot_glues;
global suc_jump;


warning off all;

global_totalpix_before =[]; 
global_totalpix_child  =[];
global_breakage_pos    =[];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
columns=9;
max_num_per_frame=500;
min_length=5;
max_jump=30;
min_scale=1e-3;
tol=1.;
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

for fig_id=1:15%15

    switch fig_id
        case 1
            first=600;last=1450;
            name_root=['\\Ifu-gwh-disk\hydromechanik3\Beat_Colloid\test_cases\Exp3b\']; 
        case 2
            first=1;last=1350;
            name_root=['\\Ifu-gwh-disk\hydromechanik3\Beat_Colloid\test_cases\Exp11b\'];
        case 3
            first=1;last=2250;
            name_root=['\\Ifu-gwh-disk\hydromechanik3\Beat_Colloid\test_cases\Exp12b\'];
        case 4
            first=1;last=750;
            name_root=['\\Ifu-gwh-disk\hydromechanik3\Beat_Colloid\test_cases\Exp18b\'];
        case 5
            first=1;last=1017;
            name_root=['\\Ifu-gwh-disk\hydromechanik3\Beat_Colloid\test_cases\Exp19b\'];
        case 6
            first=1;last=3480;
            name_root=['\\Ifu-gwh-disk\hydromechanik3\Beat_Colloid\test_cases\Exp_21b\'];
        case 7
            first=3000;last=3400;
            name_root=['\\Ifu-gwh-disk\hydromechanik3\Beat_Colloid\test_cases\Exp_22b_1\'];
        case 8
            first=3400;last=3700;
            name_root=['\\Ifu-gwh-disk\hydromechanik3\Beat_Colloid\test_cases\Exp_22b_2\'];%%-?? ok, but how to treat traj end and only one new?
        case 9
            first=1;last=3480;
            name_root=['\\Ifu-gwh-disk\hydromechanik3\Beat_Colloid\test_cases\Exp_25b\'];
        case 10
            first=1;last=1200;
            name_root=['\\Ifu-gwh-disk\hydromechanik3\Beat_Colloid\test_cases\Exp_28b\'];
        case 11
            first=4000;last=5457;
            name_root=['\\Ifu-gwh-disk\hydromechanik3\Beat_Colloid\test_cases\Exp_29b\'];
        case 12
            first=600;last=1200;
            name_root=['\\Ifu-gwh-disk\hydromechanik3\Beat_Colloid\test_cases\Exp_30b\'];%%-?? how to manage sum totalpix??? 
        case 13
            first=1000;last=2270;
            name_root=['\\Ifu-gwh-disk\hydromechanik3\Beat_Colloid\test_cases\Exp_31b\'];
        case 14
            first=1300;last=1650;
            name_root=['\\Ifu-gwh-disk\hydromechanik3\Beat_Colloid\test_cases\Exp_32b\'];
        case 15
            first=1;last=1400;
            name_root=['\\Ifu-gwh-disk\hydromechanik3\Beat_Colloid\test_cases\Exp_33b\'];
    end


    



    suc_jump=0;
    tot_glues=0;

    ptv=zeros(last-first+1,max_num_per_frame,columns+1);
    ptv(:,:,1)=-10;
    ptv(:,:,2)=-10;
    num_traj=0;
    traj_bin=[];


    for i=first:last
        if mod(i,100)==0
            i
        end
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
        [fig_id i last-first+1 tot_glues suc_jump]
        j=1;

        %%%%change to proper num per frame!!! like further below
        ind=find(ptv(i,:,1)>-10);
        num_part=length(ind);
        while j<num_part+1 %max_num_per_frame
            %find trajectories
            go=0;
            if ptv(i,j,columns+1)==0
                go=1;
                num_in_traj=init_traj(i,j);
            end
            while go==1
                old_num_in_traj=num_in_traj;
                num_in_traj=find_next(num_in_traj,columns,first,last);
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


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %plan for breakage detection
    %keep all traj and then process!
    ind_be      = find(traj_bin(:,3)==1);
    ind_en      = ind_be-1;ind_en=[ind_en; length(traj_bin)];ind_en=ind_en(2:end);
    len         = ind_en-ind_be+1;
    nu_tr       = length(ind_be);
    f_be        = traj_bin(ind_be,1);
    x_be        = traj_bin(ind_be,4);
    y_be        = traj_bin(ind_be,5);
    z_be        = traj_bin(ind_be,6);
    totalpix_be = traj_bin(ind_be,7);
    vec_be      = [x_be y_be z_be];% totalpix_be];
    %filter totalpix value

    %determine 1)large drops, 2) befor and after drop there should be a
    %difference, 3) the difference should be compensated by some next neighbour
    %trajectory.

    for j=1:length(gluedim)
        accuracy(j)=eps(gluedim(j));%[dx dy dz d_totalpix];
    end
    accuracy(4)=accuracy(4)*2;
    acc=repmat(accuracy(1:3),length(ind_be),1);

    traj_done=zeros(nu_tr,1);
    for i=1:nu_tr
        totalpix=zeros(len(i),1);
        ch_totalpix=zeros(len(i),1);
        int_ch_totalpix=zeros(len(i),1);
        for j=2:len(i)-1
            %meas. drop in totalpix
            ch_totalpix(j)=0.5*(traj_bin(ind_be(i)+j-1+1,7)-traj_bin(ind_be(i)+j-1-1,7));
            totalpix(j)=traj_bin(ind_be(i)+j-1,7);
        end
        totalpix(1)=traj_bin(ind_be(i)+1-1,7);
        totalpix(len(i))=traj_bin(ind_be(i)+len(i)-1,7);
        for j=1:len(i)
            %measure integrated drop in totalpix
            ba=j-max_jump;
            if ba<1
                ba=1;
            end
            int_ch_totalpix(j)=sum(ch_totalpix(ba:j));
        end
        %find only the local minima
        go=1;
        tmp=int_ch_totalpix;
        list=[];
        while go==1
            mini=min(tmp);
            if mini<-accuracy(4)
                ind_e=find(int_ch_totalpix==mini);
                if length(ind_e)>0 & ind_e>1
                    be=ind_e(1)-max_jump;
                    if be<1
                        be=1;
                    end
                    en=ind_e(1)+max_jump;
                    if en>length(int_ch_totalpix)
                        en=length(int_ch_totalpix);
                    end
                    int_ch_totalpix(be:ind_e(1)-1)=0;
                    int_ch_totalpix(ind_e(1)+1:en)=0;
                    tmp(be:en)=0;
                    maxi=max(totalpix(be:en));
                    ind_b=find(totalpix(be:en)==maxi)+be-1;
                    list=[list;ind_b(1) ind_e(1)]; %%list contains frames with large drops
                else
                    go=0;
                end
            else
                go=0;
            end
        end
        %preparation to determine a filtered drop measurement
        si=size(list);
        num_can=si(1,1);
        for j=1:num_can
            be=list(j,1)-max_jump;
            if be<1
                be=1;
            end
            en=list(j,2)+max_jump;
            if en>length(int_ch_totalpix)
                en=length(int_ch_totalpix);
            end
            list(j,3)=mean(totalpix(be:list(j,1)-1));
            if en>list(j,2) %% treat if drop at absolute end of traj
                list(j,4)=mean(totalpix(list(j,2)+1:en));
            else
                list(j,4)=totalpix(list(j,2));
            end
        end
        % up to now the list of potential breakage points ALONG a trajectory i is determined.

        %% here the end of every traj is added to the list
        num_can=num_can+1;
        list(num_can,1)=len(i);
        list(num_can,2)=len(i);
        be=list(num_can,1)-max_jump;
        if be<1
            be=1;
        end
        list(num_can,3)=mean(totalpix(be:list(num_can,1)-1));
        list(num_can,4)=0;

        %%check if anything is close first in time, then in space
        if num_can>0 & traj_done(i)==0 & len(i)>max_jump
            ind=find(list(:,3)-list(:,4)>accuracy(4)); %%is the diff between mean before and after large enough?
            for j=1:length(ind)
                %%here begins the story of a new potential breakage
                totalpix_before=0;
                totalpix_child=[];
                ind_time=find(traj_bin(ind_be(i)+list(ind(j),1)-1,1)-max_jump<f_be & f_be<traj_bin(ind_be(i)+list(ind(j),2)-1,1)+4*max_jump & len>max_jump & traj_done==0);
                if length(ind_time)>0
                    vec_en =[traj_bin(ind_be(i)+list(ind(j),1)-1,4) traj_bin(ind_be(i)+list(ind(j),1)-1,5) traj_bin(ind_be(i)+list(ind(j),1)-1,6)];
                    vec_en=repmat(vec_en,length(ind_time),1);
                    vec_be_proj=proj_be_to_en(vec_be,ind_time,ind_be,ind_en,max_jump,traj_bin(ind_be(i)+list(ind(j),1)-1,1));
                    dist=(vec_be_proj-vec_en)./acc(ind_time,:);
                    dist=sum(dist.^2,2).^0.5;
                    ind_begin=find(dist<5*tol*sqrt(3));
                    en_done=0;
                    for k=1:length(ind_begin)
                        scale=((traj_bin(ind_be(ind_time(ind_begin(k))),4)-traj_bin(ind_en(ind_time(ind_begin(k))),4))^2+...
                            (traj_bin(ind_be(ind_time(ind_begin(k))),5)-traj_bin(ind_en(ind_time(ind_begin(k))),5))^2+...
                            (traj_bin(ind_be(ind_time(ind_begin(k))),6)-traj_bin(ind_en(ind_time(ind_begin(k))),6))^2)^0.5;
                        if abs(ind_time(ind_begin(k))-i)>0 & scale>min_scale
                            %%%%found breakage!!!!!!!!!!!!!!!!!
                            %%%%found breakage!!!!!!!!!!!!!!!!!
                            %%%%found breakage!!!!!!!!!!!!!!!!!
                            if en_done==0
                                traj_done(i)=1;
                                en_done=1;
                                figure(fig_id);hold on;
                                scatter3(vec_en(1,1),vec_en(1,2),list(ind(j),3),10,'r');
                                totalpix_before=list(ind(j),3);
                                %%check if also begin of new child, since on the same traj
                                if(ind_be(i)+list(ind(j),2)<ind_en(i))
                                    scatter3(traj_bin(ind_be(i)+list(ind(j),2)-1,4),traj_bin(ind_be(i)+list(ind(j),2)-1,5),list(ind(j),4),10,'g');
                                    totalpix_child=[totalpix_child;list(ind(j),4)];
                                end

%                                 figure(100+0);hold on;
%                                 xx=traj_bin(ind_be(i):ind_en(i),4);
%                                 yy=traj_bin(ind_be(i):ind_en(i),5);
%                                 zz=traj_bin(ind_be(i):ind_en(i),6);
%                                 cc=traj_bin(ind_be(i):ind_en(i),7);
%                                 ii=find(yy>0.03 & yy<0.035);
%                                 scatter3(xx(ii),yy(ii),zz(ii),5,cc(ii));
                                
                                figure(100+fig_id);hold on;
                                plot3(traj_bin(ind_be(i):ind_en(i),4),...
                                    traj_bin(ind_be(i):ind_en(i),5),...
                                    traj_bin(ind_be(i):ind_en(i),7));
                                scatter3(vec_en(1,1),vec_en(1,2),list(ind(j),3),20,'r','filled');
                                if(ind_be(i)+list(ind(j),2)<ind_en(i))
                                    scatter3(traj_bin(ind_be(i)+list(ind(j),2)-1,4),traj_bin(ind_be(i)+list(ind(j),2)-1,5),list(ind(j),4),20,'g','filled');
                                    totalpix_child=[totalpix_child;list(ind(j),4)];
                                end

                            end
                            traj_done(ind_time(ind_begin(k)))=1;
                            figure(fig_id);hold on;
                            scatter3(vec_be(ind_time(ind_begin(k)),1),vec_be(ind_time(ind_begin(k)),2),totalpix_be(ind_time(ind_begin(k))),10,'g');
                            totalpix_child=[totalpix_child;totalpix_be(ind_time(ind_begin(k)))];
                            
%                             figure(100+0);hold on;
%                             xx=traj_bin(ind_be(ind_time(ind_begin(k))):ind_en(ind_time(ind_begin(k))),4);
%                             yy=traj_bin(ind_be(ind_time(ind_begin(k))):ind_en(ind_time(ind_begin(k))),5);
%                             zz=traj_bin(ind_be(ind_time(ind_begin(k))):ind_en(ind_time(ind_begin(k))),6);
%                             cc=traj_bin(ind_be(ind_time(ind_begin(k))):ind_en(ind_time(ind_begin(k))),7);
%                             ii=find(yy>0.03 & yy<0.035);
%                             scatter3(xx(ii),yy(ii),zz(ii),5,cc(ii));
                            
                            figure(100+fig_id);hold on;
                            plot3(traj_bin(ind_be(ind_time(ind_begin(k))):ind_en(ind_time(ind_begin(k))),4),...
                                traj_bin(ind_be(ind_time(ind_begin(k))):ind_en(ind_time(ind_begin(k))),5),...
                                traj_bin(ind_be(ind_time(ind_begin(k))):ind_en(ind_time(ind_begin(k))),7));
                            scatter3(vec_be(ind_time(ind_begin(k)),1),vec_be(ind_time(ind_begin(k)),2),totalpix_be(ind_time(ind_begin(k))),20,'g','filled');
                        end
                    end
                end
                %%this completes the story a particular breakage
                %here we can do some statistics
                if length(totalpix_child)>1
                    global_totalpix_before = [global_totalpix_before;totalpix_before];
                    tmp=zeros(1,10);
                    tmp(1,1:length(totalpix_child))=totalpix_child';
                    global_totalpix_child  = [global_totalpix_child;tmp];
                    global_breakage_pos    = [global_breakage_pos;...
                                              [traj_bin(ind_be(i)+list(ind(j),1)-1,4) traj_bin(ind_be(i)+list(ind(j),1)-1,5) traj_bin(ind_be(i)+list(ind(j),1)-1,6)]];
                end
            end
        end
    end
end


save stats global_totalpix_before global_totalpix_child global_breakage_pos

figure;scatter3(global_breakage_pos(:,1),global_breakage_pos(:,2),global_breakage_pos(:,3))

figure;hist(global_totalpix_before)
xlabel('totalpix before')

ind=find(global_totalpix_child>0);
figure;hist(global_totalpix_child(ind))
xlabel('totalpix children')

si=size(global_totalpix_child);
for i=1:si(1,1)
   ind=find(global_totalpix_child(i,:)>0);
   howmany(i)=length(ind);
end
figure;hist(howmany)
xlabel('# children')

figure;scatter(global_totalpix_before,howmany)
xlabel('totalpix before')
ylabel('# children')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%FUNCTIONS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function num_in_traj=init_traj(i,j)

global traj;

traj=[0 0];
traj(1,1)=i;
traj(1,2)=j;
num_in_traj=1;

function num_in_traj=find_next(num_in_traj,columns,first,last)

global traj;
global ptv;

i=traj(num_in_traj,1);
j=traj(num_in_traj,2);
if i<last-first+1
    if ptv(i,j,2)>-1 & ptv(i+1,ptv(i,j,2)+1,columns+1)==0
        num_in_traj=num_in_traj+1;
        traj(num_in_traj,1)=i+1;
        traj(num_in_traj,2)=ptv(i,j,2)+1;%% C -> Matlab
    end
end

function vec_be_proj=proj_be_to_en(vec_be,ind_time,ind_be,ind_en,max_jump,proj_frame)

global traj_bin;

vec_be_proj=[];
order=2;

for i=1:length(ind_time)
    A=zeros(max_jump+1,order+1);
    y=zeros(max_jump+1,order+1);
    count=0;
    for frame=traj_bin(ind_be(ind_time(i)),1):traj_bin(ind_be(ind_time(i)),1)+max_jump
        count=count+1;
        for or=0:order
            A(count,or+1)=frame^or;
        end
        for n=1:3
            y(count,n)=traj_bin(ind_be(ind_time(i))+count-1,3+n);
        end
    end
    X=(A'*A)\A'*y;

    %loop through jump size until max jump is reached
    proj=[];
    for n=1:n
        proj(n)=0;
        for or=0:order
            proj(n)=proj(n)+X(or+1,n)*proj_frame^or;
        end
    end
    vec_be_proj=[vec_be_proj;proj];
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
    traj_bin=[traj_bin; [frame_i' j' id' x' y' z' totalpix']];%%%%%%%%%%%%%%%%%%%%%
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






