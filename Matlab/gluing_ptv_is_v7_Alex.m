function gluing_ptv_is_v7_Alex()
% Gluing_PTV_IS_V7 is searching for trajectories that broke into pieces
% used by Debashish Saha in the aggregate breakage detection in the orifice
% experiment, see data under github.com/openptv

% global traj;
% global traj_bin;
% global ptv;

traj = [];

% global log;
log = [];

% global eps;
% global gluedim;

 %global tot_glues;
% global suc_jump;

% warning off all;

% global_totalpix_before =[];
% global_totalpix_child  =[];
% global_breakage_pos    =[];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
columns     = 32;%9;% No of variables in ptv_is files
max_num_per_frame = 100;
min_length  = 5; % Trajectories travel over at least 5 frame
max_jump    = 30;
min_scale   = 1e-3; % Trajectories' min length in metric unit
tol = 1; % spatial tol for gluing
break_tol   = 5; % spatial tol for parent child match
dx          = 0.2e-3;
dy          = 0.2e-3;
dz          = 0.4e-3;
d_totalpix  = 10;
d_xpx       = 1;
d_ypx       = 1;
d_sumgrv    = 7000;

eps         =  [0 0 dx dy dz d_totalpix d_xpx d_ypx d_sumgrv ];%S_ij,L1,L2,L3,vel
gluedim     =  [    3  4  5  6                              ];


% foll={%'3Nov2010\';
%     %'15Nov2010\';
%     %'Deb_6Nov2010\';
%     %'Exp_07112010\';
%     %'Exp_08112010\';
%     %'Exp_09112010\';
%     %'Exp_10112010\';
%     %'Exp_11112010\';
%     %'Exp_11112010_REST\';
%     %'Exp_12112010\';
%     %'Exp_14112010\';
%     % '/';
%     %'Exp3b_14112010';
%     %'missing_breakage_event\'; %%%%%% this is the showcase event
%     %'rehab_3Nov2010\';
%     %'rehab_15Nov2010\';
%     %'rehab_Exp_11112010_REST\';
%     %'rehab_Exp_14112010\';
%     %'rehab_Exp_15112010\'};
%     };
% for globi=1:length(foll) %%%%%%%
%     close all
clear dummy ar fig_id_count first_count last_count
% dummy=['/Users/alex/Dropbox/ptv_is_koni/',foll{globi}];
dummy=['/Users/alex/Dropbox/ptv_is_koni/'];


disp('looking into folder, gathering experiments')
ar=dir([dummy,'Exp*']);

br = cell(length(ar),1);
for k=1:length(ar)
%     fig_id_count(k)=k;
    br{k} = dir([dummy,ar(k).name,'/ptv_is_koni*']);
end


for k=1:length(ar)
    for l=1:length(br{k})
        cr{k,l}=regexp(br{k}(l).name,'\.','split');
        dr(l)=str2num(cr{k,l}{2});
    end
    first_count(k)=min(dr);
    last_count(k)=max(dr);
    clear dr
end

for fig_id = 1:length(ar) %[7 8 9 11 12 13 14 15]%[1 2 4 5 7]% 9 10 11 12 13 14]% 5 6 7 8 9 10 11 12 13 14]%2:2%15%[1 2 4 5 7]
    %     close all %%%%%
    te=['dealing with ',ar(fig_id).name];
    disp(te);
    
    first=first_count(fig_id);
    last=last_count(fig_id);
    name_root=[dummy,ar(fig_id).name,filesep];
    
    break_pos=[];
    unbreak_pos_par=[];
    unbreak_pos_child=[];
    start_child=[];
    pixels=[];
    
    suc_jump=0;
    tot_glues=0;
    
    ptv=zeros(last-first+1,max_num_per_frame,columns+1);
    ptv(:,:,1)=-10;
    ptv(:,:,2)=-10;
    num_traj=0;
    traj_bin = [];
    
    disp('reading ptv files')
    for i = first:last
        %         if mod(i,100)==0
        %             i
        %         end
        name=[name_root,'ptv_is_koni.',num2str(i)];
        fid = fopen(name, 'r');
        num_points = fscanf(fid, '%i', [1 1]);    % It has two rows now.
        tmp = fscanf(fid, '%i %i %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f', [columns num_points]);
        A=tmp';
        A(:,3:5)=A(:,3:5)*0.001; %%%------ in extended ptv_is_v2,x,y,z are in meter
        A(:,23:25)=A(:,23:25)*0.001; %%%------ in extended ptv_is_v2,x,y,z are in meter
        ptv(i-first+1,1:num_points,1:columns)=A;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fclose(fid);
    end
    
    disp('gluing')
    go=0;
    for i=1:last-first+1
        %         [fig_id i last-first+1 tot_glues suc_jump]
        j=1;
        
        %%%%change to proper num per frame!!! like further below
        ind=find(ptv(i,:,1)>-10);
        num_part=length(ind);
        while j<num_part+1 %max_num_per_frame
            %find trajectories
            go=0;
            if ptv(i,j,columns+1)==0
                go=1;
                [traj, num_in_traj] = init_traj(traj, i,j);
            end
            while go==1
                old_num_in_traj=num_in_traj;
                [traj, num_in_traj] = find_next(traj, ptv, num_in_traj,columns,first,last);
                si=size(traj);
                if num_in_traj>last-first-i+1 || si(1,1)>last-first-1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    go=0;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %filter, plot, etc....
                    if si(1,1) > min_length %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        num_traj = num_traj+1;
                        ptv = filter_traj(traj, ptv);%%%%%%%%%%%%%%%%%%%%%%%%
                        traj_bin = plot_traj(traj, traj_bin, ptv, fig_id,first+i-1);
                        ptv = update(traj, ptv, columns);
                    end
                elseif old_num_in_traj==num_in_traj%%%%%%%%%
                    go=0;
                    if si(1,1)>min_length
                        %HERE IS THE ACTUAL CHECK IF IT CAN BE GLUED
                        % go=glue_traj(max_jump,i,first,last,columns,tol);
                        [ptv, tot_glues, suc_jump, go] = glue_traj(traj, ptv, log, eps, gluedim, tot_glues, suc_jump,max_jump,i,first,last,columns,tol);
                    end
                    if go==0
                        %filter, plot, etc....
                        if si(1,1)>min_length
                            num_traj = num_traj+1;
                            ptv = filter_traj(traj, ptv);
                            traj_bin = plot_traj(traj, traj_bin, ptv, fig_id,first+i-1);
                            ptv = update(traj, ptv, columns);
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
    
    disp('now doing the actual stuff, the breakage detection')
    
    %plan for breakage detection
    %keep all traj and then process!
    ind_be      = find(traj_bin(:,3)==1);
    ind_en      = ind_be-1;ind_en=[ind_en; length(traj_bin)];ind_en=ind_en(2:end);
    len         = ind_en-ind_be+1;
    nu_tr       = length(ind_be);%% there will be only one 1?
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    f_be        = traj_bin(ind_be,1);
    x_be        = traj_bin(ind_be,4);
    y_be        = traj_bin(ind_be,5);
    z_be        = traj_bin(ind_be,6);
    totalpix_be = traj_bin(ind_be,7);
    strain_be = traj_bin(ind_be,20);
    lambda1_be = traj_bin(ind_be,21);
    lambda2_be = traj_bin(ind_be,22);
    lambda3_be = traj_bin(ind_be,23);
    gradmeas_be = traj_bin(ind_be,24);
    frame_be = traj_bin(ind_be,1);
    size_meas1_be=traj_bin(ind_be,27);
    size_meas2_be=traj_bin(ind_be,28);
    grv_be=traj_bin(ind_be,29);
    xpx_be=traj_bin(ind_be,25);
    ypx_be=traj_bin(ind_be,26);
    vec_be      = [x_be y_be z_be];% totalpix_be];
    
    used=zeros(length(len));
    
    %filter totalpix value
    
    %determine 1)large drops, 2) befor and after drop there should be a
    %difference, 3) the difference should be compensated by some next neighbour
    %trajectory.
    
    %     for j=1:length(gluedim)
    %         accuracy(j)=eps(gluedim(j));%[dx dy dz d_totalpix];
    %     end
    accuracy=eps(gluedim);
    accuracy(4)=accuracy(4)*2;
    acc=repmat(accuracy(1:3),length(ind_be),1);
    
    for i=1:nu_tr
        child_count(i)=0;
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
        %             figure;
        %             plot(ch_totalpix,'r')
        for j=1:len(i)
            %measure integrated drop in totalpix
            ba=j-max_jump;
            if ba<1
                ba=1;
            end
            int_ch_totalpix(j)=sum(ch_totalpix(ba:j));
        end
        
        %             figure;plot(int_ch_totalpix)
        %             hold on;plot(ch_totalpix,'r')
        %
        %             plot(280:310,ch_totalpix(280:310),'.k')
        
        %             figure;
        %             plot(totalpix,'g')
        %             hold on
        %             plot(280:310,totalpix(280:310),'.k')
        %find only the local minima
        go=1;
        tmp=int_ch_totalpix;
        list=[];
        while go==1
            mini=min(tmp);
            if mini<-accuracy(4)
                ind_e=find(tmp==mini);%%%%%%%%%%%%%%%
                if length(ind_e)>0 & ind_e>1
                    be=ind_e(1)-max_jump;
                    if be<1
                        be=1;
                    end
                    en=ind_e(1)+max_jump;%%%%%%%%%%%%%%%%%%
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
            list(j,3)=nanmean(totalpix(be:list(j,1)-1));
            if en>list(j,2) %% treat if drop at absolute end of traj
                list(j,4)=nanmean(totalpix(list(j,2)+1:en));
            else
                list(j,4)=totalpix(list(j,2));
            end
        end
        % up to now the list of potential breakage points ALONG a trajectory i is determined.
        
        % here the end of every traj is added to the list, since it could
        % also be a breakage
        num_can=num_can+1;
        list(num_can,1)=len(i);
        list(num_can,2)=len(i);
        be=list(num_can,1)-max_jump;
        if be<1
            be=1;
        end
        list(num_can,3)=mean(totalpix(be:list(num_can,1)-1));
        list(num_can,4)=0;
        
        %check if anything is close first in time, then in space
        if num_can>0 && len(i)>max_jump
            % now check if the diff between mean before and after 'large enough?'
            % remember to check the quivalent for the traj end, i.e. is the
            % check again whether drop between traj end and begin of new traj large enough?
            ind=find(list(:,3)-list(:,4)>accuracy(4));
            for j=1:length(ind) % j loops through all drops along i
                %%here begins the story of a new potential breakage
                totalpix_before=0;
                totalpix_child=[];
                % look for child in time
                ind_time=find(traj_bin(ind_be(i)+list(ind(j),1)-1,1)-max_jump<f_be & f_be<traj_bin(ind_be(i)+list(ind(j),2)-1,1)+4*max_jump & len>max_jump );
                if ~isempty(ind_time)
                    vec_en =[traj_bin(ind_be(i)+list(ind(j),1)-1,4) traj_bin(ind_be(i)+list(ind(j),1)-1,5) traj_bin(ind_be(i)+list(ind(j),1)-1,6)];
                    strain_en =traj_bin(ind_be(i)+list(ind(j),1)-1,20);
                    lambda1_en =traj_bin(ind_be(i)+list(ind(j),1)-1,21);
                    lambda2_en =traj_bin(ind_be(i)+list(ind(j),1)-1,22);
                    lambda3_en =traj_bin(ind_be(i)+list(ind(j),1)-1,23);
                    grad_meas_en =traj_bin(ind_be(i)+list(ind(j),1)-1,24);
                    size_meas1_en =traj_bin(ind_be(i)+list(ind(j),1)-1,27);
                    size_meas2_en =traj_bin(ind_be(i)+list(ind(j),1)-1,28);
                    frame_en =traj_bin(ind_be(i)+list(ind(j),1)-1,1);
                    xpx_en =traj_bin(ind_be(i)+list(ind(j),1)-1,25);
                    ypx_en =traj_bin(ind_be(i)+list(ind(j),1)-1,26);
                    vec_en=repmat(vec_en,length(ind_time),1);
                    % look for child in space, to do that the potential begins
                    % of new chlidren are projected to the time frame of vec_en
                    vec_be_proj=proj_be_to_en(traj_bin, vec_be,ind_time,ind_be,ind_en,max_jump,traj_bin(ind_be(i)+list(ind(j),1)-1,1));
                    % here the delta's are normailzed with measurement
                    % accuracy
                    dist=(vec_be_proj-vec_en)./acc(ind_time,:);
                    dist=sum(dist.^2,2).^0.5;
                    ind_begin=find(dist<break_tol*sqrt(3));
                    en_done=0;
                    for k=1:length(ind_begin)+1 % k loops through ALL potential children of i,j
                        
                        if k==length(ind_begin)+1
                            if(ind_be(i)+list(ind(j),2)<ind_en(i)) && ind_en(i)-ind_be(i)>3*max_jump
                                figure(fig_id);hold on;
                                xp=traj_bin(ind_be(i)+list(ind(j),1)-1,4);
                                yp=traj_bin(ind_be(i)+list(ind(j),1)-1,5);
                                zp=traj_bin(ind_be(i)+list(ind(j),1)-1,7);
                                zps=traj_bin(ind_be(i)+list(ind(j),1)-1,20);
                                xc=traj_bin(ind_be(i)+list(ind(j),2)-1,4);
                                yc=traj_bin(ind_be(i)+list(ind(j),2)-1,5);
                                zc=traj_bin(ind_be(i)+list(ind(j),2)-1,7);
                                zcs=traj_bin(ind_be(i)+list(ind(j),2)-1,20);
                                
                                
                                %                                     plot3(traj_bin(ind_be(i):ind_be(i)+list(ind(j),1)-1,4),...
                                %                                     traj_bin(ind_be(i):ind_be(i)+list(ind(j),1)-1,5),...
                                %                                     traj_bin(ind_be(i):ind_be(i)+list(ind(j),1)-1,7));
                                
                                
                                % draw parent end
                                scatter3(xp,yp,zp,10,'r');
                                % draw child beg
                                scatter3(xc,yc,zc,20,'g','filled');
                                % draw line parent-child
                                plot3([xp xc],[yp yc],[zp zc],'c')
                                
                                % clean totalpix figure
                                figure(100+fig_id);hold on;title('total pix')
                                
                                plot3(traj_bin(ind_be(i):ind_en(i),4),...
                                    traj_bin(ind_be(i):ind_en(i),5),...
                                    traj_bin(ind_be(i):ind_en(i),7));
                                
                                
                                break_pos=[break_pos;xp,yp,traj_bin(ind_be(i)+list(ind(j),1)-1,6),zp,zps,traj_bin(ind_be(i)+list(ind(j),1)-1,1),traj_bin(ind_be(i)+list(ind(j),1)-1,21),traj_bin(ind_be(i)+list(ind(j),1)-1,22),traj_bin(ind_be(i)+list(ind(j),1)-1,23),traj_bin(ind_be(i)+list(ind(j),1)-1,24),traj_bin(ind_be(i)+list(ind(j),1)-1,27),traj_bin(ind_be(i)+list(ind(j),1)-1,28),i,traj_bin(ind_be(i)+list(ind(j),1)-1,29),traj_bin(ind_be(i)+list(ind(j),1)-1,25),traj_bin(ind_be(i)+list(ind(j),1)-1,26)]; %#ok<AGROW>
                                unbreak_pos_par=[unbreak_pos_par;traj_bin(ind_be(i):ind_en(i),4),traj_bin(ind_be(i):ind_en(i),5),traj_bin(ind_be(i):ind_en(i),6),traj_bin(ind_be(i):ind_en(i),7),traj_bin(ind_be(i):ind_en(i),20),traj_bin(ind_be(i):ind_en(i),1),traj_bin(ind_be(i):ind_en(i),21),traj_bin(ind_be(i):ind_en(i),22),traj_bin(ind_be(i):ind_en(i),23),traj_bin(ind_be(i):ind_en(i),24),traj_bin(ind_be(i):ind_en(i),27),traj_bin(ind_be(i):ind_en(i),28),i*ones(length(traj_bin(ind_be(i):ind_en(i),4)),1),traj_bin(ind_be(i):ind_en(i),29),traj_bin(ind_be(i):ind_en(i),25),traj_bin(ind_be(i):ind_en(i),26)]; %#ok<AGROW>
                                start_child=[start_child;xc,yc,traj_bin(ind_be(i)+list(ind(j),2)-1,6),zc,zcs,traj_bin(ind_be(i)+list(ind(j),2)-1,1),traj_bin(ind_be(i)+list(ind(j),2)-1,21),traj_bin(ind_be(i)+list(ind(j),2)-1,22),traj_bin(ind_be(i)+list(ind(j),2)-1,23),traj_bin(ind_be(i)+list(ind(j),2)-1,24),traj_bin(ind_be(i)+list(ind(j),2)-1,27),traj_bin(ind_be(i)+list(ind(j),2)-1,28),i,traj_bin(ind_be(i)+list(ind(j),2)-1,29),traj_bin(ind_be(i)+list(ind(j),2)-1,25),traj_bin(ind_be(i)+list(ind(j),2)-1,26)]; %#ok<AGROW>
                                
                                
                                
                                % draw parent end
                                scatter3(xp,yp,zp,10,'r');
                                % draw child beg
                                scatter3(xc,yc,zc,20,'g','filled');
                                % draw line parent-child
                                plot3([xp xc],[yp yc],[zp zc],'c')
                                
                                
                                
                                
                                % clean strain figure
                                figure(200+fig_id);hold on;title('strain')
                                
                                plot3(traj_bin(ind_be(i):ind_en(i),4),...
                                    traj_bin(ind_be(i):ind_en(i),5),...
                                    traj_bin(ind_be(i):ind_en(i),20));
                                
                                
                                
                                % draw parent end
                                scatter3(xp,yp,zps,10,'r');
                                % draw child beg
                                scatter3(xc,yc,zcs,20,'g','filled');
                                % draw line parent-child
                                plot3([xp xc],[yp yc],[zps zcs],'c')
                                
                                
                                
                                
                            end
                        else
                            
                            scale=((traj_bin(ind_be(ind_time(ind_begin(k))),4)-traj_bin(ind_en(ind_time(ind_begin(k))),4))^2+...
                                (traj_bin(ind_be(ind_time(ind_begin(k))),5)-traj_bin(ind_en(ind_time(ind_begin(k))),5))^2+...
                                (traj_bin(ind_be(ind_time(ind_begin(k))),6)-traj_bin(ind_en(ind_time(ind_begin(k))),6))^2)^0.5;
                            %now check if this new trajectory is long enough
                            %as measured in distance, not in terms of time
                            if abs(ind_time(ind_begin(k))-i)>0 && scale>min_scale && used(ind_time(ind_begin(k)))==0
                                % and now FINALLY, again check if drop is large enough
                                % i.e. check again whether drop between traj end and begin of new traj large enough?
                                % above it was ind=find(list(:,3)-list(:,4)>accuracy(4));
                                totalpix=zeros(len(ind_time(ind_begin(k))),1);
                                for l=1:len(ind_time(ind_begin(k)))
                                    totalpix(l)=traj_bin(ind_be(ind_time(ind_begin(k)))+l-1,7);
                                end
                                list(ind(j),4)=mean(totalpix(1:max_jump)); %#ok<AGROW>
                                if list(ind(j),3)-list(ind(j),4)>accuracy(4)
                                    %%%%found breakage!!!!!!!!!!!!!!!!!
                                    %%%%found breakage!!!!!!!!!!!!!!!!!
                                    %%%%found breakage!!!!!!!!!!!!!!!!!
                                    
                                    % render stuff
                                    figure(fig_id);hold on;
                                    xp=traj_bin(ind_be(i)+list(ind(j),1)-1,4);
                                    yp=traj_bin(ind_be(i)+list(ind(j),1)-1,5);
                                    zp=traj_bin(ind_be(i)+list(ind(j),1)-1,7);
                                    zps=traj_bin(ind_be(i)+list(ind(j),1)-1,20);
                                    xc=vec_be(ind_time(ind_begin(k)),1);
                                    yc=vec_be(ind_time(ind_begin(k)),2);
                                    zc=totalpix_be(ind_time(ind_begin(k)));
                                    used(ind_time(ind_begin(k)))=1;
                                    zcs=strain_be(ind_time(ind_begin(k)));
                                    
                                    % draw parent end
                                    scatter3(xp,yp,zp,10,'r');
                                    % draw child beg
                                    scatter3(xc,yc,zc,20,'g','filled');
                                    % draw line parent-child
                                    plot3([xp xc],[yp yc],[zp zc],'c')
                                    
                                    %                                 if(ind_be(i)+list(ind(j),2)<ind_en(i))
                                    %                                     scatter3(traj_bin(ind_be(i)+list(ind(j),2)-1,4),traj_bin(ind_be(i)+list(ind(j),2)-1,5),list(ind(j),4),20,'g','filled');
                                    % %                                     totalpix_child=[totalpix_child;list(ind(j),4)];
                                    %                                 end
                                    
                                    % clean totalpix figure
                                    figure(100+fig_id);hold on;title('total pix')
                                    % draw parent traj
                                    plot3(traj_bin(ind_be(i):ind_en(i),4),...
                                        traj_bin(ind_be(i):ind_en(i),5),...
                                        traj_bin(ind_be(i):ind_en(i),7));
                                    
                                    break_pos=[break_pos;xp,yp,traj_bin(ind_be(i)+list(ind(j),1)-1,6),zp,zps,traj_bin(ind_be(i)+list(ind(j),1)-1,1),traj_bin(ind_be(i)+list(ind(j),1)-1,21),traj_bin(ind_be(i)+list(ind(j),1)-1,22),traj_bin(ind_be(i)+list(ind(j),1)-1,23),traj_bin(ind_be(i)+list(ind(j),1)-1,24),traj_bin(ind_be(i)+list(ind(j),1)-1,27),traj_bin(ind_be(i)+list(ind(j),1)-1,28),i,traj_bin(ind_be(i)+list(ind(j),1)-1,29),traj_bin(ind_be(i)+list(ind(j),1)-1,25),traj_bin(ind_be(i)+list(ind(j),1)-1,26)]; %#ok<AGROW>
                                    unbreak_pos_par=[unbreak_pos_par;traj_bin(ind_be(i):ind_en(i),4),traj_bin(ind_be(i):ind_en(i),5),traj_bin(ind_be(i):ind_en(i),6),traj_bin(ind_be(i):ind_en(i),7),traj_bin(ind_be(i):ind_en(i),20),traj_bin(ind_be(i):ind_en(i),1),traj_bin(ind_be(i):ind_en(i),21),traj_bin(ind_be(i):ind_en(i),22),traj_bin(ind_be(i):ind_en(i),23),traj_bin(ind_be(i):ind_en(i),24),traj_bin(ind_be(i):ind_en(i),27),traj_bin(ind_be(i):ind_en(i),28),i*ones(length(traj_bin(ind_be(i):ind_en(i),4)),1),traj_bin(ind_be(i):ind_en(i),29),traj_bin(ind_be(i):ind_en(i),25),traj_bin(ind_be(i):ind_en(i),26)]; %#ok<*AGROW>
                                    
                                    
                                    % draw child traj
                                    plot3(traj_bin(ind_be(ind_time(ind_begin(k))):ind_en(ind_time(ind_begin(k))),4),...
                                        traj_bin(ind_be(ind_time(ind_begin(k))):ind_en(ind_time(ind_begin(k))),5),...
                                        traj_bin(ind_be(ind_time(ind_begin(k))):ind_en(ind_time(ind_begin(k))),7));
                                    % draw parent end
                                    scatter3(xp,yp,zp,10,'r');
                                    % draw child beg
                                    scatter3(xc,yc,zc,20,'g','filled');
                                    % draw line parent-child
                                    plot3([xp xc],[yp yc],[zp zc],'c')
                                    
                                    
                                    unbreak_pos_child=[unbreak_pos_child;traj_bin(ind_be(ind_time(ind_begin(k))):ind_en(ind_time(ind_begin(k))),4),traj_bin(ind_be(ind_time(ind_begin(k))):ind_en(ind_time(ind_begin(k))),5),traj_bin(ind_be(ind_time(ind_begin(k))):ind_en(ind_time(ind_begin(k))),6),traj_bin(ind_be(ind_time(ind_begin(k))):ind_en(ind_time(ind_begin(k))),7),traj_bin(ind_be(ind_time(ind_begin(k))):ind_en(ind_time(ind_begin(k))),20),traj_bin(ind_be(ind_time(ind_begin(k))):ind_en(ind_time(ind_begin(k))),1),traj_bin(ind_be(ind_time(ind_begin(k))):ind_en(ind_time(ind_begin(k))),21),traj_bin(ind_be(ind_time(ind_begin(k))):ind_en(ind_time(ind_begin(k))),22),traj_bin(ind_be(ind_time(ind_begin(k))):ind_en(ind_time(ind_begin(k))),23),traj_bin(ind_be(ind_time(ind_begin(k))):ind_en(ind_time(ind_begin(k))),24),traj_bin(ind_be(ind_time(ind_begin(k))):ind_en(ind_time(ind_begin(k))),27),traj_bin(ind_be(ind_time(ind_begin(k))):ind_en(ind_time(ind_begin(k))),28),i*ones(length(traj_bin(ind_be(ind_time(ind_begin(k))):ind_en(ind_time(ind_begin(k))),4)),1),traj_bin(ind_be(ind_time(ind_begin(k))):ind_en(ind_time(ind_begin(k))),29),traj_bin(ind_be(ind_time(ind_begin(k))):ind_en(ind_time(ind_begin(k))),25),traj_bin(ind_be(ind_time(ind_begin(k))):ind_en(ind_time(ind_begin(k))),26)];
                                    start_child=[start_child;xc,yc,vec_be(ind_time(ind_begin(k)),3),zc,strain_be(ind_time(ind_begin(k))),frame_be(ind_time(ind_begin(k))),lambda1_be(ind_time(ind_begin(k))),lambda2_be(ind_time(ind_begin(k))),lambda3_be(ind_time(ind_begin(k))),gradmeas_be(ind_time(ind_begin(k))),size_meas1_be(ind_time(ind_begin(k))),size_meas2_be(ind_time(ind_begin(k))),i,grv_be(ind_time(ind_begin(k))),xpx_be(ind_time(ind_begin(k))),ypx_be(ind_time(ind_begin(k)))];
                                    
                                    %clean strain figure
                                    figure(200+fig_id);hold on;title('strain')
                                    % draw parent traj
                                    plot3(traj_bin(ind_be(i):ind_en(i),4),...
                                        traj_bin(ind_be(i):ind_en(i),5),...
                                        traj_bin(ind_be(i):ind_en(i),20));
                                    % draw child traj
                                    plot3(traj_bin(ind_be(ind_time(ind_begin(k))):ind_en(ind_time(ind_begin(k))),4),...
                                        traj_bin(ind_be(ind_time(ind_begin(k))):ind_en(ind_time(ind_begin(k))),5),...
                                        traj_bin(ind_be(ind_time(ind_begin(k))):ind_en(ind_time(ind_begin(k))),20));
                                    % draw parent end
                                    scatter3(xp,yp,zps,10,'r');
                                    % draw child beg
                                    scatter3(xc,yc,zcs,20,'g','filled');
                                    % draw line parent-child
                                    plot3([xp xc],[yp yc],[zps zcs],'c')
                                    
                                    % do stats
                                    % break pos, size, strain etc
                                    % num_child
                                    % child size
                                    
                                    
                                    
                                end
                                %%%%end of found breakage!!!!!!!!!!!!!!!!!
                                %%%%end of found breakage!!!!!!!!!!!!!!!!!
                                %%%%end of found breakage!!!!!!!!!!!!!!!!!
                            end
                        end
                    end
                end
            end
            %%this completes the story a particular potential breakage
            %             %here we can do some statistics
            %             if length(totalpix_child)>1
            %                 global_totalpix_before = [global_totalpix_before;totalpix_before];
            %                 tmp=zeros(1,10);
            %                 tmp(1,1:length(totalpix_child))=totalpix_child';
            %                 global_totalpix_child  = [global_totalpix_child;tmp];
            %                 global_breakage_pos    = [global_breakage_pos;...
            %                     [traj_bin(ind_be(i)+list(ind(j),1)-1,4) traj_bin(ind_be(i)+list(ind(j),1)-1,5) traj_bin(ind_be(i)+list(ind(j),1)-1,6)]];
            %             end
        end
    end
    
    
    %             if exist(['Result_unfiltered_v7\',foll{globi}])==0
    %             mkdir(['Result_unfiltered_v7\',foll{globi}]);
    %             end
    %
    %             break_pos = unique(break_pos, 'rows');
    %             unbreak_pos_par = unique(unbreak_pos_par,'rows');
    %             unbreak_pos_child = unique(unbreak_pos_child,'rows');
    %             start_child = unique(start_child,'rows');
    %             save(['Result_unfiltered_v7/',foll{globi},'break_pos_',ar(fig_id).name],'break_pos','unbreak_pos_par','unbreak_pos_child','start_child')
    %         %
    %
    %             if size(unbreak_pos_par,1)>0
    %             indd=find(unbreak_pos_par(:,10)<=0.2);
    %             unbreak_pos_par=unbreak_pos_par(indd,:);
    %             end
    %
    %             if size(unbreak_pos_child,1)>0
    %             indd_2=find(unbreak_pos_child(:,10)<=0.2);
    %             unbreak_pos_child=unbreak_pos_child(indd_2,:);
    %             end
    %
    %             if exist(['Result_filtered_0p2_v7\',foll{globi}])==0
    %             mkdir(['Result_filtered_0p2_v7\',foll{globi}]);
    %             end
    %             save(['Result_filtered_0p2_v7/',foll{globi},'break_pos_',ar(fig_id).name],'break_pos','unbreak_pos_par','unbreak_pos_child','start_child')
    
    
    
end
% here the loop through all trajectories is completed

% end

disp('finished')





% save stats global_totalpix_before global_totalpix_child global_breakage_pos

% figure;scatter3(global_breakage_pos(:,1),global_breakage_pos(:,2),global_breakage_pos(:,3))
%
% figure;hist(global_totalpix_before)
% xlabel('totalpix before')
%
% ind=find(global_totalpix_child>0);
% figure;hist(global_totalpix_child(ind))
% xlabel('totalpix children')
%
% si=size(global_totalpix_child);
% for i=1:si(1,1)
%    ind=find(global_totalpix_child(i,:)>0);
%    howmany(i)=length(ind);
% end
% figure;hist(howmany)
% xlabel('# children')
%
% figure;scatter(global_totalpix_before,howmany)
% xlabel('totalpix before')
% ylabel('# children')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%FUNCTIONS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




















