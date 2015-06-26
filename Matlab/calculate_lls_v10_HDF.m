% September 19, 2004
% The new data is of enourmous size, the usual way of saving all B 
% and relative matrices is not working anymore. we have to split it
% into separate frames and save into HDF per frame. HDF is open all the
% time and more and more lines are added during the run.
% solves only the MEMORY problem
% 


if hdf('H','ishdf', [hfile(1:end-4),'_LPQR.hdf'])
    error('File exists, delete it manually') % this option should be changed, if somebody wants to append 28.11.03
else
    access = 'create';
    % open new or existing hfile
    file_id = hdf('H','open', [hfile(1:end-4),'_LPQR.hdf'], access, 0);
    if file_id == -1
        error('HDF hopen failed');
    end
    
    % initialize the V interface
    status = hdf('V','start',file_id);
    if status == -1
        error('HDF vstart failed');
    end
    
    % initialize the SD interface
    sd_id = hdf('SD','start', hfile, 'write');
    if sd_id == -1
        error('HDF sdstart failed');
    end
    
    %     % assign optional top level SD attributes
    %     astr.NumPoints       =  int2str(numPoints);
    %     astr.numTraj         = int2str(numTraj);
    % %     astr.trajLen         = cat(1,traj.trajLen);
    %     
    %     
    %     afields = fieldnames(astr);
    %     for j = 1 : length(afields)
    %         status = hdf('SD','setattr', sd_id, afields{j}, astr.(afields{j}));
    %         if status == -1
    %             error('HDF sdsetattr failed')
    %         end
    %     end
    
    
    
end

fieldname_list = 'b,pb,qb,rb,db,vb,w,pw,qw,rw,dw,vw,t,pt,qt,rt,dt,vt,l,lt,zeta,age,coslwls,coslwlu,cosllambda1,cosllambda2,cosllambda3,d,v';
dfields = {'b','pb','qb','rb','db','vb','w','pw','qw','rw','dw','vw',...
        't','pt','qt','rt','dt','vt','l','lt','zeta','age',...
        'coslwls','coslwlu','cosllambda1','cosllambda2','cosllambda3','d','v'};
order = {9,1,1,1,3,9,9,1,1,1,3,9,9,1,1,1,3,9,N*3,N*3,N,1,N,N,N,N,N,3,9};

access  =  'w';
vdata_ref  =  -1;                                         
vdata_id = hdf('VS','attach', file_id, vdata_ref, access);
if vdata_id == -1
    error('HDF vdata attach failed')
end
sname  =  sprintf('%s', hfile(1:end-4));
status  =  hdf('VS','setname', vdata_id, sname);
if status  ==  -1
    error('HDF vsetname failed')
end

% assign the group the class "Trajectory"
sclass  =  'MaterialData';
status  =  hdf('VS','setclass', vdata_id, sclass);
if status  ==  -1
    error('HDF vsetclass failed')
end
for j  =  1 : length(dfields)
    status  =  hdf('VS','fdefine',vdata_id, dfields{j},'double',order{j});
end
status  =  hdf('VS','setfields',vdata_id,fieldname_list);
if status == -1
    error('HDF VS setfields failed')
end
status  =  hdf('VS','setinterlace',vdata_id,'full');
if status == -1
    error('HDF VS setinterlace failed')
end


% Runga-Kutta 3rd order integration of B tensor
dt=1/60/tau_eta;
% dt=1/60; % /tau_eta;
for i = 1:numTraj                                           % for all trajectories loop
    
    fprintf(2,'%6d',i);
    
    
    % numPoints = trajLen(i);
    
    
    [PB,QB,RB,PAB,QAB,RAB,...
            PT,QT,RT,PAT,QAT,RAT,...
            PW,QW,RW,PAW,QAW,RAW]       = deal(zeros(trajLen(i),1));
    [DB,DW,DT,D]                        = deal(zeros([1 3 trajLen(i)]));
    [B,AB,T,AT,VB,W,AW,VW,VT,V]         = deal(zeros([3 3 trajLen(i)]));
    
    [L,Lt,Wls,Wlu]                              = deal(zeros(N,3,trajLen(i)));
    [zeta,...
            cos_l_lambda1,...
            cos_l_lambda2,...
            cos_l_lambda3,...
        coslWls,coslWlu]                      = deal(zeros(N,1,trajLen(i)));
    
    L0 = rand(N,3)- 0.5;
    L0 =  L0./repmat(sqrt(sum(L0.^2,2)),[1 3]);
    L(1:N,1:3,1) = L0;
    B(:,:,1) = eye(3);
    
    tmp = tau_eta * ...
        [dudx(start_points(i)),dudy(start_points(i)),dudz(start_points(i));...
            dvdx(start_points(i)),dvdy(start_points(i)),dvdz(start_points(i));...
            dwdx(start_points(i)),dwdy(start_points(i)),dwdz(start_points(i))];           % !!! July 19, non-dimensional treatment.
    
    
    
    Wls(1:N,1:3,1) = [ L(1:N,1:3,1)*[dudx(start_points(i));0.5*(dudy(start_points(i))+dvdx(start_points(i)));0.5*(dudz(start_points(i))+dwdx(start_points(i)))],...
            L(:,:,1)*[0.5*(dudy(start_points(i))+dvdx(start_points(i)));dvdy(start_points(i));0.5*(dwdy(start_points(i))+dvdz(start_points(i)))],...
            L(:,:,1)*[0.5*(dudz(start_points(i))+dwdx(start_points(i)));0.5*(dwdy(start_points(i))+dvdz(start_points(i)));dwdz(start_points(i))]];
    
    Wlu(1:N,1:3,1) = ...
        [L(:,:,1)*[dudx(start_points(i));dudy(start_points(i));dudz(start_points(i))], ...
            L(:,:,1)*[dvdx(start_points(i));dvdy(start_points(i));dvdz(start_points(i))],...
            L(:,:,1)*[dwdx(start_points(i));dwdy(start_points(i));dwdz(start_points(i))]];
    
    % LLS(1:N,1,j) = dot(L(1:N,1:3,j),Wls(1:N,1:3,j),2);
    
    %     Omega = dudx(k,:) - s(k,:);
    %     Wlw(1:N,1:3) = [L(:,:,k)*Omega(1:3)', L(:,:,k)*Omega(4:6)', L(:,:,k)*Omega(7:9)'];
    
    
    coslWls(1:N,1,1)  = cosine(L(1:N,1:3,1),Wls(1:N,1:3,1));
    coslWlu(1:N,1,1)  = cosine(L(1:N,1:3,1),Wlu(1:N,1:3,1));
    %     coslWlw(1:N,k)  = cosine(L(1:N,1:3,k),Wlw(1:N,1:3));
    
    
    
    
    [v,d] = eig(0.5*(tmp + tmp')/tau_eta);
    [d,m] = sort(diag(d)); 	% ascending order
    D(:,:,1) = d(end:-1:1); 	% descending order		
    V(:,:,1) = v(:,m(end:-1:1)); 	% the same order for eigenvectors 
    
    
    
    
    for n = 1:N
        cos_l_lambda1(n,1,1) = cosine(squeeze(L(n,1:3,1)),squeeze(V(1:3,1,1))');
        cos_l_lambda2(n,1,1) = cosine(squeeze(L(n,1:3,1)),squeeze(V(1:3,2,1))'); 
        cos_l_lambda3(n,1,1) = cosine(squeeze(L(n,1:3,1)),squeeze(V(1:3,3,1))'); 
    end 
    
    
    [PB(1),QB(1),RB(1)] = invariants(B(:,:,1));
    
    % Deviator of B
    AB(:,:,1) = deviator(B(:,:,1));
    [PAB(1),QAB(1),RAB(1)] = invariants(AB(:,:,1));
    
    % T_{kn} tensor = B_{ik} B_{jk} S_{ij}
    T(:,:,1) = B(:,:,1)*B(:,:,1)*.5*(tmp + tmp');
    [PT(1),QT(1),RT(1)] = invariants(T(:,:,1));
    
    [v,d] = eig(T(:,:,1));
    [junk,m] = sort(abs(diag(d))); 	% ascending order
    DT(:,:,1) = diag(d(m(end:-1:1),m(end:-1:1))); 	% descending order		
    VT(:,:,1) = v(:,m(end:-1:1)); 	% the same order for eigenvectors 
    
    
    
    % Deviator of T
    AT(:,:,1) = deviator(T(:,:,1));
    [PAT(1),QAT(1),RAT(1)] = invariants(AT(:,:,1));
    
    
    [v,d] = eig(B(:,:,1));
    [junk,m] = sort(abs(diag(d))); 	% ascending order
    %         m = m + 1;
    DB(:,:,1) = diag(d(m(end:-1:1),m(end:-1:1))); 	% descending order		
    VB(:,:,1) = v(:,m(end:-1:1)); 	% the same order for eigenvectors 
    
    W(:,:,1) = B(:,:,1)*B(:,:,1)';
    [v,d] = eig(W(:,:,1));
    [d,m] = sort(diag(d)); 	% ascending order
    %         m = m + 1;
    DW(:,:,1) = d(end:-1:1); 	% descending order		
    VW(:,:,1) = v(:,m(end:-1:1)); 	% the same order for eigenvectors 
    [PW(1),QW(1),RW(1)] = invariants(W(:,:,1));
    
    % Deviator of W
    AW(:,:,1) = deviator(W(:,:,1));
    [PAW(1),QAW(1),RAW(1)] = invariants(AW(:,:,1));
    
    
    for j = 2:trajLen(i)                                % for all points in trajectory loop, ...
        k = start_points(i) + j - 1;
        %         tmp =  reshape(dudx(k,:),[3 3]);           % !!! July 19, non-dimensional treatment.
        %         tmp = tau_eta * reshape(dudx(k,:),[3 3])';           % !!! July 19, non-dimensional treatment.
        tmp = tau_eta * [dudx(k),dudy(k),dudz(k);...
                dvdx(k),dvdy(k),dvdz(k);...
                dwdx(k),dwdy(k),dwdz(k)];
        
        B1 = B(:,:,j-1) + dt * tmp*B(:,:,j-1);
        B2 = 3/4*B(:,:,j-1) + 1/4*B1 + 1/4* dt * tmp * B1;
        B(:,:,j) = 1/3*B(:,:,j-1) + 2/3*B2 + 2/3* dt * tmp * B2;
        
        L(1:N,1:3,j) = L(1:N,1:3,1)*B(:,:,j)';
        Lt(1:N,1:3,j) = L(1:N,1:3,j)./repmat(sqrt(sum(L(1:N,1:3,j).^2,2)),[1 3]);
        
        zeta(1:N,1,j) = diag((0.5*(tmp + tmp')*Lt(1:N,1:3,j)')'*Lt(1:N,1:3,j)');
        
        
        Wls(1:N,1:3,j) = [L(1:N,1:3,j)*[dudx(k);0.5*(dudy(k)+dvdx(k));0.5*(dudz(k)+dwdx(k))],...
                L(:,:,j)*[0.5*(dudy(k)+dvdx(k));dvdy(k);0.5*(dwdy(k)+dvdz(k))],...
                L(:,:,j)*[0.5*(dudz(k)+dwdx(k));0.5*(dwdy(k)+dvdz(k));dwdz(k)]];
        
        Wlu(1:N,1:3,j) = [L(:,:,j)*[dudx(k);dudy(k);dudz(k)], ...
                L(:,:,j)*[dvdx(k);dvdy(k);dvdz(k)],...
                L(:,:,j)*[dwdx(k);dwdy(k);dwdz(k)]];
        
        % LLS(1:N,1,j) = dot(L(1:N,1:3,j),Wls(1:N,1:3,j),2);
        
        %     Omega = dudx(k,:) - s(k,:);
        %     Wlw(1:N,1:3) = [L(:,:,k)*Omega(1:3)', L(:,:,k)*Omega(4:6)', L(:,:,k)*Omega(7:9)'];
        
        
        coslWls(1:N,1,j)  = cosine(L(1:N,1:3,j),Wls(1:N,1:3,j));
        coslWlu(1:N,1,j)  = cosine(L(1:N,1:3,j),Wlu(1:N,1:3,j));
        %     coslWlw(1:N,k)  = cosine(L(1:N,1:3,k),Wlw(1:N,1:3));
        
        
        
        
        [v,d] = eig(0.5*(tmp + tmp')/tau_eta);
        [d,m] = sort(diag(d)); 	% ascending order
        D(:,:,j) = d(end:-1:1); 	% descending order		
        V(:,:,j) = v(:,m(end:-1:1)); 	% the same order for eigenvectors 
        
        
        
        
        for n = 1:N
            cos_l_lambda1(n,1,j) = cosine(squeeze(L(n,1:3,j)),squeeze(V(1:3,1,j))');
            cos_l_lambda2(n,1,j) = cosine(squeeze(L(n,1:3,j)),squeeze(V(1:3,2,j))'); 
            cos_l_lambda3(n,1,j) = cosine(squeeze(L(n,1:3,j)),squeeze(V(1:3,3,j))'); 
        end 
        
        
        
        
        [PB(j),QB(j),RB(j)] = invariants(B(:,:,j));
        
        % Deviator of B
        AB(:,:,j) = deviator(B(:,:,j));
        [PAB(j),QAB(j),RAB(j)] = invariants(AB(:,:,j));
        
        % T_{jn} tensor = B_{ij} B_{jj} S_{ij}
        T(:,:,j) = B(:,:,j)'*(0.5*(tmp + tmp')*B(:,:,j));
        [PT(j),QT(j),RT(j)] = invariants(T(:,:,j));
        [v,d] = eig(T(:,:,j));
        [junj,m] = sort(abs(diag(d))); 	% ascending order
        DT(:,:,j) = diag(d(m(end:-1:1),m(end:-1:1))); 	% descending order		
        VT(:,:,j) = v(:,m(end:-1:1)); 	% the same order for eigenvectors 
        
        
        % Deviator of T
        AT(:,:,j) = deviator(T(:,:,j));
        [PAT(j),QAT(j),RAT(j)] = invariants(AT(:,:,j));
        
        [v,d] = eig(B(:,:,j));
        [junj,m] = sort(abs(diag(d))); 	% ascending order
        %         m = m + 1;
        DB(:,:,j) = diag(d(m(end:-1:1),m(end:-1:1))); 	% descending order		
        VB(:,:,j) = v(:,m(end:-1:1)); 	% the same order for eigenvectors 
        
        W(:,:,j) = B(:,:,j)*B(:,:,j)';
        [v,d] = eig(W(:,:,j));
        [d,m] = sort(diag(d)); 	% ascending order
        %         m = m + 1;
        DW(:,:,j) = d(end:-1:1); 	% descending order		
        VW(:,:,j) = v(:,m(end:-1:1)); 	% the same order for eigenvectors 
        [PW(j),QW(j),RW(j)] = invariants(W(:,:,j));
        
        % Deviator of W
        AW(:,:,j) = deviator(W(:,:,j));
        [PAW(j),QAW(j),RAW(j)] = invariants(AW(:,:,j));
        
    end % for j
    
    
    num_of_records = hdf('VS','write',vdata_id,...
        {reshape(B,9,[]),reshape(PB,1,[]),reshape(QB,1,[]),reshape(RB,1,[]),reshape(DB,3,[]),reshape(VB,9,[]),...
            reshape(W,9,[]),reshape(PW,1,[]),reshape(QW,1,[]),reshape(RW,1,[]),reshape(DW,3,[]),reshape(VW,9,[]),...
            reshape(T,9,[]),reshape(PT,1,[]),reshape(QT,1,[]),reshape(RT,1,[]),reshape(DT,3,[]),reshape(VT,9,[]),...
            reshape(L,N*3,[]),reshape(Lt,N*3,[]),reshape(zeta,N,[]),(0:trajLen(i)-1),...
            reshape(coslWls,N,[]),reshape(coslWlu,N,[]),...
            reshape(cos_l_lambda1,N,[]),reshape(cos_l_lambda2,N,[]),reshape(cos_l_lambda3,N,[]),...
            reshape(D,3,[]),reshape(V,9,[])});
    
    
    fprintf(2,'\b\b\b\b\b\b');
end % for i
status  =  hdf('VS','detach',vdata_id);
if status == -1
    error('HDF vsdetach failed')
end
% end vgroup interface access
status = hdf('V','end',file_id);
if status == -1
    error('HDF v end failed')
end

% end vgroup interface access
status = hdf('SD','end',sd_id); 
if status == -1
    error('HDF sd end failed')
end

% close the HDF file
status = hdf('H','close',file_id);
if status == -1
    error('HDF hclose failed')
end
