% function statisticsFromGridDiv2HDF_v2(basename,first,last,dataNum)
% First attempt to use HDF5 infrastructure
% 
% To write, we keep in memory or write down the hdf5.h5compound objects.

% we keep it as a script in the debugging version
dataNum = 8;

if strcmp(getenv('computername'),'IBM-LIBERZON')
    basename{8} = 'D:\0104_w\gridFilt.';
    basename{9} = 'D:\0104_p\gridFilt.';
    basename{10} = 'D:\02Nov_w\gridFilt.';
    basename{11} = 'D:\02Nov_20\gridFilt.';
end
        
first = 10020;
last  = 12979;

load(['meanField_',int2str(dataNum)]);
len = length(mea.x);

% fieldname_list = ['u,v,w,ux,uy,uz,vx,vy,vz,wx,wy,wz,lambda1,lambda2,lambda3,wws,w2',...
%     'cosss,cosww,coswl1,coswl2,coswl3,reldiv,ii,cosulambda1m,cosulambda2m,cosulambda3m,t'];
dfields = {'x','y','z','u','v','w','ux','uy','uz','vx','vy','vz','wx','wy','wz','lambda1','lambda2',...
    'lambda3','wws','w2','cosss','cosww','coswl1','coswl2','coswl3','reldiv',...
    'cosulambda1m','cosulambda2m','cosulambda3m','t'};

gridPoint = hdf5.h5compound(dfields{:});
gridData = repmat(gridPoint,len,1);


[wws, sss, coswW, ...
    coswl1,coswl2, coswl3,...
    lambda1,lambda2,lambda3, Ssq, Energy, cosS,...
    reldiv, ...
    cosuLam1, cosuLam2, cosuLam3] = deal(zeros(len,1));


hfile = ['data_',int2str(dataNum),'.h5'];
% access = 'create';
% file_id = hdf('H','open', hfile, access, 0);
% if file_id == -1
%     error('HDF hopen failed');
% end

% initialize the V interface
% status = hdf('V','start',file_id);
% if status == -1
%     error('HDF vstart failed');
% end
% end

% x,y,z are inserted once
for m = 1:len
            gridData(m).Data(1,1) = {mea.x(m)};
            gridData(m).Data(2,1) = {mea.y(m)};
            gridData(m).Data(3,1) = {mea.z(m)};
%            gridData(m).Data(find(strcmp(dfields,'ii'))) = m; % stupid,
%            the number of gridData itself is a number of the gridPoint
%
end

t = 0;
for i = first:last

    t = t + 1;

    clear f;
    fid = fopen([basename{dataNum},int2str(i)],'r');
    f = textscan(fid,'%f');
    if isempty(f), continue, end
    f = reshape(f{:},32,len).';
    fclose(fid)

    %% Fluctuating quantities
%     u = f(:,4) - mea.u;
%     v = f(:,5) - mea.v;
%     w = f(:,6) - mea.w;
%     Energy = 0.5*(u.^2 + v.^2 + w.^2);

%     ux = f(:,7) - mea.ux;
%     uy = f(:,8) - mea.uy;
%     uz = f(:,9) - mea.uz;
%     vx = f(:,10) - mea.vx;
%     vy = f(:,11) - mea.vy;
%     vz = f(:,12) - mea.vz;
%     wx = f(:,13) - mea.wx;
%     wy = f(:,14) - mea.wy;
%     wz = f(:,15) - mea.wz;
    
    
    meaFields = fieldnames(mea);
    
    for m = 1:len
        for n = 4:6
            gridData(m).Data(n,t) = {f(m,n) - mea.(meaFields{n})(m)}; 
        end
        for n = 7:15
            gridData(m).Data(n,t) = {f(m,n) - mea.(meaFields{n+1})(m)}; 
        end
    end
    

    w1 = wy - vz;
    w2 = uz - wx;
    w3 = vx - uy;
    absw = (w1.^2 + w2.^2 + w3.^2).^0.5;

    % enstrophy = w1.^2+w2.^2+w3.^2;
    s11 = 0.5*(ux + ux);
    s12 = 0.5*(uy + vx);
    s13 = 0.5*(uz + wx);
    s22 = 0.5*(vy + vy);
    s23 = 0.5*(vz + wy);
    s33 = 0.5*(wz + wz);

    wS1 = w1.*s11 + w2.*s12 + w3.*s13;
    wS2 = w1.*s12 + w2.*s22 + w3.*s23;
    wS3 = w1.*s13 + w2.*s23 + w3.*s33;
    absws = (wS1.^2 + wS2.^2 + wS3.^2).^0.5;

    reldiv = abs(ux+vy+wz)./(abs(ux) + abs(vy) + abs(wz));





    % ------------- HDF part -------------------
    % create a new vs data
    access  =  'w';
    vdata_ref  =  -1;                                         % flag to create
    vdata_id  =  hdf('VS','attach', file_id, vdata_ref, access);
    if vdata_id  ==  -1, error('HDF vdata attach failed'), end

    % sname  =  sprintf('%d', runIndx);
    status  =  hdf('VS','setname', vdata_id, 'gridPoint');
    if status  ==  -1, error('HDF vsetname failed'), end

    % assign the group the class "Trajectory"
    sclass  =  'Grid';
    status  =  hdf('VS','setclass', vdata_id, sclass);
    if status  ==  -1, error('HDF vsetclass failed'),  end

    for j  =  1 : length(dfields)                 % 32 fields
        hdf('VS','fdefine',vdata_id, dfields{j},'double',1);
    end

    status  =  hdf('VS','setfields',vdata_id,fieldname_list);
    if status  ==  -1, error('HDF VS setfields failed'), end


    status  =  hdf('VS','setinterlace',vdata_id,'full');
    if status  ==  -1, error('HDF VS setinterlace failed'), end


    for ii = 1:len
        A = [s11(ii) s12(ii) s13(ii);
            s12(ii) s22(ii) s23(ii);
            s13(ii) s23(ii) s33(ii)];

        [V,D] = eig(A);
        [D,k]  =  sort(diag(D)); 	% ascending order
        D  =  diag(D(end:-1:1)); 			% descending order
        V  =  V(:,k(end:-1:1)); 	% the same order for eigenvectors

        lambda1(ii) = D(1,1);
        lambda2(ii) = D(2,2);
        lambda3(ii) = D(3,3);

        wws(ii)  =  ([w1(ii) w2(ii) w3(ii)]*[wS1(ii) wS2(ii) wS3(ii)]');

        strain = [s11(ii), s12(ii), s13(ii) ;
            s12(ii), s22(ii), s23(ii);
            s13(ii) s23(ii) s33(ii)];
        S = -reshape((strain^2)',[1 9]);
        cosS(ii) = cosine(strain(:),S);
        % sum(strain(:)'.*S,2)./norm(strain(:))./norm(S);


        cosuLam1(ii) = cosine([u(ii) v(ii) w(ii)],[mea.V1x(ii) mea.V1y(ii) mea.V1z(ii)]);
        cosuLam2(ii) = cosine([u(ii) v(ii) w(ii)],[mea.V2x(ii) mea.V2y(ii) mea.V2z(ii)]);
        cosuLam3(ii) = cosine([u(ii) v(ii) w(ii)],[mea.V3x(ii) mea.V3y(ii) mea.V3z(ii)]);


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        coswW(ii) = cosine([w1(ii) w2(ii) w3(ii)],[wS1(ii) wS2(ii) wS3(ii)]);
        %')/(absw(ii)*absws(ii));

        coswl1(ii) = ([w1(ii) w2(ii) w3(ii)]*V(:,1))/absw(ii);
        coswl2(ii) = ([w1(ii) w2(ii) w3(ii)]*V(:,2))/absw(ii);
        coswl3(ii) = ([w1(ii) w2(ii) w3(ii)]*V(:,3))/absw(ii);

    end
end

% sss  =  -lambda1.*lambda2.*lambda3;
% Ssq  =  lambda1.^2 + lambda2.^2 + lambda3.^2;


% write the buffered data into the first vdata with full interlace mode
hdf('VS','write',vdata_id,{u',v',w',ux',uy',uz',vx',vy',vz',wx',wy',wz',...
lambda1',lambda2',lambda3',wws',Energy(1:ii,1)',...
    absw(1:ii,1)'.^2,cosS(1:ii,1)',coswWArray(1:ii,1)',coswl1(1:ii,1)',...
    coswl2(1:ii,1)',coswl3(1:ii,1)',reldiv(1:ii,1)',ii(1:ii,1)',...
    cosulambdam(1:ii,1)',cosulambdam(1:ii,2)',cosulambdam(1:ii,3)'});

status  =  hdf('VS','detach',vdata_id);
if status  ==  -1, error('HDF vsdetach failed'), end

end % for s(1,1)
end % first:last

status  =  hdf('V','end',file_id);
if status  ==  -1
    error('HDF v end failed')
end

% close the HDF file
status  =  hdf('H','close',file_id);
if status  ==  -1, error('HDF hclose failed'), end


