function [] = productionFromHDF(dataNumVec)

if nargin == 0,
    dataNumVec = [4];
end

fastcosine = @(a,b)sum(a.*b,2)./sqrt(sum(a.^2,2))./sqrt(sum(b.^2,2));

%% Statistics
%     fields = {'u','v','w',...1-3
...'ux','uy','uz','vx','vy','vz','wx','wy','wz',...4-12,
    ...'lambda1','lambda2','lambda3',...13-15,
    ... 'wws','enst',...16-17
    ...'cosss','cosww',...18-19
    ...'coswl1','coswl2','coswl3',...20-22
    ...'cosul1','cosul2','cosul3',...23-35
    ...'reldiv'};...26


%% The main loop, through various flow cases
for dataNum = dataNumVec


    hfile = ['statistics_',int2str(dataNum),'.hdf'];
    finfo = hdfinfo(hfile);
    len = [finfo.SDS.Dims.Size];
    disp('Production: ...')

    gridData = hdfread(hfile,'gridData','Index', {[1.0 1.0 1.0 ],[1.0 1.0 1.0 ],[len(1) 3.0 len(3)]});
    load(['medField_',int2str(dataNum)],'med');

    disp('done')

    disp(['Analyzing ... ',datestr(now)])
    %% Production initialization, see below the statistics for the rest
    prodfields = {'uu','uv','uw','vv','vw','ww',... % reynolds stresses
        'p1','p2','p3','p',... % 3 components and production
        };

    for j = 1:length(prodfields)
        meaCor.(prodfields{j}) = zeros(len(3),1);
    end

    %     meaCor.uu = squeeze(trimmean(gridData(:,1,:).^2,10,1));
    %     meaCor.vv = squeeze(trimmean(gridData(:,2,:).^2,10,1));
    %     meaCor.ww = squeeze(trimmean(gridData(:,3,:).^2,10,1));
    %
    %     meaCor.uv = squeeze(trimmean(gridData(:,1,:).*gridData(:,2,:),10,1));
    %     meaCor.uw = squeeze(trimmean(gridData(:,1,:).*gridData(:,3,:),10,1));
    %     meaCor.vw = squeeze(trimmean(gridData(:,2,:).*gridData(:,3,:),10,1));

    meaCor.uu = squeeze(median(gridData(:,1,:).^2,1));
    meaCor.vv = squeeze(median(gridData(:,2,:).^2,1));
    meaCor.ww = squeeze(median(gridData(:,3,:).^2,1));

    meaCor.uv = squeeze(median(gridData(:,1,:).*gridData(:,2,:),1));
    meaCor.uw = squeeze(median(gridData(:,1,:).*gridData(:,3,:),1));
    meaCor.vw = squeeze(median(gridData(:,2,:).*gridData(:,3,:),1));

    %% Production contributions :-)

    for i = 1:len(1) % for all frames


        meaCor.p1 = meaCor.p1 + squeeze(gridData(i,1,:).^2 + gridData(i,2,:).^2 + gridData(i,1,:).^2).*...
            med.lambda1(:).*(fastcosine([squeeze(gridData(i,1,:)) squeeze(gridData(i,2,:)) squeeze(gridData(i,3,:))],[med.v1x.' med.v1y.' med.v1z.']).^2);
        meaCor.p2 = meaCor.p2 + squeeze(gridData(i,1,:).^2 + gridData(i,2,:).^2 + gridData(i,1,:).^2).*...
            med.lambda2(:).*(fastcosine([squeeze(gridData(i,1,:)) squeeze(gridData(i,2,:)) squeeze(gridData(i,3,:))],[med.v2x.' med.v2y.' med.v2z.']).^2);
        meaCor.p3 = meaCor.p3 + squeeze(gridData(i,1,:).^2 + gridData(i,2,:).^2 + gridData(i,1,:).^2).*...
            med.lambda3(:).*(fastcosine([squeeze(gridData(i,1,:)) squeeze(gridData(i,2,:)) squeeze(gridData(i,3,:))],[med.v3x.' med.v3y.' med.v3z.']).^2);
    end
    % %% Averaging
    for j = 7:9 % all besides the coordinates
        meaCor.(prodfields{j}) = meaCor.(prodfields{j})/(len(1));
    end

     s11m = med.ux;
    s12m = 0.5*(med.uy + med.vx);
    s13m = 0.5*(med.uz + med.wx);
    s22m = med.vy;
    s23m = 0.5*(med.vz + med.wy);
    s33m = med.wz;
    
    meaCor.p = -1*(...
        meaCor.uu.*s11m + meaCor.uv.*s12m + meaCor.uw.*s13m + ...
        meaCor.uv.*s12m + meaCor.vv.*s22m + meaCor.vw.*s23m + ...
        meaCor.uw.*s13m + meaCor.vw.*s23m + meaCor.ww.*s33m);


    save(['meanCorField_',int2str(dataNum)],'meaCor')
    disp(['Done ... ',datestr(now)])

    % clear dataset mea meaCor

end % files loop
