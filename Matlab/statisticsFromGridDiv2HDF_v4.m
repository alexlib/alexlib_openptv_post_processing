function statisticsFromGridDiv2HDF_v4(basename,first,last,dataNum)
% Modified at: 14-Feb-2005
% $Revision: 1.0 $  $Date: 14-Feb-2005 09:57:55$
% $Revision: 12.0 $  $Date: 14-Feb-2005 21:57:55$
% - let's try the SD interface

if nargin == 0
    if strcmp(getenv('computername'),'IHW-LIBERZON')
        basename = {'I:\water\gridFilt.',...
            'I:\20ppm\gridFilt.',...
            'I:\50ppm\gridFilt.',...
            'I:\100ppm\gridFilt.',...
            'I:\watercenter\gridFilt.',...
            'I:\20ppmcenter\gridFilt.',...
            'G:\PTV\Sep2004\290904\res_2950\gridFilt.',...
            'G:\PTV\April2004\010404\res_newsoftware\gridFilt.',...
            'G:\PTV\April2004\010404\res_hompol\gridFilt.',...
            'I:\res_02Novw\gridFilt.',...
            'I:\res_02Nov20\gridFilt.',...
            'I:\res_02Nov50\gridFilt.',...
            'I:\res_02Nov100\gridFilt.'};
    elseif strcmp(getenv('computername'),'IHW-LIBERZON-2')
        basename = {'E:\TKEproduction\water\gridFilt.',...
            'E:\TKEproduction\20ppm\gridFilt.',...
            'E:\TKEproduction\50ppm\gridFilt.',...
            'E:\TKEproduction\100ppm\gridFilt.',...
            'E:\TKEproduction\watercenter\gridFilt.',...
            'E:\TKEproduction\20ppmcenter\gridFilt.',...
            'G:\PTV\Sep2004\290904\res_2950\gridFilt.',...
            'E:\TKEproduction\res_newsoftware\gridFilt.',...
            'E:\TKEproduction\res_hompol\gridFilt.',...
            'E:\TKEproduction\res_02Novw\gridFilt.',...
            'E:\TKEproduction\res_02Nov20\gridFilt.',...
            'E:\TKEproduction\res_02Nov50\gridFilt.',...
            'E:\TKEproduction\res_02Nov100\gridFilt.'};
    elseif strcmp(getenv('computername'),'IBM-LIBERZON')
        basename{8} = 'D:\0104_w\gridFilt.';
        basename{9} = 'D:\0104_p\gridFilt.';
        basename{10} = 'D:\02Nov_w\gridFilt.';
        basename{11} = 'D:\02Nov_20\gridFilt.';
    end

    % dataNumVec = (8:11);
    first = 10020;
    last  = 10025;

    dataNum = 11;
end
%%
load(['meanField_',int2str(dataNum)]);
len = length(mea.x);

%%
% % Define Field names
fields = {'u','v','w','ux','uy','uz','vx','vy','vz','wx','wy','wz',...
    'lambda1','lambda2','lambda3','wws','enst',...
    'cosss','cosww','coswl1','coswl2','coswl3',...
    'cosul1','cosul2','cosul3','reldiv'};

for j = 1:length(fields) % initialize the structure
    stat.(fields{j}) = zeros(len,1);
end


fieldname_list = sprintf('%s,',fields{:});
fieldname_list = fieldname_list(1:end-1);



%% Initialize HDF file and contents
hfile = ['statistics_',int2str(dataNum),'.hdf'];


sd_id = hdf('SD','start', hfile, 'create');
if sd_id == -1, error('HDF sdstart failed'); end

hdf('SD','setattr',sd_id,'fields',fieldname_list);

dimsizes = [(last-first+1) length(fields) len ]; 
rank = length(dimsizes);
start = zeros(1,rank);
stride = ones(1,rank);

sds_id = hdf('SD','create',sd_id,'gridData','float64',rank,dimsizes);

%% Main loop through the files (frames)
t = 0;
for i = first:last
    t = t + 1;

    clear f;

    fid = fopen([basename{dataNum},int2str(i)],'r');
    f = textscan(fid,'%f');
    fclose(fid);
    f = reshape(f{:},32,[]).';
    if isempty(f), continue, end

%% Fluctuating quantities

    for j = 1:12 % all the rest
        stat.(fields{j}) = f(:,j+3) - mea.(fields{j});
    end


    % stat.Energy(CounterN,1) = 0.5*(u(ii)^2+v(ii)^2+w(ii)^2);

    stat.reldiv = abs(stat.ux + stat.vy + stat.wz)./(abs(stat.ux) + abs(stat.vy) + abs(stat.wz));

    w1 = stat.wy - stat.vz;
    w2 = stat.uz - stat.wx;
    w3 = stat.vx - stat.uy;
    stat.enst = (w1.^2 + w2.^2 + w3.^2);

    s11 = stat.ux;
    s12 = 0.5*(stat.uy + stat.vx);
    s13 = 0.5*(stat.uz + stat.wx);
    s22 = stat.vy;
    s23 = 0.5*(stat.vz + stat.wy);
    s33 = stat.wz;

    wS1 = w1.*s11 + w2.*s12 + w3.*s13;
    wS2 = w1.*s12 + w2.*s22 + w3.*s23;
    wS3 = w1.*s13 + w2.*s23 + w3.*s33;
    


    stat.cosul1 = abs(cosine([stat.u, stat.v, stat.w],[mea.v1x, mea.v1y mea.v1z]));
    stat.cosul2 = abs(cosine([stat.u, stat.v, stat.w],[mea.v2x, mea.v2y mea.v2z]));
    stat.cosul3 = abs(cosine([stat.u, stat.v, stat.w],[mea.v3x, mea.v3y mea.v3z]));

    stat.cosww = cosine([w1 w2 w3],[wS1 wS2 wS3]);

%     velSq = stat.u.^2 + stat.v.^2 + stat.w.^2;
%     stat.p = -1* (velSq.*mea.lambda1.*stat.cosul1.^2 + velSq.*mea.lambda2.*stat.cosul2.^2 + velSq.*mea.lambda3.*stat.cosul3.^2);
%%
    for ii = 1:len
        A = [s11(ii) s12(ii) s13(ii);
            s12(ii) s22(ii) s23(ii);
            s13(ii) s23(ii) s33(ii)];

        [V,D] = eig(A);
        [D,k]  =  sort(diag(D)); 	% ascending order
        D  =  diag(D(end:-1:1)); 			% descending order
        V  =  V(:,k(end:-1:1)); 	% the same order for eigenvectors


        stat.lambda1(ii) = D(1,1);
        stat.lambda2(ii) = D(2,2);
        stat.lambda3(ii) = D(3,3);


        % stat.sss(ii)  =  -D(1,1)*D(2,2)*D(3,3);
        % stat.sSq(ii)  =  D(1,1)^2+D(2,2)^2+D(3,3)^2;

        stat.wws(ii)  =  ([w1(ii) w2(ii) w3(ii)]*[wS1(ii) wS2(ii) wS3(ii)]');

        S = -reshape((A^2)',[1 9]);
        stat.cosss(ii) = cosine(A(:)',S);


        stat.coswl1(ii) = abs(cosine([w1(ii) w2(ii) w3(ii)],V(:,1).'));
        stat.coswl2(ii) = abs(cosine([w1(ii) w2(ii) w3(ii)],V(:,2).'));
        stat.coswl3(ii) = abs(cosine([w1(ii) w2(ii) w3(ii)],V(:,3).'));

    end % end of for ii = 1:len

%% Write each file
start(1) = t - 1;    
edges = [1 dimsizes(2:3)];
data = [];
for k = 1:length(fields), data = cat(2,data,stat.(fields{k})); end

status = hdf('SD','writedata', sds_id, start, stride, edges, data);
if status == -1
  error('HDF SDS write failed');
end 

end % i = first:last


%% end sd interface access
status  =  hdf('SD','endaccess',sds_id);
if status  ==  -1, error('HDF sds endacess failed'), end

%% close the HDF file
status  =  hdf('SD','end',sd_id);
if status  ==  -1, error('HDF SD end failed'), end




