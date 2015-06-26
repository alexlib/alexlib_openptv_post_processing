cd('/Volumes/MY PASSPORT/Documents/MATLAB/HDF/')
% hfile =  '0104_water.hdf'
hfile =  '0104_hompol.hdf'


%% reading from HDF4 file
[traj,attr] = readTrajHDF_v9(hfile,'ax',[],'ay',[],'az',[],'w1',[],'w2',[],'w3',[],...
    's11',[],'s12',[],'s13',[],'s22',[],'s23',[],'s33',[],'t',[],'minlength',20);


% % Rearrange the derivatives into tensor per each point
% dudx  =   cat(1,traj.s11);
% dudy  =   cat(1,traj.s12) - 0.5 * cat(1,traj.w3);
% dudz  =   cat(1,traj.s13) + 0.5 * cat(1,traj.w2);
% dvdx  =   cat(1,traj.s12) + 0.5 * cat(1,traj.w3);
% dvdy  =   cat(1,traj.s22);
% dvdz  =   cat(1,traj.s23) - 0.5 * cat(1,traj.w1);
% dwdx  =   cat(1,traj.s13) - 0.5 * cat(1,traj.w2);
% dwdy  =   cat(1,traj.s23) + 0.5 * cat(1,traj.w1);
% dwdz  =   cat(1,traj.s33);


numTraj = length(traj)
numPoints = length(cat(1,traj.s11))
trajLen = cat(1,traj.trajLen);




k = 0;
good = zeros(numTraj,1);
fn = fieldnames(traj);
for i=1:numTraj
    absdiv = abs(traj(i).s11 +  traj(i).s22 + traj(i).s33);
    reldiv = absdiv./(abs(traj(i).s11) +  abs(traj(i).s22) + abs(traj(i).s33));
    ind = find(reldiv < .1 & absdiv < .1);
    ind1 = diff([ind;Inf]);
    ind2 = [0; find(ind1 > 1)];
    ind3 = diff(ind2);
    [r,c] = max(ind3);
    ind = ind(ind2(c)+1:ind2(c+1));
    if length(ind) > 50
        k = k + 1;
        good(k) = i;
        for j = 1:length(fn)- 1
            %             if length(traj(i).(fn{j}) > length(ind))
            traj(i).(fn{j}) = traj(i).(fn{j})(ind);
            %             end
        end
        traj(i).trajLen = length(ind);
    end
end

good = good(1:k);

traj = traj(good);

% % Rearrange the derivatives into tensor per each point
% dudx  =   cat(1,traj.s11);
% dudy  =   cat(1,traj.s12) - 0.5 * cat(1,traj.w3);
% dudz  =   cat(1,traj.s13) + 0.5 * cat(1,traj.w2);
% dvdx  =   cat(1,traj.s12) + 0.5 * cat(1,traj.w3);
% dvdy  =   cat(1,traj.s22);
% dvdz  =   cat(1,traj.s23) - 0.5 * cat(1,traj.w1);
% dwdx  =   cat(1,traj.s13) - 0.5 * cat(1,traj.w2);
% dwdy  =   cat(1,traj.s23) + 0.5 * cat(1,traj.w1);
% dwdz  =   cat(1,traj.s33);


numTraj = length(traj)
numPoints = length(cat(1,traj.s11))
trajLen = cat(1,traj.trajLen);
trajnum = unique(cat(1,traj.trajnum));

save([hfile(1:end-4),'_good'],'good','numTraj','numPoints','trajnum','trajLen');

% clear traj

maxLen    = max(trajLen)
minLen    = min(trajLen)
meanLen    = mean(trajLen)

start_points = [1; cumsum(trajLen(1:end-1))+1];     % index of the first trajectory points


D = [];
V = [];

for i = 1:length(traj)
    for j = 1:length(traj(i).s11)         % for all time points
        %     j = inddiv(
        s = [traj(i).s11(j), traj(i).s12(j), traj(i).s13(j);
            traj(i).s12(j),traj(i).s22(j),traj(i).s23(j);
            traj(i).s13(j),traj(i).s23(j),traj(i).s33(j)];
        [v,d] = eig(s);
        [d,k] = sort(diag(d)); 	% ascending order
        %     m = m + 1;
        D = cat(3,D,d(end:-1:1)); 	% descending order
        V = cat(3,V,v(:,k(end:-1:1))); 	% the same order for eigenvectors
    end
end


% -------------------------------------------------------------------------
    omega = [cat(1,traj.w1), cat(1,traj.w2), cat(1,traj.w3)];


% -------------------------------------------------------------------------
% cos(w,\lambda_i)
cos_wlambda1 = cosine(omega(:,:),squeeze(V(1:3,1,:))');
cos_wlambda2 = cosine(omega(:,:),squeeze(V(1:3,2,:))');
cos_wlambda3 = cosine(omega(:,:),squeeze(V(1:3,3,:))');


[n1,x1] = nhist(cos_wlambda1,20);
[n2,x2] = nhist(cos_wlambda2,20);
[n3,x3] = nhist(cos_wlambda3,20);
% hf = figure, hold on
% plot(x1,n1,'ro',x2,n2,'rs',x3,n3,'r<','linewidth',1);
% legend('cos(\omega,\lambda_1)','cos(\omega,\lambda_2)','cos(\omega,\lambda_3)')
% for i=1:3
%     eval(sprintf('fnplt(csapi(x%d,n%d),''%s:'',[],0.5)',i,i,'r'))
% end

% xlabel('cos(\omega,\lambda_i)')
% saveas(hf,['pdf_cos_w_lambda_all_times',hfile(1:end-4)],'fig')


dt = 1/60;
agebins = [-0.1,1,2,3,4,6,8]*0.2/dt

age = [];
for i = 1:length(traj)
    age = cat(2,age,0:length(traj(i).s11)-1);
end

ind = bindex(age,agebins,1);




colors = {'m','r','g','b','k'};
linestyles = {':','-','--','-.'};
markers = {'v','','o','s','x','+','<','^','>','d','p','*'};

for i = 1:length(agebins)
    style{i} = [colors{mod(i,length(colors))+1},linestyles{mod(i,length(linestyles))+1},markers{mod(i,length(markers))+1}];
end

hf = figure,
hold on
plot(x1,n1,'ro','linewidth',2);

[n,x] = deal(zeros(20,length(agebins)));
for i = 1:max(ind)
    tmp = cos_wlambda1(ind == i);
    if ~isempty(tmp)
        nhist(tmp(:),21,style{i});
    end
end
hold off
% legend('\tau_{\eta} \leq 1','1 \leq \tau_{\eta} \leq 2','2 \leq \tau_{\eta} \leq 4','4 \leq \tau_{\eta} \leq 6','6 \leq \tau_{\eta}')
xlabel('cos(\omega,\lambda)')
legend('all','1','2','3','4','6','8')
saveas(hf,['pdf_coswlambda1_vs_t',hfile(1:end-4)],'jpg')


hf = figure,
hold on
plot(x2,n2,'rs','linewidth',2);

[n,x] = deal(zeros(20,length(agebins)));
for i = 1:max(ind)
    tmp = cos_wlambda2(ind == i);
    if ~isempty(tmp)
        nhist(tmp(:),21,style{i});
    end
end
hold off
% legend('\tau_{\eta} \leq 1','1 \leq \tau_{\eta} \leq 2','2 \leq \tau_{\eta} \leq 4','4 \leq \tau_{\eta} \leq 6','6 \leq \tau_{\eta}')
xlabel('cos(\omega,\lambda)')
legend('all','1','2','3','4','6','8')
saveas(hf,['pdf_coswlambda2_vs_t',hfile(1:end-4)],'jpg')

hf = figure,
hold on
plot(x3,n3,'r<','linewidth',2);

[n,x] = deal(zeros(20,length(agebins)));
for i = 1:max(ind)
    tmp = cos_wlambda3(ind == i);
    if ~isempty(tmp)
        nhist(tmp(:),21,style{i});
    end
end
hold off
% legend('\tau_{\eta} \leq 1','1 \leq \tau_{\eta} \leq 2','2 \leq \tau_{\eta} \leq 4','4 \leq \tau_{\eta} \leq 6','6 \leq \tau_{\eta}')
xlabel('cos(\omega,\lambda)')
legend('all','1','2','3','4','6','8')
saveas(hf,['pdf_coswlambda3_vs_t',hfile(1:end-4)],'jpg')


cos_wlambda11 = [];
cos_wlambda21 = [];
cos_wlambda31 = [];

for i = 1:length(traj)
    nPoints = length(traj(i).s11);
    
    omega = [traj(i).w1, traj(i).w2, traj(i).w3];
    
    D = repmat(NaN,[1 3 nPoints]);
    V = repmat(NaN,[3 3 nPoints]);
    
    for j = 1:nPoints         % for all time points
        %     j = inddiv(
        s = [traj(i).s11(j), traj(i).s12(j), traj(i).s13(j);
            traj(i).s12(j),traj(i).s22(j),traj(i).s23(j);
            traj(i).s13(j),traj(i).s23(j),traj(i).s33(j)];
        
        [v,d] = eig(s);
        [d,k] = sort(diag(d)); 	% ascending order
        %     m = m + 1;
        D(:,:,j) = d(end:-1:1); 	% descending order
        V(:,:,j) = v(:,k(end:-1:1)); 	% the same order for eigenvectors
        
 
    end
    cos_wlambda11 = cat(1,cos_wlambda11,cosine(omega(:,:),repmat(V(1:3,1,1)',nPoints,1)));
    cos_wlambda21 = cat(1,cos_wlambda21,cosine(omega(:,:),repmat(V(1:3,2,1)',nPoints,1)));
    cos_wlambda31 = cat(1,cos_wlambda31,cosine(omega(:,:),repmat(V(1:3,3,1)',nPoints,1)));
    
end
    

hf = figure,
hold on
plot(x1,n1,'ro','linewidth',2);

[n,x] = deal(zeros(20,length(agebins)));
for i = 1:max(ind)
    tmp = cos_wlambda11(ind == i);
    if ~isempty(tmp)
        nhist(tmp(:),21,style{i});
    end
end
hold off
% legend('\tau_{\eta} \leq 1','1 \leq \tau_{\eta} \leq 2','2 \leq \tau_{\eta} \leq 4','4 \leq \tau_{\eta} \leq 6','6 \leq \tau_{\eta}')
xlabel('cos(\omega,\lambda)')
legend('all','1','2','3','4','6','8')
saveas(hf,['pdf_coswlambda11',hfile(1:end-4)],'jpg')


hf = figure,
hold on
plot(x2,n2,'rs','linewidth',2);

[n,x] = deal(zeros(20,length(agebins)));
for i = 1:max(ind)
    tmp = cos_wlambda21(ind == i);
    if ~isempty(tmp)
        nhist(tmp(:),21,style{i});
    end
end
hold off
% legend('\tau_{\eta} \leq 1','1 \leq \tau_{\eta} \leq 2','2 \leq \tau_{\eta} \leq 4','4 \leq \tau_{\eta} \leq 6','6 \leq \tau_{\eta}')
xlabel('cos(\omega,\lambda)')
legend('all','1','2','3','4','6','8')
saveas(hf,['pdf_coswlambda21',hfile(1:end-4)],'jpg')

hf = figure,
hold on
plot(x3,n3,'r<','linewidth',2);

[n,x] = deal(zeros(20,length(agebins)));
for i = 1:max(ind)
    tmp = cos_wlambda31(ind == i);
    if ~isempty(tmp)
        nhist(tmp(:),21,style{i});
    end
end
hold off
% legend('\tau_{\eta} \leq 1','1 \leq \tau_{\eta} \leq 2','2 \leq \tau_{\eta} \leq 4','4 \leq \tau_{\eta} \leq 6','6 \leq \tau_{\eta}')
xlabel('cos(\omega,\lambda)')
legend('all','1','2','3','4','6','8')
saveas(hf,['pdf_coswlambda3',hfile(1:end-4)],'jpg')
