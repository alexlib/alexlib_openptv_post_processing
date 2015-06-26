% hfile =  '2604_pol.hdf'
% hfile =  '2604_pol30min.hdf'
% hfile =  '2604_water.hdf'
 hfile =  './HDF/0104_water.hdf' % 01.03.06, for wws/lls for Arkady
% hfile =  '0104_hompol.hdf'
% hfile =  'mai10.hdf'
% hfile =  '3006_water.hdf'
% hfile =  '0207_water.hdf'
% hfile =  '0207_pol.hdf'
% hfile = 'Oct03_75.hdf'
% hfile = 'Oct03_75_pol.hdf'
% hfile = '1203_water'
% hfile = '1203_pol'
% hfile = '2308.hdf'
% hfile = '100ppm.hdf'
% hfile = '0809_water.hdf'
% hfile = '0809_50.hdf'
% hfile = '0809_100.hdf'
% hfile = '2909_water.hdf'
% hfile = '2909_50.hdf'
% hfile = '2909_20.hdf'
% hfile = '2909_100.hdf'
%  hfile = '0211_w.hdf'
% hfile = '0211_20.hdf'
% hfile = '0211_50.hdf'
% hfile = '0211_100.hdf'

runtime = datestr(now)
diary(['analysisHDF_ver0804_',strrep(strrep(runtime,' ','_'),':','_'),'.txt'])

[pathstr,namestr] = fileparts(hfile);

% hfile = '0809_water.hdf', analysisHDF_ver0804, hfile = '100ppm.hdf', analysisHDF_ver0804, clear all, hfile = '50ppm.hdf', analysisHDF_ver0804

% clear all
%[traj,attr] = readTrajAccHDF(hfile,'x',[],'y',[],'z',[],'u',[],'v',[],'w',[],...
%    'dudx',[],'dudy',[],'dudz',[],...
%    'dvdx',[],'dvdy',[],'dvdz',[],...
%    'dwdx',[],'dwdy',[],'dwdz',[],...
%    'minlength',15,'frames',[1000 4000]); %,...
    %'reldiv',[0 0.1]);
    
    
    
% [traj,attr] = readTrajHDF_v9(hfile,...
%     'x',[],'y',[],'z',[],...
%     'u',[],'v',[],'w',[],...
%     'dudx',[],'dudy',[],'dudz',[],...
%     'dvdx',[],'dvdy',[],'dvdz',[],...
%     'dwdx',[],'dwdy',[],'dwdz',[],...
%     'minlength',20,'frames',[1 3000],'trajnum',[1 Inf]);

[traj,attr] = readTrajHDF_v9(hfile,...
    'x',[],'y',[],'z',[],...
    'u',[],'v',[],'w',[],...
    's11',[],'s12',[],'s13',[],...
    's22',[],'s23',[],'s33',[],...
    'w1',[],'w2',[],'w3',[],...
    'minlength',20,'frames',[1 3000]);

numTraj = length(traj)

numPoints = length(cat(1,traj.x))
trajLen = cat(1,traj.trajLen);

% 
% for one type of HDFs


for i = 1:numTraj
    traj(i).dudx = traj(i).s11;
    traj(i).dudy = traj(i).s12 - .5*traj(i).w3;
    traj(i).dudz = traj(i).s13 + .5*traj(i).w2;
    traj(i).dvdx = traj(i).s12 + .5*traj(i).w3;
    traj(i).dvdy = traj(i).s22;
    traj(i).dvdz = traj(i).s23 - .5*traj(i).w1;
    traj(i).dwdx = traj(i).s13 - .5*traj(i).w2;
    traj(i).dwdy = traj(i).s23 + .5*traj(i).w1;
    traj(i).dwdz = traj(i).s33;
end
% rmfield(traj,{'s11','s12','s13','s22','s23','s33','w1','w2','w3'})

% % 
% dudx = cat(1,traj.dudx);
% dudy = cat(1,traj.dudy);
% dudz = cat(1,traj.dudz);
% 
% dvdx = cat(1,traj.dvdx);
% dvdy = cat(1,traj.dvdy);
% dvdz = cat(1,traj.dvdz);
% 
% dwdx = cat(1,traj.dwdx);
% dwdy = cat(1,traj.dwdy);
% dwdz = cat(1,traj.dwdz);


% reldiv = abs(dudx + dvdy + dwdz)./(abs(dudx)+abs(dvdy)+abs(dwdz));
% absdiv = abs(dudx + dvdy + dwdz);
% inddiv = absdiv < 0.1 & reldiv < 0.1;

k = 0;
good = zeros(numTraj,1);
fn = fieldnames(traj); fn1 =   find(strcmp(fn,'trajLen'));
for i=1:numTraj
    absdiv = abs(traj(i).dudx +  traj(i).dvdy + traj(i).dwdz);
    reldiv = absdiv./(abs(traj(i).dudx) +  abs(traj(i).dvdy) + abs(traj(i).dwdz));
%    ind = find(reldiv < .1); %  & absdiv < .1);
    ind = find(reldiv < .1 & absdiv < .1);
    ind1 = diff([ind;Inf]);
    ind2 = [0; find(ind1 > 1)];
    ind3 = diff(ind2);
    [r,c] = max(ind3);
    ind = ind(ind2(c)+1:ind2(c+1));
    if length(ind) > 15
        k = k + 1;
        good(k) = i;
        for j = [1:fn1-1,fn1+1:length(fn)]
%             if length(traj(i).(fn{j}) > length(ind))
                traj(i).(fn{j}) = traj(i).(fn{j})(ind);
%             end
        end
        traj(i).trajLen = length(ind);
    end
%     
%     
%     if all(reldiv < .1 & absdiv < .1)
%         k = k + 1;
%         good(k) = i;
%     end
end

good = good(1:k);
traj = traj(good);

numTraj = length(traj)
numPoints = length(cat(1,traj.dudx))
trajLen = cat(1,traj.trajLen);
trajnum = cat(1,traj.trajnum);


% return



dudx = cat(1,traj.dudx);
dudy = cat(1,traj.dudy);
dudz = cat(1,traj.dudz);

dvdx = cat(1,traj.dvdx);
dvdy = cat(1,traj.dvdy);
dvdz = cat(1,traj.dvdz);

dwdx = cat(1,traj.dwdx);
dwdy = cat(1,traj.dwdy);
dwdz = cat(1,traj.dwdz);



maxLen    = max(trajLen)
minLen    = min(trajLen)
meanLen    = mean(trajLen)

save([namestr,'_good'],'good','numTraj','numPoints','trajnum','trajLen','minLen','runtime');
diary off % diary(['analysisHDF_ver0804_',strrep(strrep(datestr(now),' ','_'),':','_'),'.txt'])



start_points = [1; cumsum(trajLen(1:end-1))+1];     % index of the first trajectory points



% reldiv = abs(dudx + dvdy + dwdz)./(abs(dudx)+abs(dvdy)+abs(dwdz));
% absdiv = abs(dudx + dvdy + dwdz);
% inddiv = absdiv < 0.1 & reldiv < 0.1;




sSq = dudx.^2 + dvdy.^2 + dwdz.^2 + 2*((0.5*(dudy + dvdx)).^2 + (0.5*(dudz + dwdx)).^2 + (0.5*(dvdz + dwdy)).^2);
diss = mean(sSq)*2*1e-6
% diss = mean(sSq(inddiv))*2*1e-6

tau_eta = (1.0e-6/diss)^.5
eta = ((1.0e-6)^3/diss)^0.25

u = cat(1,traj.u);
v = cat(1,traj.v);
w = cat(1,traj.w);

% x = cat(1,traj.x);
% y = cat(1,traj.y);
% z = cat(1,traj.z);


uCube = mean(1/3*((u - mean(u)).^3 + (v-mean(v)).^3 + (w-mean(w)).^3))


% Taylor micro-scale, \lambda
Taylor_lambda = 1/3 * ( sqrt( mean(u.^2)./mean(dudx.^2)) + ...
    sqrt(mean(v.^2)./mean(dvdy.^2)) + ...
    sqrt(mean(w.^2)./mean(dwdz.^2)) )




uSq = 1/3*((u - mean(u)).^2 + (v-mean(v)).^2 + (w-mean(w)).^2);

Re_lambda =   sqrt(mean(uSq)) * Taylor_lambda / 1.e-6

% return

% -----------------------------------------------------------------------
 
 clear traj

 
 N = 5;
 calculate_lls_v10_hdf;
% 
% return


% data = hdfread([namestr,'_lpqr.hdf'],namestr,'Fields','d,v','firstrecord',1);
data = hdfread([namestr,'_lpqr.hdf'],namestr,'Fields','d','firstrecord',1);
% D = zeros(1,3,numPoints);
D = shiftdim([data{:}],-1);

mean(squeeze(D(:,1,:)))
mean(squeeze(D(:,2,:)))
mean(squeeze(D(:,3,:)))

sSq = squeeze(sum(D.^2,2));
hf = figure,
nhist(sSq(sSq < 50),200)
xlabel('\Lambda_1^2+\Lambda_2^2+\Lambda_3^2')
ylabel('PDF')
set(gca,'yscale','log')
saveas(hf,['pdf_s2_',namestr],'fig')% Plot histogram of the eigenvalues

hf = figure,
nhist(sSq(sSq < 50)./mean(sSq(sSq < 50)),200)
xlabel('\Lambda_1^2+\Lambda_2^2+\Lambda_3^2/\langle \Lambda_1^2+\Lambda_2^2+\Lambda_3^2 \rangle')
ylabel('PDF')
set(gca,'yscale','log')
saveas(hf,['pdf_s2_normalized_',namestr],'fig')% Plot histogram of the eigenvalues

% feval('save',[namestr,'_eigens.mat'],'D','V' ,'trajLen','numPoints','numTraj','inddiv');



% histogram of eigenvalues
% [n1,x1] = nhist(D(:,1,:),200);
% [n2,x2] = nhist(D(:,2,:),200);
% [n3,x3] = nhist(D(:,3,:),200);
hf = figure
nhist(D(:,1,:),201);
hold on
nhist(D(:,2,:),201);
nhist(D(:,3,:),201);
% semilogy(x1,n1,'b-',x2,n2,'b-s',x3,n3,'b-x','linewidth',1,'markersize',4);
set(gca,'yscale','log')
saveas(hf,['pdf_lambda_',namestr],'fig')% Plot histogram of the eigenvalues
% saveas(hf,['pdf_lambda_',namestr],'psc2')


% % histogram of eigenvalues
% [n1,x1] = nhist(D(:,1,inddiv)/sqrt(diss/2e-6),200);
% [n2,x2] = nhist(D(:,2,inddiv)/sqrt(diss/2e-6),200);
% [n3,x3] = nhist(D(:,3,inddiv)/sqrt(diss/2e-6),200);
% hf = figure
% semilogy(x1,n1,'r-',x2,n2,'r-s',x3,n3,'r-x','linewidth',1,'markersize',4);
% hold on



hf = figure
hold on
% histogram of eigenvalues
[n1,x1] = nhist(D(:,1,:)/sqrt(diss/2e-6),200);
[n2,x2] = nhist(D(:,2,:)/sqrt(diss/2e-6),200);
[n3,x3] = nhist(D(:,3,:)/sqrt(diss/2e-6),200);
% hf = figure
semilogy(x1,n1,'b-',x2,n2,'b-s',x3,n3,'b-x','linewidth',1,'markersize',4);



% % histogram of eigenvalues
% [n1,x1] = nhist(D(:,1,reldiv < .1)/sqrt(diss/2e-6),200);
% [n2,x2] = nhist(D(:,2,reldiv < .1)/sqrt(diss/2e-6),200);
% [n3,x3] = nhist(D(:,3,reldiv < .1)/sqrt(diss/2e-6),200);
% % hf = figure
% semilogy(x1,n1,'g-',x2,n2,'g-s',x3,n3,'g-x','linewidth',1,'markersize',4);
saveas(hf,['pdf_lambda_sumlambda2_',namestr],'fig')% Plot histogram of the eigenvalues
% saveas(hf,['pdf_lambda_sumlambda2_',namestr],'fig')% Plot histogram of the eigenvalues


% --------------------------------------------------------------------
inds = sSq < prctile(sSq,99);
mean_ssq = mean(sSq(inds))
strong_s = (sSq > 2*mean_ssq) & inds;
weak_s = (sSq < 1/3*mean_ssq) & inds;

% -------------------------------------------------------------------------
% omega = cat(2,cat(1,traj.w1), cat(1,traj.w2),cat(1,traj.w3));
% omega = [w1(ind) w2(ind) w3(ind)];
omega = [dwdy - dvdz, dudz - dwdx, dvdx - dudy];
%w1 = dwdy - dvdz;
%w2 = dudz - dwdx;
%w3 = dvdx - dudy;

enstrophy = sum(omega.^2,2);
indw = enstrophy < prctile(enstrophy,99);



hf = figure;
nhist(enstrophy(indw),100);
hold on
nhist(sSq(inds),100);
saveas(hf,['pdf_s2_enstrophy_',namestr],'fig')


mean_omegasq = mean(enstrophy(indw))
strong_omega = (enstrophy > 3*mean_omegasq) & indw;
weak_omega = (enstrophy < 1/3*mean_omegasq) & indw;




cos_sS = zeros(numPoints,1);
% Big S is the s^2 (not s.^2)

for j = 1:numPoints
    % tmp = reshape(s(i,1:9),[3 3])';
    s = [0.5*(dudx(j) + dudx(j)), 0.5*(dudy(j) + dvdx(j)), 0.5*(dudz(j) + dwdx(j));
    0.5*(dudy(j) + dvdx(j)),0.5*(dvdy(j) + dvdy(j)),0.5*(dvdz(j) + dwdy(j));
    0.5*(dudz(j) + dwdx(j)),0.5*(dvdz(j) + dwdy(j)),0.5*(dwdz(j) + dwdz(j))];
    S = -reshape((s^2)',[1 9]);
%     cos_sS(j) = cosine(s(:)',S);
    % cos_sS(j) = sum(s(:)'.*S,2)./sqrt(sum(s(:)'.^2,2))./sqrt(sum(S.^2,2));
    cos_sS(j) = sum(s(:)'.*S,2)./norm(s(:))./norm(S); 
end

% cos_sS = cosine(s,S);
% PDF of the cos(s,S)
hf = figure
nhist(cos_sS,50,'b:o');
xlabel('cos(s,S), unconditioned')
saveas(hf,['pdf_cos_sS_uncond',namestr],'fig')


hf = figure
nhist(cos_sS(strong_s),25,'b:o');
xlabel('cos(s,S), > 3 \langles^2\rangle')
saveas(hf,['pdf_cos_sS_strongs_',namestr],'fig')

hf = figure
nhist(cos_sS(weak_s),25,'g-.x');
xlabel('cos(s,S), < 1/3 \langles^2\rangle')
saveas(hf,['pdf_cos_sS_',namestr],'fig')

% -------------------------------------------------------------------------
% Joint PDF  of cos_sS versus s^2
hf = figure
[n,x,bins] = histmulti5([sSq(inds),cos_sS(inds)]);
contour(x(:,1),x(:,2),filter2(fspecial('gaussian',[10 10]),n'))
axis tight
xlabel('s^2')
ylabel('cos(s,S)')
grid on
axis square
saveas(hf,['jointpdf_s2_cos_sS_',namestr],'fig')


% Joint PDF  of cos_sS versus enstrophy
hf = figure
[n,x,bins] = histmulti5([enstrophy(inds),cos_sS(inds)]);
contour(x(:,1),x(:,2),n')
axis tight
xlabel('\omega^2')
ylabel('cos(s,S)')
grid on
axis square
saveas(hf,['jointpdf_enstrophy_cos_sS_',namestr],'fig')




% New, combined figure of all options
hf = figure
[n1,x1] = nhist(cos_sS,15);
[n2,x2] = nhist(cos_sS(strong_s),15);
[n3,x3] = nhist(cos_sS(weak_s),15);
[n4,x4] = nhist(cos_sS(strong_omega),15);
[n5,x5] = nhist(cos_sS(weak_omega),15);
[n6,x6] = nhist(cos_sS(weak_omega & strong_s),15);
[n7,x7] = nhist(cos_sS(weak_s & strong_omega),15);
[n8,x8] = nhist(cos_sS(weak_s & weak_omega),15);
[n9,x9] = nhist(cos_sS(strong_s & strong_omega),15);

plot(x1,n1,'ro',x2,n2,'bs',x3,n3,'k>',x4,n4,'md',x5,n5,'gx',x6,n6,'cv',x7,n7,'r+',x8,n8,'bp',x9,n9,'k^');
legend('no cond.','strong s','weak s','strong \omega','weak \omega','weak \omega & strong s','strong \omega & weak s','weak both','strong both')
hold on
for i=1:9
	eval(sprintf('fnplt(csapi(x%d,n%d),%s,[],0.5)',i,i,''':'''))
end
hold off
xlabel('cos(s,S)')
saveas(hf,['pdf_cos_sS_conditions_',namestr],'fig')


% -------------------------------------------------------------------------
% w_iS_{ij}
W = repmat(0,[numPoints 3]);
for j = 1:numPoints
    % tmp = reshape(s(i,1:9),[3 3]);
    W(j,1:3) = omega(j,1:3)*[0.5*(dudx(j) + dudx(j)), 0.5*(dudy(j) + dvdx(j)), 0.5*(dudz(j) + dwdx(j));
    0.5*(dudy(j) + dvdx(j)),0.5*(dvdy(j) + dvdy(j)),0.5*(dvdz(j) + dwdy(j));
    0.5*(dudz(j) + dwdx(j)),0.5*(dvdz(j) + dwdy(j)),0.5*(dwdz(j) + dwdz(j))]';
end

% feval('save',[namestr,'_wws_sS.mat'],'W','cos_sS','trajLen','numPoints','numTraj');

% coswW = dot(omega,W,2)./sqrt(sum(W.^2,2))./sqrt(sum(omega.^2,2));
coswW = cosine(omega,W);

% PDF of the cos(w,W)
hf = figure
nhist(coswW);
xlabel('cos(\omega,W), unconditioned')
saveas(hf,['pdf_cos_wW_',namestr],'fig')


hf = figure;
nhist(coswW,50,'r-');
hold on
nhist(coswW(strong_s),50,'b--');
nhist(coswW(weak_s),50,'k-.');
hold off
xlabel('cos(\omega,W)');
ylabel('PDF')
legend('uncond','strong','weak',2)
saveas(hf,['pdf_coswW_',namestr],'fig')

% New, combined figure of all options
hf = figure
[n1,x1] = nhist(coswW,15);
[n2,x2] = nhist(coswW(strong_s),15);
[n3,x3] = nhist(coswW(weak_s),15);
[n4,x4] = nhist(coswW(strong_omega),15);
[n5,x5] = nhist(coswW(weak_omega),15);
[n6,x6] = nhist(coswW(weak_omega & strong_s),15);
[n7,x7] = nhist(coswW(weak_s & strong_omega),15);
[n8,x8] = nhist(coswW(weak_s & weak_omega),15);
[n9,x9] = nhist(coswW(strong_s & strong_omega),15);
plot(x1,n1,'ro',x2,n2,'bs',x3,n3,'k>',x4,n4,'md',x5,n5,'gx',x6,n6,'cv',x7,n7,'r+',x8,n8,'bp',x9,n9,'k^');
legend('no cond.','strong s','weak s','strong \omega','weak \omega','weak \omega & strong s','strong \omega & weak s','weak both','strong both')
xlabel('cos(\omega,W)');
colors={'r','b','k','m','g','c','r','b','k'};
hold on
for i=1:9
	eval(sprintf('fnplt(csapi(x%d,n%d),''%s:'',[],0.5)',i,i,colors{i}))
end
hold off
saveas(hf,['pdf_coswW_conditions_',namestr],'fig')




% -------------------------------------------------------------------------


data = hdfread([namestr,'_lpqr.hdf'],namestr,'Fields','v','firstrecord',1);
V = reshape(shiftdim([data{:}],-1),3,3,[]);

% cos(w,\lambda_i)
cos_wlambda1 = cosine(omega(:,:),squeeze(V(1:3,1,:))');
cos_wlambda2 = cosine(omega(:,:),squeeze(V(1:3,2,:))');
cos_wlambda3 = cosine(omega(:,:),squeeze(V(1:3,3,:))');


[n1,x1] = nhist(cos_wlambda1,21);
[n2,x2] = nhist(cos_wlambda2,21);
[n3,x3] = nhist(cos_wlambda3,21);
hf = figure, hold on
plot(x1,n1,'ro',x2,n2,'rs',x3,n3,'r<','linewidth',1);
legend('cos(\omega,\lambda_1)','cos(\omega,\lambda_2)','cos(\omega,\lambda_3)')
for i=1:3
	eval(sprintf('fnplt(csapi(x%d,n%d),''%s:'',[],0.5)',i,i,'r'))
end

xlabel('cos(\omega,\lambda_i)')
saveas(hf,['pdf_cos_w_lambda_',namestr],'fig')

[n1,x1] = nhist(cos_wlambda1(strong_s));
[n2,x2] = nhist(cos_wlambda2(strong_s));
[n3,x3] = nhist(cos_wlambda3(strong_s));
hf = figure
plot(x1,n1,'r-',x2,n2,'r-s',x3,n3,'r-x','linewidth',1);
xlabel('cos(\omega,\lambda_i) > 3 \langles^2\rangle')
saveas(hf,['pdf_cos_w_lambda_strongs_',namestr],'fig')

[n1,x1] = nhist(cos_wlambda1(weak_s));
[n2,x2] = nhist(cos_wlambda2(weak_s));
[n3,x3] = nhist(cos_wlambda3(weak_s));
hf = figure
plot(x1,n1,'r-',x2,n2,'r-s',x3,n3,'r-x','linewidth',1);
xlabel('cos(\omega,\lambda_i) < 1/3 \langles^2\rangle')
saveas(hf,['pdf_cos_w_lambda_weaks_',namestr],'fig')

% -------------------------------------------------------------
% New, combined figure of all options
hf = figure
[n1,x1] = nhist(cos_wlambda1,15);
[n2,x2] = nhist(cos_wlambda1(strong_s),15);
[n3,x3] = nhist(cos_wlambda1(weak_s),15);
[n4,x4] = nhist(cos_wlambda1(strong_omega),15);
[n5,x5] = nhist(cos_wlambda1(weak_omega),15);
[n6,x6] = nhist(cos_wlambda1(weak_omega & strong_s),15);
[n7,x7] = nhist(cos_wlambda1(weak_s & strong_omega),15);
[n8,x8] = nhist(cos_wlambda1(weak_s & weak_omega),15);
[n9,x9] = nhist(cos_wlambda1(strong_s & strong_omega),15);
plot(x1,n1,'ro',x2,n2,'bs',x3,n3,'k>',x4,n4,'md',x5,n5,'gx',x6,n6,'cv',x7,n7,'r+',x8,n8,'bp',x9,n9,'k^');
legend('no cond.','strong s','weak s','strong \omega','weak \omega','weak \omega & strong s','strong \omega & weak s','weak both','strong both')
xlabel('cos(\omega,\lambda_1)');
colors={'r','b','k','m','g','c','r','b','k'};
hold on
for i=1:9
	eval(sprintf('fnplt(csapi(x%d,n%d),''%s:'',[],0.5)',i,i,colors{i}))
end
hold off
saveas(hf,['pdf_cosw_lambda1_conditions_',namestr],'fig')

% -------------------------------------------------------------
% New, combined figure of all options
hf = figure
[n1,x1] = nhist(cos_wlambda2,15);
[n2,x2] = nhist(cos_wlambda2(strong_s),15);
[n3,x3] = nhist(cos_wlambda2(weak_s),15);
[n4,x4] = nhist(cos_wlambda2(strong_omega),15);
[n5,x5] = nhist(cos_wlambda2(weak_omega),15);
[n6,x6] = nhist(cos_wlambda2(weak_omega & strong_s),15);
[n7,x7] = nhist(cos_wlambda2(weak_s & strong_omega),15);
[n8,x8] = nhist(cos_wlambda2(weak_s & weak_omega),15);
[n9,x9] = nhist(cos_wlambda2(strong_s & strong_omega),15);
plot(x1,n1,'ro',x2,n2,'bs',x3,n3,'k>',x4,n4,'md',x5,n5,'gx',x6,n6,'cv',x7,n7,'r+',x8,n8,'bp',x9,n9,'k^');
legend('no cond.','strong s','weak s','strong \omega','weak \omega','weak \omega & strong s','strong \omega & weak s','weak both','strong both')
xlabel('cos(\omega,\lambda_2)');
colors={'r','b','k','m','g','c','r','b','k'};
hold on
for i=1:9
	eval(sprintf('fnplt(csapi(x%d,n%d),''%s:'',[],0.5)',i,i,colors{i}))
end
hold off
saveas(hf,['pdf_cosw_lambda2_conditions_',namestr],'fig')

% -------------------------------------------------------------
% New, combined figure of all options
hf = figure
[n1,x1] = nhist(cos_wlambda3,15);
[n2,x2] = nhist(cos_wlambda3(strong_s),15);
[n3,x3] = nhist(cos_wlambda3(weak_s),15);
[n4,x4] = nhist(cos_wlambda3(strong_omega),15);
[n5,x5] = nhist(cos_wlambda3(weak_omega),15);
[n6,x6] = nhist(cos_wlambda3(weak_omega & strong_s),15);
[n7,x7] = nhist(cos_wlambda3(weak_s & strong_omega),15);
[n8,x8] = nhist(cos_wlambda3(weak_s & weak_omega),15);
[n9,x9] = nhist(cos_wlambda3(strong_s & strong_omega),15);
plot(x1,n1,'ro',x2,n2,'bs',x3,n3,'k>',x4,n4,'md',x5,n5,'gx',x6,n6,'cv',x7,n7,'r+',x8,n8,'bp',x9,n9,'k^');
legend('no cond.','strong s','weak s','strong \omega','weak \omega','weak \omega & strong s','strong \omega & weak s','weak both','strong both')
xlabel('cos(\omega,\lambda_3)');
colors={'r','b','k','m','g','c','r','b','k'};
hold on
for i=1:9
	eval(sprintf('fnplt(csapi(x%d,n%d),''%s:'',[],0.5)',i,i,colors{i}))
end
hold off
saveas(hf,['pdf_cosw_lambda3_conditions_',namestr],'fig')
% -------------------------------------------------------------------------
% New figures, 23.11.04, for Michele
% different constraints
%
% No. 1-3: strong strain, weak enstrophy, 2 times

% Make it all in once, by using COMBPLOT, thanks to Toby Driscoll
close all

ind1 = (sSq > 1*mean_ssq) & inds & (enstrophy < 1*mean_omegasq) & indw;
[n1,x1] = nhist(cos_wlambda1(ind1));
[n2,x2] = nhist(cos_wlambda2(ind1));
[n3,x3] = nhist(cos_wlambda3(ind1));
hf = figure
plot(x1,n1,'r--o',x2,n2,'g--s',x3,n3,'b--^','linewidth',1);
xlabel(['$cos(\omega,\lambda_i), s^2 > \langle s^2 \rangle \&\ \omega^2 < \langle \omega^2 \rangle$, No.= ',int2str(sum(ind1))],'Interpreter','Latex')

% saveas(hf,['pdf_cos_w_lambda_Michele1_',namestr],'fig')

ind1 = (sSq > 2*mean_ssq) & inds & (enstrophy < 0.5*mean_omegasq) & indw;
[n1,x1] = nhist(cos_wlambda1(ind1));
[n2,x2] = nhist(cos_wlambda2(ind1));
[n3,x3] = nhist(cos_wlambda3(ind1));
hf = figure
plot(x1,n1,'r--o',x2,n2,'g--s',x3,n3,'b--^','linewidth',1);
xlabel(['$cos(\omega,\lambda_i), s^2 > 2\langle s^2 \rangle \&\ \omega^2 < 1/2\langle \omega^2 \rangle$, No.= ',int2str(sum(ind1))],'Interpreter','Latex')

% saveas(hf,['pdf_cos_w_lambda_Michele2_',namestr],'fig')

ind1 = (sSq > 3*mean_ssq) & inds & (enstrophy < 0.33*mean_omegasq) & indw;
[n1,x1] = nhist(cos_wlambda1(ind1));
[n2,x2] = nhist(cos_wlambda2(ind1));
[n3,x3] = nhist(cos_wlambda3(ind1));
hf = figure
plot(x1,n1,'r--o',x2,n2,'g--s',x3,n3,'b--^','linewidth',1);
xlabel(['$cos(\omega,\lambda_i), s^2 > 3\langle s^2 \rangle \&\ \omega^2 < 1/3\langle \omega^2 \rangle$, No.= ',int2str(sum(ind1))],'Interpreter','Latex')

% saveas(hf,['pdf_cos_w_lambda_Michele3_',namestr],'fig')

% No. 4-6: strong enstrophy, weak strain

ind1 = (sSq < 1*mean_ssq) & inds & (enstrophy > 1*mean_omegasq) & indw;
[n1,x1] = nhist(cos_wlambda1(ind1));
[n2,x2] = nhist(cos_wlambda2(ind1));
[n3,x3] = nhist(cos_wlambda3(ind1));
hf = figure
plot(x1,n1,'r--o',x2,n2,'g--s',x3,n3,'b--^','linewidth',1);
xlabel(['$cos(\omega,\lambda_i), s^2 < \langle s^2 \rangle \&\ \omega^2 > \langle \omega^2 \rangle$, No.= ',int2str(sum(ind1))],'Interpreter','Latex')

% saveas(hf,['pdf_cos_w_lambda_Michele4_',namestr],'fig')

ind1 = (sSq < 0.5*mean_ssq) & inds & (enstrophy > 2*mean_omegasq) & indw;
[n1,x1] = nhist(cos_wlambda1(ind1));
[n2,x2] = nhist(cos_wlambda2(ind1));
[n3,x3] = nhist(cos_wlambda3(ind1));
hf = figure
plot(x1,n1,'r--o',x2,n2,'g--s',x3,n3,'b--^','linewidth',1);
xlabel(['$cos(\omega,\lambda_i), s^2 < 1/2\langle s^2 \rangle \&\ \omega^2 > 2\langle \omega^2 \rangle$, No.= ',int2str(sum(ind1))],'Interpreter','Latex')

% saveas(hf,['pdf_cos_w_lambda_Michele5_',namestr],'fig')

ind1 = (sSq < 0.33*mean_ssq) & inds & (enstrophy > 3*mean_omegasq) & indw;
[n1,x1] = nhist(cos_wlambda1(ind1));
[n2,x2] = nhist(cos_wlambda2(ind1));
[n3,x3] = nhist(cos_wlambda3(ind1));
hf = figure
plot(x1,n1,'r--o',x2,n2,'g--s',x3,n3,'b--^','linewidth',1);
xlabel(['$cos(\omega,\lambda_i), s^2 < 1/3\langle s^2 \rangle \&\ \omega^2 > 3\langle \omega^2 \rangle$, No.= ',int2str(sum(ind1))],'Interpreter','Latex')

combplot(hf-5:hf,[2 3])
saveas(gcf,['pdf_cos_w_lambda_Michele_',namestr],'fig')




% -------------------------------------------------------------------------
% wws
wws = repmat(0,[numPoints 1]);
for i = 1:numPoints
    wws(i) = omega(i,1:3)*W(i,1:3)';
end
wws_w2 = wws./sum(omega.^2,2);


% -------------------------------------------------------------
% New, combined figure of all options
hf = figure
[n1,x1] = nhist(wws,15);
[n2,x2] = nhist(wws(strong_s),15);
[n3,x3] = nhist(wws(weak_s),15);
[n4,x4] = nhist(wws(strong_omega),15);
[n5,x5] = nhist(wws(weak_omega),15);
[n6,x6] = nhist(wws(weak_omega & strong_s),15);
[n7,x7] = nhist(wws(weak_s & strong_omega),15);
[n8,x8] = nhist(wws(weak_s & weak_omega),15);
[n9,x9] = nhist(wws(strong_s & strong_omega),15);
plot(x1,n1,'ro',x2,n2,'bs',x3,n3,'k>',x4,n4,'md',x5,n5,'gx',x6,n6,'cv',x7,n7,'r+',x8,n8,'bp',x9,n9,'k^');
set(gca,'yscale','log')
legend('no cond.','strong s','weak s','strong \omega','weak \omega','weak \omega & strong s','strong \omega & weak s','weak both','strong both')
xlabel('\omega_{i}\omega_{j}s_{ij}');
colors={'r','b','k','m','g','c','r','b','k'};
hold on
for i=1:9
	eval(sprintf('fnplt(csapi(x%d,n%d),''%s:'',[],0.5)',i,i,colors{i}))
end
hold off
saveas(hf,['pdf_wws_conditions_',namestr],'fig')

% -------------------------------------------------------------
% New, combined figure of all options
hf = figure
[n1,x1] = nhist(wws_w2,15);
[n2,x2] = nhist(wws_w2(strong_s),15);
[n3,x3] = nhist(wws_w2(weak_s),15);
[n4,x4] = nhist(wws_w2(strong_omega),15);
[n5,x5] = nhist(wws_w2(weak_omega),15);
[n6,x6] = nhist(wws_w2(weak_omega & strong_s),15);
[n7,x7] = nhist(wws_w2(weak_s & strong_omega),15);
[n8,x8] = nhist(wws_w2(weak_s & weak_omega),15);
[n9,x9] = nhist(wws_w2(strong_s & strong_omega),15);
plot(x1,n1,'ro',x2,n2,'bs',x3,n3,'k>',x4,n4,'md',x5,n5,'gx',x6,n6,'cv',x7,n7,'r+',x8,n8,'bp',x9,n9,'k^');
legend('no cond.','strong s','weak s','strong \omega','weak \omega','weak \omega & strong s','strong \omega & weak s','weak both','strong both')
xlabel('\omega_{i}\omega_{j}s_{ij}/\omega^2');
colors={'r','b','k','m','g','c','r','b','k'};
hold on
for i=1:9
	eval(sprintf('fnplt(csapi(x%d,n%d),''%s:'',[],0.5)',i,i,colors{i}))
end
hold off
saveas(hf,['pdf_wws_w2_conditions_',namestr],'fig')


%hf = figure,
%nhist(wws_w2);
%set(gca,'yscale','log')
%xlabel('\omega_{i}\omega_{j}s_{ij}/\omega^2','fontname','euclid','fontangle','italic')
%saveas(hf,['pdf_wws_w2_',namestr],'fig')
%
%% conditioned wws_w2, strong omega
%hf = figure,
%nhist(wws_w2(strong_omega));
%set(gca,'yscale','log')
%xlabel('\omega_{i}\omega_{j}s_{ij}/\omega^2, > 3 \langle\omega^2\rangle','fontname','euclid','fontangle','italic')
%saveas(hf,['pdf_wws_w2_strong_w',namestr],'fig')
%
%% conditioned wws_w2, weak omega
%hf = figure,
%nhist(wws_w2(weak_omega));
%set(gca,'yscale','log')
%xlabel('\omega_{i}\omega_{j}s_{ij}/\omega^2, < 1/3 \langle\omega^2\rangle','fontname','euclid','fontangle','italic')
%saveas(hf,['pdf_wws_w2_weak_w',namestr],'fig')
%% -------------------------------------------------------------------------
%% One special figure is low-pressure  = strong enstrophy and weak
%% dissipation
%hf = figure,
%nhist(wws_w2(strong_omega & weak_s));
%set(gca,'yscale','log')
%xlabel('\omega_{i}\omega_{j}s_{ij}/\omega^2, > 3\langle\omega^2\rangle & < 1/3 \langles^2\rangle','fontname','euclid','fontangle','italic')
%saveas(hf,['pdf_wws_w2_strong_w_weak_s',namestr],'fig')
%
%% -------------------------------------------------------------------------
%% Second special figure is low-pressure  = strong dissipation and weak
%% enstrophy
%hf = figure,
%nhist(wws_w2(strong_s & weak_omega));
%set(gca,'yscale','log')
%xlabel('\omega_{i}\omega_{j}s_{ij}/\omega^2, < 1/3 \langle\omega^2\rangle & > 3 \langles^2\rangle','fontname','euclid','fontangle','italic')
%saveas(hf,['pdf_wws_w2_strong_s_weak_w_',namestr],'fig')


% -------------------------------------------------------------------------
[sss,sss_s2] = deal(repmat(0,[numPoints, 1]));
i = 1:numPoints;
%     s11 = 0.5*(dudx + dudx);
%     s12 = 0.5*(dudy + dvdx);
%     s13 = 0.5*(dudz(i) + dwdx(i));
%     s22 = 0.5*(dvdy + dvdy);
%     s23 = 0.5*(dvdz(i) + dwdy(i));
%     s33 = 0.5*(dwdz + dwdz);

s111 = dudx(i).^3;
s222 = dvdy(i).^3;
s333 = dwdz(i).^3;
s112 = dudx(i).*((0.5*(dudy(i) + dvdx(i))).^2);
s113 = dudx(i).*((0.5*(dudz(i) + dwdx(i))).^2);
s221 = dvdy(i).*((0.5*(dudy(i) + dvdx(i))).^2);
s223 = dvdy(i).*((0.5*(dvdz(i) + dwdy(i))).^2);
s331 = dwdz(i).*((0.5*(dudz(i) + dwdx(i))).^2);
s332 = dwdz(i).*((0.5*(dvdz(i) + dwdy(i))).^2);
s123 = (0.5*(dudy(i) + dvdx(i))).*(0.5*(dvdz(i) + dwdy(i))).*(0.5*(dudz(i) + dwdx(i)));
sss(i) = -1 * (s111 + s222 + s333 +...
    3.*(s112 + s113 + s221 + s223 + s331 + s332) +...
    6.*s123);

clear s111 s222 s333 s112 s113 s221 s223 s331 s332 s123    
sss_s2(i) = sss(i)./sSq(i);

feval('save',[namestr,'_wws_sss.mat'],'wws','wws_w2','sss','sss_s2');


% -------------------------------------------------------------
% New, combined figure of all options
hf = figure
[n1,x1] = nhist(sss,15);
[n2,x2] = nhist(sss(strong_s),15);
[n3,x3] = nhist(sss(weak_s),15);
[n4,x4] = nhist(sss(strong_omega),15);
[n5,x5] = nhist(sss(weak_omega),15);
[n6,x6] = nhist(sss(weak_omega & strong_s),15);
[n7,x7] = nhist(sss(weak_s & strong_omega),15);
[n8,x8] = nhist(sss(weak_s & weak_omega),15);
[n9,x9] = nhist(sss(strong_s & strong_omega),15);
plot(x1,n1,'ro',x2,n2,'bs',x3,n3,'k>',x4,n4,'md',x5,n5,'gx',x6,n6,'cv',x7,n7,'r+',x8,n8,'bp',x9,n9,'k^');
legend('no cond.','strong s','weak s','strong \omega','weak \omega','weak \omega & strong s','strong \omega & weak s','weak both','strong both')
xlabel('s_{ij}s_{jk}s_{ki}');
colors={'r','b','k','m','g','c','r','b','k'};
hold on
for i=1:9
	eval(sprintf('fnplt(csapi(x%d,n%d),''%s:'',[],0.5)',i,i,colors{i}))
end
hold off
saveas(hf,['pdf_sss_conditions_',namestr],'fig')

% -------------------------------------------------------------
% New, combined figure of all options
hf = figure
[n1,x1] = nhist(sss_s2,15);
[n2,x2] = nhist(sss_s2(strong_s),15);
[n3,x3] = nhist(sss_s2(weak_s),15);
[n4,x4] = nhist(sss_s2(strong_omega),15);
[n5,x5] = nhist(sss_s2(weak_omega),15);
[n6,x6] = nhist(sss_s2(weak_omega & strong_s),15);
[n7,x7] = nhist(sss_s2(weak_s & strong_omega),15);
[n8,x8] = nhist(sss_s2(weak_s & weak_omega),15);
[n9,x9] = nhist(sss_s2(strong_s & strong_omega),15);
semilogy(x1,n1,'ro',x2,n2,'bs',x3,n3,'k>',x4,n4,'md',x5,n5,'gx',x6,n6,'cv',x7,n7,'r+',x8,n8,'bp',x9,n9,'k^');
legend('no cond.','strong s','weak s','strong \omega','weak \omega','weak \omega & strong s','strong \omega & weak s','weak both','strong both')
xlabel('s_{ij}s_{jk}s_{ki}/s^2');
colors={'r','b','k','m','g','c','r','b','k'};
hold on
for i=1:9
	eval(sprintf('semilogy(x%d,ppval(x%d,csapi(x%d,n%d)),''%s:'')',i,i,i,i,colors{i}))
end
hold off
saveas(hf,['pdf_sss_s2_conditions_',namestr],'fig')



hf = figure
nhist(wws(inds),500,'r')
hold on
nhist(4/3*sss(inds),500,'b')
set(gca,'yscale','log')
saveas(hf,['pdf_sss_wws_',namestr],'fig')



hf = figure,
nhist(sss_s2(inds),200);
set(gca,'yscale','log')
xlabel('-s_{ij}s_{ji}s_{ki}/s^2','fontname','euclid','fontangle','italic')
saveas(hf,['pdf_sss_s2_',namestr],'fig')
%
%% conditioned sss_s2
%hf = figure,
%nhist(sss_s2(strong_s));
%set(gca,'yscale','log')
%xlabel('s_{ij}s_{ij}s_{ij}/s^2, > 3 \langles^2\rangle','fontname','euclid','fontangle','italic')
%saveas(hf,['pdf_sss_s2_strongs',namestr],'fig')
%
%% conditioned sss_s2
%hf = figure,
%nhist(sss_s2(weak_s));
%set(gca,'yscale','log')
%xlabel('s_{ij}s_{jk}s_{ki}/s^2,  < 1/3 \langles^2\rangle','fontname','euclid','fontangle','italic')
%saveas(hf,['pdf_sss_s2_weaks',namestr],'fig')


% -------------------------------------------------------------------------
%calculate_lls(traj,10);
%feval('load',[namestr,'_lines'],'lB','lls','Wls','Wlu','coslWlu','coslWls','lls_l2');


 % New piece: plot of the values from the HDF files:
 data = hdfread([namestr,'_lpqr.hdf'],namestr,'Fields','coslwls,coslwlu,age','firstrecord',1);
 coslWls = data{1};
 coslWlu = data{2};
 age = data{3};
 
dt = 1/60;
% agebins = [-0.5 0.5 1.5 2.5 3.5 4.5 5.5 6.5]*tau_eta/dt
agebins = [0:6]*tau_eta/dt

% agebins = [15 30 60 90 120];  %*tau_eta/dt;
% for i = 1:length(traj)
%     traj(i).age = 0:traj(i).trajLen-1;
% end
% age = cat(2,traj.age);
ind = bindex(age,agebins,1);   




colors = {'m','r','g','b','k'};
linestyles = {':','-','--','-.'};
markers = {'v','','o','s','x','+','<','^','>','d','p','*'};

for i = 1:length(agebins)
    style{i} = [colors{mod(i,length(colors))+1},linestyles{mod(i,length(linestyles))+1},markers{mod(i,length(markers))+1}];
end

hf = figure,
hold on
[n,x] = deal(zeros(20,length(agebins)));
for i = 1:length(agebins)
    tmp = coslWls(:,(ind == i));
    agebins(i),length(tmp)
    if ~isempty(tmp)
        nhist(tmp(:),21,style{i});
    end
end
hold off
legend('\tau_{\eta} \leq 1','1 \leq \tau_{\eta} \leq 2','2 \leq \tau_{\eta} \leq 4','4 \leq \tau_{\eta} \leq 6','6 \leq \tau_{\eta}')
xlabel('cos(l,W^{l,s})')
saveas(hf,['pdf_coslWls_',namestr],'fig')


hf = figure,
hold on
[n,x] = deal(zeros(20,length(agebins)));
for i = 1:length(agebins)
    tmp = coslWlu(:,(ind == i));
    agebins(i),length(tmp)
    if ~isempty(tmp)
        nhist(tmp(:),21,style{i});
    end
end
hold off
legend('\tau_{\eta} \leq 1','1 \leq \tau_{\eta} \leq 2','2 \leq \tau_{\eta} \leq 4','4 \leq \tau_{\eta} \leq 6','6 \leq \tau_{\eta}')
xlabel('cos(l,W^{l,u})')
saveas(hf,['pdf_coslWlu_',namestr],'fig')


data = hdfread([namestr,'_lpqr.hdf'],namestr,'Fields','cosllambda1,cosllambda2,cosllambda3','firstrecord',1);
[cosllambda1,cosllambda2,cosllambda3] = deal(data{1:3});


hf = figure,
hold on
[n,x] = deal(zeros(20,length(agebins)));
for i = 1:length(agebins)
    tmp = cosllambda1(:,(ind == i));
    agebins(i),length(tmp)
    if ~isempty(tmp)
        nhist(tmp(:),21,style{i});
    end
end
hold off
legend('\tau_{\eta} \leq 1','1 \leq \tau_{\eta} \leq 2','2 \leq \tau_{\eta} \leq 4','4 \leq \tau_{\eta} \leq 6','6 \leq \tau_{\eta}')
xlabel('cos(l,\lambda_1)')
saveas(hf,['pdf_cosllambda1_',namestr],'fig')

hf = figure,
hold on
[n,x] = deal(zeros(20,length(agebins)));
for i = 1:length(agebins)
    tmp = cosllambda2(:,(ind == i));
    agebins(i),length(tmp)
    if ~isempty(tmp)
        nhist(tmp(:),21,style{i});
    end
end
hold off
legend('\tau_{\eta} \leq 1','1 \leq \tau_{\eta} \leq 2','2 \leq \tau_{\eta} \leq 4','4 \leq \tau_{\eta} \leq 6','6 \leq \tau_{\eta}')
xlabel('cos(l,\lambda_2)')
saveas(hf,['pdf_cosllambda2_',namestr],'fig')


hf = figure,
hold on
[n,x] = deal(zeros(20,length(agebins)));
for i = 1:length(agebins)
    tmp = cosllambda3(:,(ind == i));
    agebins(i),length(tmp)
    if ~isempty(tmp)
        nhist(tmp(:),21,style{i});
    end
end
hold off
legend('\tau_{\eta} \leq 1','1 \leq \tau_{\eta} \leq 2','2 \leq \tau_{\eta} \leq 4','4 \leq \tau_{\eta} \leq 6','6 \leq \tau_{\eta}')
xlabel('cos(l,\lambda_3)')
saveas(hf,['pdf_cosllambda3_',namestr],'fig')


% Plot of material line stretching rate
data = hdfread([namestr,'_lpqr.hdf'],namestr,'Fields','zeta,age','firstrecord',1);

[zeta,age2] = deal(data{1:2});
% [zeta] = deal(data{1});
% Age is arleady loaded, if not load it again
agebins = [0:.1:6]*tau_eta/dt
% ind = bindex(age,agebins,1);   
ind = bindex(age2,agebins,1);   
hf = figure,
tmp = zeros(size(agebins));
for i = 1:length(agebins)
%    tmp(i) = mean(mean(dlnLdt(1:N,1,(ind == i)),1),3);
        tmp(i) = mean(mean(zeta(1:N,(ind == i)),1),2);
        %tmp2(i) = mean(mean(zeta(2,(ind == i)),1),2);
        %tmp3(i) = mean(mean(zeta(3,(ind == i)),1),2);
end
plot(agebins(1:end-1)/tau_eta*dt,tmp(1:end-1))
ylabel('dln(l)/dt','fontname','euclid','fontangle','italic')
xlabel('\tau_{\eta}')
% saveas(hf,['mean_zeta_tau_',namestr],'fig')

% 13-12-2004, Alex.
% Comparative figure, if you need:
% Plot of material line stretching rate
namestr2 = '0211_20';
data = hdfread([namestr2,'_lpqr.hdf'],namestr2,'Fields','zeta,age','firstrecord',1);

[zeta2,age2] = deal(data{1:2});
ind = bindex(age2,agebins,1);   
tmp = zeros(size(agebins));
for i = 1:length(agebins)
        tmp(i) = mean(mean(zeta2(1:N,(ind == i)),1),2);
end
figure(hf); hold on
plot(agebins(1:end-1)/tau_eta*dt,tmp(1:end-1),'r--')


% 14-12-2004. Plot of ln(w) from the HDF
% 

data = hdfread([namestr,'_lpqr.hdf'],namestr,'Fields','dw,age','firstrecord',1);
% [zeta,age2] = deal(data{1:2});
DW = shiftdim([data{1}],-1);
age = data{2};

agebins = (0:.1:6)*tau_eta/dt
 ind = bindex(age,agebins,0); 
 hf = figure
 % [n,x] = deal(zeros(20,length(agebins)));
 for i = 1:length(agebins)
     tmp1(i) = mean(DW(1,1,(ind == i)),3);
     tmp2(i) = mean(DW(1,2,(ind == i)),3);
     tmp3(i) = mean(DW(1,3,(ind == i)),3);
     tmp4(i) = mean(DW(1,1,(ind == i)).*DW(1,2,(ind == i)).*DW(1,3,(ind == i)));
 end
 plot(agebins*dt/tau_eta,log(tmp1),agebins*dt/tau_eta,log(tmp2),agebins*dt/tau_eta,log(tmp3))
 % legend('\beta^2_1','\beta^2_2','\beta^2_3')
 % legend('\tau_{\eta} \leq 1','1 \leq \tau_{\eta} \leq 2','2 \leq \tau_{\eta} \leq 4','4 \leq \tau_{\eta} \leq 6')
 xlabel('\tau_{\eta}')
 saveas(hf,['mean_lnw_',hfile(1:end-4)],'fig')


% PDF(Q(W)) for different time moments
%

style = {'r-','g-x','b--s','k-.o','m:<','r-.','b-o','g--<','k->','m:'};
data = hdfread([namestr,'_lpqr.hdf'],namestr,'Fields','qw,age','firstrecord',1);
% [zeta,age2] = deal(data{1:2});
QW = data{1}';
age = data{2};
agebins = (0:1:6)*tau_eta/dt
ind = bindex(age,agebins,0); 
hf = figure,
hold on
[n,x] = deal(zeros(20,length(agebins)));
for i = 1:length(agebins)
    tmp = QW((ind == i),1);
    if ~isempty(tmp)
        [n(:,i),x(:,i)] = nhist(tmp(:),20);
        plot(x(:,i),n(:,i),style{i});
    end
end
hold off
legend('\tau_{\eta} \leq 1','1 \leq \tau_{\eta} \leq 2','2 \leq \tau_{\eta} \leq 4','4 \leq \tau_{\eta} \leq 6','6 \leq \tau_{\eta}')
xlabel('Q(W)')
saveas(hf,['pdf_QW_',namestr],'fig')


data = hdfread([namestr,'_lpqr.hdf'],namestr,'Fields','dt,age','firstrecord',1);
DT = shiftdim([data{1}],-1);
age = data{2};
agebins = (0:1:6)*tau_eta/dt;
ind = bindex(age,agebins,0);


% PDF of DT_1 at different time moments
hf = figure,
hold on
style = {'r-','g-x','b--s','k-.o','m:<','r-.','b-o','g--<','k->','m:'};
[n,x] = deal(zeros(20,length(agebins)));
for i = 1:length(agebins)
    tmp = DT(:,1,(ind == i));
    if ~isempty(tmp)
        [n(:,i),x(:,i)] = nhist(tmp(:),20);
        plot(x(:,i),n(:,i),style{i});
    end
end
hold off
legend('\tau_{\eta} \leq 1','1 \leq \tau_{\eta} \leq 2','2 \leq \tau_{\eta} \leq 4','4 \leq \tau_{\eta} \leq 6','6 \leq \tau_{\eta}')
xlabel('T_1')
saveas(hf,['pdf_lambdaT_1_',namestr],'fig')


return


% % [Wls,Wlu,Wlw] = deal(repmat(0,[N,3,numPoints]));
% [coslWls,coslWlu,coslWlw,LLS] = deal(repmat(0,[N,numPoints]));
% for i = 1:numPoints
% %     Wls(1:N,1:3,i) = [L(:,:,i)*s(i,1:3)', L(:,:,i)*s(i,4:6)', L(:,:,i)*s(i,7:9)'];
% %     Wlu(1:N,1:3,i) = [L(:,:,i)*dudx(i,1:3)', L(:,:,i)*dudx(i,4:6)', L(:,:,i)*dudx(i,7:9)'];
%     Wls(1:N,1:3) = [L(:,:,i)*s(i,1:3)', L(:,:,i)*s(i,4:6)', L(:,:,i)*s(i,7:9)'];
%     LLS(1:N,i) = dot(L(1:N,1:3,i),Wls(1:N,1:3),2);
%     Wlu(1:N,1:3) = [L(:,:,i)*dudx(i,1:3)', L(:,:,i)*dudx(i,4:6)', L(:,:,i)*dudx(i,7:9)'];
%     Omega = dudx(i,:) - s(i,:);
%     Wlw(1:N,1:3) = [L(:,:,i)*Omega(1:3)', L(:,:,i)*Omega(4:6)', L(:,:,i)*Omega(7:9)'];
%     
% 
%     coslWls(1:N,i)  = cosine(L(1:N,1:3,i),Wls(1:N,1:3));
%     coslWlu(1:N,i)  = cosine(L(1:N,1:3,i),Wlu(1:N,1:3));
%     coslWlw(1:N,i)  = cosine(L(1:N,1:3,i),Wlw(1:N,1:3));
% end
% 
% LLS_L2 = LLS./squeeze(sum(L.^2,2));


% -------------------------------------------------------------------------
% PDF of omega^2 and s^2
% w2 = cat(1,traj.w1).^2+cat(1,traj.w2).^2+cat(1,traj.w3).^2;
% 
% 
% s2 = (cat(1,traj.s11).^2+2*cat(1,traj.s12).^2+2*cat(1,traj.s13).^2 + ...
%     cat(1,traj.s22).^2 + 2*cat(1,traj.s23).^2  + cat(1,traj.s33).^2);
% 
% [n1,x1] = nhist(w2);
% [n2,x2] = nhist(s2);    
% 
% hf = figure
% semilogy(x1,n1,'b-',x2,n2,'r--','linewidth',1,'markersize',4);
% legend('\langle\omega^2\rangle','\langles^2\rangle');
% saveas(hf,['pdf_w2s2_',namestr],'fig')






% -------------------------------------------------------------------------
% Different time moments, Binning:
% % 
% dt = 1/60;
% % agebins = [-0.5 0.5 1.5 2.5 3.5 4.5 5.5 6.5]*tau_eta/dt
% agebins = [0 1 2 4 6]*tau_eta/dt
% 
% % agebins = [15 30 60 90 120];  %*tau_eta/dt;
% for i = 1:length(traj)
%     traj(i).age = 0:traj(i).trajLen-1;
% end
% age = cat(2,traj.age);
% ind = bindex(age,agebins,1);    

% -------------------------------------------------------------------------
% cos(l, \omega) for the whole field
coslw = zeros(numPoints,N);
for i = 1:N,
    coslw(1:numPoints,i) = cosine(squeeze(L(i,1:3,1:numPoints))',omega(1:numPoints,1:3));
end


% -------------------------------------------------------------------------
% PDF of omega^2/s^2 for different time moments, one of the Beat's checks
hf = figure,
hold on
[n,x] = deal(zeros(length(agebins),20));
clear tmp
for i = 1:length(agebins)
    tmp(i) = mean(sum(omega((ind == i)).^2,2))./mean(sum(s(ind == i).^2,2));
    agebins(i),length(tmp)
end
plot(agebins/tau_eta*dt,tmp,'o')
hold off
%legend('\tau_{\eta} \leq 1','1 \leq \tau_{\eta} \leq 2','2 \leq \tau_{\eta} \leq 4','4 \leq \tau_{\eta} \leq 6')
xlabel('\tau_{\eta}')
ylabel('\langle\omega^2\rangle/\langles^2\rangle')
saveas(hf,['mean_w2overs2_',namestr],'fig')
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% PDF of lls at different time moments
hf = figure,
hold on
style = {'r-','b-x','g--s','k-.o','y.'}
[n,x] = deal(zeros(20,length(agebins)));
for i = 1:length(agebins)-1
    tmp = LLS(:,(ind == i));
%     agebins(i),length(tmp)
    if ~isempty(tmp)
        [n(:,i),x(:,i)] = nhist(tmp(tmp < 500),20);
        semilogy(x(:,i),n(:,i),style{i});
    end
end
hold off
legend('\tau_{\eta} \leq 1','1 \leq \tau_{\eta} \leq 2','2 \leq \tau_{\eta} \leq 4','4 \leq \tau_{\eta} \leq 6')
xlabel('l_il_js_{ij}')
saveas(hf,['pdf_lls_',namestr],'fig')


% -------------------------------------------------------------------------
% PDF of lls/l_2 at different time moments
hf = figure,
hold on
style = {'r-','b-x','g--s','k-.o','m.'};
[n,x] = deal(zeros(20,length(agebins)));
for i = 1:length(agebins)
    tmp = LLS_L2(:,(ind == i));
    if ~isempty(tmp)
        [n(:,i),x(:,i)] = nhist(tmp(:),20);
        plot(x(:,i),n(:,i),style{i});
    end
end
hold off
legend('\tau_{\eta} \leq 1','1 \leq \tau_{\eta} \leq 2','2 \leq \tau_{\eta} \leq 4','4 \leq \tau_{\eta} \leq 6')
xlabel('l_il_js_{ij}/l^2')
saveas(hf,['pdf_lls_l2_',namestr],'fig')


    




% --------------------------------------------------------------------------
% PDF of cos(l,\omega) at different time moments
hf = figure,
hold on
style = {'r-','b-x','g--s','k-.o','m.','gx'}
[n,x] = deal(zeros(20,length(agebins)));
for i = 1:length(agebins)
    tmp = coslw((ind == i),:);
    agebins(i),length(tmp)
    if ~isempty(tmp)
        [n(:,i),x(:,i)] = nhist(tmp(:),20);
        plot(x(:,i),n(:,i),style{i});
    end
end
hold off
legend('\tau_{\eta} \leq 1','1 \leq \tau_{\eta} \leq 2','2 \leq \tau_{\eta} \leq 4','4 \leq \tau_{\eta} \leq 6','6 \leq \tau_{\eta}')
xlabel('cos(l,\omega)')
saveas(hf,['pdf_coslomega_',namestr],'fig')


% -------------------------------------------------------------------------
% PDF of cos(l,Wls) at different time moments

hf = figure,
hold on
style = {'r-','b-x','g--s','k-.o','y.'}
[n,x] = deal(zeros(20,length(agebins)));
for i = 1:length(agebins)
    tmp = coslWls(:,(ind == i));
    agebins(i),length(tmp)
    if ~isempty(tmp)
        [n(:,i),x(:,i)] = nhist(tmp(:),20);
        plot(x(:,i),n(:,i),style{i});
    end
end
hold off
legend('\tau_{\eta} \leq 1','1 \leq \tau_{\eta} \leq 2','2 \leq \tau_{\eta} \leq 4','4 \leq \tau_{\eta} \leq 6','6 \leq \tau_{\eta}')
xlabel('cos(l,W^{l,s})')
saveas(hf,['pdf_coslWls_',namestr],'fig')

% -------------------------------------------------------------------------
% PDF of cos(l,Wls) at different time moments, strong s

hf = figure,
hold on
style = {'r-','b-x','g--s','k-.o','y.'}
[n,x] = deal(zeros(20,length(agebins)));
for i = 1:length(agebins)
    tmp = coslWls(:,[(ind == i) & strong_s']);
    agebins(i),length(tmp)
    if ~isempty(tmp)
        [n(:,i),x(:,i)] = nhist(tmp(:),20);
        plot(x(:,i),n(:,i),style{i});
    end
end
hold off
legend('\tau_{\eta} \leq 1','1 \leq \tau_{\eta} \leq 2','2 \leq \tau_{\eta} \leq 4','4 \leq \tau_{\eta} \leq 6','6 \leq \tau_{\eta}')
xlabel('cos(l,W^{l,s}), s^2 > 3 \langle s^2 \rangle')
saveas(hf,['pdf_coslWls_strongs_',namestr],'fig')

% -------------------------------------------------------------------------
% PDF of cos(l,Wls) at different time moments, weak s

hf = figure,
hold on
style = {'r-','b-x','g--s','k-.o','m.'}
[n,x] = deal(zeros(20,length(agebins)));
for i = 1:length(agebins)
    tmp = coslWls(:,[(ind == i) & weak_s']);
    agebins(i),length(tmp)
    if ~isempty(tmp)
        [n(:,i),x(:,i)] = nhist(tmp(:),20);
        plot(x(:,i),n(:,i),style{i});
    end
end
hold off
legend('\tau_{\eta} \leq 1','1 \leq \tau_{\eta} \leq 2','2 \leq \tau_{\eta} \leq 4','4 \leq \tau_{\eta} \leq 6','6 \leq \tau_{\eta}')
xlabel('cos(l,W^{l,s}), s^2 < 1/3 \langle s^2 \rangle')
saveas(hf,['pdf_coslWls_weaks_',namestr],'fig')


% -------------------------------------------------------------------------
% PDF of cos(l,Wlu) at different time moments

hf = figure,
hold on
style = {'r-','b-x','g--s','k-.o','y.'}
[n,x] = deal(zeros(20,length(agebins)));
for i = 1:length(agebins)
    tmp = coslWlw(:,(ind == i));
    agebins(i),length(tmp)
    if ~isempty(tmp)
        [n(:,i),x(:,i)] = nhist(tmp(:),20);
        plot(x(:,i),n(:,i),style{i});
    end
end
hold off
legend('\tau_{\eta} \leq 1','1 \leq \tau_{\eta} \leq 2','2 \leq \tau_{\eta} \leq 4','4 \leq \tau_{\eta} \leq 6','6 \leq \tau_{\eta}')
xlabel('cos(l,W^{l,\Omega})')
saveas(hf,['pdf_coslWlw_',namestr],'fig')



% -------------------------------------------------------------------------
% PDF of cos(l,Wlu) at different time moments

hf = figure,
hold on
style = {'r-','b-x','g--s','k-.o','m.'}
[n,x] = deal(zeros(20,length(agebins)));
for i = 1:length(agebins)
    tmp = coslWlu(:,(ind == i));
    agebins(i),length(tmp)
    if ~isempty(tmp)
        [n(:,i),x(:,i)] = nhist(tmp(:),20);
        plot(x(:,i),n(:,i),style{i});
    end
end
hold off
legend('\tau_{\eta} \leq 1','1 \leq \tau_{\eta} \leq 2','2 \leq \tau_{\eta} \leq 4','4 \leq \tau_{\eta} \leq 6','6 \leq \tau_{\eta}')
xlabel('cos(l,W^{l,\nabla u})')
saveas(hf,['pdf_coslWlu_',namestr],'fig')

% -------------------------------------------------------------------------
% PDF of cos(l,Wlu) at different time moments, strong S condition

hf = figure,
hold on
style = {'r-','b-x','g--s','k-.o','m--.'}
[n,x] = deal(zeros(20,length(agebins)));
for i = 1:length(agebins)
    tmp = coslWlu(:,[(ind == i) & strong_s']);
    agebins(i),length(tmp)
    if ~isempty(tmp)
        [n(:,i),x(:,i)] = nhist(tmp(:),20);
        plot(x(:,i),n(:,i),style{i});
    end
end
hold off
legend('\tau_{\eta} \leq 1','1 \leq \tau_{\eta} \leq 2','2 \leq \tau_{\eta} \leq 4','4 \leq \tau_{\eta} \leq 6','6 \leq \tau_{\eta}')
xlabel('cos(l,W^{l,\nabla u}), s^2 > 3 \langle s^2 \rangle')
saveas(hf,['pdf_coslWlu_strongs_',namestr],'fig')

% -------------------------------------------------------------------------
% PDF of cos(l,Wlu) at different time moments, weak S condition

hf = figure,
hold on
style = {'r-','b-x','g--s','k-.o','m--.'}
[n,x] = deal(zeros(20,length(agebins)));
for i = 1:length(agebins)
    tmp = coslWlu(:,[(ind == i) & weak_s']);
    agebins(i),length(tmp)
    if ~isempty(tmp)
        [n,x] = nhist(tmp(:),20);
        plot(x,n,style{i});
    end
end
hold off
legend('\tau_{\eta} \leq 1','1 \leq \tau_{\eta} \leq 2','2 \leq \tau_{\eta} \leq 4','4 \leq \tau_{\eta} \leq 6','6 \leq \tau_{\eta}')
xlabel('cos(l,W^{l,\nabla u}), s^2 < 1/3 \langle s^2 \rangle')
saveas(hf,['pdf_coslWlu_weaks_',namestr],'fig')



% ----------------------------------------------------------------------
% cos(l,\lambda_i)
[cos_l_lambda1,cos_l_lambda2,cos_l_lambda3] = deal(zeros(numPoints,10));
for i = 1:10
    cos_l_lambda1(1:numPoints,i) = cosine(squeeze(L(i,:,:))',squeeze(V(1:3,1,:))');
    cos_l_lambda2(1:numPoints,i) = cosine(squeeze(L(i,:,:))',squeeze(V(1:3,2,:))'); 
    cos_l_lambda3(1:numPoints,i) = cosine(squeeze(L(i,:,:))',squeeze(V(1:3,3,:))'); 
end

% -------------------------------------------------------------------------
% cos(l, lambda_1)
hf = figure,
hold on
style = {'r-','b-x','g--s','k-.o','m--.'}
[n,x] = deal(zeros(20,N));
for i = 1:length(agebins)
    for j = 1:N
        tmp = cos_l_lambda1([(ind == i)],j);
%         agebins(i),length(tmp)
        if ~isempty(tmp)
            [n(:,j),x(:,j)] = nhist(tmp(:),20);
            
        end
    end
    plot(x(:,1),mean(n,2),style{i});
end
hold off
legend('\tau_{\eta} \leq 1','1 \leq \tau_{\eta} \leq 2','2 \leq \tau_{\eta} \leq 4','4 \leq \tau_{\eta} \leq 6','6 \leq \tau_{\eta}')
xlabel('cos(l,\lambda_1)')
saveas(hf,['pdf_cos_l_lambda1_',namestr],'fig')

% -------------------------------------------------------------------------
% cos(l, lambda_2)
hf = figure,
hold on
style = {'r-','b-x','g--s','k-.o','m--.'}
[n,x] = deal(zeros(20,N));
for i = 1:length(agebins)
    for j = 1:N
        tmp = cos_l_lambda2([(ind == i)],j);
%         agebins(i),length(tmp)
        if ~isempty(tmp)
            [n(:,j),x(:,j)] = nhist(tmp(:),20);
            
        end
    end
    plot(x(:,1),mean(n,2),style{i});
end
hold off
legend('\tau_{\eta} \leq 1','1 \leq \tau_{\eta} \leq 2','2 \leq \tau_{\eta} \leq 4','4 \leq \tau_{\eta} \leq 6','6 \leq \tau_{\eta}')
xlabel('cos(l,\lambda_2)')
saveas(hf,['pdf_cos_l_lambda2_',namestr],'fig')
% -------------------------------------------------------------------------
% cos(l, lambda_3)
hf = figure,
hold on
style = {'r-','b-x','g--s','k-.o','m--.'}
[n,x] = deal(zeros(20,N));
for i = 1:length(agebins)
    for j = 1:N
        tmp = cos_l_lambda3([(ind == i)],j);
%         agebins(i),length(tmp)
        if ~isempty(tmp)
            [n(:,j),x(:,j)] = nhist(tmp(:),20);
            
        end
    end
    plot(x(:,1),mean(n,2),style{i});
end
hold off
legend('\tau_{\eta} \leq 1','1 \leq \tau_{\eta} \leq 2','2 \leq \tau_{\eta} \leq 4','4 \leq \tau_{\eta} \leq 6','6 \leq \tau_{\eta}')
xlabel('cos(l,\lambda_3)')
saveas(hf,['pdf_cos_l_lambda3_',namestr],'fig')




% -------------------------------------------------------------------------
% Joint PDF  of s^2 versus omega^2
hf = figure
[n,x,bins] = histmulti5([s2(:),w2(:)],[(0:.5:70)',(0:1:140)']);
contour(x(:,1),x(:,2),n')
axis tight
xlabel('s^2')
ylabel('\omega^2')
grid on
axis square
saveas(hf,['jointpdf_s2w2_',namestr],'fig')


% -------------------------------------------------------------------------
% Joint PDF  of lls/l^2 versus s^2 for different time moments
s2 = sum(s.^2,2);

tau = [1,2,4,6];
for i = 1:length(agebins)
    tmp1 = LLS_L2(:,ind == i);
    tmp2 = s2(ind == i);
    if ~isempty(tmp1) & ~isempty(tmp2)
            [n,x,bins] = histmulti5([tmp2(:),tmp1(1,:)']);
        for j = 2:N
            [n1,x,bins] = histmulti5([tmp2(:),tmp1(j,:)']);
            n = n + n1; 
        end    
        hf = figure
        contour(x(:,1),x(:,2),n')
        axis tight
        xlabel('s^2')
        ylabel(sprintf('lls/l^2, \\tau = %d',tau(i)))
        grid on
        axis square
        colorbar vert
        saveas(hf,sprintf('jointpdf_llsl2_s2_tau_%d_%s',tau(i),namestr),'fig')    
    end
end
% legend('\tau_{\eta} \leq 1','1 \leq \tau_{\eta} \leq 2','2 \leq \tau_{\eta} \leq 4','4 \leq \tau_{\eta} \leq 6')
% xlabel('cos(l,W^{l,\nabla u}), s^2 < 1/3 \langle s^2 \rangle')




% -----------------------------------------------------
% Mean stretching rates - plot of quantities against time
% first binning into 20 bins from 0 to 6 Kolmogorov times
% i.e. from 0 to 90 (for water) frames
% We have above this binning, use 'ind' variable
dt = 1/60;
agebins = 0:2:10*tau_eta/dt;
% agebins = 0:20:maxLen;
%
for i=1:numTraj, traj(i).age = 0:trajLen(i)-1; end
ind = bindex(cat(2,traj.age),agebins,1);    


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -----------------------------------------------------------------------
% \langle \omega^2 \rangle / \langle s^2 \rangle for different time moments
hf = figure,
style = {'r-','b-x','g--s','k-.o','y.'}
tmp = zeros(size(agebins));
for i = 1:length(agebins)
    tmp(i) = mean(sum(omega(ind == i,:).^2,2))/mean(sum(s(ind == i,:).^2,2));
end
plot(agebins(1:end-1)/tau_eta*dt,tmp(1:end-1),'o')
ylabel('\langle \omega^2 \rangle / \langle s^2 \rangle')
xlabel('\tau_{\eta}')
saveas(hf,['mean_w2overs2_tau_',namestr],'fig')




% ---------------------------------------------------------------------
hf = figure,
style = {'r-','b-x','g--s','k-.o','y.'}
tmp = zeros(size(agebins));
for i = 1:length(agebins)
    tmp(i) = mean(sss((ind == i)));
end
plot(agebins(1:end-1)/tau_eta*dt,tmp(1:end-1))
ylabel('s_{ij}s_{jk}s_{ki}','fontname','euclid','fontangle','italic')
xlabel('\tau_{\eta}')
saveas(hf,['mean_sss_tau_',namestr],'fig')



hf = figure,
style = {'r-','b-x','g--s','k-.o','y.'}
tmp = zeros(size(agebins));
for i = 1:length(agebins)
    tmp(i) = mean(wws(find(ind == i)));
end
plot(agebins(1:end-1)/tau_eta*dt,tmp(1:end-1))
ylabel('\omega_{i} \omega_j s_{ij}','fontname','euclid','fontangle','italic')
xlabel('\tau_{\eta}')
saveas(hf,['mean_wws_tau_',namestr],'fig')


hf = figure,
style = {'r-','b-x','g--s','k-.o','y.'}
tmp = zeros(size(agebins));
for i = 1:length(agebins)
    tmp(i) = mean(mean(LLS(:,(ind == i))));
end
plot(agebins(1:end-1)/tau_eta*dt,tmp(1:end-1))
ylabel('l_{i} l_j s_{ij}','fontname','euclid','fontangle','italic')
xlabel('\tau_{\eta}')
saveas(hf,['mean_lls_tau_',namestr],'fig')

% ----------------------------------------------------------------------
% dlnLdt(:,:,end+1) = dlnLdt(:,:,end);

hf = figure,
tmp = zeros(size(agebins));
for i = 1:length(agebins)
%    tmp(i) = mean(mean(dlnLdt(1:N,1,(ind == i)),1),3);
        tmp(i) = mean(mean(zeta(1:N,(ind == i & inddiv')),1),2);
        %tmp2(i) = mean(mean(zeta(2,(ind == i)),1),2);
        %tmp3(i) = mean(mean(zeta(3,(ind == i)),1),2);


end
plot(agebins(1:end-1)/tau_eta*dt,tmp(1:end-1))
ylabel('dln(l)/dt','fontname','euclid','fontangle','italic')
xlabel('\tau_{\eta}')
saveas(hf,['mean_zeta_tau_',namestr],'fig')




% sss/s2


hf = figure,
style = {'r-','b-x','g--s','k-.o','y.'}
tmp = zeros(size(agebins));
for i = 1:length(agebins)
    tmp(i) = mean(sss_s2(ind == i));
end
plot(agebins(1:end-1)/tau_eta*dt, -tmp(1:end-1),style{1})
ylabel('s_{ij}s_{jk}s_{ki}/s^2','fontname','euclid','fontangle','italic')
xlabel('\tau_{\eta}')
saveas(hf,['mean_sss_s2_tau_',namestr],'fig')

hf = figure,
tmp = zeros(size(agebins));
for i = 1:length(agebins)
    tmp(i) = mean(wws_w2((ind == i)));
end
plot(agebins(1:end-1)/tau_eta*dt,tmp(1:end-1),style{2})
ylabel('\omega_{i} \omega_j s_{ij}/\omega^2','fontname','euclid','fontangle','italic')
xlabel('\tau_{\eta}')
saveas(hf,['mean_wws_w2_tau_',namestr],'fig')


% hf = figure,
% style = {'r-','b-x','g--s','k-.o','y.'}
% tmp = zeros(size(agebins));
% for i = 1:length(agebins)
%     tmp(i) = mean(mean(lls_l2(:,(ind == i))));
% end
% plot(agebins(1:end-1)/tau_eta*dt,tmp(1:end-1),style{3})
% ylabel('l_{i} l_j s_{ij}/l^2','fontname','euclid','fontangle','italic')
% xlabel('\tau_{\eta}')
% saveas(hf,['mean_lls_l2_tau_',namestr],'fig')



