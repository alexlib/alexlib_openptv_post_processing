curdir = cd;
hfile = [curdir(1),':\Matlab\Alex\HDF\0104_water.hdf'];

%% reading from HDF4 file
[traj,attr] = readTrajHDF_v9(hfile,'ax',[],'ay',[],'az',[],'w1',[],'w2',[],'w3',[],...
    's11',[],'s12',[],'s13',[],'s22',[],'s23',[],'s33',[],'t',[],'minlength',20);


%% Rearrange the derivatives into tensor per each point
dudx  =   cat(1,traj.s11);
dudy  =   cat(1,traj.s12) - 0.5 * cat(1,traj.w3);
dudz  =   cat(1,traj.s13) + 0.5 * cat(1,traj.w2);
dvdx  =   cat(1,traj.s12) + 0.5 * cat(1,traj.w3);
dvdy  =   cat(1,traj.s22);
dvdz  =   cat(1,traj.s23) - 0.5 * cat(1,traj.w1);
dwdx  =   cat(1,traj.s13) - 0.5 * cat(1,traj.w2);
dwdy  =   cat(1,traj.s23) + 0.5 * cat(1,traj.w1);
dwdz  =   cat(1,traj.s33);


numTraj = length(traj)
numPoints = length(dudx)
trajLen = cat(1,traj.trajLen);


%% selecting only 'best' results
%{
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
    if length(ind) > 15
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

% Rearrange the derivatives into tensor per each point
dudx  =   cat(1,traj.s11);
dudy  =   cat(1,traj.s12) - 0.5 * cat(1,traj.w3);
dudz  =   cat(1,traj.s13) + 0.5 * cat(1,traj.w2);
dvdx  =   cat(1,traj.s12) + 0.5 * cat(1,traj.w3);
dvdy  =   cat(1,traj.s22);
dvdz  =   cat(1,traj.s23) - 0.5 * cat(1,traj.w1);
dwdx  =   cat(1,traj.s13) - 0.5 * cat(1,traj.w2);
dwdy  =   cat(1,traj.s23) + 0.5 * cat(1,traj.w1);
dwdz  =   cat(1,traj.s33);


numTraj = length(traj)
numPoints = length(dudx)
trajLen = cat(1,traj.trajLen);
trajnum = unique(cat(1,traj.trajnum));

save([hfile(1:end-4),'_good'],'good','numTraj','numPoints','trajnum','trajLen');

% clear traj

maxLen    = max(trajLen)
minLen    = min(trajLen)
meanLen    = mean(trajLen)

start_points = [1; cumsum(trajLen(1:end-1))+1];     % index of the first trajectory points

% reldiv = abs(dudx + dvdy + dwdz)./(abs(dudx)+abs(dvdy)+abs(dwdz));
% absdiv = abs(dudx + dvdy + dwdz);
% inddiv = absdiv < 0.1 & reldiv < 0.1;

%}

%%
sSq = dudx.^2 + dvdy.^2 + dwdz.^2 + 2*(0.5*(dudy + dvdx).^2 + 0.5*(dudz + dwdx).^2 + 0.5*(dvdz + dwdy).^2);
diss = mean(sSq)*2*1e-6
% diss = mean(sSq(inddiv))*2*1e-6

tau_eta = (1.0e-6/diss)^.5
eta = ((1.0e-6)^3/diss)^0.25
% --------------------------------------------------------------------
inds = sSq < prctile(sSq,99);
mean_ssq = mean(sSq(inds))
strong_s = (sSq > 3*mean_ssq) & inds;
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
% saveas(hf,['pdf_s2_enstrophy_',hfile(1:end-4)],'fig')


mean_omegasq = mean(enstrophy(:))
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
saveas(hf,['pdf_cos_sS_',hfile(1:end-4)],'fig')


hf = figure
nhist(cos_sS(strong_s),25,'b:o');
xlabel('cos(s,S), > 3 \langles^2\rangle')
saveas(hf,['pdf_cos_sS_strongs_',hfile(1:end-4)],'fig')

hf = figure
nhist(cos_sS(weak_s),25,'g-.x');
xlabel('cos(s,S), < 1/3 \langles^2\rangle')
saveas(hf,['pdf_cos_sS_omega_',hfile(1:end-4)],'fig')

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
saveas(hf,['pdf_cos_sS_conditions_',hfile(1:end-4)],'fig')


% -------------------------------------------------------------------------
% w_iS_{ij}
W = repmat(0,[numPoints 3]);
for j = 1:numPoints
    % tmp = reshape(s(i,1:9),[3 3]);
    W(j,1:3) = omega(j,1:3)*[0.5*(dudx(j) + dudx(j)), 0.5*(dudy(j) + dvdx(j)), 0.5*(dudz(j) + dwdx(j));
        0.5*(dudy(j) + dvdx(j)),0.5*(dvdy(j) + dvdy(j)),0.5*(dvdz(j) + dwdy(j));
        0.5*(dudz(j) + dwdx(j)),0.5*(dvdz(j) + dwdy(j)),0.5*(dwdz(j) + dwdz(j))]';
end

feval('save',[hfile(1:end-4),'_wws_sS.mat'],'W','cos_sS','trajLen','numPoints','numTraj');

% coswW = dot(omega,W,2)./sqrt(sum(W.^2,2))./sqrt(sum(omega.^2,2));
coswW = cosine(omega,W);

% PDF of the cos(w,W)
hf = figure
[n,x] = nhist(coswW);
plot(x,n,'bo');
hold on
fnplt(csapi(x,n),'b:',[],0.5)
hold off
xlabel('cos(\omega,W), unconditioned')
saveas(hf,['pdf_cos_wW_',hfile(1:end-4)],'fig')


hf = figure;
[n,x] = nhist(coswW,50);
plot(x,n,'r-');
hold on
[n,x] = nhist(coswW(strong_s),50);
plot(x,n,'b--');
[n,x] = nhist(coswW(weak_s),50);
plot(x,n,'k-.');
hold off
xlabel('cos(\omega,W)');
ylabel('PDF')
legend('uncond','strong','weak',2)
saveas(hf,['pdf_coswW_',hfile(1:end-4)],'fig')

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
saveas(hf,['pdf_coswW_conditions_',hfile(1:end-4)],'fig')




% -------------------------------------------------------------------------
% cos(w,\lambda_i)
cos_wlambda1 = cosine(omega(:,:),squeeze(V(1:3,1,:))');
cos_wlambda2 = cosine(omega(:,:),squeeze(V(1:3,2,:))');
cos_wlambda3 = cosine(omega(:,:),squeeze(V(1:3,3,:))');


[n1,x1] = nhist(cos_wlambda1(inddiv),20);
[n2,x2] = nhist(cos_wlambda2(inddiv),20);
[n3,x3] = nhist(cos_wlambda3(inddiv),20);
hf = figure, hold on
plot(x1,n1,'ro',x2,n2,'rs',x3,n3,'r<','linewidth',1);
legend('cos(\omega,\lambda_1)','cos(\omega,\lambda_2)','cos(\omega,\lambda_3)')
for i=1:3
    eval(sprintf('fnplt(csapi(x%d,n%d),''%s:'',[],0.5)',i,i,'r'))
end

xlabel('cos(\omega,\lambda_i)')
saveas(hf,['pdf_cos_w_lambda_',hfile(1:end-4)],'fig')

[n1,x1] = nhist(cos_wlambda1(strong_s));
[n2,x2] = nhist(cos_wlambda2(strong_s));
[n3,x3] = nhist(cos_wlambda3(strong_s));
hf = figure
plot(x1,n1,'r-',x2,n2,'r-s',x3,n3,'r-x','linewidth',1);
xlabel('cos(\omega,\lambda_i) > 3 \langles^2\rangle')
saveas(hf,['pdf_cos_w_lambda_strongs_',hfile(1:end-4)],'fig')

[n1,x1] = nhist(cos_wlambda1(weak_s));
[n2,x2] = nhist(cos_wlambda2(weak_s));
[n3,x3] = nhist(cos_wlambda3(weak_s));
hf = figure
plot(x1,n1,'r-',x2,n2,'r-s',x3,n3,'r-x','linewidth',1);
xlabel('cos(\omega,\lambda_i) < 1/3 \langles^2\rangle')
saveas(hf,['pdf_cos_w_lambda_weaks_',hfile(1:end-4)],'fig')

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
saveas(hf,['pdf_cosw_lambda1_conditions_',hfile(1:end-4)],'fig')

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
saveas(hf,['pdf_cosw_lambda2_conditions_',hfile(1:end-4)],'fig')

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
saveas(hf,['pdf_cosw_lambda3_conditions_',hfile(1:end-4)],'fig')

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
saveas(hf,['pdf_wws_conditions_',hfile(1:end-4)],'fig')

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
saveas(hf,['pdf_wws_w2_conditions_',hfile(1:end-4)],'fig')


%hf = figure,
%nhist(wws_w2);
%set(gca,'yscale','log')
%xlabel('\omega_{i}\omega_{j}s_{ij}/\omega^2','fontname','euclid','fontangle','italic')
%saveas(hf,['pdf_wws_w2_',hfile(1:end-4)],'fig')
%
%% conditioned wws_w2, strong omega
%hf = figure,
%nhist(wws_w2(strong_omega));
%set(gca,'yscale','log')
%xlabel('\omega_{i}\omega_{j}s_{ij}/\omega^2, > 3 \langle\omega^2\rangle','fontname','euclid','fontangle','italic')
%saveas(hf,['pdf_wws_w2_strong_w',hfile(1:end-4)],'fig')
%
%% conditioned wws_w2, weak omega
%hf = figure,
%nhist(wws_w2(weak_omega));
%set(gca,'yscale','log')
%xlabel('\omega_{i}\omega_{j}s_{ij}/\omega^2, < 1/3 \langle\omega^2\rangle','fontname','euclid','fontangle','italic')
%saveas(hf,['pdf_wws_w2_weak_w',hfile(1:end-4)],'fig')
%% -------------------------------------------------------------------------
%% One special figure is low-pressure  = strong enstrophy and weak
%% dissipation
%hf = figure,
%nhist(wws_w2(strong_omega & weak_s));
%set(gca,'yscale','log')
%xlabel('\omega_{i}\omega_{j}s_{ij}/\omega^2, > 3\langle\omega^2\rangle & < 1/3 \langles^2\rangle','fontname','euclid','fontangle','italic')
%saveas(hf,['pdf_wws_w2_strong_w_weak_s',hfile(1:end-4)],'fig')
%
%% -------------------------------------------------------------------------
%% Second special figure is low-pressure  = strong dissipation and weak
%% enstrophy
%hf = figure,
%nhist(wws_w2(strong_s & weak_omega));
%set(gca,'yscale','log')
%xlabel('\omega_{i}\omega_{j}s_{ij}/\omega^2, < 1/3 \langle\omega^2\rangle & > 3 \langles^2\rangle','fontname','euclid','fontangle','italic')
%saveas(hf,['pdf_wws_w2_strong_s_weak_w_',hfile(1:end-4)],'fig')


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
sss_s2(i) = - sss(i)./sSq(i);

feval('save',[hfile(1:end-4),'_wws_sss.mat'],'wws','wws_w2','sss','sss_s2');


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
saveas(hf,['pdf_sss_conditions_',hfile(1:end-4)],'fig')

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
saveas(hf,['pdf_sss_s2_conditions_',hfile(1:end-4)],'fig')



%hf = figure,
%nhist(sss_s2);
%set(gca,'yscale','log')
%xlabel('-s_{ij}s_{ji}s_{ki}/s^2','fontname','euclid','fontangle','italic')
%saveas(hf,['pdf_sss_s2_',hfile(1:end-4)],'fig')
%
%% conditioned sss_s2
%hf = figure,
%nhist(sss_s2(strong_s));
%set(gca,'yscale','log')
%xlabel('s_{ij}s_{ij}s_{ij}/s^2, > 3 \langles^2\rangle','fontname','euclid','fontangle','italic')
%saveas(hf,['pdf_sss_s2_strongs',hfile(1:end-4)],'fig')
%
%% conditioned sss_s2
%hf = figure,
%nhist(sss_s2(weak_s));
%set(gca,'yscale','log')
%xlabel('s_{ij}s_{jk}s_{ki}/s^2,  < 1/3 \langles^2\rangle','fontname','euclid','fontangle','italic')
%saveas(hf,['pdf_sss_s2_weaks',hfile(1:end-4)],'fig')


% -------------------------------------------------------------------------
%calculate_lls(traj,10);
%feval('load',[hfile(1:end-4),'_lines'],'lB','lls','Wls','Wlu','coslWlu','coslWls','lls_l2');


% New piece: plot of the values from the HDF files:
data = hdfread([hfile(1:end-4),'_lpqr.hdf'],hfile(1:end-4),'Fields','coslwls,coslwlu,zeta,age','firstrecord',1);
coslWls = data{1};
coslWlu = data{2};
zeta = data{3};
age = data{4};


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
saveas(hf,['pdf_coslWls_',hfile(1:end-4)],'fig')


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
saveas(hf,['pdf_coslWlu_',hfile(1:end-4)],'fig')





data = hdfread('0809_water_lpqr.hdf','0809_water','Fields','cosllambda1,cosllambda2,cosllambda3','firstrecord',1);
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
saveas(hf,['pdf_cosllambda1_',hfile(1:end-4)],'fig')

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
saveas(hf,['pdf_cosllambda2_',hfile(1:end-4)],'fig')


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
saveas(hf,['pdf_cosllambda3_',hfile(1:end-4)],'fig')


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
w2 = cat(1,traj.w1).^2+cat(1,traj.w2).^2+cat(1,traj.w3).^2;


s2 = (cat(1,traj.s11).^2+2*cat(1,traj.s12).^2+2*cat(1,traj.s13).^2 + ...
    cat(1,traj.s22).^2 + 2*cat(1,traj.s23).^2  + cat(1,traj.s33).^2);

[n1,x1] = nhist(w2);
[n2,x2] = nhist(s2);    

hf = figure
semilogy(x1,n1,'b-',x2,n2,'r--','linewidth',1,'markersize',4);
legend('\langle\omega^2\rangle','\langles^2\rangle');
saveas(hf,['pdf_w2s2_',hfile(1:end-4)],'fig')






% -------------------------------------------------------------------------
% Different time moments, Binning:
% 
dt = 1/60;
% agebins = [-0.5 0.5 1.5 2.5 3.5 4.5 5.5 6.5]*tau_eta/dt
agebins = [0 1 2 4 6]*tau_eta/dt

% agebins = [15 30 60 90 120];  %*tau_eta/dt;
for i = 1:length(traj)
    traj(i).age = 0:traj(i).trajLen-1;
end
age = cat(2,traj.age);
ind = bindex(age,agebins,1);    

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
saveas(hf,['mean_w2overs2_',hfile(1:end-4)],'fig')
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
saveas(hf,['pdf_lls_',hfile(1:end-4)],'fig')


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
saveas(hf,['pdf_lls_l2_',hfile(1:end-4)],'fig')







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
saveas(hf,['pdf_coslomega_',hfile(1:end-4)],'fig')


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
saveas(hf,['pdf_coslWls_',hfile(1:end-4)],'fig')

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
saveas(hf,['pdf_coslWls_strongs_',hfile(1:end-4)],'fig')

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
saveas(hf,['pdf_coslWls_weaks_',hfile(1:end-4)],'fig')


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
saveas(hf,['pdf_coslWlw_',hfile(1:end-4)],'fig')



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
saveas(hf,['pdf_coslWlu_',hfile(1:end-4)],'fig')

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
saveas(hf,['pdf_coslWlu_strongs_',hfile(1:end-4)],'fig')

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
saveas(hf,['pdf_coslWlu_weaks_',hfile(1:end-4)],'fig')



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
saveas(hf,['pdf_cos_l_lambda1_',hfile(1:end-4)],'fig')

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
saveas(hf,['pdf_cos_l_lambda2_',hfile(1:end-4)],'fig')
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
saveas(hf,['pdf_cos_l_lambda3_',hfile(1:end-4)],'fig')




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
saveas(hf,['jointpdf_s2w2_',hfile(1:end-4)],'fig')


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
        saveas(hf,sprintf('jointpdf_llsl2_s2_tau_%d_%s',tau(i),hfile(1:end-4)),'fig')    
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
% for i=1:numTraj, traj(i).age = 0:trajLen(i)-1; end
ind = bindex(age,agebins,1);    


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
saveas(hf,['mean_w2overs2_tau_',hfile(1:end-4)],'fig')




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
saveas(hf,['mean_sss_tau_',hfile(1:end-4)],'fig')



hf = figure,
style = {'r-','b-x','g--s','k-.o','y.'}
tmp = zeros(size(agebins));
for i = 1:length(agebins)
    tmp(i) = mean(wws(find(ind == i)));
end
plot(agebins(1:end-1)/tau_eta*dt,tmp(1:end-1))
ylabel('\omega_{i} \omega_j s_{ij}','fontname','euclid','fontangle','italic')
xlabel('\tau_{\eta}')
saveas(hf,['mean_wws_tau_',hfile(1:end-4)],'fig')


hf = figure,
style = {'r-','b-x','g--s','k-.o','y.'}
tmp = zeros(size(agebins));
for i = 1:length(agebins)
    tmp(i) = mean(mean(LLS(:,(ind == i))));
end
plot(agebins(1:end-1)/tau_eta*dt,tmp(1:end-1))
ylabel('l_{i} l_j s_{ij}','fontname','euclid','fontangle','italic')
xlabel('\tau_{\eta}')
saveas(hf,['mean_lls_tau_',hfile(1:end-4)],'fig')

% ----------------------------------------------------------------------
% dlnLdt(:,:,end+1) = dlnLdt(:,:,end);

hf = figure,
tmp = zeros(size(agebins));
for i = 1:length(agebins)
    %    tmp(i) = mean(mean(dlnLdt(1:N,1,(ind == i)),1),3);
    tmp(i) = mean(mean(zeta(:,(ind == i)),1),2);
    %tmp2(i) = mean(mean(zeta(2,(ind == i)),1),2);
    %tmp3(i) = mean(mean(zeta(3,(ind == i)),1),2);
    
    
end
plot(agebins(:)/tau_eta*dt,tmp(:))
ylabel('dln(l)/dt','fontname','euclid','fontangle','italic')
xlabel('\tau_{\eta}')
saveas(hf,['mean_zeta_tau_',hfile(1:end-4)],'fig')




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
saveas(hf,['mean_sss_s2_tau_',hfile(1:end-4)],'fig')

hf = figure,
tmp = zeros(size(agebins));
for i = 1:length(agebins)
    tmp(i) = mean(wws_w2((ind == i)));
end
plot(agebins(1:end-1)/tau_eta*dt,tmp(1:end-1),style{2})
ylabel('\omega_{i} \omega_j s_{ij}/\omega^2','fontname','euclid','fontangle','italic')
xlabel('\tau_{\eta}')
saveas(hf,['mean_wws_w2_tau_',hfile(1:end-4)],'fig')


% hf = figure,
% style = {'r-','b-x','g--s','k-.o','y.'}
% tmp = zeros(size(agebins));
% for i = 1:length(agebins)
%     tmp(i) = mean(mean(lls_l2(:,(ind == i))));
% end
% plot(agebins(1:end-1)/tau_eta*dt,tmp(1:end-1),style{3})
% ylabel('l_{i} l_j s_{ij}/l^2','fontname','euclid','fontangle','italic')
% xlabel('\tau_{\eta}')
% saveas(hf,['mean_lls_l2_tau_',hfile(1:end-4)],'fig')



