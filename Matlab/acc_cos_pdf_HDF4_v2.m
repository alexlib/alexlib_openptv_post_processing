function acc_cos_pdf_HDF4_v2(hdfFile,firstFrame,lastFrame,stepFrames,threshold)

%% Prepare necessary quantities:
if nargin == 0
    nBins = 51;
    first = 1;
    last = 2990; % 2990
    stepFrames = 1;
    threshold = 1.37;
    curdir = cd; % can work on external disk on any computer, assigning different drive letters
    hdfFile = [curdir(1),':\Matlab\Alex\HDF\0104_hompol.hdf'];
elseif nargin ~= 5
    error('Usage: acc_cos_pdf_HDF4_v2(hdfFile,firstFrame,lastFrame,stepFrames,threshold)')
end

% or
% hfile = [curdir(1),':\Matlab\Alex\HDF\0104_hompol.hdf'];
addpath([curdir(1),':\Matlab_stuff']);

%% Loop through files (comparisons)

%% reading from HDF4 file
[traj,attr] = readTrajHDF_v9(hdfFile,'u',[],'v',[],'w',[],'ax',[],'ay',[],'az',[],'w1',[],'w2',[],'w3',[],...
    's11',[],'s12',[],'s13',[],'s22',[],'s23',[],'s33',[],'t',[],'alx',[],'aly',[],...
    'alz',[],'minlength',20,'frames',[first last]);

%% Call (nested) functions
% parse_inputs
fields = fieldnames(traj);
for i = 1:length(fields)-2
    eval([fields{i},'=cat(1,traj.',fields{i},');']);
end
clear traj

%% Nested functions
%% Parsing inputs
%    function parse_inputs

%% Velocity

absu = (u.^2 + v.^2 + w.^2).^0.5;

%% $$u_i s_{ij}$$
usx = u.*s11 + v.*s12 + w.*s13;
usy = u.*s12 + v.*s22 + w.*s23;
usz = u.*s13 + v.*s23 + w.*s33;
absus = (usx.^2+usy.^2+usz.^2).^0.5;


absa = (ax.^2 + ay.^2 + az.^2).^0.5;
absal=(alx.^2 + aly.^2 + alz.^2).^0.5;
enstro = w1.^2+w2.^2+w3.^2;

%% Velocity gradient $$\partial u_i/\partial x_j$$
% ux  =   s11;
% uy  =   s12 - 0.5 * w3;
% uz  =   s13 + 0.5 * w2;
% vx  =   s12 + 0.5 * w3;
% vy  =   s22;
% vz  =   s23 - 0.5 * w1;
% wx  =   s13 - 0.5 * w2;
% wy  =   s23 + 0.5 * w1;
% wz  =   s33;

acx = u.*s11              + v.*(s12 - 0.5 * w3) + w.*(s13 + 0.5 * w2);
acy = u.*(s12 + 0.5 * w3) + v.*s22              + w.*(s23 - 0.5 * w1);
acz = u.*(s13 - 0.5 * w2) + v.*(s23 + 0.5 * w1) +   w.*s33;
absac = (acx.^2+acy.^2+acz.^2).^0.5;

ua = u.*acx + v.*acy + w.*acz;


aLamx = w2.*w - w3.*v;%
aLamy = w3.*u - w1.*w;%
aLamz = w1.*v - w2.*u;
% absaLam=(aLamx.^2+aLamy.^2+aLamz.^2).^0.5;

aBx = acx-aLamx;%
aBy = acy-aLamy;%
aBz = acz-aLamz;
absaB=(aBx.^2+aBy.^2+aBz.^2).^0.5;%%

unx = u./absu;
uny = v./absu;
unz = w./absu;
aun = ax.*unx + ay.*uny + az.*unz;

aPx=aun.*unx;%%
aPy=aun.*uny;
aPz=aun.*unz;%
absaP = (aPx.^2+aPy.^2+aPz.^2).^0.5;%

aOx = ax - aPx;
aOy = ay - aPy;%%
aOz = az - aPz;
absaO = (aOx.^2+aOy.^2+aOz.^2).^0.5;%%

cux=(u.*s11+v.*(s12 - 0.5 * w3)+w.*(s13 + 0.5 * w2))./(absu.^2);
cuy=(u.*(s12 + 0.5 * w3)+v.*s22+w.*(s23 - 0.5 * w1))./(absu.^2);
cuz=(u.*(s13 - 0.5 * w2)+v.*(s23 + 0.5 * w1)+w.*s33)./(absu.^2);
curvGrad=(cux.^2+cuy.^2+cuz.^2).^0.5;
curv=(((v.*az-w.*ay).^2+(w.*ax-u.*az).^2+(u.*ay-v.*ax).^2).^0.5)./(absu.^3);
curvac=(((v.*acz-w.*acy).^2+(w.*acx-u.*acz).^2+(u.*acy-v.*acx).^2).^0.5)./(absu.^3);
strain = s11.^2 + s22.^2 + s33.^2 + 2*(s12.^2 + s23.^2 + s13.^2);


% o11 = 0;
% o12 = - 0.5 * w3;
% o13 = 0.5 * w2;
% o21 = -o12;
% o22 = 0;
% o23 = - 0.5 * w1;
% o31 = - 0.5 * w2;
% o32 = 0.5 * w1;
% o33 = 0;

uox = -0.5*v.*w3 + 0.5*w.*w2;
uoy = 0.5*u.*w3  -0.5*w.*w1;
uoz = -0.5*u.*w2 + 0.5*v.*w1;
absuo = ((-0.5*v.*w3 + 0.5*w.*w2).^2 +...
    (0.5*u.*w3 - 0.5* w.*w1).^2 +...
    (-0.5*u.*w2 + 0.5*v.*w1).^2).^0.5;



%{
        daxdx=f(:,22);
        daxdy=f(:,23);
        daxdz=f(:,24);
        daydx=f(:,25);
        daydy=f(:,26);
        daydz=f(:,27);
        dazdx=f(:,28);
        dazdy=f(:,29);
        dazdz=f(:,30);
%}
%         reldiv=f(:,31);
%         age=f(:,32);
%    end

%% Main calculations

V = zeros(length(s11),3,3);
D = zeros(length(s11),3);

%% Heavy loop - calculating eigenvalues/eigenvectors
for ii = 1:length(s11)
    %         A = [...
    %             s11(ii) s12(ii) s13(ii);
    %             s12(ii) s22(ii) s23(ii);
    %             s13(ii) s23(ii) s33(ii)...
    %             ];
    S = [...
        s11(ii) s12(ii) s13(ii);
        s12(ii) s22(ii) s23(ii);
        s13(ii) s23(ii) s33(ii);...
        ];
    [VV,DD] = eig(S);
    [DD,k] = sort(diag(DD));        % ascending order
    D(ii,:) = DD(end:-1:1); 		% descending order
    V(ii,:,:) = VV(:,k(end:-1:1)); 	% the same order for eigenvectors
end

%% Prepare relevant indices: all, low, strain, enstro, etc.

%{
% example
negQ = Q < -1;
pos_ua = (u.*ax+v.*ay+w.*az)./(absa.*absu)>0.7;
neg_ua = (u.*ax+v.*ay+w.*az)./(absa.*absu)< -0.7;
posQweak = Q > 1.5 & Q < 5;
posQstrong = Q > 5;


statistic(:,1) = cosine([ax,ay,az],V(:,:,1));
statistic(:,2) = cosine([ax,ay,az],V(:,:,2));
statistic(:,3) = cosine([ax,ay,az],V(:,:,3));
statistic(:,4) = cosine([ax,ay,az],[alx,aly,alz]);
statistic(:,5) = cosine([ax,ay,az],ac);
statistic(:,6) = cosine([alx,aly,alz],ac);

%}

% indAll =  reldiv < 0.2 & absa.^2/8e-4 < 100 & curv < 3000;
indLow = enstro < 10^threshold & 2*strain < 10^threshold;
indEnstro = 2*strain<enstro & enstro>10^threshold;
indStrain =  2*strain>enstro & 2*strain>10^threshold;


cos_a_w  = cosine([ax, ay, az ],[w1,w2,w3]);
cos_ao_w = cosine([aOx,aOy,aOz],[w1,w2,w3]);
cos_ap_w = cosine([aPx,aPy,aPz],[w1,w2,w3]);
cos_ac_w = cosine([acx,acy,acz],[w1,w2,w3]);
cos_al_w = cosine([alx,aly,alz],[w1,w2,w3]);
cos_us_w = cosine([usx,usy,usz],[w1,w2,w3]);
cos_uo_w = cosine([uox,uoy,uoz],[w1,w2,w3]);

cos_a_l1  = cosine([ax, ay, az ],V(:,:,1));
cos_ao_l1 = cosine([aOx,aOy,aOz],V(:,:,1));
cos_ap_l1 = cosine([aPx,aPy,aPz],V(:,:,1));
cos_ac_l1 = cosine([acx,acy,acz],V(:,:,1));
cos_al_l1 = cosine([alx,aly,alz],V(:,:,1));
cos_us_l1 = cosine([usx,usy,usz],V(:,:,1));
cos_uo_l1 = cosine([uox,uoy,uoz],V(:,:,1));

cos_a_l2  = cosine([ax, ay, az ],V(:,:,2));
cos_ao_l2 = cosine([aOx,aOy,aOz],V(:,:,2));
cos_ap_l2 = cosine([aPx,aPy,aPz],V(:,:,2));
cos_ac_l2 = cosine([acx,acy,acz],V(:,:,2));
cos_al_l2 = cosine([alx,aly,alz],V(:,:,2));
cos_us_l2 = cosine([usx,usy,usz],V(:,:,2));
cos_uo_l2 = cosine([uox,uoy,uoz],V(:,:,2));

cos_a_l3  = cosine([ax, ay, az ],V(:,:,3));
cos_ao_l3 = cosine([aOx,aOy,aOz],V(:,:,3));
cos_ap_l3 = cosine([aPx,aPy,aPz],V(:,:,3));
cos_ac_l3 = cosine([acx,acy,acz],V(:,:,3));
cos_al_l3 = cosine([alx,aly,alz],V(:,:,3));
cos_us_l3 = cosine([usx,usy,usz],V(:,:,3));
cos_uo_l3 = cosine([uox,uoy,uoz],V(:,:,3));

cos_a_u  = cosine([ax, ay, az ],[u,v,w]);
cos_ao_u = cosine([aOx,aOy,aOz],[u,v,w]);
cos_ap_u = cosine([aPx,aPy,aPz],[u,v,w]);
cos_ac_u = cosine([acx,acy,acz],[u,v,w]);
cos_al_u = cosine([alx,aly,alz],[u,v,w]);
cos_us_u = cosine([usx,usy,usz],[u,v,w]);
cos_uo_u = cosine([uox,uoy,uoz],[u,v,w]);

cos_alac = cosine([alx,aly,alz],[acx,acy,acz]);
cos_usac = cosine([usx,usy,usz],[acx,acy,acz]);
cos_uoac = cosine([uox uoy uoz],[acx acy acz]);
cos_usuo = cosine([usx usy usz],[uox uoy uoz]);
cos_apac = cosine([aPx aPy aPz],[acx acy acz]);
%% Prepare graphics

%% typical figure construction: 4 conditions

styles = {'k','r','g','b'};
ind = {':','indLow','indEnstro','indStrain'};
quantity = {'absa','absaO','absaP','absal','absac','curv','curvac','absus','absuo'};
xlabels = {'|a|','|a_{\perp}|','|a_{||}|','|a_l|','|a_c|','1/r','1/r acc','|us|','|uo|'};

for k = 1:length(quantity)
    figure, hold on, grid on, box on
    for m = 1:length(ind)
        [n,x] = nhist(eval([quantity{k},'(',ind{m},')']),nBins);
        plot(x,n,styles{m});
        set(gca,'XScale','log','YScale','log')
    end
end



%% cosines figure construction:
figure;
hold on;
grid on;
box on;
[n,x]=nhist(cos_a_w,nBins);
plot(x,n,'k');
[n,x]=nhist(cos_ao_w,nBins);
plot(x,n,'r');
[n,x]=nhist(cos_ap_w,nBins);
plot(x,n,'g');
[n,x]=nhist(cos_al_w,nBins);
plot(x,n,'b');
[n,x]=nhist(cos_ac_w,nBins);
plot(x,n,'c');
[n,x]=nhist(cos_us_w,nBins);
plot(x,n,'m');
[n,x]=nhist(cos_uo_w,nBins);
plot(x,n,'y');
xlabel('cos(a,\omega) all','FontSize',12,'FontName','Times New Roman')
ylabel('pdf','FontSize',12,'FontName','Times New Roman')
hleg = legend('$\cos(a,\omega)$','$\cos(a_{\perp},\omega)$','$\cos(a_{||},\omega)$','$\cos(a_L,\omega)$',...
    '$\cos(a_c,\omega)$','$\cos(u_i s_{ij},\omega)$','$\cos(u_i r_{ij},\omega)$')
set(hleg,'interpreter','latex')
axis([-1 1 0 2])

figure;
hold on;
grid on;
box on;
[n,x]=nhist(cos_a_w(indLow),nBins);
plot(x,n,'k');
[n,x]=nhist(cos_ao_w(indLow),nBins);
plot(x,n,'r');
[n,x]=nhist(cos_ap_w(indLow),nBins);
plot(x,n,'g');
[n,x]=nhist(cos_al_w(indLow),nBins);
plot(x,n,'b');
[n,x]=nhist(cos_ac_w(indLow),nBins);
plot(x,n,'c');
[n,x]=nhist(cos_us_w(indLow),nBins);
plot(x,n,'m');
[n,x]=nhist(cos_uo_w(indLow),nBins);
plot(x,n,'y');
xlabel('cos(a,\omega) low','FontSize',12,'FontName','Times New Roman')
ylabel('pdf','FontSize',12,'FontName','Times New Roman')
axis([-1 1 0 2])

figure;
hold on;
grid on;
box on;
[n,x]=nhist(cos_a_w(indEnstro),nBins);
plot(x,n,'k');
[n,x]=nhist(cos_ao_w(indEnstro),nBins);
plot(x,n,'r');
[n,x]=nhist(cos_ap_w(indEnstro),nBins);
plot(x,n,'g');
[n,x]=nhist(cos_al_w(indEnstro),nBins);
plot(x,n,'b');
[n,x]=nhist(cos_ac_w(indEnstro),nBins);
plot(x,n,'c');
[n,x]=nhist(cos_us_w(indEnstro),nBins);
plot(x,n,'m');
[n,x]=nhist(cos_uo_w(indEnstro),nBins);
plot(x,n,'y');
xlabel('cos(a,\omega) enstro','FontSize',12,'FontName','Times New Roman')
ylabel('pdf','FontSize',12,'FontName','Times New Roman')
axis([-1 1 0 2])

figure;
hold on;
grid on;
box on;
[n,x]=nhist(cos_a_w(indStrain),nBins);
plot(x,n,'k');
[n,x]=nhist(cos_ao_w(indStrain),nBins);
plot(x,n,'r');
[n,x]=nhist(cos_ap_w(indStrain),nBins);
plot(x,n,'g');
[n,x]=nhist(cos_al_w(indStrain),nBins);
plot(x,n,'b');
[n,x]=nhist(cos_ac_w(indStrain),nBins);
plot(x,n,'c');
[n,x]=nhist(cos_us_w(indStrain),nBins);
plot(x,n,'m');
[n,x]=nhist(cos_uo_w(indStrain),nBins);
plot(x,n,'y');
xlabel('cos(a,\omega) strain','FontSize',12,'FontName','Times New Roman')
ylabel('pdf','FontSize',12,'FontName','Times New Roman')
axis([-1 1 0 2])

figure;
hold on;
grid on;
box on;
[n,x]=nhist(cos_a_u,nBins);
plot(x,n,'k');
[n,x]=nhist(cos_al_u,nBins);
plot(x,n,'b');
[n,x]=nhist(cos_ac_u,nBins);
plot(x,n,'c');
[n,x]=nhist(cos_us_u,nBins);
plot(x,n,'m');
[n,x]=nhist(cos_uo_u,nBins);
plot(x,n,'y');
xlabel('cos(u,a) all','FontSize',12,'FontName','Times New Roman')
ylabel('pdf','FontSize',12,'FontName','Times New Roman')
axis([-1 1 0 2])

figure;
hold on;
grid on;
box on;
[n,x]=nhist(cos_a_u(indLow),nBins);
plot(x,n,'k');
[n,x]=nhist(cos_al_u(indLow),nBins);
plot(x,n,'b');
[n,x]=nhist(cos_ac_u(indLow),nBins);
plot(x,n,'c');
[n,x]=nhist(cos_us_u(indLow),nBins);
plot(x,n,'m');
[n,x]=nhist(cos_uo_u(indLow),nBins);
plot(x,n,'y');
xlabel('cos(u,a) low','FontSize',12,'FontName','Times New Roman')
ylabel('pdf','FontSize',12,'FontName','Times New Roman')
axis([-1 1 0 2])

figure;
hold on;
grid on;
box on;
[n,x]=nhist(cos_a_u(indEnstro),nBins);
plot(x,n,'k');
[n,x]=nhist(cos_al_u(indEnstro),nBins);
plot(x,n,'b');
[n,x]=nhist(cos_ac_u(indEnstro),nBins);
plot(x,n,'c');
[n,x]=nhist(cos_us_u(indEnstro),nBins);
plot(x,n,'m');
[n,x]=nhist(cos_uo_u(indEnstro),nBins);
plot(x,n,'y');
xlabel('cos(u,a) enstro','FontSize',12,'FontName','Times New Roman')
ylabel('pdf','FontSize',12,'FontName','Times New Roman')
axis([-1 1 0 2])

figure;
hold on;
grid on;
box on;
[n,x]=nhist(cos_a_u(indStrain),nBins);
plot(x,n,'k');
[n,x]=nhist(cos_al_u(indStrain),nBins);
plot(x,n,'b');
[n,x]=nhist(cos_ac_u(indStrain),nBins);
plot(x,n,'c');
[n,x]=nhist(cos_us_u(indStrain),nBins);
plot(x,n,'m');
[n,x]=nhist(cos_uo_u(indStrain),nBins);
plot(x,n,'y');
xlabel('cos(u,a) strain','FontSize',12,'FontName','Times New Roman')
ylabel('pdf','FontSize',12,'FontName','Times New Roman')
axis([-1 1 0 2])

figure;
subplot(1,3,1);
hold on;
grid on;
box on;
[n,x]=nhist(cos_a_l1,nBins);
plot(x,n,'k');
[n,x]=nhist(cos_ao_l1,nBins);
plot(x,n,'r');
[n,x]=nhist(cos_ap_l1,nBins);
plot(x,n,'g');
[n,x]=nhist(cos_al_l1,nBins);
plot(x,n,'b');
[n,x]=nhist(cos_ac_l1,nBins);
plot(x,n,'c');
[n,x]=nhist(cos_us_l1,nBins);
plot(x,n,'m');
[n,x]=nhist(cos_uo_l1,nBins);
plot(x,n,'y');
xlabel('cos(a,\lambda_1) all','FontSize',12,'FontName','Times New Roman')
ylabel('pdf','FontSize',12,'FontName','Times New Roman')
axis([-1 1 0 2])

subplot(1,3,2);
hold on;
grid on;
box on;
[n,x]=nhist(cos_a_l2,nBins);
plot(x,n,'k');
[n,x]=nhist(cos_ao_l2,nBins);
plot(x,n,'r');
[n,x]=nhist(cos_ap_l2,nBins);
plot(x,n,'g');
[n,x]=nhist(cos_al_l2,nBins);
plot(x,n,'b');
[n,x]=nhist(cos_ac_l2,nBins);
plot(x,n,'c');
[n,x]=nhist(cos_us_l2,nBins);
plot(x,n,'m');
[n,x]=nhist(cos_uo_l2,nBins);
plot(x,n,'y');
xlabel('cos(a,\lambda_2) all','FontSize',12,'FontName','Times New Roman')
ylabel('pdf','FontSize',12,'FontName','Times New Roman')
axis([-1 1 0 2])

subplot(1,3,3);
hold on;
grid on;
box on;
[n,x]=nhist(cos_a_l3,nBins);
plot(x,n,'k');
[n,x]=nhist(cos_ao_l3,nBins);
plot(x,n,'r');
[n,x]=nhist(cos_ap_l3,nBins);
plot(x,n,'g');
[n,x]=nhist(cos_al_l3,nBins);
plot(x,n,'b');
[n,x]=nhist(cos_ac_l3,nBins);
plot(x,n,'c');
[n,x]=nhist(cos_us_l3,nBins);
plot(x,n,'m');
[n,x]=nhist(cos_uo_l3,nBins);
plot(x,n,'y');
xlabel('cos(a,\lambda_3) all','FontSize',12,'FontName','Times New Roman')
ylabel('pdf','FontSize',12,'FontName','Times New Roman')
axis([-1 1 0 2])



figure;
subplot(1,3,1);
hold on;
grid on;
box on;
[n,x]=nhist(cos_a_l1(indLow),nBins);
plot(x,n,'k');
[n,x]=nhist(cos_ao_l1(indLow),nBins);
plot(x,n,'r');
[n,x]=nhist(cos_ap_l1(indLow),nBins);
plot(x,n,'g');
[n,x]=nhist(cos_al_l1(indLow),nBins);
plot(x,n,'b');
[n,x]=nhist(cos_ac_l1(indLow),nBins);
plot(x,n,'c');
[n,x]=nhist(cos_us_l1(indLow),nBins);
plot(x,n,'m');
[n,x]=nhist(cos_uo_l1(indLow),nBins);
plot(x,n,'y');
xlabel('cos(a,\lambda_1) low','FontSize',12,'FontName','Times New Roman')
ylabel('pdf','FontSize',12,'FontName','Times New Roman')
axis([-1 1 0 2])

subplot(1,3,2);
hold on;
grid on;
box on;
[n,x]=nhist(cos_a_l2(indLow),nBins);
plot(x,n,'k');
[n,x]=nhist(cos_ao_l2(indLow),nBins);
plot(x,n,'r');
[n,x]=nhist(cos_ap_l2(indLow),nBins);
plot(x,n,'g');
[n,x]=nhist(cos_al_l2(indLow),nBins);
plot(x,n,'b');
[n,x]=nhist(cos_ac_l2(indLow),nBins);
plot(x,n,'c');
[n,x]=nhist(cos_us_l2(indLow),nBins);
plot(x,n,'m');
[n,x]=nhist(cos_uo_l2(indLow),nBins);
plot(x,n,'y');
xlabel('cos(a,\lambda_2) low','FontSize',12,'FontName','Times New Roman')
ylabel('pdf','FontSize',12,'FontName','Times New Roman')
axis([-1 1 0 2])

subplot(1,3,3);
hold on;
grid on;
box on;
[n,x]=nhist(cos_a_l3(indLow),nBins);
plot(x,n,'k');
[n,x]=nhist(cos_ao_l3(indLow),nBins);
plot(x,n,'r');
[n,x]=nhist(cos_ap_l3(indLow),nBins);
plot(x,n,'g');
[n,x]=nhist(cos_al_l3(indLow),nBins);
plot(x,n,'b');
[n,x]=nhist(cos_ac_l3(indLow),nBins);
plot(x,n,'c');
[n,x]=nhist(cos_us_l3(indLow),nBins);
plot(x,n,'m');
[n,x]=nhist(cos_uo_l3(indLow),nBins);
plot(x,n,'y');
xlabel('cos(a,\lambda_3) low','FontSize',12,'FontName','Times New Roman')
ylabel('pdf','FontSize',12,'FontName','Times New Roman')
axis([-1 1 0 2])

figure;
subplot(1,3,1);
hold on;
grid on;
box on;
[n,x]=nhist(cos_a_l1(indEnstro),nBins);
plot(x,n,'k');
[n,x]=nhist(cos_ao_l1(indEnstro),nBins);
plot(x,n,'r');
[n,x]=nhist(cos_ap_l1(indEnstro),nBins);
plot(x,n,'g');
[n,x]=nhist(cos_al_l1(indEnstro),nBins);
plot(x,n,'b');
[n,x]=nhist(cos_ac_l1(indEnstro),nBins);
plot(x,n,'c');
[n,x]=nhist(cos_us_l1(indEnstro),nBins);
plot(x,n,'m');
[n,x]=nhist(cos_uo_l1(indEnstro),nBins);
plot(x,n,'y');
xlabel('cos(a,\lambda_1) enstro','FontSize',12,'FontName','Times New Roman')
ylabel('pdf','FontSize',12,'FontName','Times New Roman')
axis([-1 1 0 2])

subplot(1,3,2);
hold on;
grid on;
box on;
[n,x]=nhist(cos_a_l2(indEnstro),nBins);
plot(x,n,'k');
[n,x]=nhist(cos_ao_l2(indEnstro),nBins);
plot(x,n,'r');
[n,x]=nhist(cos_ap_l2(indEnstro),nBins);
plot(x,n,'g');
[n,x]=nhist(cos_al_l2(indEnstro),nBins);
plot(x,n,'b');
[n,x]=nhist(cos_ac_l2(indEnstro),nBins);
plot(x,n,'c');
[n,x]=nhist(cos_us_l2(indEnstro),nBins);
plot(x,n,'m');
[n,x]=nhist(cos_uo_l2(indEnstro),nBins);
plot(x,n,'y');
xlabel('cos(a,\lambda_2) enstro','FontSize',12,'FontName','Times New Roman')
ylabel('pdf','FontSize',12,'FontName','Times New Roman')
axis([-1 1 0 2])

subplot(1,3,3);
hold on;
grid on;
box on;
[n,x]=nhist(cos_a_l3(indEnstro),nBins);
plot(x,n,'k');
[n,x]=nhist(cos_ao_l3(indEnstro),nBins);
plot(x,n,'r');
[n,x]=nhist(cos_ap_l3(indEnstro),nBins);
plot(x,n,'g');
[n,x]=nhist(cos_al_l3(indEnstro),nBins);
plot(x,n,'b');
[n,x]=nhist(cos_ac_l3(indEnstro),nBins);
plot(x,n,'c');
[n,x]=nhist(cos_us_l3(indEnstro),nBins);
plot(x,n,'m');
[n,x]=nhist(cos_uo_l3(indEnstro),nBins);
plot(x,n,'y');
xlabel('cos(a,\lambda_3) enstro','FontSize',12,'FontName','Times New Roman')
ylabel('pdf','FontSize',12,'FontName','Times New Roman')
axis([-1 1 0 2])

figure; % 21
subplot(1,3,1);
hold on;
grid on;
box on;
[n,x]=nhist(cos_a_l1(indStrain),nBins);
plot(x,n,'k');
[n,x]=nhist(cos_ao_l1(indStrain),nBins);
plot(x,n,'r');
[n,x]=nhist(cos_ap_l1(indStrain),nBins);
plot(x,n,'g');
[n,x]=nhist(cos_al_l1(indStrain),nBins);
plot(x,n,'b');
[n,x]=nhist(cos_ac_l1(indStrain),nBins);
plot(x,n,'c');
[n,x]=nhist(cos_us_l1(indStrain),nBins);
plot(x,n,'m');
[n,x]=nhist(cos_uo_l1(indStrain),nBins);
plot(x,n,'y');
xlabel('cos(a,\lambda_1) strain','FontSize',12,'FontName','Times New Roman')
ylabel('pdf','FontSize',12,'FontName','Times New Roman')
axis([-1 1 0 2])

subplot(1,3,2);
hold on;
grid on;
box on;
[n,x]=nhist(cos_a_l2(indStrain),nBins);
plot(x,n,'k');
[n,x]=nhist(cos_ao_l2(indStrain),nBins);
plot(x,n,'r');
[n,x]=nhist(cos_ap_l2(indStrain),nBins);
plot(x,n,'g');
[n,x]=nhist(cos_al_l2(indStrain),nBins);
plot(x,n,'b');
[n,x]=nhist(cos_ac_l2(indStrain),nBins);
plot(x,n,'c');
[n,x]=nhist(cos_us_l2(indStrain),nBins);
plot(x,n,'m');
[n,x]=nhist(cos_uo_l2(indStrain),nBins);
plot(x,n,'y');
xlabel('cos(a,\lambda_2) strain','FontSize',12,'FontName','Times New Roman')
ylabel('pdf','FontSize',12,'FontName','Times New Roman')
axis([-1 1 0 2])

subplot(1,3,3);
hold on;
grid on;
box on;
[n,x]=nhist(cos_a_l3(indStrain),nBins);
plot(x,n,'k');
[n,x]=nhist(cos_ao_l3(indStrain),nBins);
plot(x,n,'r');
[n,x]=nhist(cos_ap_l3(indStrain),nBins);
plot(x,n,'g');
[n,x]=nhist(cos_al_l3(indStrain),nBins);
plot(x,n,'b');
[n,x]=nhist(cos_ac_l3(indStrain),nBins);
plot(x,n,'c');
[n,x]=nhist(cos_us_l3(indStrain),nBins);
plot(x,n,'m');
[n,x]=nhist(cos_uo_l3(indStrain),nBins);
plot(x,n,'y');
xlabel('cos(a,\lambda_3) strain','FontSize',12,'FontName','Times New Roman')
ylabel('pdf','FontSize',12,'FontName','Times New Roman')
axis([-1 1 0 2])


figure; %22
hold on;
grid on;
box on;
[n,x]=nhist(cos_alac,nBins);
plot(x,n,'k');
[n,x]=nhist(cos_alac(indLow),nBins);
plot(x,n,'r');
[n,x]=nhist(cos_alac(indEnstro),nBins);
plot(x,n,'g');
[n,x]=nhist(cos_alac(indStrain),nBins);
plot(x,n,'b');
xlabel('cos(a_l,a_c) ','FontSize',12,'FontName','Times New Roman')
ylabel('pdf','FontSize',12,'FontName','Times New Roman')
axis([-1 1 0 10])

figure; % 23
hold on;
grid on;
box on;
[n,x]=nhist(cos_usac,nBins);
plot(x,n,'k');
[n,x]=nhist(cos_usac(indLow),nBins);
plot(x,n,'r');
[n,x]=nhist(cos_usac(indEnstro),nBins);
plot(x,n,'g');
[n,x]=nhist(cos_usac(indStrain),nBins);
plot(x,n,'b');
xlabel('cos(us,a_c) ','FontSize',12,'FontName','Times New Roman')
ylabel('pdf','FontSize',12,'FontName','Times New Roman')
axis([-1 1 0 8])

figure; % 24
hold on;
grid on;
box on;
[n,x]=nhist(cos_uoac,nBins);
plot(x,n,'k');
[n,x]=nhist(cos_uoac(indLow),nBins);
plot(x,n,'r');
[n,x]=nhist(cos_uoac(indEnstro),nBins);
plot(x,n,'g');
[n,x]=nhist(cos_uoac(indStrain),nBins);
plot(x,n,'b');
xlabel('cos(uo,a_c) ','FontSize',12,'FontName','Times New Roman')
ylabel('pdf','FontSize',12,'FontName','Times New Roman')
axis([-1 1 0 8])

figure; %25
hold on;
grid on;
box on;
[n,x]=nhist(cos_usuo,nBins);
plot(x,n,'k');
[n,x]=nhist(cos_usuo(indLow),nBins);
plot(x,n,'r');
[n,x]=nhist(cos_usuo(indEnstro),nBins);
plot(x,n,'g');
[n,x]=nhist(cos_usuo(indStrain),nBins);
plot(x,n,'b');
xlabel('cos(us,uo)')
ylabel('pdf')
axis([-1 1 0 2])

figure; % 26
hold on;
grid on;
box on;
[n,x]=nhist(cos_apac,nBins);
plot(x,n,'k');
[n,x]=nhist(cos_apac(indLow),nBins);
plot(x,n,'r');
[n,x]=nhist(cos_apac(indEnstro),nBins);
plot(x,n,'g');
[n,x]=nhist(cos_apac(indStrain),nBins);
plot(x,n,'b');
xlabel('cos(a_{||},a_c) ','FontSize',12,'FontName','Times New Roman')
ylabel('pdf','FontSize',12,'FontName','Times New Roman')
axis([-1 1 0 2])

%% Save figures
[fpath,fname] = fileparts(hdfFile);
save_open_figures(fname);

function save_open_figures(baseName,fmt)
if nargin == 0
    baseName = 'fig';
    fmt = 'fig';
elseif nargin == 1
    fmt = 'fig';
end

hf = sort(findobj(0,'type','figure'));
for i = 1:length(hf)
    saveas(hf(i),[baseName,'_',int2str(i)],fmt);
    close(hf(i));
end




';
end

hf = sort(findobj(0,'type','figure'));
for i = 1:length(hf)
    saveas(hf(i),[baseName,'_',int2str(i)],fmt);
    close(hf(i));
end




