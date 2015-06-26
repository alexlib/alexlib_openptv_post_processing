% function a_ratios_cond_diva_ua_alex2(first,last)

clear all

%% Prepare necessary quantities:
first = 1;
last = 2990;


%% Assign the working data
curdir = cd;
% hfile = [curdir(1),':\Matlab\Alex\HDF\0104_water.hdf'];
% or
hfile = [curdir(1),':\Matlab\Alex\HDF\0104_hompol.hdf'];

%% Main loop%
% for hfile = hfiles

clear statistics

%% reading from HDF4 file
traj = readTrajHDF_v9(hfile,'u',[],'v',[],'w',[],'ax',[],'ay',[],'az',[],'w1',[],'w2',[],'w3',[],...
    's11',[],'s12',[],'s13',[],'s22',[],'s23',[],'s33',[],'t',[],'alx',[],'aly',[],...
    'alz',[],'minlength',20,'frames',[first last]);

ux  =   cat(1,traj.s11);
uy  =   cat(1,traj.s12) - 0.5 * cat(1,traj.w3);
uz  =   cat(1,traj.s13) + 0.5 * cat(1,traj.w2);
vx  =   cat(1,traj.s12) + 0.5 * cat(1,traj.w3);
vy  =   cat(1,traj.s22);
vz  =   cat(1,traj.s23) - 0.5 * cat(1,traj.w1);
wx  =   cat(1,traj.s13) - 0.5 * cat(1,traj.w2);
wy  =   cat(1,traj.s23) + 0.5 * cat(1,traj.w1);
wz  =   cat(1,traj.s33);

u = cat(1,traj.u);
v = cat(1,traj.v);
w = cat(1,traj.w);
absu=(u.^2+v.^2+w.^2).^0.5;

ax = cat(1,traj.ax); ay = cat(1,traj.ay); az = cat(1,traj.az);

absa=(ax.^2 + ay.^2 + az.^2).^0.5;

% absal=(cat(1,traj.alx).^2 + cat(1,traj.aly).^2 + cat(1,traj.alz).^2).^0.5;%%%%

w1 = cat(1,traj.w1);
w2 = cat(1,traj.w2);
w3 = cat(1,traj.w3);

alx = cat(1,traj.alx);
aly = cat(1,traj.aly);
alz = cat(1,traj.alz);

clear traj

ac = [...
    u.*ux+v.*uy+w.*uz,...
    u.*vx+v.*vy+w.*vz,...
    u.*wx+v.*wy+w.*wz...
    ];


% absac = (acx.^2+acy.^2+acz.^2).^0.5;

strain = ux.^2 + vy.^2 + wz.^2 + 2*(0.5*(uy + vx).^2 + 0.5*(uz + wx).^2 + 0.5*(vz + wy).^2);

enstrophy = sum([w1, w2, w3].^2,2);




Q = 0.25*(enstrophy-2*strain);

D = zeros(length(Q),3);
V = zeros(length(Q),3,3);

for ii = 1:length(Q)
    %         A = [...
    %             s11(ii) s12(ii) s13(ii);
    %             s12(ii) s22(ii) s23(ii);
    %             s13(ii) s23(ii) s33(ii)...
    %             ];
    A = [...
        ux(ii) uy(ii) uz(ii);
        vx(ii) vy(ii) vz(ii);
        wx(ii) wy(ii) wz(ii);...
        ];
    [VV,DD] = eig(0.5*(A + A'));
    [DD,k] = sort(diag(DD)); 	% ascending order
    D(ii,:) = DD(end:-1:1); 			% descending order
    V(ii,:,:) = VV(:,k(end:-1:1)); 	% the same order for eigenvectors
end


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


%% Figures
%% a
% cosines a,\lambda_i - all
[N1,X1] = nhist(statistic(:,1),1000); %%%cos(a,\lambda_1)
[N2,X2] = nhist(statistic(:,2),1000); %%%cos(a,\lambda_2)
[N3,X3] = nhist(statistic(:,3),1000); %%%cos(a,\lambda_3)
hf1 = createfigurecosines(X1,N1,X2,N2,X3,N3);

% cosines a,al, a,ac, al,ac - all
[N1,X1] = nhist(statistic(:,4),1000); %%%cos(a,\lambda_1)
[N2,X2] = nhist(statistic(:,5),1000); %%%cos(a,\lambda_2)
[N3,X3] = nhist(statistic(:,6),1000); %%%cos(a,\lambda_3)
hf2 = createfigurecosines(X1,N1,X2,N2,X3,N3);

% combine:
figs2subplots([hf1 hf2],[2 1])
close(hf1,hf2)


% conditions;

% negQ and pos_div_a
ind = negQ & pos_ua;
[N1,X1] = nhist(statistic(ind,1),1000); %%%cos(a,\lambda_1)
[N2,X2] = nhist(statistic(ind,2),1000); %%%cos(a,\lambda_2)
[N3,X3] = nhist(statistic(ind,3),1000); %%%cos(a,\lambda_3)
hf1 = createfigurecosines(X1,N1,X2,N2,X3,N3);

% cosines a,al, a,ac, al,ac - all
[N1,X1] = nhist(statistic(ind,4),1000); %%%cos(a,\lambda_1)
[N2,X2] = nhist(statistic(ind,5),1000); %%%cos(a,\lambda_2)
[N3,X3] = nhist(statistic(ind,6),1000); %%%cos(a,\lambda_3)
hf2 = createfigurecosines(X1,N1,X2,N2,X3,N3);

% combine:
figs2subplots([hf1 hf2],[2 1])
title('$Q < -1\ \& \ u \cdot a > 0.7$','interpreter','latex')
close(hf1,hf2)

% negQ and neg_div_a
ind = negQ & neg_ua;
[N1,X1] = nhist(statistic(ind,1),1000); %%%cos(a,\lambda_1)
[N2,X2] = nhist(statistic(ind,2),1000); %%%cos(a,\lambda_2)
[N3,X3] = nhist(statistic(ind,3),1000); %%%cos(a,\lambda_3)
hf1 = createfigurecosines(X1,N1,X2,N2,X3,N3);

% cosines a,al, a,ac, al,ac - all
[N1,X1] = nhist(statistic(ind,4),1000); %%%cos(a,\lambda_1)
[N2,X2] = nhist(statistic(ind,5),1000); %%%cos(a,\lambda_2)
[N3,X3] = nhist(statistic(ind,6),1000); %%%cos(a,\lambda_3)
hf2 = createfigurecosines(X1,N1,X2,N2,X3,N3);

% combine:
figs2subplots([hf1 hf2],[2 1])
title('$Q < -1\ \& \ u \cdot a < 0.7$','interpreter','latex')
close(hf1,hf2)

% negQ and pos_div_a
ind = posQweak;
[N1,X1] = nhist(statistic(ind,1),1000); %%%cos(a,\lambda_1)
[N2,X2] = nhist(statistic(ind,2),1000); %%%cos(a,\lambda_2)
[N3,X3] = nhist(statistic(ind,3),1000); %%%cos(a,\lambda_3)
hf1 = createfigurecosines(X1,N1,X2,N2,X3,N3);

% cosines a,al, a,ac, al,ac - all
[N1,X1] = nhist(statistic(ind,4),1000); %%%cos(a,\lambda_1)
[N2,X2] = nhist(statistic(ind,5),1000); %%%cos(a,\lambda_2)
[N3,X3] = nhist(statistic(ind,6),1000); %%%cos(a,\lambda_3)
hf2 = createfigurecosines(X1,N1,X2,N2,X3,N3);

% combine:
figs2subplots([hf1 hf2],[2 1])
title('$Q > 1.5 \ \& \ Q < 5$, weak','interpreter','latex')
close(hf1,hf2)

% negQ and pos_div_a
ind = posQstrong;
[N1,X1] = nhist(statistic(ind,1),1000); %%%cos(a,\lambda_1)
[N2,X2] = nhist(statistic(ind,2),1000); %%%cos(a,\lambda_2)
[N3,X3] = nhist(statistic(ind,3),1000); %%%cos(a,\lambda_3)
hf1 = createfigurecosines(X1,N1,X2,N2,X3,N3);

% cosines a,al, a,ac, al,ac - all
[N1,X1] = nhist(statistic(ind,4),1000); %%%cos(a,\lambda_1)
[N2,X2] = nhist(statistic(ind,5),1000); %%%cos(a,\lambda_2)
[N3,X3] = nhist(statistic(ind,6),1000); %%%cos(a,\lambda_3)
hf2 = createfigurecosines(X1,N1,X2,N2,X3,N3);

% combine:
figs2subplots([hf1 hf2],[2 1])
title('$Q > 5 \ $, strong','interpreter','latex')
close(hf1,hf2)


