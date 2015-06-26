% function a_ratios_cond_diva_ua_alex2(first,last)


%% Prepare necessary quantities:
first = 1;
last = 2990; %2990;
tail = 3;
nBins = 100;
edges = logspace(-3,3,nBins); % for binning of enstrophy and 2*strain
%% One time calculations for later graphics
xedges = 0.5*(edges(1:(nBins-1))+edges(2:nBins)); % xTickLabels
dA = diff(edges)'*diff(edges); % normalization of JPDF


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

    absal=(cat(1,traj.alx).^2 + cat(1,traj.aly).^2 + cat(1,traj.alz).^2).^0.5;%%%%

    w1 = cat(1,traj.w1);
    w2 = cat(1,traj.w2);
    w3 = cat(1,traj.w3);

    clear traj

    acx = u.*ux+v.*uy+w.*uz;
    acy = u.*vx+v.*vy+w.*vz;
    acz = u.*wx+v.*wy+w.*wz;

    absac = (acx.^2+acy.^2+acz.^2).^0.5;

    %     aLamx=w2.*w-w3.*v;%
    %     aLamy=w3.*u-w1.*w;%
    %     aLamz=w1.*v-w2.*u;
    %     absaLam=(aLamx.^2+aLamy.^2+aLamz.^2).^0.5;

    %     aBx=acx-aLamx;%
    %     aBy=acy-aLamy;%
    %     aBz=acz-aLamz;
    %     absaB=(aBx.^2+aBy.^2+aBz.^2).^0.5;%%

    %     unx=u./absu;
    %     uny=v./absu;
    %     unz=w./absu;
    %     aun=ax.*unx+ay.*uny+az.*unz;

    %     aPx=aun.*unx;%%
    %     aPy=aun.*uny;
    %     aPz=aun.*unz;%
    %     absaP=(aPx.^2+aPy.^2+aPz.^2).^0.5;%

    %     aOx=ax-aPx;
    %     aOy=ay-aPy;
    %     aOz=az-aPz;
    %     absaO=(aOx.^2+aOy.^2+aOz.^2).^0.5;


    strain = ux.^2 + vy.^2 + wz.^2 + 2*(0.5*(uy + vx).^2 + 0.5*(uz + wx).^2 + 0.5*(vz + wy).^2);

    enstrophy = sum([w1, w2, w3].^2,2);




    Q = 0.25*(enstrophy-2*strain);

    negQ = Q < -1;
    pos_ua = (u.*ax+v.*ay+w.*az)./(absa.*absu)>0.7;
    neg_ua = (u.*ax+v.*ay+w.*az)./(absa.*absu)< -0.7;
    posQweak = Q > 1.5 & Q < 5;
    posQstrong = Q > 5;

    statistic(:,1) = absa./absal;
    statistic(:,2) = absa./absac;
    statistic(:,3) = absal./absac;


    [N1,X1] = nhist(statistic(:,1),1000); %%%cos(a,\lambda_1)
    [N2,X2] = nhist(statistic(:,2),1000); %%%cos(a,\lambda_2)
    [N3,X3] = nhist(statistic(:,3),1000); %%%cos(a,\lambda_3)
    createfigure1(X1,N1,X2,N2,X3,N3);

    % negQ and pos_div_a
    ind = negQ & pos_ua;
    [N1,X1] = nhist(statistic(ind,1),1000); %%%cos(a,\lambda_1)
    [N2,X2] = nhist(statistic(ind,2),1000); %%%cos(a,\lambda_2)
    [N3,X3] = nhist(statistic(ind,3),1000); %%%cos(a,\lambda_3)
    createfigure1(X1,N1,X2,N2,X3,N3);
    title('$Q < -1\ \& \ u \cdot a > 0.7$','interpreter','latex')

    % negQ and neg_div_a
    ind = negQ & neg_ua;
    [N1,X1] = nhist(statistic(ind,1),1000); %%%cos(a,\lambda_1)
    [N2,X2] = nhist(statistic(ind,2),1000); %%%cos(a,\lambda_2)
    [N3,X3] = nhist(statistic(ind,3),1000); %%%cos(a,\lambda_3)
    createfigure1(X1,N1,X2,N2,X3,N3);
    title('$Q < -1\ \& \ u \cdot a < 0.7$','interpreter','latex')

    % negQ and pos_div_a
    ind = posQweak;
    [N1,X1] = nhist(statistic(ind,1),1000); %%%cos(a,\lambda_1)
    [N2,X2] = nhist(statistic(ind,2),1000); %%%cos(a,\lambda_2)
    [N3,X3] = nhist(statistic(ind,3),1000); %%%cos(a,\lambda_3)
    createfigure1(X1,N1,X2,N2,X3,N3);
    title('$Q > 1.5 \ \& \ Q < 5$, weak','interpreter','latex')


    % negQ and pos_div_a
    ind = posQstrong;
    [N1,X1] = nhist(statistic(ind,1),1000); %%%cos(a,\lambda_1)
    [N2,X2] = nhist(statistic(ind,2),1000); %%%cos(a,\lambda_2)
    [N3,X3] = nhist(statistic(ind,3),1000); %%%cos(a,\lambda_3)
    createfigure1(X1,N1,X2,N2,X3,N3);
    title('$Q > 5 \ $, strong','interpreter','latex')
    
    

% end % hfiles

