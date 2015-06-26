function test
%% ANALYSIS_15-Nov-2007.M
% for acceleration of water/polymers
% Alex, 15.11.07    

%% Prepare necessary quantities:
    nBins = 100;
    edges = logspace(-2,3,nBins); % for binning of enstrophy and 2*strain
    xedges = 0.5*(edges(1:(nBins-1))+edges(2:nBins));

%% Assign the working data
curdir = cd;
hfiles{1} = [curdir(1),':\Matlab\Alex\HDF\0104_water.hdf'];
% or
hfiles{2} = [curdir(1),':\Matlab\Alex\HDF\0104_hompol.hdf'];

%% Main loop
for hfile = hfiles

    %% reading from HDF4 file
    [traj,attr] = readTrajHDF_v9(hfile,'u',[],'v',[],'w',[],'ax',[],'ay',[],'az',[],'w1',[],'w2',[],'w3',[],...
        's11',[],'s12',[],'s13',[],'s22',[],'s23',[],'s33',[],'t',[],'alx',[],'aly',[],...
        'alz',[],'minlength',20);

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
    
    acx=u.*ux+v.*uy+w.*uz;
    acy=u.*vx+v.*vy+w.*vz;
    acz=u.*wx+v.*wy+w.*wz;

    absac=(acx.^2+acy.^2+acz.^2).^0.5;

    aLamx=w2.*w-w3.*v;%
    aLamy=w3.*u-w1.*w;%
    aLamz=w1.*v-w2.*u;
    absaLam=(aLamx.^2+aLamy.^2+aLamz.^2).^0.5;

    aBx=acx-aLamx;%
    aBy=acy-aLamy;%
    aBz=acz-aLamz;
    absaB=(aBx.^2+aBy.^2+aBz.^2).^0.5;%%

    unx=u./absu;
    uny=v./absu;
    unz=w./absu;
    aun=ax.*unx+ay.*uny+az.*unz;

    aPx=aun.*unx;%%
    aPy=aun.*uny;
    aPz=aun.*unz;%
    absaP=(aPx.^2+aPy.^2+aPz.^2).^0.5;%

    aOx=ax-aPx;
    aOy=ay-aPy;%%
    aOz=az-aPz;
    absaO=(aOx.^2+aOy.^2+aOz.^2).^0.5;%%

    cux=(u.*ux+v.*uy+w.*uz)./(absu.^2);
    cuy=(u.*vx+v.*vy+w.*vz)./(absu.^2);
    cuz=(u.*wx+v.*wy+w.*wz)./(absu.^2);
    curvGrad=(cux.^2+cuy.^2+cuz.^2).^0.5;
    curv=(((v.*az-w.*ay).^2+(w.*ax-u.*az).^2+(u.*ay-v.*ax).^2).^0.5)./(absu.^3);

    strain = ux.^2 + vy.^2 + wz.^2 + 2*(0.5*(uy + vx).^2 + 0.5*(uz + wx).^2 + 0.5*(vz + wy).^2);

    enstrophy = sum([w1, w2, w3].^2,2);

    %% Joint PDFs
    % 1) enstrophy vs strain
    mhist2d = hist2d([enstrophy(:),strain(:)],edges,edges);
    figure,
    contourf(log10(xedges),log10(xedges),log10(mhist2d))
    xlabel('$\omega^2$','interpreter','latex')
    ylabel('$2s^2$','interpreter','latex')
    colorbar('vert')

    % 2) conditional JPDF of absa
    [mhist2d,cmhist2d] = condhist2d([enstrophy(:),strain(:),absa(:).^2],edges,edges);
    figure,
    contourf(log10(xedges),log10(xedges),log10(cmhist2d))
    xlabel('$\omega^2$','interpreter','latex')
    ylabel('$2s^2$','interpreter','latex')
    title('$|a|$','interpreter','latex')
    colorbar('vert')

    % 3) conditional JPDF of absaP
    [mhist2d,cmhist2d] = condhist2d([enstrophy(:),strain(:),absaP(:).^2],edges,edges);
    figure,
    contourf(log10(xedges),log10(xedges),log10(cmhist2d))
    xlabel('$\omega^2$','interpreter','latex')
    ylabel('$2s^2$','interpreter','latex')
    title('$|a_p|$','interpreter','latex')
    colorbar('vert')

    % 4) conditional JPDF of absaO
    [mhist2d,cmhist2d] = condhist2d([enstrophy(:),strain(:),absaO(:).^2],edges,edges);
    figure,
    contourf(log10(xedges),log10(xedges),log10(cmhist2d))
    xlabel('$\omega^2$','interpreter','latex')
    ylabel('$2s^2$','interpreter','latex')
    title('$|a_O|$','interpreter','latex')
    colorbar('vert')

    % 5) conditional JPDF of absa
    [mhist2d,cmhist2d] = condhist2d([enstrophy(:),strain(:),absaO(:).^2./(absaP(:).^2 + absaO(:).^2)],edges,edges);
    figure,
    contourf(log10(xedges),log10(xedges),log10(cmhist2d))
    xlabel('$\omega^2$','interpreter','latex')
    ylabel('$2s^2$','interpreter','latex')
    title('$|a_O/(a_P + a_O)|$','interpreter','latex')
    colorbar('vert')

    % 6) conditional JPDF of curvature
    [mhist2d,cmhist2d] = condhist2d([enstrophy(:),strain(:),curv(:)],edges,edges);
    figure,
    contourf(log10(xedges),log10(xedges),log10(cmhist2d))
    xlabel('$\omega^2$','interpreter','latex')
    ylabel('$2s^2$','interpreter','latex')
    title('$r$','interpreter','latex')
    colorbar('vert')

    % 7) conditional JPDF of curvature
    [mhist2d,cmhist2d] = condhist2d([enstrophy(:),strain(:),1./curv(:)],edges,edges);
    figure,
    contourf(log10(xedges),log10(xedges),log10(cmhist2d))
    xlabel('$\omega^2$','interpreter','latex')
    ylabel('$2s^2$','interpreter','latex')
    title('$1/r$','interpreter','latex')
    colorbar('vert')
end

%% end %%
