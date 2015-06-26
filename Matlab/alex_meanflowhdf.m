
% first = 1
% last = 2048
% loadResult=0
% hfile = 'F:\SPTV\Sourcecode\ptv_RP\exp_6_18\res\exp_6_16.hdf'
[hfile,hpath] = uigetfile('*.hdf','Pick HDF file');
[traj,attr]=readtrajhdf(fullfile(hpath,hfile),'u',[],'v',[],'w',[],'x',[],'y',[],'z',[]);




u = cat(1,traj.u);
v = cat(1,traj.v);
w = cat(1,traj.w);
x=cat(1,traj.x);
y=cat(1,traj.y);
z=cat(1,traj.z);

binX = linspace(min(x), max(x), 30);
binY = linspace(min(y), max(y), 25);
binZ = linspace(min(z), max(z), 20);

[gridY,gridX,gridZ] = meshgrid(binY,binX,binZ);

indX = bindex(x,binX,0);
indY = bindex(y,binY,0);
indZ = bindex(z,binZ,0);

% nhist(sqrt(u.^2+v.^2+w.^2),100)

[U,V,W] = deal(zeros(size(gridX)));

for i = 1:length(binX)
    tmpjk = (indX == i);
    for j = 1:length(binY)
        tmpk = ((indY == j) & tmpjk);
        for k = 1:length(binZ)
            tmp = (tmpk & (indZ == k));
            if any(tmp)
%                 [i,j,k]
                U(i,j,k) = mean(u(tmp));
                V(i,j,k) = mean(v(tmp));
                W(i,j,k) = mean(w(tmp));
            end
        end
    end
end

% U(isnan(U)) = 0;
% V(isnan(V)) = 0;
% W(isnan(W)) = 0;

figure
mesh(squeeze(mean(U,1)))
zlabel('U, [m/s]')

figure
mesh(squeeze(mean(U,2)))
zlabel('U, [m/s]')

figure
mesh(squeeze(mean(U,3)))
zlabel('U, [m/s]')


figure
mesh(squeeze(mean(V,1)))
zlabel('V, [m/s]')

figure
mesh(squeeze(mean(V,2)))
zlabel('V, [m/s]')

figure
mesh(squeeze(mean(V,3)))
zlabel('V, [m/s]')

figure
mesh(squeeze(mean(W,1)))
zlabel('W, [m/s]')

figure
mesh(squeeze(mean(W,2)))
zlabel('W, [m/s]')

figure
mesh(squeeze(mean(W,3)))
zlabel('W, [m/s]')

feval('save',[hpath,filesep,hfile,'.mat']);
