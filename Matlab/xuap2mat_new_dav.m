clear all

% Initialization

directoryName = 'res_scene32'; %update the file name we run

dx =-0.05:0.005:0.03;
dz =-0.04:0.005:0.04;
dy=-0.045:0.005:0.045;

% Read the xuap files.
tmp = readXUAPFiles(directoryName);
% tmp = building_trajectories(tmp);
[XI,YI,ZI] = meshgrid(dx,dy,dz);
for i=1:length(tmp) 
    ind = (find(abs(tmp(i).uf) > 10e-5 & abs(tmp(i).vf) > 10e-5 & abs(tmp(i).wf) > 10e-5));
    if length(ind) > 10
        tm(i).u= griddata3(tmp(i).xf(ind),tmp(i).yf(ind),tmp(i).zf(ind),tmp(i).uf(ind),XI,YI,ZI);
        tm(i).v= griddata3(tmp(i).xf(ind),tmp(i).yf(ind),tmp(i).zf(ind),tmp(i).vf(ind),XI,YI,ZI);
        tm(i).w= griddata3(tmp(i).xf(ind),tmp(i).yf(ind),tmp(i).zf(ind),tmp(i).wf(ind),XI,YI,ZI);
    end
end