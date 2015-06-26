clear all

% Initialization

directoryName = 'res_scene32'; %update the file name we run

dx =-0.05:0.005:0.05;
dz =-0.05:0.005:0.05;
% Creating the particle counting matrix it must be the same size as data

% Read the xuap files.
tmp = readXUAPFiles(directoryName);
% tm = struct('par_count',zeros(21,21,21),'u',zeros(21,21,21),'v',zeros(21,21,21),...
%     'w',zeros(21,21,21));
tm = struct('par_count',zeros(21,21,21));
tm = repmat(tm,[1,length(tmp)]);
% tmp = building_trajectories(tmp);
[XI,YI,ZI] = meshgrid(dx,dx,dz);
for i=1:length(tmp)
    ind = (find(abs(tmp(i).uf) > 10e-5 & abs(tmp(i).vf) > 10e-5 & abs(tmp(i).wf) > 10e-5));
    if length(ind) > 10
        for a=1:21
            for b=1:21
                for c=1:21
                    ra=-0.05+0.005*(a-1);
                    ra1=-0.05+0.005*a;
                    rb=-0.05+0.005*(b-1);
                    rb1=-0.05+0.005*b;
                    rc=-0.05+0.005*(c-1);
                    rc1=-0.05+0.005*c;
                    tmp2 = tmp(i).xf(ind)>=ra & tmp(i).xf(ind)<ra1 &  tmp(i).yf(ind)>=rb & tmp(i).yf(ind)<rb1 & tmp(i).zf(ind)>=rc & tmp(i).zf(ind)<rc1;
                    tm(i).count(a,b,c)=tm(i).par_count(a,b,c)+sum(tmp2);
                end
            end
        end
    end
    %      tm(i).u= griddata3(tmp(i).xf(ind),tmp(i).yf(ind),tmp(i).zf(ind),tmp(i).uf(ind),XI,YI,ZI);
    %      tm(i).v= griddata3(tmp(i).xf(ind),tmp(i).yf(ind),tmp(i).zf(ind),tmp(i).vf(ind),XI,YI,ZI);
    %      tm(i).w= griddata3(tmp(i).xf(ind),tmp(i).yf(ind),tmp(i).zf(ind),tmp(i).wf(ind),XI,YI,ZI);
end

