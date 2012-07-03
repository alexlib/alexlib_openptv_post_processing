clear all 

% Initialization

directoryName = 'test'; %update the file name we run 

dx = [-0.05:0.005:0.05];
dz = [-0.05:0.005:0.05];

% Read the xuap files.
tmp = readXUAPFiles(directoryName);
% tmp = building_trajectories(tmp);
[XI,YI,ZI] = meshgrid(dx,dx,dz);
[ui_sum,vi_sum,wi_sum] = deal(zeros(size(XI)));
count = 0;
for i=1:length(tmp)
ind = (find(abs(tmp(i).uf) > 10e-5 & abs(tmp(i).vf) > 10e-5 & abs(tmp(i).wf) > 10e-5));
if length(ind) > 10 
ui = griddata3(tmp(i).xf(ind),tmp(i).yf(ind),tmp(i).zf(ind),tmp(i).uf(ind),XI,YI,ZI); 
%or should i run on it through the frame? ui(i), ui=u interpolation 

vi = griddata3(tmp(i).xf(ind),tmp(i).yf(ind),tmp(i).zf(ind),tmp(i).vf(ind),XI,YI,ZI);
wi = griddata3(tmp(i).xf(ind),tmp(i).yf(ind),tmp(i).zf(ind),tmp(i).wf(ind),XI,YI,ZI);
count = count + 1;
%for the average calculation 
ui_sum = ui_sum + ui; 
vi_sum = vi_sum + vi;
wi_sum = wi_sum + wi;
end 
end
% Average:
% ui_sum = ui_sum/count;
% vi_sum = vi_sum/count;
% wi_sum = wi_sum/count;
% % Plot in 2D relates to the average values
% quiver3(XI,YI,ZI,ui_sum,vi_sum,wi_sum)
% xlabel('x'),ylabel('y'),zlabel('z') 
% view(2)
% Plot in 3D relates to the average values
% figure
% quiver3(XI,YI,ZI,ui_sum,vi_sum,wi_sum)
% xlabel('x'),ylabel('y'),zlabel('z') 
% view(3)