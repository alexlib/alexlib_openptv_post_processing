function track_2D_Sobel
% track_2D_Sobel - uses TXT files from MC1324 Sobel camera, 
% use subfunction txt2poslist to load the data in the 'trackable' format
% and uses 'track.m' from http://physics.georgetown.edu/matlab to track the
% particles.

pos_list = txt2poslist;
param.mem = 1;
param.good = 5;
param.dim = 2;
param.quiet = 0;
tr = track(pos_list, 10, param);
traj = struct('x',[],'y',[],'u',[],'v',[],'ax',[],'ay',[],'id',[]);
% figure, hold on
for i = 1:length(tr)
    ind = find(tr(:,4) == i); 
    if ~isempty(ind)
    % plot(tr(ind,1),tr(ind,2));
    traj(i).x = tr(ind,1); 
    traj(i).y = 1024 - tr(ind,2);
    traj(i).t = tr(ind,3);
    traj(i).u = diff([traj(i).x;traj(i).x(length(ind))]);
    traj(i).v = diff([traj(i).y;traj(i).y(length(ind))]);
    traj(i).ax = diff([traj(i).u;traj(i).u(length(ind))]);
    traj(i).ay = diff([traj(i).v;traj(i).v(length(ind))]);
    traj(i).id = tr(ind(1),4);
    end
    % drawnow 
end

figure, quiver(cat(1,traj.x),cat(1,traj.y),cat(1,traj.u),cat(1,traj.v))
axis([00 1100 50 1000])
% axis tight
set(gca,'ydir','reverse')

save traj.mat

function pos_list = txt2poslist
% TXT2POSLIST converts TXT files from Sobel-filtered camera MC1324 into
% a single position list file for track.m (see
% http://physics.georgetown.edu/matlab/tutorial.html for more details)
% run this command in Matlab in the same directory where the TXT files are.
% Usage: 
% >> cd /myTXTresults/
% >> txt2rtis
% 

% Author: Alex Liberzon
% Date: 26.May.2008

d = dir('*.txt');
% fid = fopen('pos_list.txt','w');
tmp = load(d(1).name);
pos_list = cat(2,tmp,repmat(0,length(tmp),1)); % zeros(100000,3);

for i = 2:length(d)
    tmp = load(d(i).name);
    
    pos_list = cat(1,pos_list,cat(2,tmp,repmat(i-1,length(tmp),1)));
%     for ii = 1:size(tmp,1)
%         fprintf(fid,'%6.4f %6.4f %d\n',  tmp(ii,1), tmp(ii,2), i-1);
%     end  
end
% fclose(fid);

    