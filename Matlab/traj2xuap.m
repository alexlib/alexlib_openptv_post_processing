function xuap = traj2xuap(traj)
% XUAP = TRAJ2XUAP(TRAJ) converts the structure TRAJ to the structure
% XUAP which is a frame by frame structure of X, U, A
% last column is ID


first_frame = traj(1).t(1);
last_frame = traj(end).t(end);

x = cat(1,traj.xf);
y = cat(1,traj.yf);
z = cat(1,traj.zf);
u = cat(1,traj.uf);
v = cat(1,traj.vf);
w = cat(1,traj.wf);
ax = cat(1,traj.axf);
ay = cat(1,traj.ayf);
az = cat(1,traj.azf);
id = cat(1,traj.trajid);
t = cat(1,traj.t);

xuap = struct('x',[],'y',[],'z',[],...
    'u',[],'v',[],'w',[],...
    'ax',[],'ay',[],'az',[],'id',[]);
k = 0;

for i = first_frame:last_frame
    ind = t == i;
    k = k + 1;
    xuap(k).x = x(ind);
    xuap(k).y = y(ind);
    xuap(k).z = z(ind);
    xuap(k).u = u(ind);
    xuap(k).v = v(ind);
    xuap(k).w = w(ind);
    xuap(k).ax = ax(ind);
    xuap(k).ay = ay(ind);
    xuap(k).az = az(ind);
    xuap(k).id = id(ind);
end


    
    


