function [newtraj] = test_rbf_interp(traj,smoothness)
% should work for scattered data
x = cat(1,traj.xf);
y = cat(1,traj.yf);
z = cat(1,traj.zf);
t = cat(1,traj.t);

x = x(:)';
y = y(:)';
z = z(:)';
t = t(:)';

ti = linspace(t(1),t(end),length(t)*2);


figure
plot(t,x,'o',t,y,'s',t,z,'^');
hold on

if nargin < 2
    smoothness = 0.1;
end


coeffx = rbfcreate(t,x,'RBFsmooth',smoothness);
coeffy = rbfcreate(t,y,'RBFsmooth',smoothness);
coeffz = rbfcreate(t,z,'RBFsmooth',smoothness);



xi = rbfinterp(ti,coeffx);
yi = rbfinterp(ti,coeffy);
zi = rbfinterp(ti,coeffz);

% plot(ti,xi,'r--',ti,yi,'g--',ti,zi,'b--')
plot(ti,xi,'r--',ti,yi,'g--',ti,zi,'b--')

newtraj.t = ti;
newtraj.xf = xi;
newtraj.yf = yi;
newtraj.zf = zi;

newtraj.uf = gradient(newtraj.xf,newtraj.t/50);
newtraj.vf = gradient(newtraj.yf,newtraj.t/50);
newtraj.wf = gradient(newtraj.zf,newtraj.t/50);

newtraj.axf = gradient(newtraj.uf,newtraj.t/50);
newtraj.ayf = gradient(newtraj.vf,newtraj.t/50);
newtraj.azf = gradient(newtraj.wf,newtraj.t/50);

newtraj.trajid = ones(length(xi),1);

