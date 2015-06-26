curdir = cd;
hfile = [curdir(1),':\Matlab\Alex\HDF\0104_water.hdf'];
% or 
 % hfile = [curdir(1),':\Matlab\Alex\HDF\0104_hompol.hdf'];

    
[traj,attr] = readTrajHDF_v9(hfile,'u',[],'v',[],'w',[],'ax',[],'ay',[],'az',[],'w1',[],'w2',[],'w3',[],...
's11',[],'s12',[],'s13',[],'s22',[],'s23',[],'s33',[],'t',[],'minlength',20);


figure, hold on
nhist(cat(1,traj.ax)./std(cat(1,traj.ax)),1000,'r')
nhist(cat(1,traj.ay)./std(cat(1,traj.ay)),1000,'g')
nhist(cat(1,traj.az)./std(cat(1,traj.az)),1000,'b')

dudx  =   cat(1,traj.s11);
dudy  =   cat(1,traj.s12) - 0.5 * cat(1,traj.w3);
dudz  =   cat(1,traj.s13) + 0.5 * cat(1,traj.w2);
dvdx  =   cat(1,traj.s12) + 0.5 * cat(1,traj.w3);
dvdy  =   cat(1,traj.s22);
dvdz  =   cat(1,traj.s23) - 0.5 * cat(1,traj.w1);
dwdx  =   cat(1,traj.s13) - 0.5 * cat(1,traj.w2);
dwdy  =   cat(1,traj.s23) + 0.5 * cat(1,traj.w1);
dwdz  =   cat(1,traj.s33);

sSq = dudx.^2 + dvdy.^2 + dwdz.^2 + 2*(0.5*(dudy + dvdx).^2 + 0.5*(dudz + dwdx).^2 + 0.5*(dvdz + dwdy).^2);

aSq = cat(1,traj.ax).^2 + cat(1,traj.ay).^2 + cat(1,traj.az).^2;

omega = [dwdy - dvdz, dudz - dwdx, dvdx - dudy];

enstrophy = sum(omega.^2,2);

urms = (cat(1,traj.u).^2 + cat(1,traj.v).^2 + cat(1,traj.w).^2).^0.5;

% -------------------------------------------------------------------------
% plot()

% Joint PDF  of s^2 versus a^2
figure;
% [n,x,bins] = histmulti5([sSq(:),aSq(:)],[(0:.5:70)',(0:1:140)']);
[n,x,bins] = histmulti5([aSq(:),sSq(:)],[linspace(0,0.012,100)',linspace(0,200,100)']);
contour(x(:,1),x(:,2),n')
axis tight
ylabel('2s^2')
xlabel('a^2')
grid on
axis square


inda = bindex(aSq(:),x(:,1),1); 
% inds = bindex(1.74*sSq(:),x(:,2),1); 
inds = bindex(sSq(:),x(:,2),1); 

tmp1 = zeros(length(x(:,1)),1);
tmp2 =  zeros(length(x(:,2)),1);


hfline = figure;
for i = 1:length(x(:,1))
    tmp1(i) = mean(sSq(inda == i));
    %     tmp1(i) = mean(1.74*sSq(inda == i));
    tmp2(i) = mean(aSq(inds == i));
end

hold on
plot(x(:,1),tmp1,'b',tmp2,x(:,2),'r');




% Joint PDF  of s^2 versus w^2
figure;
% [n,x,bins] = histmulti5([sSq(:),aSq(:)],[(0:.5:70)',(0:1:140)']);
[n,x,bins] = histmulti5([aSq(:),enstrophy(:)],[linspace(0,0.012,100)',linspace(0,200,100)']);
contour(x(:,1),x(:,2),n')
axis tight
ylabel('w^2')
xlabel('a^2')
grid on
axis square


figure
% inda = bindex(aSq(:),x(:,1),1); 
indw = bindex(enstrophy(:),x(:,2),1); 

tmp1 = zeros(length(x(:,1)),1);
tmp2 =  zeros(length(x(:,2)),1);

for i = 1:length(x(:,1))
    tmp1(i) = mean(enstrophy(inda == i));
    tmp2(i) = mean(aSq(indw == i));
end

figure(hfline)
hold on
plot(x(:,1),tmp1,'b--',tmp2,x(:,2),'g');


figure(hfline)
hl = findobj(gca,'type','line');
figure, hold on,
for i = 1:length(hl),
xd = get(hl(i),'xdata');
yd = get(hl(i),'ydata');
 plot(xd(1:end-1)/mean(aSq(:)), yd(1:end-1)/mean(enstrophy(:)));
end


%% aSq vs u_rms
% Joint PDF  of s^2 versus a^2
figure;
% [n,x,bins] = histmulti5([sSq(:),aSq(:)],[(0:.5:70)',(0:1:140)']);
[n,x,bins] = histmulti5([aSq(:),urms(:)],[linspace(0,0.012,100)',linspace(0,0.03,100)']);
contour(x(:,1),x(:,2),n')
axis tight
ylabel('u_{rms}')
xlabel('a^2')
grid on
axis square


inda = bindex(aSq(:),x(:,1),1); 
% inds = bindex(1.74*sSq(:),x(:,2),1); 
inds = bindex(urms(:),x(:,2),1); 

tmp1 = zeros(length(x(:,1)),1);
tmp2 =  zeros(length(x(:,2)),1);


hfline = figure;
for i = 1:length(x(:,1))
    tmp1(i) = mean(urms(inda == i));
    %     tmp1(i) = mean(1.74*sSq(inda == i));
    tmp2(i) = mean(aSq(inds == i));
end

hold on
ind1 = ~isnan(tmp1);
ind2 = ~isnan(tmp2);
plot(x(ind1,1),tmp1(ind1),'b',tmp2(ind2),x(ind2,2),'r');
xlabel('u_{rms}')
ylabel('a^2')

%% Save figures
[fpath,fname] = fileparts(hfile);
saveAllFigures(fname)



return

% --------------------
% cos(u,a)
cosua = cosine([cat(1,traj.u),cat(1,traj.v),cat(1,traj.w)], [cat(1,traj.ax),cat(1,traj.ay),cat(1,traj.az)]);

% Joint PDF  of s^2 versus w^2
% hf = figure
% % [n,x,bins] = histmulti5([sSq(:),aSq(:)],[(0:.5:70)',(0:1:140)']);
% [n,x,bins] = histmulti5([aSq(:),cosua(:)],[linspace(0,0.025,200)',linspace(-1,1,200)']);
% contour(x(:,1),x(:,2),n')
% axis tight
% ylabel('cosua')
% xlabel('a^2')
% grid on
% axis square

indpos = find(cosua > .7);
indneg = find(cosua < -.7);


% Joint PDF  of s^2 versus a^2
hf = figure
% [n,x,bins] = histmulti5([sSq(:),aSq(:)],[(0:.5:70)',(0:1:140)']);
[n,x,bins] = histmulti5([aSq(indpos),1.74*sSq(indpos)],[linspace(0,0.025,200)',linspace(0,400,200)']);
contour(x(:,1),x(:,2),n')
axis tight
ylabel('2s^2')
xlabel('a^2')
grid on
axis square


tmpa = aSq(indpos);
tmps = 1.74*sSq(indpos);

inda = bindex(tmpa,x(:,1),1); 
inds = bindex(tmps,x(:,2),1); 

tmp1 = zeros(length(x(:,1)),1);
tmp2 =  zeros(length(x(:,2)),1);


hfline = figure;
for i = 1:length(x(:,1))
    tmp1(i) = mean(tmps(inda == i));
    tmp2(i) = mean(tmpa(inds == i));
end

hfline = figure
plot(x(:,1)/mean(tmpa),tmp1/mean(tmps),'r-',tmp2/mean(tmpa),x(:,2)/mean(tmps),'r--');

% -------------------------
% negative
% Joint PDF  of s^2 versus a^2
hf = figure
% [n,x,bins] = histmulti5([sSq(:),aSq(:)],[(0:.5:70)',(0:1:140)']);
[n,x,bins] = histmulti5([aSq(indneg),1.74*sSq(indneg)],[linspace(0,0.025,200)',linspace(0,400,200)']);
contour(x(:,1),x(:,2),n')
axis tight
ylabel('2s^2')
xlabel('a^2')
grid on
axis square


tmpa = aSq(indneg);
tmps = 1.74*sSq(indneg);

inda = bindex(tmpa,x(:,1),1); 
inds = bindex(tmps,x(:,2),1); 

tmp1 = zeros(length(x(:,1)),1);
tmp2 =  zeros(length(x(:,2)),1);


for i = 1:length(x(:,1))
    tmp1(i) = mean(tmps(inda == i));
    tmp2(i) = mean(tmpa(inds == i));
end

figure(hfline)
hold on
plot(x(:,1)/mean(tmpa),tmp1/mean(tmps),'b-',tmp2/mean(tmpa),x(:,2)/mean(tmps),'b--');

% -----------------
% enstrophy, conditioned on cos(u,a)


% Joint PDF  of w^2 versus a^2
hf = figure
% [n,x,bins] = histmulti5([sSq(:),aSq(:)],[(0:.5:70)',(0:1:140)']);
[n,x,bins] = histmulti5([aSq(indpos),enstrophy(indpos)],[linspace(0,0.025,200)',linspace(0,400,200)']);
contour(x(:,1),x(:,2),n')
axis tight
ylabel('2s^2')
xlabel('a^2')
grid on
axis square


tmpa = aSq(indpos);
tmps = enstrophy(indpos);

inda = bindex(tmpa,x(:,1),1); 
inds = bindex(tmps,x(:,2),1); 

tmp1 = zeros(length(x(:,1)),1);
tmp2 =  zeros(length(x(:,2)),1);


for i = 1:length(x(:,1))
    tmp1(i) = mean(tmps(inda == i));
    tmp2(i) = mean(tmpa(inds == i));
end


hfline = figure
hold on
plot(x(:,1)/mean(tmpa),tmp1/mean(tmps),'r-',tmp2/mean(tmpa),x(:,2)/mean(tmps),'r--');

% -------------------------
% negative
% Joint PDF  of w^2 versus a^2
hf = figure
% [n,x,bins] = histmulti5([sSq(:),aSq(:)],[(0:.5:70)',(0:1:140)']);
[n,x,bins] = histmulti5([aSq(indneg),enstrophy(indneg)],[linspace(0,0.025,200)',linspace(0,400,200)']);
contour(x(:,1),x(:,2),n')
axis tight
ylabel('2s^2')
xlabel('a^2')
grid on
axis square


tmpa = aSq(indneg);
tmps = enstrophy(indneg);

inda = bindex(tmpa,x(:,1),1); 
inds = bindex(tmps,x(:,2),1); 

tmp1 = zeros(length(x(:,1)),1);
tmp2 =  zeros(length(x(:,2)),1);



for i = 1:length(x(:,1))
    tmp1(i) = mean(tmps(inda == i));
    tmp2(i) = mean(tmpa(inds == i));
end

figure(hfline)
hold on
plot(x(:,1)/mean(tmpa),tmp1/mean(tmps),'b-',tmp2/mean(tmpa),x(:,2)/mean(tmps),'b--');
