hfile = 'a.hdf'

[traj,attr] = readTrajAccHDF(hfile,...
    'u',[ ],'v',[ ],'w',[ ],...
    'dudx',[ ],'dudy',[ ],'dudz',[ ],...
    'dvdx',[ ],'dvdy',[ ],'dvdz',[ ],...
    'dwdx',[ ],'dwdy',[ ],'dwdz',[ ],...
    'ax',[ ],'ay',[ ],'az',[ ]); %,...
    %'daxdx',[ ],'daxdy',[ ],'daxdz',[ ],...
    %'daydx',[ ],'daydy',[ ],'daydz',[ ],...
    %'dazdx',[ ],'dazdy',[ ],'dazdz',[ ]);

numPoints = length(cat(1,traj.dudx))
numTraj = length(traj)

% Velocity vector
vel = cat(2,cat(1,traj.u),cat(1,traj.v),cat(1,traj.w));


dudx = cat(1,traj.dudx);
dudy = cat(1,traj.dudy);
dudz = cat(1,traj.dudz);

dvdx = cat(1,traj.dvdx);
dvdy = cat(1,traj.dvdy);
dvdz = cat(1,traj.dvdz);

dwdx = cat(1,traj.dwdx);
dwdy = cat(1,traj.dwdy);
dwdz = cat(1,traj.dwdz);

reldiv = abs(dudx + dvdy + dwdz)./(abs(dudx)+abs(dvdy)+abs(dwdz));

inddiv = reldiv < 0.1;


 % 1
 % Divergence check with derivatives
 % Joint PDF's
 figure
 [n,x,bins] = histmulti5([-squeeze(dudx(:,1)),squeeze(dudx(:,5)+dudx(:,9))],500);
 ax(1) = subplot(131);
 contour(x(:,1),x(:,2),n')
 axis tight
 xlabel('-\partialu/\partialx')
 ylabel('\partialv/\partialy+\partialw/\partialz')
 grid on
 axis square
 set(gca,'CLim',[0 1])
 % saveas(gcf,'dudx','fig');
 
 [n,x,bins] = histmulti5([-squeeze(dudx(:,5)),squeeze(dudx(:,1)+dudx(:,9))],500);
 ax(2) = subplot(132);
 contour(x(:,1),x(:,2),n')
 axis tight
 xlabel('-\partialv/\partialy')
 ylabel('\partialu/\partialx+\partialw/\partialz')
 grid on
 axis square
 set(gca,'CLim',[0 1])
 % saveas(gcf,'dvdy','fig');
 
 [n,x,bins] = histmulti5([-squeeze(dudx(:,9)),squeeze(dudx(:,1)+dudx(:,5))],500);
 ax(3) = subplot(133);
 contour(x(:,1),x(:,2),n')
 axis tight
 xlabel('-\partialw/\partialz')
 ylabel('\partialu/\partialx+\partialv/\partialy')
 grid on
 axis square
 set(gca,'CLim',[0 1])
 pos = get(ax(3), 'Position');
 
 h = colorbar('horiz');
 dp = get(ax(3),'position') - pos;
 
 for i=1:2
 set(ax(i),'position', get(ax(i), 'Position') + dp) ;
 end
 set(h,'position',[.11 .31 .81 .02])
 saveas(gcf,'divcheck','fig');


% 2
% Acceleration check (previously al_TrajecPointAAALAC, al_trajPointAAALAC)
% Joint PDF's
figure
[n,x] = histmulti5([cat(1,traj.ax),al(:,1) + ac(:,1)],1000);
ax(1) = subplot(131);
contour(x(:,1),x(:,2),n')
axis tight
xlabel('a_x [m^2/s]')
ylabel('a_{l,x}+a_{c,x}, [m^2/s]')
grid on
axis square
% set(gca,'CLim',[0 1])
% saveas(gcf,'dudx','fig');




[n,x,bins] = histmulti5([cat(1,traj.ay),al(:,2) + ac(:,2)],1000);
ax(2) = subplot(132);
contour(x(:,1),x(:,2),n')
axis tight
xlabel('a_y [m^2/s]')
ylabel('a_{l,y}+a_{c,y}, [m^2/s]')
grid on
axis square
% set(gca,'CLim',[0 1])


[n,x,bins] = histmulti5([cat(1,traj.az),al(:,3) + ac(:,3)],1000);
ax(3) = subplot(133);
contour(x(:,1),x(:,2),n')
axis tight
xlabel('a_z [m^2/s]')
ylabel('a_{l,z}+a_{c,z}, [m^2/s]')
grid on
axis square
% set(gca,'CLim',[0 1])


% One colorbar for subplots
pos = get(ax(3), 'Position');
h = colorbar('horiz');
dp = get(ax(3),'position') - pos;

for i=1:2
    set(ax(i),'position', get(ax(i), 'Position') + dp) ;
end
set(h,'position',[.11 .31 .81 .02])

saveas(gcf,'acceleration_check','fig');



% Additional stuff

D = repmat(NaN,[1 3 numPoints]);
V = repmat(NaN,[3 3 numPoints]);
    for j = 1:numPoints         % for all time points
    s = [0.5*(dudx(j) + dudx(j)), 0.5*(dudy(j) + dvdx(j)), 0.5*(dudz(j) + dwdx(j));
    0.5*(dudy(j) + dvdx(j)),0.5*(dvdy(j) + dvdy(j)),0.5*(dvdz(j) + dwdy(j));
    0.5*(dudz(j) + dwdx(j)),0.5*(dvdz(j) + dwdy(j)),0.5*(dwdz(j) + dwdz(j))];
        [v,d] = eig(s);
        [d,k] = sort(diag(d)); 	% ascending order
        %     m = m + 1;
        D(:,:,j) = d(end:-1:1); 	% descending order		
        V(:,:,j) = v(:,k(end:-1:1)); 	% the same order for eigenvectors 
    end



mean(squeeze(D(:,1,:)))
mean(squeeze(D(:,2,:)))
mean(squeeze(D(:,3,:)))


% histogram of eigenvalues
[n1,x1] = nhist(D(:,1,:),200);
[n2,x2] = nhist(D(:,2,:),200);
[n3,x3] = nhist(D(:,3,:),200);
hf = figure
semilogy(x1,n1,'r-',x2,n2,'r-s',x3,n3,'r-x','linewidth',1,'markersize',4);
saveas(hf,['pdf_lambda_',hfile(1:end-4)],'fig')% Plot histogram of the eigenvalues





%% Convective acceleration
%ac = cat(2,...
%    dot(vel,dudx(:,1:3),2), ... % acx = u.*ux + v.*uy + w.*uz;
%    dot(vel,dudx(:,4:6),2), ... % acy = u.*vx + v.*vy + w.*vz;
%    dot(vel,dudx(:,7:9),2));    % acz = u.*wx + v.*wy + w.*wz;

% absac = (acx.^2+acy.^2+acz.^2).^0.5;  


%% cos(al,ac)
%cos_al_ac = dot(al,ac,2)./sqrt(sum(al.^2,2))./sqrt(sum(ac.^2,2));
%
%% PDF
%figure;
%[n,x] = nhist(cos_al_ac,100);
%plot(x,n);
%xlabel('cos(a_l,a_c)','FontName','Euclid','FontAngle','italic');
%ylabel('Probability density');
%saveas(gcf,['pdf_accel_',hfile(1:end-4)],'fig')