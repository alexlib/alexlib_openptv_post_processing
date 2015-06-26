% Analysis of lalb HDF data
% Reads HDF data from [pairsFile,'_L.hdf'] and prepares figures
%

% Created: 02-May-2006
% Author: Alex Liberzon
% E-Mail : liberzon@ifu.baug.ethz.ch
% Phone : +41 (0)44 633 3754
% Copyright (c) 2006 IFU, ETH Zurich
%
%
% $Revision: 1.0 $  $Date: 03-May-2006 11:20:21$


pairsFile = input('HDF file name (e.g. mai10) ','s')
% pairsFile = 'April26'
N = 1;
dt = 1/60; tau = [0:.1:10];
tau_eta = 0.23; % sec
agebins = tau*tau_eta/dt;

colors = {'m','r','g','b','k'};
linestyles = {':','-','--','-.'};
markers = {'v','','o','s','x','+','<','^','>','d','p','*'};

for i = 1:length(agebins)-1
    style{i} = [colors{mod(i,length(colors))+1},linestyles{mod(i,length(linestyles))+1},markers{mod(i,length(markers))+1}];
end


% Load data:
data = hdfread([pairsFile,'_L.hdf'],pairsFile,'Fields','zetaa,zetab,age','firstrecord',1);
zetaa = data{1};
zetab = data{2};
age = data{3};

ind = bindex(age,agebins,1);

% Stretching rates, zetaA, zetaB
hf = figure,
[tmpa,tmpb] = deal(zeros(size(agebins)));
for i = 1:length(agebins)
%    tmp(i) = mean(mean(dlnLdt(1:N,1,(ind == i)),1),3);
        tmpa(i) = mean(zetaa(1:N,ind == i),2);
        tmpb(i) = mean(zetab(1:N,ind == i),2);
        %tmp2(i) = mean(mean(zeta(2,(ind == i)),1),2);
        %tmp3(i) = mean(mean(zeta(3,(ind == i)),1),2);
end
plot(agebins(1:end-1)/tau_eta*dt,tmpa(1:end-1))
hold on
plot(agebins(1:end-1)/tau_eta*dt,tmpb(1:end-1))
ylabel('dln(l)/dt','fontname','euclid','fontangle','italic')
xlabel('\tau_{\eta}')
% saveas(hf,['mean_zeta_tau_',namestr],'fig')


data = hdfread([pairsFile,'_L.hdf'],pairsFile,'Fields','llsa,llsb','firstrecord',1);
llsa = data{1};
llsb = data{2};


% Stretching, llsA, llsB
hf = figure,
[tmpa,tmpb] = deal(zeros(size(agebins)));
for i = 1:length(agebins)
%    tmp(i) = mean(mean(dlnLdt(1:N,1,(ind == i)),1),3);
        tmpa(i) = mean(llsa(1:N,ind == i),2);
        tmpb(i) = mean(llsb(1:N,ind == i),2);
        %tmp2(i) = mean(mean(zeta(2,(ind == i)),1),2);
        %tmp3(i) = mean(mean(zeta(3,(ind == i)),1),2);
end
plot(agebins(1:end-1)/tau_eta*dt,tmpa(1:end-1))
hold on
plot(agebins(1:end-1)/tau_eta*dt,tmpb(1:end-1))
ylabel('l_il_js_{ij}','fontname','euclid','fontangle','italic')
xlabel('\tau_{\eta}')


% Mean 
k = 1;
figure, hold on
for i = [1 2 4]
    tmp1 = llsa(1,ind == i);
    tmp2 = llsb(1,ind == i);
    if ~isempty(tmp1) & ~isempty(tmp2)
        [n,x] = hist3([tmp2(:),tmp1(1,:)'],{-5:.5:20,-5:.5:20});
        for j = 2:N-1
            [n1,x] = hist3([tmp2(:),tmp1(j,:)'],x);
            n = n + n1;
        end

        subplot(3,2,k); k = k + 2;
        contour(x{1},x{2},n',30)
        axis tight
        xlabel('$l^A_i l^A_j s_{ij}$','interpreter','latex')
        ylabel(sprintf('$l^B_i l^B_j s_{ij}, \\tau = %d \\tau_{\\eta}$',tau(i)),'interpreter','latex')
        grid on
        axis square
        colorbar vert
        %         saveas(gcf,sprintf('jointpdf_lls_wws_rand_tau_%d_0104',tau(i)),'fig')
    end
end


% Test of the la,lb consistency.
k = 1;
figure, hold on
for i = [1 2 4]
    tmp1 = llsa(1,ind == i);
    tmp2 = llsb(1,ind == i);
    if ~isempty(tmp1) & ~isempty(tmp2)
        [n,x] = hist3([tmp2(:),tmp1(1,:)'],{-5:.5:20,-5:.5:20});
        for j = 2:N-1
            [n1,x] = hist3([tmp2(:),tmp1(j,:)'],x);
            n = n + n1;
        end

        subplot(3,2,k); k = k + 2;
        contour(x{1},x{2},n',30)
        axis tight
        xlabel('$l^A_i l^A_j s_{ij}$','interpreter','latex')
        ylabel(sprintf('$l^B_i l^B_j s_{ij}, \\tau = %d \\tau_{\\eta}$',tau(i)),'interpreter','latex')
        grid on
        axis square
        colorbar vert
        %         saveas(gcf,sprintf('jointpdf_lls_wws_rand_tau_%d_0104',tau(i)),'fig')
    end
end




omega = [dwdy - dvdz, dudz - dwdx, dvdx - dudy];

W = repmat(0,[numPoints 3]);
for j = 1:numPoints
    W(j,1:3) = omega(j,1:3)*...
        [0.5*(dudx(j) + dudx(j)), 0.5*(dudy(j) + dvdx(j)), 0.5*(dudz(j) + dwdx(j));
        0.5*(dudy(j) + dvdx(j)),0.5*(dvdy(j) + dvdy(j)),0.5*(dvdz(j) + dwdy(j));
        0.5*(dudz(j) + dwdx(j)),0.5*(dvdz(j) + dwdy(j)),0.5*(dwdz(j) + dwdz(j))]';
end

wws = repmat(0,[numPoints 1]);
for i = 1:numPoints
    wws(i) = omega(i,1:3)*W(i,1:3)';
end
% wws_w2 = wws./sum(omega.^2,2);

data = hdfread([hfile(1:end-4),'_L.hdf'],namestr,'Fields','lls,age','firstrecord',1);
lls = data{1};
age = data{2};



save 0104_wws_lls wws lls age

dt = 1/60; tau = [1:5];
tau_eta = 0.23; % sec
agebins = tau*tau_eta/dt;
ind = bindex(age,agebins,1);

colors = {'m','r','g','b','k'};
linestyles = {':','-','--','-.'};
markers = {'v','','o','s','x','+','<','^','>','d','p','*'};

for i = 1:length(agebins)-1
    style{i} = [colors{mod(i,length(colors))+1},linestyles{mod(i,length(linestyles))+1},markers{mod(i,length(markers))+1}];
end

k = 1; figure, hold on
for i = [1,2,4]
    tmp1 = lls(1:N-1,ind == i);
    tmp2 = wws(ind == i);
    if ~isempty(tmp1) & ~isempty(tmp2)
        [n,x] = hist3([tmp2(:),tmp1(1,:)'],{-30:1:30,-5:.5:20});
        for j = 2:N-1
            [n1,x] = hist3([tmp2(:),tmp1(j,:)'],x);
            n = n + n1;
        end

        subplot(3,2,k); k = k + 2;
        contour(x{1},x{2},n',30)
        axis tight
        xlabel('$\wws$','interpreter','latex')
        ylabel(sprintf('$l_i l_j s_{ij}, \\tau = %d \\tau_{\\eta}$',tau(i)),'interpreter','latex')
        grid on
        axis square
        colorbar vert
        %         saveas(gcf,sprintf('jointpdf_lls_wws_rand_tau_%d_0104',tau(i)),'fig')
    end
end

k = 2;
for i = [1,2,4]
    tmp1 = lls(1:N,ind == i);
    tmp2 = wws(ind == i);
    if ~isempty(tmp1) & ~isempty(tmp2)
        [n,x] = hist3([tmp2(:),tmp1(1,:)'],{-30:1:30,-5:.5:20});
        %         for j = 2:N-1
        %             [n1,x] = hist3([tmp2(:),tmp1(j,:)'],x);
        %             n = n + n1;
        %         end
        subplot(3,2,k); k = k + 2;
        contour(x{1},x{2},n',30)
        axis tight
        xlabel('$\wws$','interpreter','latex')
        ylabel(sprintf('$l_i l_j s_{ij}, \\tau = %d \\tau_{\\eta}$',tau(i)),'interpreter','latex')
        grid on
        axis square
        colorbar vert

    end
end
saveas(gcf,sprintf('jointpdf_lls_wws_rand_spec_tau_0104'),'fig')

