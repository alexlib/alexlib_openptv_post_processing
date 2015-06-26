% hfile = '0809_water.hdf'
% hfile = '100ppm.hdf'
% hfile = '50ppm.hdf'
 hfile = '0104_water.hdf';
% hfile = '0104_hompol.hdf';

directory = '/Users/alex/Documents/Projects/Polymers/';

hdffile = fullfile(directory,[hfile(1:end-4),'_lpqr.hdf']);

tau_eta = 0.2;



% New piece: plot of the values from the HDF files:
data = hdfread(hdffile,hfile(1:end-4),'Fields','coslwls,coslwlu,zeta,age','firstrecord',1);
coslWls = data{1};
coslWlu = data{2};
zeta = data{3};
age = data{4};

dt = 1/60;
agebins = [0:6]*tau_eta/dt
ind = bindex(age,agebins,1);

colors = {'m','r','g','b','k'};
linestyles = {':','-','--','-.'};
markers = {'v','','o','s','x','+','<','^','>','d','p','*'};
for i = 1:length(agebins)
style{i} = [colors{mod(i,length(colors))+1},linestyles{mod(i,length(linestyles))+1},markers{mod(i,length(markers))+1}];
end

hf = figure,
hold on
[n,x] = deal(zeros(20,length(agebins)));
for i = 1:length(agebins)
tmp = coslWls(:,(ind == i));
agebins(i),length(tmp)
if ~isempty(tmp)
nhist(tmp(:),21,style{i});
end
end
hold off
legend('\tau_{\eta} \leq 1','1 \leq \tau_{\eta} \leq 2','2 \leq \tau_{\eta} \leq 4','4 \leq \tau_{\eta} \leq 6','6 \leq \tau_{\eta}')
xlabel('cos(l,W^{l,s})')
% saveas(hf,['pdf_coslWls_',hfile(1:end-4)],'fig')

hf = figure,
hold on
[n,x] = deal(zeros(20,length(agebins)));
for i = 1:length(agebins)
tmp = coslWlu(:,(ind == i));
agebins(i),length(tmp)
if ~isempty(tmp)
    nhist(tmp(:),21,style{i});
end
end
hold off
legend('\tau_{\eta} \leq 1','1 \leq \tau_{\eta} \leq 2','2 \leq \tau_{\eta} \leq 4','4 \leq \tau_{\eta} \leq 6','6 \leq \tau_{\eta}')
xlabel('cos(l,W^{l,u})')
% saveas(hf,['pdf_coslWlu_',hfile(1:end-4)],'fig')

agebins = 0:2:10*tau_eta/dt;
ind = bindex(age,agebins,1);

hf = figure,
tmp = zeros(size(agebins));
for i = 1:length(agebins)
    tmp(i) = mean(mean(zeta(:,(ind == i)),1),2);
end
plot(agebins(:)/tau_eta*dt,tmp(:))
ylabel('dln(l)/dt','fontname','euclid','fontangle','italic')
xlabel('\tau_{\eta}')
% saveas(hf,['mean_zeta_tau_',hfile(1:end-4)],'fig')