% hfile = '100ppm.hdf'
% hfile = '0809_water.hdf'
% hfile = '50ppm.hdf'


data = hdfread([hfile(1:end-4),'_lpqr.hdf'],hfile(1:end-4),'Fields','d','firstrecord',1);
% D = zeros(1,3,numPoints);
D(:,:,:) = shiftdim([data{:}],-1);

mean(squeeze(D(:,1,:)))
mean(squeeze(D(:,2,:)))
mean(squeeze(D(:,3,:)))

sSq = squeeze(sum(D.^2,2));
hf = figure,
nhist(sSq,200)
xlabel('\Lambda_1^2+\Lambda_2^2+\Lambda_3^2')
ylabel('PDF')
set(gca,'yscale','log')
saveas(hf,['pdf_s2_',hfile(1:end-4)],'fig')% Plot histogram of the eigenvalues

hf = figure,
nhist(sSq./mean(sSq),200)
xlabel('\Lambda_1^2+\Lambda_2^2+\Lambda_3^2/\langle \Lambda_1^2+\Lambda_2^2+\Lambda_3^2 \rangle')
ylabel('PDF')
set(gca,'yscale','log')
saveas(hf,['pdf_s2_normalized_',hfile(1:end-4)],'fig')% Plot histogram of the eigenvalues