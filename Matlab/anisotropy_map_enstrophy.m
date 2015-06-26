% function anisotropy_map_enstrophy(dataNumVec)
% % extended version of plot_anisotropy_map_convection
% % that adds the anisotropy number to the plot
% % and makes also anisotropy maps of w_iw_j and s_{ik}s_{kj}
% %
%
% if ~nargin,
dataNumVec = [8,9,10,11];
quantityName = {'\omega','s'};
for dataNum = dataNumVec

    for i = 1:2
        hf(i) = figure;
        axis([-.04 .08 0 .3]);
        hold on, grid on, box on
        plot_Lumley_triangle(hf(i));
        xlabel('III')
        ylabel('II')
    end

    hfile = ['statistics_',int2str(dataNum),'.hdf'];
    finfo = hdfinfo(hfile);
    len = [finfo.SDS.Dims.Size];
    disp('reading ...')
    gridData = hdfread(hfile,'gridData','Index', {[1.0 4.0 1.0 ],[1.0 1.0 1.0 ],[len(1) 9.0 len(3)]});
    disp('done')



    ux = squeeze(gridData(:,1,:));
    uy = squeeze(gridData(:,2,:));
    uz = squeeze(gridData(:,3,:));

    vx = squeeze(gridData(:,4,:));
    vy = squeeze(gridData(:,5,:));
    vz = squeeze(gridData(:,6,:));

    wx = squeeze(gridData(:,7,:));
    wy = squeeze(gridData(:,8,:));
    wz = squeeze(gridData(:,9,:));

    clear gridData

    w1 = wy - vz;
    w2 = uz - wx;
    w3 = vx - uy;

    s11 = ux;
    s12 = 0.5*(uy + vx);
    s13 = 0.5*(uz + wx);
    s22 = vy;
    s23 = 0.5*(vz + wy);
    s33 = wz;

    [ww,ss] = deal(zeros(3,3,len(1),len(3)));

    for i = 1:len(1)
        for j = 1:len(3)
            ww(:,:,i,j) = [w1(i,j) w2(i,j) w3(i,j)]'*[w1(i,j) w2(i,j) w3(i,j)];
            ss(:,:,i,j) = [...
                s11(i,j) s12(i,j) s13(i,j);...
                s12(i,j) s22(i,j) s23(i,j);...
                s13(i,j) s23(i,j) s33(i,j)]^2;
        end
    end
    ww = squeeze(mean(ww,3));
    ss = squeeze(mean(ss,3));


   
    [ii,iii] = deal(cell(2,len(3)));

    for i = 1:len(3)
        [ii{1,i},iii{1,i}] = anisotropy(ww(:,:,i));
        %
        [ii{2,i},iii{2,i}] = anisotropy(ss(:,:,i));
    end

    for i = 1:2
        figure(hf(i));
        scatter([iii{i,:}],-[ii{i,:}],'r.');
        AIM = mean(...
            ([iii{i,:}]/2./(-[ii{i,:}]/3).^(3/2)));
        title(sprintf('$\\langle A_{%s} \\rangle = $ %3.1f',quantityName{i},AIM),'interpreter','latex');
        saveas(hf(i),['anisotropy_smallscale_',int2str(dataNum),'_',int2str(i),'.fig'],'fig');
    end
    
    clear ww ss

end
