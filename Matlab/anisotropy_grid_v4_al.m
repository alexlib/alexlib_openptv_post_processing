% x_min=0.005;
% x_max=0.053;
% x_step=0.003;
% y_min=0.005;
% y_max=0.053;
% y_step=0.003;
% z_min=-0.015;
% z_max=0.012;
% z_step=0.003;

%path='Y:\rotating\080811\run1';
% path='..\run1';
% cd(path)

%mov1 = avifile('vel_quiver.avi','fps',2,'compression','None');
%mov2 = avifile('w_contour.avi','fps',2,'compression','None');

for i = 1:1
    hf(i) = figure(3);
    axis([-.04 .08 0 .35]);
    hold on, grid on, box on
    plot_Lumley_triangle(hf(i));
    xlabel('III')
    ylabel('II')
end

for timestep=800:800 %:5:1500 %timesteps

    filename_in=strcat('X:/rotating/PTV08/080811/run1/grid_filt.',num2str(timestep));   
    D=load(filename_in);

    last_pos=length(D(:,1));

    %         x_min=D(1,1);
    %         x_max=D(last_pos,1);
    %         y_min=D(1,2);
    %         y_max=D(last_pos,2);
    %         z_min=D(1,3);
    %         z_max=D(last_pos,3);
    %         z_step=D(2,3)-D(1,3);
    %         x_step=z_step;
    %         y_step=z_step;
    %
    %         s1=(x_max-x_min)/x_step+1;
    %         s2=(y_max-y_min)/y_step+1;
    %         s3=(z_max-z_min)/z_step+1;
    %         j=1;
    %         for i=1:s1
    %             for ii=1:s2
    %                 for iii=1:s3
    %                     x(i,ii,iii)=D(j,1);
    %                     y(i,ii,iii)=D(j,2);
    %                     z(i,ii,iii)=D(j,3);
    %                     u(i,ii,iii)=D(j,4);
    %                     v(i,ii,iii)=D(j,5);
    %                     w(i,ii,iii)=D(j,6);
    %     %                 w1(i,ii,iii)=D(j,10);
    %     %                 w2(i,ii,iii)=D(j,11);
    %     %                 w3(i,ii,iii)=D(j,12);
    %
    %                     j=j+1;
    %                 end
    %             end
    %         end

    x = permute(reshape(D(:,1),10,17,17),[3 2 1]);
    y = permute(reshape(D(:,2),10,17,17),[3 2 1]);
    z = permute(reshape(D(:,3),10,17,17),[3 2 1]);

    u = permute(reshape(D(:,4),10,17,17),[3 2 1]);
    v = permute(reshape(D(:,5),10,17,17),[3 2 1]);
    w = permute(reshape(D(:,6),10,17,17),[3 2 1]);


    for i = 1:5 % size(y,2) % 17 different levels of y
        
        
        uu = nanmean(nanmean(squeeze(u(:,i,:).*u(:,i,:))));
        vv = nanmean(nanmean(squeeze(v(:,i,:).*v(:,i,:))));
        ww = nanmean(nanmean(squeeze(w(:,i,:).*w(:,i,:))));
        uv = nanmean(nanmean(squeeze(u(:,i,:).*v(:,i,:))));
        uw = nanmean(nanmean(squeeze(u(:,i,:).*w(:,i,:))));
        vw = nanmean(nanmean(squeeze(v(:,i,:).*w(:,i,:))));

        Rs = [uu,uv,uw;uv,vv,vw;uw,vw,ww];

        [II,III] = anisotropy(Rs);
        figure(hf(1));
        scatter(III,-II,'.','DisplayName',num2str(y(1,i,1)));drawnow;
    end
end
for i=1:5
    [i y(1,i,1) ]
end

%     timestep
%
%
%     for i = 1:length(D)
%         u = D(i,4);
%         v = D(i,5);
%         w = D(i,6);
% %         i, u,v,w
%         if ~isnan(u) && ~isnan(v) && ~isnan(w)
%             Rs = [u*u,u*v,u*w;
%                 u*v,v*v,v*w;
%                 u*w,v*w,w*w];
%             [II,III] = anisotropy(Rs);
%             figure(hf(1));
%             plot(III,-II,'r.'); drawnow;
%         end
%
%     end
% end

%
%
%     for k = 2:3:size(x2,2)-2
% %         [ux,uz] = gradient(squeeze(u1(:,k,:)));
% %         [wx,wz] = gradient(squeeze(u3(:,k,:)));
% %         pcolor(uz - wx);
%
%         h = streamslice(squeeze(u1(:,k,:)),squeeze(u3(:,k,:)),[],'noarrows',[],[],10);
%         for i=1:length(h);
%             %         zi = interp2(unique(x2(:,k,:)),get(h(i),'xdata'),get(h(i),'ydata'));
%             set(h(i),'zdata',repmat(unique(x2(:,k,:)),size([get(h(i),'xdata')])));
%         end
%     end
%     view(36,14); axis square, axis off
%     if framestep == 1, camzoom(1.35), end
%     M(framestep) = getframe;
%     framestep = framestep + 1;
%     cla
% end
%
%
% save tmp.mat M

