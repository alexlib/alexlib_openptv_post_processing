% analysis of pairs

tau_eta = 1; % tau_eta = 15;
cosine = @(A,B) (sum(A.*B,2)./sqrt(sum(A.^2,2))./sqrt(sum(B.^2,2)));
R = @(x) (sum(x.^2,2).^.5);
angle = @(A,B) (atan2(norm(cross(A,B)),dot(A,B)));

basename = 'mai10_pairs_';
combine_figures = false;

matfiles = dir([basename,'*.mat']);
for m = 1:1 % length(matfiles)
    % matfile = 'mai10_pairs_0.000100_0.000120.mat';
    matfile = matfiles(m).name;
    [~,name] = fileparts(matfile);
    
    figtitle = strrep(strrep(name,basename,'R0 = '),'_','-');
    
    load(matfile)
    
    % 3d plot of pairs motion
    %     figure; hold on
    %     for i = 1:numel(pair)
    %         plot3(pair(i).x(:,1),pair(i).x(:,2),pair(i).x(:,3));
    %         plot3(pair(i).x(1,1),pair(i).x(1,2),pair(i).x(1,3),'ro');
    %         plot3(pair(i).x(end,1),pair(i).x(end,2),pair(i).x(end,3),'k.');
    %     end
    %     axis on; grid on; view(3)
    %     title(figtitle);
    
    
%     hf = figure;
%     hold on
    hf1 = figure;
    hold on
    for i = 1:length(pair)
        t = 0:size(pair(i).r,1)-1;
        t = t/tau_eta;
        tmp = log(R(pair(i).r)./R(pair(i).r(1,:)));
        tmp2 = R(pair(i).r).^(-3)./R(pair(i).r(1,:)).^(-3); 
        tmp3 = angle(pair(i).r.',pair(i).u.');
        % figure(hf);
        % plot(diff(tmp,2))
        if all(abs(diff(tmp,2)) < 1)
            plot(t,tmp3);
            % plot(t,tmp3);
        end
        
        
        
        %         if max(abs(diff(tmp))) > .5
        %             plot(t,tmp,'r');
        %         else
%         if sum(diff(tmp)<0) < 5
%             figure(hf)
%             plot(t, tmp, 'b-');
%         else
%             figure(hf1);
%             % plot(t,tmp,'c:');
%                 quiver3(pair(i).x(:,1),pair(i).x(:,2),pair(i).x(:,3),...
%         pair(i).u(:,1),pair(i).u(:,2),pair(i).u(:,3),'b');
%     quiver3(pair(i).x(:,1),pair(i).x(:,2),pair(i).x(:,3),...
%         pair(i).r(:,1),pair(i).r(:,2),pair(i).r(:,3),'r');
%     plot3(pair(i).x(1,1),pair(i).x(1,2),pair(i).x(1,3),'bo');
%         end
        %         end
    end
    grid on; box on;
    xlabel('$t/t_k $','Interpreter','latex')
    ylabel('$\ln(R(t)/R(0))$','Interpreter','latex')
    title(figtitle);
    
    
    %     figure;
    %     hold on
    %     for i = 1:length(pair)
    %         t = 1:size(pair(i).r,1);
    %         t = t/tau_eta;
    %         plot(t, ( R(pair(i).r)./R(pair(i).r(1,:)) ).^(-3) );
    %     end
    %     grid on; box on;
    %     xlabel('$t/t_k $','Interpreter','latex')
    %     ylabel('$(R(t)/R(0))^{-3})$','Interpreter','latex')
    %     title(figtitle);
    
    
    
    % cosines = nan(length(pair),1);
        angles = nan(length(pair),1);
        for i = 1:length(pair)
            if length(pair(i).u(:,1)) >= tau_eta
                % cosines(i) = cosine(pair(i).r(1,:),pair(i).u(tau_eta,:));
                angles(i) = angle(pair(i).r(1,:),pair(i).u(tau_eta,:));
            end
        end
    %
    
    %     figure,
    %     nhist(cosines,51)
    %     % smoothaxis
    %     set(gca,'xlim',[-1,1])
    %     xlabel('$\cos(R,u)$','Interpreter','latex')
    %     ylabel('PDF')
    %     title(figtitle);
    
%         figure,
%         nhist(angles,51)
%         % smoothaxis
%         set(gca,'xlim',[0,pi])
%         xlabel('$\theta(R,u) [rad]$','Interpreter','latex')
%         ylabel('PDF')
%         title(figtitle);
%     
    
    % dlmwrite('cosines.txt',length(cosines),'-append');
    % dlmwrite('cosines.txt',cosines,'-append');
    % dlmwrite('cosines.txt',-999,'-append');
    
    
    %
    % close all
end

if combine_figures
    saveAllFigures(basename);
    close all
    combining_figures(basename);
end
