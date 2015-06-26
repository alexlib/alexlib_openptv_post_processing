function [B,L,N,BW] = track_enstrophy_boundary(v,plot_result)
% [B,L,N,BW] = track_enstrophy_boundary(v) shows the evolution of the
% boundary of the matrix v in time
%
% % Usage:
% load ../../water_10_53Hz_run3.mat
% track_enstrophy_boundary(Vsmooth{1})
%
% % another example:
% [B,L,N,BW] = track_enstrophy_boundary(Vsmooth{1}(500));
% [th,r] = cart2pol(abs(boundary(:,2)-85), boundary(:,1)-50);
% polar(th,r);
% view([180,90])

if nargin < 2
    plot_result = true;
end

% assuming v is a pivmat matrix of length N
for nMap = 1:length(v)
    ens = vec2scal(v(nMap),'enstrophy');
    I = ens.w'; % think of enstrophy map as an image
    BW = im2bw(I,0.1*graythresh(I)); % adaptive to the image graythreshold
    
    % [B,L,N] = bwboundaries(BW);
    % figure, imshow(BW); hold on;
    % for k=1:length(B),
    %     boundary = B{k};
    %     if(k > N)
    %         plot(boundary(:,2), boundary(:,1), 'g','LineWidth',2);
    %     else
    %         plot(boundary(:,2), boundary(:,1), 'r','LineWidth',2);
    %     end
    % end
    
    % % no holes
    % [B,L,N] = bwboundaries(BW,'noholes');
    % figure, imshow(BW); hold on;
    % for k=1:length(B),
    %     boundary = B{k};
    %     if(k > N)
    %         plot(boundary(:,2), boundary(:,1), 'g','LineWidth',2);
    %     else
    %         plot(boundary(:,2), boundary(:,1), 'r','LineWidth',2);
    %     end
    % end
    
    % remove the small objects
    BW2 = bwareaopen(BW, 10);
    % fill the holes first
    BW_filled = imfill(BW2,'holes');
    [B,L,N] = bwboundaries(BW_filled);%,'noholes');
    if plot_result
        imshow(BW_filled); hold on;
        for k=1:length(B),
            boundary = B{k};
            if(k > N)
                plot(boundary(:,2), boundary(:,1), 'g','LineWidth',2);
            else
                plot(boundary(:,2), boundary(:,1), 'r','LineWidth',2);
            end
        end
        hold off
        drawnow
    end
end


