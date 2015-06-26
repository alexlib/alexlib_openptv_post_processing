% function [P,R,Q,hf] = calculate_PRQ(traj,varargin)
% CALCULATE_PRQ  calculates and saves R,Q for the later use.
%
% [P,R,Q] = CALCULATE_RQ(TRAJ)
%
% Inputs:
%		TRAJ - structure with at least 6 fields: s11,s12,...
% Outputs
%       P,R,Q - 3 invariants of the 3x3 gradient tensor du/dx

% Author: Alex Liberzon
% Copyright (c) 2004, IHW, ETHZ
% Last modified: July 1, 2004

dataNumVec = 10:11;
%
for dataNum = dataNumVec

    hfile = ['statistics_',int2str(dataNum),'.hdf'];
    finfo = hdfinfo(hfile);
    len = [finfo.SDS.Dims.Size];
    disp('reading ...')
    gridData = hdfread(hfile,'gridData','Index', {[1.0 4.0 1.0 ],[1.0 1.0 1.0 ],[len(1) 9.0 len(3)]});
    load(['meanField_',int2str(dataNum),'.mat']);
    load(['meanCorField_',int2str(dataNum),'.mat']);
    disp('done')
    
    numPoints = len(1)*len(3);
    
    tmp = reshape(permute(gridData,[2 1 3]),9,[]);
    
    reldiv = hdfread(hfile,'gridData','Index', {[1.0 len(2) 1.0 ],[1.0 1.0 1.0 ],[len(1) 1.0 len(3)]});
    
    reldiv = reshape(reldiv(:,1,:),1,[]);
    ind = find(reldiv <= 0.1);
    [P,R,Q] = deal(repmat(0,[length(ind) 1]));
    
    m = 0;
    for i = 1:length(ind)
        du = reshape(tmp(1:9,ind(i)),[3 3])';
        s = 0.5*(du + du');
        r = 0.5*(du - du');
        P(i) = trace(du);
        P(i) = 0;
        Q(i) = 0.5*(P(i)^2 - trace(du*du));
        R(i) = 1/3*(-P(i)^3 + 3*P(i)*Q(i) - trace(du^3));
%         if ~mod(i,500), fprintf(1,'.'); end
%         if ~mod(i,50000), fprintf(1,'\n'); end
    end

    feval('save',['RQ_',hfile,'.mat'],'P','R','Q');

    limit = 10;

    ind = R > -limit & R < limit & Q > -limit & Q < limit;


    hf = figure
    hold on
    [n,x,bins] = histmulti5([R(ind),Q(ind)]);
    % contour(x(:,1),x(:,2),n')
    plot(R(ind(1:50:end)),Q(ind(1:50:end)),'.','markersize',0.5)
    [cs,h] = contour(x(:,1),x(:,2),filter2(fspecial('average',7),n','same'),logspace(-3,1,10));
    clabel(cs,h,'color','k','rotation',0,'labelspacing',300)
    grid on
    set(gca,'xlim',[-limit limit], 'ylim',[-limit limit]);
    daspect([1 1 1]);

    Qbound = [-limit:-1,0,-1:-1:-limit];
    Rbound = sqrt(-4/27*Qbound.^3);
    Rbound(1:(length(Qbound)+1)/2) = -1 * Rbound(1:(length(Qbound)+1)/2);
    plot(Rbound, Qbound, 'k-');

    xlabel('R'),ylabel('Q')
    hold off
    saveas(hf,['RQ_',hfile(1:end-4)],'fig')

    % [EOF] calculate_eigens


    
    hf = figure
    hold on
    [n,x,bins] = histmulti5([R(ind),Q(ind)]);
    % bins = linspace(min(x(1,:)),max(x(1,:)),50);
    % [n,x,bins] = histmulti5([R(ind),Q(ind)],[linspace(min(x(:,1)),max(x(:,1)),50)',linspace(min(x(:,2)),max(x(:,2)),50)']);
    n = filter2(fspecial('average',7),n','same');
    n(n < .0001) = NaN;
    % contour(x(:,1),x(:,2),n')
    % plot(R(ind(1:50:end)),Q(ind(1:50:end)),'.','markersize',0.5)
    [cs,h] = contourf(x(:,1),x(:,2),n,logspace(-3,-1,10)); % ,logspace(-3,1,10));
    % clabel(cs,h,'color','k','rotation',0,'labelspacing',300)
    grid on
    set(gca,'xlim',[-limit limit], 'ylim',[-limit limit]);
    colormap(1-gray)
    daspect([1 1 1]);

    Qbound = [-limit:-1,0,-1:-1:-limit];
    Rbound = sqrt(-4/27*Qbound.^3);
    Rbound(1:(length(Qbound)+1)/2) = -1 * Rbound(1:(length(Qbound)+1)/2);
    plot(Rbound, Qbound, 'k-','linewidth',2);

    xlabel('R'),ylabel('Q')
    hold off
    saveas(hf,['RQ_',hfile(1:end-4)],'fig')

end % dataNum