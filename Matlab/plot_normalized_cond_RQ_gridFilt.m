dataNumVec = 8:11;
%
for dataNum = dataNumVec

    hfile = ['statistics_',int2str(dataNum),'.hdf'];
    
    feval('load',['RQ_',hfile,'.mat'],'P','R','Q');
    
    
    
    Qnorm = median(Q);
    R = R./abs(Qnorm^(3/2));
    Q = Q./abs(Qnorm);
        
    limit = 10;

    ind = R > -limit & R < limit & Q > -limit & Q < limit;
    
    
    hf = figure
    hold on
    
    [n,x] = hist3([R(ind),Q(ind)],[500 500]);
    n = filter2(fspecial('average',3),n,'same');
    n(n < 0.01) = NaN;
    [nr,xr] = hist(R,x{1});
    [nq,xq] = hist(Q,x{2});
    n = n./(nr'*nq);

%     [cs,h] = contourf(x{1},x{2},n,logspace(-4,-1,10)); % ,logspace(-3,1,10));
    [cs,h] = contourf(x{1},x{2},n'); 
    % clabel(cs,h,'color','k','rotation',0,'labelspacing',300)
    grid on
    set(gca,'xlim',[-limit/2 limit/2], 'ylim',[-limit/2 limit/2]);
    colormap(1-gray)
    daspect([1 1 1]);

    Qbound = [-limit/2:-1,0,-1:-1:-limit/2];
    Rbound = sqrt(-4/27*Qbound.^3);
    Rbound(1:(length(Qbound)+1)/2) = -1 * Rbound(1:(length(Qbound)+1)/2);
    plot(Rbound, Qbound, 'k-','linewidth',2);

    xlabel('$R/\langle Q \rangle ^{3/2}$','Interpreter','Latex');
    ylabel('$Q/\langle Q \rangle$','Interpreter','Latex')
    hold off
    saveas(hf,['RQ_normalized_conditional',hfile(1:end-4)],'fig')
    
end