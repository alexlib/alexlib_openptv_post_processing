function AlexFromMicheleAfterBeat(sixinone,first,last)

if ~nargin
    sixinone = 0;
    first = 1
    last = 6
end

if sixinone == 1
    load(['meanField']);
    load(['meanCorField']);
    hfile = ['meanCorField','.hdf']    
else
    for loop = first:last     
        run = 1;
        load(['meanField_',int2str(run),'_',int2str(loop)]);
        load(['meanCorField_',int2str(run),'_',int2str(loop)]);
        hfile = ['meanCorField_',int2str(run),'_',int2str(loop),'.hdf']
    
    
    
    
    [data] = readStatisticsHDF(hfile,'lambda1',[],'lambda2',[],'lambda3',[],'ii',[],'sss',[],'wws',[],'s2',[],'w2',[],...
        'cosss',[],'cosww',[],'coswl1',[],'coswl2',[],'coswl3',[]);
    
    
    threshold = 1e-6;
    indexp = find(meaCor.production > threshold);
    indexn = find(meaCor.production < -threshold);
    
    %%%%%%%% ISMEMBER IS THE ANSWER %%%%%%
    ipos = ismember(cat(1,data.ii),indexp);
    ineg = ismember(cat(1,data.ii),indexn);
    
    lambda_1 = cat(1,data.lambda1);
    lambda_2 = cat(1,data.lambda2);
    lambda_3 = cat(1,data.lambda3);
    coswW = cat(1,data.cosww);
    cossS = cat(1,data.cosss);
    ssq = cat(1,data.s2);
    wsq = cat(1,data.w2);
    sssok = cat(1,data.sss);
    wwsok = cat(1,data.wws);
    coswlam1 = cat(1,data.coswl1);
    coswlam2 = cat(1,data.coswl2);
    coswlam3 = cat(1,data.coswl3);
    
    ratesss = -sssok./ssq;
    ratewws = wwsok./wsq;
    lam1cube = lambda_1.^3;
    lam2cube = lambda_2.^3;
    lam3cube = lambda_3.^3;
    
    
    
     %%%%%%%% ISMEMBER IS THE ANSWER %%%%%%
    valuable_ssq = ssq < prctile(ssq,99);
    ipos = ipos & valuable_ssq;
    ineg = ineg & valuable_ssq;

    
    
    %figures%%%%%%%%%%%%%%%%
    nbins = 20
    
    
    figure
    nhist(lambda_1(ipos),nbins,'k');
    hold on
    nhist(lambda_2(ipos),nbins,'k.-');
    nhist(lambda_3(ipos),nbins,'k--');
    nhist(lambda_1(ineg),nbins,'b');
    nhist(lambda_2(ineg),nbins,'b.-');
    nhist(lambda_3(ineg),nbins,'b--');
    grid on
    axis([-2 2 0 3])
    legend('\Lambda_1 +','\Lambda_2 +','\Lambda_3 +','\Lambda_1 -','\Lambda_2 -','\Lambda_3 -');
    xlabel('[s^{-1}]','FontSize',12);
    ylabel('pdf','FontSize',12);
    saveas(gcf,['Lambda_',int2str(run),'_',int2str(loop)]);
    
    
    
    figure
    nhist(lam1cube(ipos),nbins,'k');
    hold on
    nhist(lam2cube(ipos),nbins,'k.-');
    nhist(lam3cube(ipos),nbins,'k--');
    nhist(lam1cube(ineg),nbins,'b');
    nhist(lam2cube(ineg),nbins,'b.-');
    nhist(lam3cube(ineg),nbins,'b--');
    grid on
    %axis([-2 2 0 3])
    legend('\Lambda_1^3 +','\Lambda_2^3 +','\Lambda_3^3 +','\Lambda_1^3  -','\Lambda_2^3 -','\Lambda_3^3 -');
    xlabel('[s^{-3}]','FontSize',12);
    ylabel('pdf','FontSize',12)
    saveas(gcf,['LambdaCube_',int2str(run),'_',int2str(loop)]);
    
    
    figure
    grid on
    nhist(coswW(ipos),nbins,'k');
    hold on
    nhist(coswW(ineg),nbins,'b');
    legend('cos(W,\omega)+','cos(W,\omega) -');
    ylabel('pdf','FontSize',12);
    saveas(gcf,['coswW_',int2str(run),'_',int2str(loop)]);

    
    figure
    nhist(cossS(ipos),nbins,'k.-');
    hold on
    nhist(cossS(ineg),nbins,'b.-');
    grid on
    %axis([-1 1 0 1])
    legend('cos(s,S)+','cos(s,S) -');
    %xlabel('[s^{-1}]','FontSize',12);
    ylabel('pdf','FontSize',12);
   saveas(gcf,['cossS_',int2str(run),'_',int2str(loop)]);

    
    figure
    nhist(-sssok(ipos),nbins,'k');
    hold on
    nhist(wwsok(ipos),nbins,'k.-');
    nhist(-sssok(ineg),nbins,'b');
    nhist(wwsok(ineg),nbins,'b.-');
    grid on
    %axis([-1 1 0 1])
    legend('-s_{ij}s_{jk}s_{ki} +','\omega_{i} \omega_{j}s_{ij} +','-s_{ij}s_{jk}s_{ki} -','\omega_{i} \omega_{j}s_{ij} -');
    xlabel('[s^{-3}]','FontSize',12);
    ylabel('pdf','FontSize',12);
    saveas(gcf,['sss_',int2str(run),'_',int2str(loop)]);

    
    figure
    nhist(ratesss(ipos),nbins,'k');
    hold on
    nhist(ratewws(ipos),nbins,'k.-');
    nhist(ratesss(ineg),nbins,'b');
    nhist(ratewws(ineg),nbins,'b.-');
    grid on
    %axis([-1 1 0 1])
    legend('-s_{ij}s_{jk}s_{ki}/s^2 +','\omega_{i} \omega_{j}s_{ij}/\omega^2 +','-s_{ij}s_{jk}s_{ki}/s^2 -','\omega_{i} \omega_{j}s_{ij}/\omega^2 -');
    xlabel('[s^{-1}]','FontSize',12);
    ylabel('pdf','FontSize',12);
    saveas(gcf,['sss_s2_',int2str(run),'_',int2str(loop)]);
    
    
    figure
    nhist(ssq(ipos),nbins,'k');
    hold on
    nhist(wsq(ipos),nbins,'k.-');
    nhist(ssq(ineg),nbins,'b');
    nhist(wsq(ineg),nbins,'b.-');
    grid on
    %axis([-1 1 0 1])
    legend('s^2 +','\omega^2 +','s^2 -','\omega^2 -');
    xlabel('[s^{-2}]','FontSize',12);
    ylabel('pdf','FontSize',12);
    saveas(gcf,['s2_',int2str(run),'_',int2str(loop)]);
    
    figure
    nhist(coswlam1(ipos),nbins,'k');
    hold on
    nhist(coswlam2(ipos),nbins,'k.-');
    nhist(coswlam3(ipos),nbins,'k--');
    nhist(coswlam1(ineg),nbins,'b');
    nhist(coswlam2(ineg),nbins,'b.-');
    nhist(coswlam3(ineg),nbins,'b--');
    grid on
    %axis([-1 1 0 3])
    legend('cos(\omega,\lambda_1) +','cos(\omega,\lambda_2) +','cos(\omega,\lambda_3) +','cos(\omega,\lambda_1) -','cos(\omega,\lambda_2) -','cos(\omega,\lambda_3) -');
    xlabel('[s^{-1}]','FontSize',12);
    ylabel('pdf','FontSize',12);
    saveas(gcf,['coswLambda_',int2str(run),'_',int2str(loop)]);
    
%     close all
    end
end