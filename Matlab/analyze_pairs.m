% analysis of pairs

tau_eta = 15;
cosine = @(A,B) (sum(A.*B,2)./sqrt(sum(A.^2,2))./sqrt(sum(B.^2,2)));
R = @(x) (sum(x.^2,2).^.5);


matfiles = dir('mai10_pairs*.mat');
for m = 1:length(matfiles)
    % matfile = 'mai10_pairs_0.000100_0.000120.mat';
    matfile = matfiles(m).name;
    [~,name] = fileparts(matfile);
    
    figtitle = strrep(strrep(name,'mai10_pairs_','R0 = '),'_','-');
    
    load(matfile)
    

    % find maximum length
    
    lenPairs = zeros(length(pair),1);
    
    
    for i = 1:numel(pair)
        lenPairs(i) = length(pair(i).x);
    end
    
    [maxLen,longest] = max(lenPairs);


    % averaging
    
    t = (0:maxLen-1)/tau_eta;
    counter = zeros(maxLen,1);
    mean_lnR = zeros(maxLen,1);
    r3 = mean_lnR;
    
    for i = 1:length(pair)
        tmp = log(R(pair(i).r)./R(pair(i).r(1,:)));
        ind = 1:length(tmp);
        mean_lnR(ind) = mean_lnR(ind) + tmp;
        counter(ind) = counter(ind) + 1;
        tmp = (R(pair(i).r)./R(pair(i).r(1,:))).^(-3);
        r3(ind) = r3(ind) + tmp;
    end
    
    figure(5),
    plot(t,mean_lnR./counter)
    grid on; box on;
    xlabel('$t/t_k $','Interpreter','latex')
    ylabel('$<\ln(R(t)/R(0))>$','Interpreter','latex')
    title(figtitle);
    
    figure(6),
    plot(t,r3./counter)
    grid on; box on;
    xlabel('$t/t_k $','Interpreter','latex')
    ylabel('$(R(t)/R(0))^{-3})$','Interpreter','latex')
    title(figtitle);
    drawnow
    
    
    saveAllFigures(name,[5,6]);
    close all
end
