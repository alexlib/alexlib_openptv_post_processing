% Analysis of Lamb, helicity and additional quantities
% version 2, modified on
%-- 5/9/05  9:59 PM --%

dataNumVec = [6,7];

names = {'1','2', '3','4',...
    '5','6', '7', ...
    };

style = {{'-o','MarkerFaceColor','m','Color','m'},{'-s','MarkerFaceColor','w','Color','m'},... % 1-2
    {'-^','MarkerFaceColor',[.8 .8 .8],'Color','m'},{'-d','MarkerFaceColor',[.5 .5 .5],'Color','m'},... % 3-4
    {'--o','MarkerFaceColor','b','Color','k'},{'--s','MarkerFaceColor','w','Color','k'},... %5-6
    {'--^','MarkerFaceColor',[.8 .8 .8],'Color','k'},...% 7
    {'-o','MarkerFaceColor','r','Color','r'},{'-o','MarkerFaceColor','w','Color',[1 0 0]},... % 8 - 9
    {'--s','MarkerFaceColor','b','Color','b'},{'--s','MarkerFaceColor','w','Color','b'},... % 10 - 11
    {':^','MarkerFaceColor',[.8 .8 .8],'Color','m'},{':d','MarkerFaceColor',[0.5 0.5 0.5],'Color','m'},...
    {'-.s','MarkerFaceColor','m','Color','m'},{'-.s','MarkerFaceColor','w','Color','m'},...% 14-15 - 10/11 A
    {'-.s','MarkerFaceColor','w','Color','m'}}; % 16 = 11B

nBins = 50;

titles = {'','$P>0$','$P<0$'};
for j = 1:4
    hf(j) = figure;
    for i = 1:3
        subplot(3,1,i); hold on; grid on;title(titles{i});
    end
end

% style = {'r-o','b--s','k-.^'};

% table = zeros(length(dataNumVec),9);
table2 = zeros(length(dataNumVec),3);
table3 = zeros(length(dataNumVec),3);
table4 = zeros(length(dataNumVec),3);

m = 0;
for dataNum = dataNumVec
    m = m + 1;
    load(['LH_',int2str(dataNum)]);
    load(['production_',int2str(dataNum)],'cos_Rs_S');

    [posIndx,negIndx] = deal(zeros(size(uNu)));
    posIndx(:,cos_Rs_S > .1) = 1;
    posIndx = find(posIndx);
    negIndx = zeros(size(uNu));
    negIndx(:,cos_Rs_S < -.1) = 1;
    negIndx = find(negIndx);

%     for i = 1:3
%         figure(hf(i))
%         subplot(3,1,1);
%         nhist(L(:,i),20,style{:,dataNum}{:},'DisplayName',names{dataNum});
%         subplot(3,1,2);
%         nhist(L(posIndx,i),20,style{:,dataNum}{:},'DisplayName',names{dataNum});
%         subplot(3,1,3);
%         nhist(L(negIndx,i),20,style{:,dataNum}{:},'DisplayName',names{dataNum});
%     end


%     table(m,1:3)  = trimmean(L,5,1);
%     table(m,4:6) = trimmean(L(posIndx,:),5,1);
%     table(m,7:9) = trimmean(L(negIndx,:),5,1);
    
    table3(m,1) = trimmean(LL,5,1);
    table3(m,2) = trimmean(LL(posIndx),5,1);
    table3(m,3) = trimmean(LL(negIndx),5,1);

    %% helicity density
    figure(hf(1))
    subplot(3,1,1)
%     nhist(1-hh.^2,nBins,style{:,dataNum}{:},'DisplayName',names{dataNum})
    nhist(sqrt(LL.^2),nBins,style{:,dataNum}{:},'DisplayName',names{dataNum})
    subplot(3,1,2)
%     nhist(1-hh(posIndx).^2,nBins,style{:,dataNum}{:},'DisplayName',names{dataNum})
    nhist(sqrt(LL(posIndx).^2),nBins,style{:,dataNum}{:},'DisplayName',names{dataNum})

    subplot(3,1,3)
%     nhist(1-hh(negIndx).^2,nBins,style{:,dataNum}{:},'DisplayName',names{dataNum})
    nhist(sqrt(LL(negIndx).^2),nBins,style{:,dataNum}{:},'DisplayName',names{dataNum})

    
    table2(m,1) = trimmean(abs(hh),5,1);
    table2(m,2) = trimmean(abs(hh(posIndx,:)),5,1);
    table2(m,3) = trimmean(abs(hh(negIndx,:)),5,1);
    
    table4(m,1) = trimmean(hh,5,1);
    table4(m,2) = trimmean(hh(posIndx,:),5,1);
    table4(m,3) = trimmean(hh(negIndx,:),5,1);


    %% u_i s_{ij}
%     figure(hf(5))
%     subplot(3,1,1)
%     nhist(us1(:),nBins,style{:,dataNum}{:},'DisplayName',names{dataNum})
%     subplot(3,1,2)
%     nhist(us1(posIndx),nBins,style{:,dataNum}{:},'DisplayName',names{dataNum})
%     subplot(3,1,3)
%     nhist(us1(negIndx),nBins,style{:,dataNum}{:},'DisplayName',names{dataNum})
%     
%     figure(hf(6))
%     subplot(3,1,1)
%     nhist(us2(:),nBins,style{:,dataNum}{:},'DisplayName',names{dataNum})
%     subplot(3,1,2)
%     nhist(us2(posIndx),nBins,style{:,dataNum}{:},'DisplayName',names{dataNum})
%     subplot(3,1,3)
%     nhist(us2(negIndx),nBins,style{:,dataNum}{:},'DisplayName',names{dataNum})
%     
%     figure(hf(7))
%     subplot(3,1,1)
%     nhist(us3(:),nBins,style{:,dataNum}{:},'DisplayName',names{dataNum})
%     subplot(3,1,2)
%     nhist(us3(posIndx),nBins,style{:,dataNum}{:},'DisplayName',names{dataNum})
%     subplot(3,1,3)
%     nhist(us3(negIndx),nBins,style{:,dataNum}{:},'DisplayName',names{dataNum})
    
    
    figure(hf(2))
    subplot(3,1,1)
    nhist(LL(:),nBins,style{:,dataNum}{:},'DisplayName',names{dataNum})
    subplot(3,1,2)
    nhist(LL(posIndx),nBins,style{:,dataNum}{:},'DisplayName',names{dataNum})
    subplot(3,1,3)
    nhist(LL(negIndx),nBins,style{:,dataNum}{:},'DisplayName',names{dataNum})
    
    figure(hf(3))
    subplot(3,1,1)
    nhist(hh(:),nBins,style{:,dataNum}{:},'DisplayName',names{dataNum})
    subplot(3,1,2)
    nhist(hh(posIndx),nBins,style{:,dataNum}{:},'DisplayName',names{dataNum})
    subplot(3,1,3)
    nhist(hh(negIndx),nBins,style{:,dataNum}{:},'DisplayName',names{dataNum})
    
    figure(hf(4))
    subplot(3,1,1)
    nhist(abs(hh(:)),nBins,style{:,dataNum}{:},'DisplayName',names{dataNum})
    subplot(3,1,2)
    nhist(abs(hh(posIndx)),nBins,style{:,dataNum}{:},'DisplayName',names{dataNum})
    subplot(3,1,3)
    nhist(abs(hh(negIndx)),nBins,style{:,dataNum}{:},'DisplayName',names{dataNum})
    
end

% columnLabels = {'$\lange L_1 \rangle$', '$\lange L_2 \rangle$','$\lange L_3 \rangle$',...
%     '$\lange L_1 \rangle, P > 0$', '$\lange L_2 \rangle, P > 0$','$\lange L_3 \rangle, P > 0$',...
%     '$\lange L_1 \rangle, P < 0$', '$\lange L_2 \rangle, P < 0$','$\lange L_3 \rangle, P < 0$',...
%     };
rowLabels = {'1', '2','3','4','5','6','7'};
% matrix2latex(table*1e4, 'mean_lamb.tex', 'rowLabels', rowLabels, 'columnLabels', columnLabels, 'alignment', 'c', 'format', '%-6.2f', 'size', 'small');
columnLabels = {'$\lange h \rangle$', '$\lange h  \rangle , P > 0$','$\lange h \rangle, P < 0$',...
    };
matrix2latex(table4, 'mean_abs_helicity_density.tex', 'rowLabels', rowLabels, 'columnLabels', columnLabels, 'alignment', 'c', 'format', '%-6.2f', 'size', 'small');

columnLabels = {'$|h|$', '$|h | , P > 0$','$|h|, P < 0$',...
    };

matrix2latex(table2, 'mean_helicity_density.tex', 'rowLabels', rowLabels, 'columnLabels', columnLabels, 'alignment', 'c', 'format', '%-6.2f', 'size', 'small');

columnLabels = {'$\lange L \rangle$', '$\lange L  \rangle , P > 0$','$\lange L \rangle, P < 0$',...
    };
matrix2latex(table3, 'mean_Lamb_density.tex', 'rowLabels', rowLabels, 'columnLabels', columnLabels, 'alignment', 'c', 'format', '%-6.2f', 'size', 'small');


for i = 1:4
    saveas(hf(i),['analysis_lamb_helicity_',int2str(i),'.fig'],'fig')
end


