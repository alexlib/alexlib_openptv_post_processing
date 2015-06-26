function out = calculate_integral_overgrid(dataNumVec,quantityName,varargin)
%CALCULATE_INTEGRAL_OVERGRID calculates the integral over volume
%of the quantity QUANTITYNAME
% quantity. 
%   CALCULATE_INTEGRAL_OVERGRID(DATANUMVEC,QUANTITYNAME)
%
%   inputs:
%      dataNumVec    - Vector of data set numbers:
%         Experimental cases:
%         1 - water, disks, 29-Sep-04
%         2 - 20 ppm, disks, 29-Sep-04
%         3 - 50 ppm, disks, 29-Sep-04
%         4 - 100 ppm, disks, 29-Sep-04, 
%         5 - water, center, 29-Sep-04
%         6 - 20 ppm, center, 29-Sep-04
%         7 - 50 ppm, center, 29-Sep-04
%         8 - water, center, 01-Apr-04
%         9 - 20ppm, center ,01-Apr-04
%         10 - water, center, baffles, 02-Nov-04
%         11 - 20ppm, center, baffles, 02-Nov-04
%         12 - 50ppm, center, baffles, 02-Nov-04
%         13 - 100ppm, center, baffles, 02-Nov-04

%
%
%
%      quantityName  - Name of the quantity to plot, e.g., 'mea.D1', or
%      'mea.absu', or 'meaCor.prod'
%
%
%   outputs:
%           No outputs
%
%   example:
%
%       prepare_PDF_comparison([1:3,5:12],'mea.D1');
%
%   notes:
%
%   See also HELP 
%

% Created: 16-Jan-2005
% Author: Alex Liberzon 
% E-Mail : liberzon@ihw.baug.ethz.ch 
% Phone : +41 (0)1 633 3754 
% Copyright (c) 2005 IHW, ETH Zurich 
%
%
% $Revision: 1.0 $  $Date: 16-Jan-2005 21:31$ 
%

% Static list of quantities
% names = {'water, disks','20ppm, disks', '50ppm, disks','100ppm, disks',...
%     'water, 29.09','20ppm, 29.09', '50ppm, 29.09', ...
%     'water, 01.04','20ppm, 01.04',...
%     'water, baffles','20ppm, baffles','50ppm, baffles'};
% style = {{'-o','MarkerFaceColor','r','Color','r'},{'-s','MarkerFaceColor','w','Color','r'},... % 1-2
%     {'-^','MarkerFaceColor',[.8 .8 .8],'Color','r'},{'-d','MarkerFaceColor',[.5 .5 .5],'Color','r'},... % 3-4
%     {'--o','MarkerFaceColor','b','Color','b'},{'--s','MarkerFaceColor','w','Color','b'},... %5-6
%     {'--^','MarkerFaceColor',[.8 .8 .8],'Color','b'},...% 7
%     {'-.o','MarkerFaceColor',[0 .5 0],'Color',[0 0.5 0]},{'-.s','MarkerFaceColor','w','Color',[0 .5 0]},... % 8 -9
%     {':o','MarkerFaceColor','k','Color','k'},{':s','MarkerFaceColor','w','Color','k'},... % 10 - 11
%     {':^','MarkerFaceColor',[.8 .8 .8],'Color','k'}}; % 12
% 
% nBins = 21;


% hf = figure;
% hold on
out = zeros(dataNumVec,1);

for dataNum = dataNumVec  % usually it is [1:3,5:12]
    load(['meanField_',int2str(dataNum),'.mat']);
    load(['meanCorField_',int2str(dataNum),'.mat']);

%     figure(hf);
%     feval('nhist',eval(quantityName),nBins,style{:,dataNum}{:},'DisplayName',names{dataNum})
x1 = unique(mea.x);
y1 = unique(mea.y);
z1 = unique(mea.z);
tmp = feval('reshape',eval(quantityName),[length(z1) length(y1) length(x1)]);
tmp(isnan(tmp)) = 0;

    out(dataNum,1) = trapz(x1,trapz(y1,trapz(z1,tmp)))/((max(x1)-min(x1))*(max(y1)-min(y1))*(max(z1)-min(z1)));
end

% figure(hf)
% xlabel(quantityName,'Interpreter','none')
% smoothaxis
% hl = legend(gca,'toggle'); 
% set(gcf,'Name',[quantityName,' ',mat2str(dataNumVec)])
% saveas(gcf,['pdf_',quantityName,' ',mat2str(dataNumVec),'.fig'])

if nargin > 2
    disp(strrep(mat2str(out'*varargin{1},3),' ',' & '))
end


