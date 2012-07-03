uiopen('/Users/alex/Downloads/drag_force_first.fig',1)
% select the blue line
hl = sort(findobj(gcf,'type','line','Color',[0 0 1]));

% manually set up the threshold
thresh = 4e-10;


% select the data of the line 
xd = get(hl,'xdata'); yd = get(hl,'ydata');

% select the circle that belongs to that line:
xde = get(hl+1,'xdata'); yde = get(hl+1,'ydata');

% plot again to be sure it's the data you want
figure, plot(xd,yd,xde,yde,'o')
% choose threshold

figure, plot(xd,yd,xde,yde,'o',[xd(1), xd(end)],[thresh thresh],'r--')

% select the data for the integral:
% a) above the line
% b) before the event
ind = yd > thresh & xd <= xde;
figure, plot(xd,yd,xde,yde,'o',[xd(1), xd(end)],[thresh thresh],'r--',xd(ind),yd(ind),'m.')

% subtract the threshold and get the area under the curve
yd = yd - thresh;
area = trapz(xd,yd);
figure, plot(xd(ind),yd(ind),xde,yde,'o')
legend(sprintf('Area = %6.4f',area))


