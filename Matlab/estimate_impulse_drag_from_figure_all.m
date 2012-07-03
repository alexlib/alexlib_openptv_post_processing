uiopen('drag_110_25.fig',1);
% uiopen('/Users/alex/Downloads/drag_force_first.fig',1)
% select the blue line
hl = sort(findobj(1,'type','line','Color',[0 0 1])); % all the blue lines
ho = sort(findobj(gcf,'type','line','Color',[1 0 0])); % all the red dots and the red line

% find out the red line
for i = 1:length(ho)
    yd = get(ho(i),'YData');
    if length(yd) > 1
        thresh = yd(1);
        lineInd = i;
    end
end

ho(lineInd) = [];


%
%
% if there is no line, then comment out the above portion
% manually set up the threshold
thresh = 0.4e-6;


% sometimes it is best to test the way we read the file by
% plotting it again:

%figure, hold on;
 %for i = 1:length(hl), plot(get(hl(i),'xdata'), get(hl(i),'ydata'),'b-'); plot(get(ho(i),'xdata'),get(ho(i),'ydata'),'ro'); end


[usquare,area] = deal(zeros(size(hl)));

for i = 1:length(hl) % 1 is the red line 
    
% select the data of the line 
xd = get(hl(i),'xdata'); yd = get(hl(i),'ydata');

% select the circle that belongs to that line:
xde = get(ho(i),'xdata'); yde = get(ho(i),'ydata');

% let's store the events:
usquare(i) = yde;

% plot again to be sure it's the data you want
% figure, plot(xd,yd,xde,yde,'o')

% figure, plot(xd,yd,xde,yde,'o',[xd(1), xd(end)],[thresh thresh],'r--')

% select the data for the integral:
% a) above the line
% b) before the event
ind = yd > thresh & xd <= xde;
% figure, plot(xd,yd,xde,yde,'o',[xd(1), xd(end)],[thresh thresh],'r--',xd(ind),yd(ind),'m.')

% subtract the threshold and get the area under the curve
yd = yd - thresh;
area(i) = trapz(xd,yd)/160;
set(hl(i),'DisplayName',sprintf('%6.4f',area(i)));

% figure, plot(xd(ind),yd(ind),xde,yde,'o')
% legend(sprintf('Area = %6.4f',area))
end

figure, plot(usquare,area,'.');
xlabel('drag force');
ylabel('impulse');




